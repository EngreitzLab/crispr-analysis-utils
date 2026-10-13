"""`summarize`, the output tables, and the coordinates they report."""

import pandas as pd
import pytest
import yaml
from guide_alignment_fixtures import (
    SPACER_A,
    SPACER_C,
    SPACER_G,
    SPACER_T,
    Contig,
    write_fasta,
)

from crispr_analysis_utils.guide_alignment import (
    Site,
    fan_out_sites,
    guide_columns,
    read_guides,
    scan_guides,
    summarize,
    write_columns_yaml,
    write_guides_bed,
    write_guides_tsv,
    write_qc_tsv,
)
from crispr_analysis_utils.guide_alignment.iupac import reverse_complement
from crispr_analysis_utils.guide_alignment.summary import (
    BED_COLUMNS,
    COLUMNS,
    GROUPS,
)

SPACER_FLIPPED = "CAGTTCGAACTGGTACCATG"
SPACER_NOTHING = "GGTACCTTAGCAATCGCTAG"
SPACERS = {
    "uniq": SPACER_A,
    "uniq_minus": SPACER_C,
    "multi": SPACER_G,
    "flipped": SPACER_FLIPPED,
    "no_pam": SPACER_T,
    "nothing": SPACER_NOTHING,
    "dup": SPACER_A,
}
# One verdict per source row: 'dup' is row 6 but takes row 0's verdict.
VERDICTS = {
    0: "forward",
    1: "forward",
    2: "forward",
    3: "reverse_complement",
    4: "no_pam",
    5: "no_match",
}
EVIDENCE = {3: "ambiguous"}

COLUMN_NAMES = [
    "guide_id",
    "alias",
    "spacer",
    "site_class",
    "n_matches",
    "n_cut_sites",
    "orientation",
    "orientation_evidence",
    "chr",
    "start",
    "end",
    "strand",
    "cut_site",
    "pam_matches",
    "mapped_sequence",
]
SITE_COLUMN_NAMES = COLUMN_NAMES[8:]


def site(spacer, start, strand="+", mismatches=0, pam_class="NGG", pam="TGG"):
    """One site of the guide whose sequence is `spacer`."""
    return Site(
        guide_id=spacer + "NGG",
        chr="chr1",
        start=start,
        end=start + 20,
        strand=strand,
        pam=pam,
        pam_class=pam_class,
        n_mismatches=mismatches,
        mismatch_positions=tuple(range(1, mismatches + 1)),
        mapped_sequence=spacer + pam,
        genomic_protospacer=spacer,
        genomic_5p_base=None,
    )


@pytest.fixture
def guides():
    return read_guides(pd.DataFrame(list(SPACERS.items())))


@pytest.fixture
def sites():
    """Exact sites, as the scan reports them: every match, PAM or no PAM."""
    return {
        # Two matches of the same protospacer, only one of them with a PAM.
        0: [site(SPACER_A, 100), site(SPACER_A, 500, pam_class="none", pam="TTT")],
        1: [site(SPACER_C, 200, strand="-", pam="AGG")],
        2: [site(SPACER_G, 600), site(SPACER_G, 700, strand="-", pam="CGG")],
        3: [
            site(SPACER_FLIPPED, 1000, pam_class="none", pam="TAG"),
            site(SPACER_FLIPPED, 1100, pam_class="none", pam="TCC"),
        ],
        4: [site(SPACER_T, 800, pam_class="none", pam="TTT")],
    }


def summarized(guides, sites, orientation=VERDICTS, **kwargs):
    return summarize(
        guides,
        sites,
        orientation,
        evidence=EVIDENCE if orientation is VERDICTS else None,
        **kwargs,
    )


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def rows_by_alias(table):
    """The table's rows as dicts, grouped by alias: a guide can have several."""
    grouped = {}
    for row in table.to_dict("records"):
        grouped.setdefault(row["alias"], []).append(row)
    return grouped


def write_tables(path, guides, summaries, sites):
    """Write every table `run` writes into `path`, and read the two of them back."""
    write_guides_tsv(path / "guides_alignment_metadata.tsv", guides, summaries, sites)
    write_guides_bed(path / "guides_mapping.bed", guides, sites)
    write_columns_yaml(path / "columns.yaml")
    table = read_table(path / "guides_alignment_metadata.tsv")
    bed = [
        line.split("\t")
        for line in (path / "guides_mapping.bed")
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    return table, bed


def planted(tmp_path):
    """Four guides carrying an added G, planted at known coordinates on chr1."""
    guides = read_guides(
        pd.DataFrame(
            [
                ("plus", SPACER_G),
                ("minus", SPACER_C),
                ("plus_g", SPACER_A),
                ("kept_g", SPACER_T),
            ]
        ),
        add_leading_g=True,
    )
    contig = Contig()
    contig.add("A" * 100)
    assert contig.add(SPACER_G + "TGG") == 100
    contig.add("A" * 27)
    assert contig.add_minus("T" + SPACER_C + "AGG") == 150
    contig.add("A" * 26)
    assert contig.add("C" + SPACER_A + "GGG") == 200
    contig.add("A" * 10)
    assert contig.add("G" + SPACER_T + "TGG") == 234
    contig.add("A" * 10)
    fasta = write_fasta(tmp_path / "ref.fa", {"chr1": contig.sequence})
    return guides, contig, fasta


def test_classes_counts_and_coordinates(guides, sites):
    summaries = summarized(guides, sites)

    by_alias = {summary.alias: summary for summary in summaries}
    assert [summary.alias for summary in summaries] == list(SPACERS)
    assert [summary.guide_id for summary in summaries] == [
        spacer + "NGG" for spacer in SPACERS.values()
    ]
    assert {name: s.site_class for name, s in by_alias.items()} == {
        "uniq": "unique",
        "uniq_minus": "unique",
        "multi": "multi",
        "flipped": "no_site",
        "no_pam": "no_site",
        "nothing": "no_site",
        "dup": "unique",
    }
    uniq = by_alias["uniq"]
    # The PAM-less match counts as a match, but neither as a cut site nor
    # towards multi.
    assert (uniq.n_matches, uniq.n_sites) == (2, (1, 0, 0, 0))
    assert (uniq.guide_chr, uniq.guide_start, uniq.guide_end, uniq.strand) == (
        "chr1",
        100,
        120,
        "+",
    )
    assert (uniq.site_pam, uniq.cut_site, uniq.has_perfect_site) == ("TGG", 117, True)
    assert (uniq.mapped_sequence, uniq.input_sequence) == (SPACER_A + "TGG", SPACER_A)
    assert by_alias["uniq_minus"].cut_site == 203
    assert by_alias["uniq_minus"].site_pam == "AGG"
    multi = by_alias["multi"]
    assert (
        multi.has_perfect_site,
        multi.guide_chr,
        multi.cut_site,
        multi.site_pam,
        multi.mapped_sequence,
    ) == (True, None, None, None, None)
    flipped = by_alias["flipped"]
    assert (flipped.n_matches, flipped.n_sites, flipped.min_mismatches) == (
        2,
        (0, 0, 0, 0),
        None,
    )
    assert by_alias["no_pam"].n_matches == 1
    assert by_alias["nothing"].n_matches == 0
    assert by_alias["nothing"].n_sites == (0, 0, 0, 0)


def test_mismatched_sites_are_counted_by_how_many_mismatches_they_have(guides):
    perfect = summarize(
        guides, {0: [site(SPACER_A, 100), site(SPACER_A, 300, mismatches=2)]}
    )[0]
    assert (perfect.n_sites, perfect.min_mismatches, perfect.site_class) == (
        (1, 0, 1, 0),
        0,
        "unique",
    )

    offt = summarize(
        guides,
        {0: [site(SPACER_A, 800, mismatches=1), site(SPACER_A, 900, mismatches=3)]},
    )[0]
    assert (offt.has_perfect_site, offt.min_mismatches, offt.n_sites) == (
        False,
        1,
        (0, 1, 0, 1),
    )
    assert offt.site_class == "no_site"


def test_orientation_and_its_evidence_reach_the_summaries(guides, sites):
    summaries = summarized(guides, sites)

    assert [s.orientation for s in summaries] == [
        "forward",
        "forward",
        "forward",
        "reverse_complement",
        "no_pam",
        "no_match",
        "forward",
    ]
    # Evidence is only recorded for a reverse_complement entry.
    assert [s.orientation_evidence for s in summaries] == [
        None,
        None,
        None,
        "ambiguous",
        None,
        None,
        None,
    ]


def test_a_duplicate_gets_the_first_guide_s_sites_and_verdict(guides, sites):
    summaries = summarized(guides, sites)

    uniq, dup = summaries[0], summaries[6]
    assert guides[6].duplicate_of == "uniq"
    assert guides[0].duplicate_of is None
    assert dup.orientation == uniq.orientation == "forward"
    for field in ("site_class", "guide_start", "cut_site", "n_sites", "site_pam"):
        assert getattr(dup, field) == getattr(uniq, field)


def test_entries_sharing_a_sequence_have_their_sites_fanned_out_once(guides, sites):
    fanned = fan_out_sites(guides, sites)

    # 'uniq' and 'dup' are one guide: one guide_id, one set of site rows.
    assert guides[6].guide_id == guides[0].guide_id
    assert [s.start for s in fanned if s.guide_id == SPACER_A + "NGG"] == [100, 500]
    assert len(fanned) == 2 + 1 + 2 + 2 + 1


def test_without_an_orientation_check_every_guide_is_not_checked(guides, sites):
    summaries = summarize(guides, sites)

    assert {summary.orientation for summary in summaries} == {"not_checked"}
    assert {summary.orientation_evidence for summary in summaries} == {None}


def test_max_mismatches_sets_the_count_columns(guides):
    summaries = summarize(guides, {0: [site(SPACER_A, 100)]}, max_mismatches=0)

    assert summaries[0].n_sites == (1,)


def test_the_columns_are_the_fifteen_the_data_dictionary_names(guides):
    assert guide_columns() == COLUMN_NAMES
    assert guide_columns() == list(COLUMNS)
    # Every column belongs to exactly one documented group.
    assert sum(n for _, n in GROUPS) == len(COLUMNS)


@pytest.mark.parametrize(
    ("sites", "kwargs", "message"),
    [
        ({0: [site(SPACER_A, 1, mismatches=4)]}, {}, "more than max_mismatches=3"),
        ({0: [site(SPACER_A, 1, pam_class="NAG")]}, {}, "not one of"),
        ({}, {"orientation": {0: "forward"}}, "No orientation verdict"),
        ({}, {"max_mismatches": -1}, "non-negative integer"),
    ],
)
def test_summarize_rejects_inconsistent_inputs(guides, sites, kwargs, message):
    with pytest.raises(ValueError, match=message):
        summarize(guides, sites, **kwargs)


def test_guides_tsv(tmp_path, guides, sites):
    summaries = summarized(guides, sites)

    table, _ = write_tables(tmp_path, guides, summaries, sites)

    lines = (
        (tmp_path / "guides_alignment_metadata.tsv")
        .read_text(encoding="utf-8")
        .splitlines()
    )
    assert lines[0].split("\t") == COLUMN_NAMES
    assert lines[1].split("\t") == (
        [SPACER_A + "NGG", "uniq", SPACER_A, "unique", "2", "1", "forward", ""]
        + ["chr1", "100", "120", "+", "117", "TRUE", SPACER_A + "TGG"]
    )
    # The guide's own columns repeat on its second row; only the site changes.
    assert lines[2].split("\t")[:8] == lines[1].split("\t")[:8]
    assert lines[2].split("\t")[8:] == (
        ["chr1", "500", "520", "+", "517", "FALSE", SPACER_A + "TTT"]
    )
    # One row per entry per site, so the 11 sites of the library are 11 rows.
    assert len(lines) == 1 + 11
    rows = rows_by_alias(table)
    assert [row["start"] for row in rows["uniq"]] == ["100", "500"]
    assert [row["start"] for row in rows["dup"]] == ["100", "500"]
    assert [row["guide_id"] for row in rows["dup"]] == [SPACER_A + "NGG"] * 2
    assert [row["start"] for row in rows["flipped"]] == ["1000", "1100"]
    assert {row["pam_matches"] for row in rows["flipped"]} == {"FALSE"}
    multi = rows["multi"]
    assert [(row["start"], row["strand"]) for row in multi] == [
        ("600", "+"),
        ("700", "-"),
    ]
    assert {row["site_class"] for row in multi} == {"multi"}
    assert {row["n_cut_sites"] for row in multi} == {"2"}


def test_a_protospacer_matching_twice_with_one_pam_is_unique(tmp_path, guides, sites):
    summaries = summarized(guides, sites)

    table, _ = write_tables(tmp_path, guides, summaries, sites)

    uniq = summaries[0]
    assert (uniq.n_matches, uniq.n_sites[0], uniq.site_class) == (2, 1, "unique")
    rows = rows_by_alias(table)["uniq"]
    assert [
        (row["n_matches"], row["n_cut_sites"], row["site_class"]) for row in rows
    ] == [("2", "1", "unique")] * 2
    assert [row["pam_matches"] for row in rows] == ["TRUE", "FALSE"]


def test_two_pam_bearing_matches_are_multi_with_a_bed_row_each(tmp_path, guides, sites):
    summaries = summarized(guides, sites)

    table, bed = write_tables(tmp_path, guides, summaries, sites)

    multi = summaries[2]
    assert (multi.site_class, multi.n_sites[0], multi.n_matches) == ("multi", 2, 2)
    rows = rows_by_alias(table)["multi"]
    assert [row["n_cut_sites"] for row in rows] == ["2", "2"]
    named = [row for row in bed if row[3] == SPACER_G]
    assert named == [
        ["chr1", "600", "620", SPACER_G, "0", "+"],
        ["chr1", "700", "720", SPACER_G, "0", "-"],
    ]


def test_a_guide_matching_nowhere_gets_one_row_with_empty_site_columns(
    tmp_path, guides, sites
):
    summaries = summarized(guides, sites)

    table, bed = write_tables(tmp_path, guides, summaries, sites)

    (row,) = rows_by_alias(table)["nothing"]
    assert [row[name] for name in SITE_COLUMN_NAMES] == [""] * 7
    assert (row["n_matches"], row["n_cut_sites"]) == ("0", "0")
    assert (row["site_class"], row["orientation"]) == ("no_site", "no_match")
    assert row["spacer"] == SPACER_NOTHING
    assert not [line for line in bed if line[3] == SPACER_NOTHING]


def test_a_match_without_a_pam_is_left_out_of_the_bed_but_kept_in_the_tsv(
    tmp_path, guides, sites
):
    summaries = summarized(guides, sites)

    table, bed = write_tables(tmp_path, guides, summaries, sites)

    (row,) = rows_by_alias(table)["no_pam"]
    assert (row["chr"], row["start"], row["pam_matches"]) == ("chr1", "800", "FALSE")
    assert row["mapped_sequence"] == SPACER_T + "TTT"
    assert row["site_class"] == "no_site"
    # Nothing can be cut there, so the BED leaves it out.
    assert not [line for line in bed if line[3] == SPACER_T]
    assert "800" not in {line[1] for line in bed}


def test_the_bed_name_is_the_spacer_and_its_span_is_the_protospacer(
    tmp_path, guides, sites
):
    summaries = summarized(guides, sites)

    _, bed = write_tables(tmp_path, guides, summaries, sites)

    # One row per distinct guide, not per entry: 'uniq' and 'dup' are the same
    # sequence at the same locus, so bedtools must not see it twice.
    assert bed == [
        ["chr1", "100", "120", SPACER_A, "0", "+"],
        ["chr1", "200", "220", SPACER_C, "0", "-"],
        ["chr1", "600", "620", SPACER_G, "0", "+"],
        ["chr1", "700", "720", SPACER_G, "0", "-"],
    ]
    assert len(bed) == len({tuple(row) for row in bed})
    for chrom, start, end, name, score, strand in bed:
        assert name in SPACERS.values()
        # The span is the protospacer: no PAM, and no base 5' of it.
        assert int(end) - int(start) == len(name) == 20
        assert (chrom, score) == ("chr1", "0")
        assert strand in ("+", "-")


def test_columns_yaml_documents_exactly_what_is_written(tmp_path, guides, sites):
    summaries = summarized(guides, sites)

    table, bed = write_tables(tmp_path, guides, summaries, sites)

    document = yaml.safe_load((tmp_path / "columns.yaml").read_text(encoding="utf-8"))
    assert list(document) == ["guides_alignment_metadata.tsv", "guides_mapping.bed"]
    # The data dictionary names every column of the table, in the table's order.
    assert list(document["guides_alignment_metadata.tsv"]) == list(table.columns)
    assert list(document["guides_alignment_metadata.tsv"]) == guide_columns()
    assert list(document["guides_mapping.bed"]) == list(BED_COLUMNS)
    assert len(bed[0]) == len(BED_COLUMNS)
    assert sum(n for _, n in GROUPS) == len(COLUMNS)
    assert all(document["guides_alignment_metadata.tsv"][name] for name in COLUMNS)


def test_guides_tsv_rejects_summaries_that_are_not_one_per_guide(
    tmp_path, guides, sites
):
    summaries = summarized(guides, sites)

    with pytest.raises(ValueError, match="6 summaries for 7 guides"):
        write_guides_tsv(tmp_path / "guides.tsv", guides, summaries[:-1], sites)


def test_the_cut_site_column_is_three_bases_from_the_pam(tmp_path, guides, sites):
    # The boundary is end - 3 on + and start + 3 on -, counted from the PAM
    # either way.
    summaries = summarized(guides, sites)
    by_alias = {s.alias: s for s in summaries}
    table, _ = write_tables(tmp_path, guides, summaries, sites)
    rows = rows_by_alias(table)

    assert by_alias["uniq"].strand == "+"
    assert by_alias["uniq"].guide_end - by_alias["uniq"].cut_site == 3
    assert by_alias["uniq_minus"].strand == "-"
    assert by_alias["uniq_minus"].cut_site - by_alias["uniq_minus"].guide_start == 3
    assert (by_alias["uniq"].cut_site, by_alias["uniq_minus"].cut_site) == (117, 203)
    assert rows["uniq"][0]["cut_site"] == "117"
    assert rows["uniq_minus"][0]["cut_site"] == "203"


def test_qc_tsv_counts_sum_to_the_entries(tmp_path, guides, sites):
    summaries = summarized(guides, sites)

    write_qc_tsv(
        tmp_path / "qc.tsv",
        summaries,
        counts={"scan": {"n_matches": 12, "n_no_pam": 2}},
    )

    table = read_table(tmp_path / "qc.tsv")
    assert list(table.columns) == ["section", "name", "value"]
    rows = {(r.section, r.name): r.value for r in table.itertuples()}
    assert rows[("guides", "n_entries")] == "7"
    assert rows[("guides", "n_guides")] == "6"
    assert rows[("guides", "n_shared_spacers")] == "1"
    classes = {name: int(v) for (s, name), v in rows.items() if s == "class"}
    assert classes == {"unique": 3, "multi": 1, "no_site": 3}
    verdicts = {name: int(v) for (s, name), v in rows.items() if s == "orientation"}
    assert verdicts == {
        "forward": 4,
        "reverse_complement": 1,
        "no_pam": 1,
        "no_match": 1,
        "not_checked": 0,
    }
    assert sum(classes.values()) == sum(verdicts.values()) == len(guides)
    assert rows[("scan", "n_no_pam")] == "2"


def test_qc_tsv_without_an_orientation_check(tmp_path, guides, sites):
    write_qc_tsv(tmp_path / "qc.tsv", summarize(guides, sites))

    table = read_table(tmp_path / "qc.tsv")
    verdicts = table[table.section == "orientation"].set_index("name")["value"]
    assert verdicts.to_dict() == {
        "forward": "0",
        "reverse_complement": "0",
        "no_pam": "0",
        "no_match": "0",
        "not_checked": "7",
    }


def test_qc_tsv_rejects_a_clashing_section(tmp_path, guides, sites):
    with pytest.raises(ValueError, match="built-in section names: class"):
        write_qc_tsv(
            tmp_path / "qc.tsv",
            summarize(guides, sites),
            counts={"class": {"n": 1}},
        )


def test_coordinates_and_cut_sites_from_a_hand_built_reference(tmp_path):
    guides, contig, fasta = planted(tmp_path)
    assert [guide.extra_5p for guide in guides] == [True] * 4

    scan = scan_guides(guides, fasta)
    summaries = summarize(guides, scan.sites, scan.verdicts, max_mismatches=0)

    # + : protospacer [100, 120), cut between 116 and 117 (b = 117).
    # - : PAM [150, 153), protospacer [153, 173), b = 156.
    # + with the added G on a C: protospacer [201, 221), b = 218.
    # + with the added G on a G: protospacer [235, 255), b = 252.
    assert [(s.guide_start, s.guide_end, s.cut_site) for s in summaries] == [
        (100, 120, 117),
        (153, 173, 156),
        (201, 221, 218),
        (235, 255, 252),
    ]
    assert [s.site_pam for s in summaries] == ["TGG", "AGG", "GGG", "TGG"]
    assert [s.strand for s in summaries] == ["+", "-", "+", "+"]
    assert [scan.sites[row][0].genomic_5p_base for row in range(4)] == [
        "A",
        "T",
        "C",
        "G",
    ]
    # The 3 PAM-proximal spacer bases lie on the PAM side of each cut.
    assert contig.sequence[117:120] == SPACER_G[-3:]
    assert reverse_complement(contig.sequence[153:156]) == SPACER_C[-3:]
    assert contig.sequence[218:221] == SPACER_A[-3:]
    assert contig.sequence[252:255] == SPACER_T[-3:]

    table, bed = write_tables(tmp_path, guides, summaries, scan.sites)
    rows = rows_by_alias(table)
    assert [rows[alias][0]["cut_site"] for alias in ("plus", "minus")] == ["117", "156"]
    assert sorted(line[1] for line in bed) == ["100", "153", "201", "235"]


def test_guide_id_ends_in_the_pam_pattern_and_mapped_sequence_in_the_real_pam(tmp_path):
    guides, _, fasta = planted(tmp_path)

    scan = scan_guides(guides, fasta)
    summaries = summarize(guides, scan.sites, scan.verdicts, max_mismatches=0)
    table, _ = write_tables(tmp_path, guides, summaries, scan.sites)

    rows = {alias: found[0] for alias, found in rows_by_alias(table).items()}
    # Every guide_id ends in the pattern; every mapped_sequence in the real PAM.
    assert [row["guide_id"][-3:] for row in rows.values()] == ["NGG"] * 4
    assert [row["mapped_sequence"][-3:] for row in rows.values()] == [
        "TGG",
        "AGG",
        "GGG",
        "TGG",
    ]
    # 'kept_g' sits on a genomic G, so its added G pairs and is kept: the two
    # sequences line up base for base.
    assert rows["kept_g"]["guide_id"] == "G" + SPACER_T + "NGG"
    assert rows["kept_g"]["mapped_sequence"] == "G" + SPACER_T + "TGG"
    assert len(rows["kept_g"]["mapped_sequence"]) == len(rows["kept_g"]["guide_id"])
    # 'plus' sits on an A, so the G does not pair and is left out.
    assert rows["plus"]["guide_id"] == "G" + SPACER_G + "NGG"
    assert rows["plus"]["mapped_sequence"] == SPACER_G + "TGG"
    assert len(rows["plus"]["mapped_sequence"]) == len(rows["plus"]["guide_id"]) - 1
    assert [s.mapped_sequence for s in summaries] == [
        SPACER_G + "TGG",
        SPACER_C + "AGG",
        SPACER_A + "GGG",
        "G" + SPACER_T + "TGG",
    ]

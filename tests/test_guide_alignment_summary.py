"""`summarize`, the output tables, and the coordinates they report."""

import pandas as pd
import pytest
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
    write_guides_bed,
    write_guides_tsv,
    write_qc_tsv,
    write_sites_tsv,
)
from crispr_analysis_utils.guide_alignment.iupac import reverse_complement

SPACERS = {
    "uniq": SPACER_A,
    "uniq_minus": SPACER_C,
    "multi": SPACER_G,
    "offt": SPACER_T,
    "alt_only": "CAGTTCGAACTGGTACCATG",
    "nothing": "GGTACCTTAGCAATCGCTAG",
    "dup": SPACER_A,
}
ALT_PAMS = ["NAG", "NGA"]
# One verdict per source row: 'dup' is row 6 but takes row 0's verdict.
VERDICTS = {
    0: "forward",
    1: "forward",
    2: "reverse_complement",
    3: "no_pam",
    4: "no_perfect_hit",
    5: "forward",
}
EVIDENCE = {2: "ambiguous"}


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
        mapped_sequence=spacer,
        genomic_protospacer="A" * 20,
        added_g_base=None,
    )


@pytest.fixture
def guides():
    return read_guides(pd.DataFrame(list(SPACERS.items())))


@pytest.fixture
def sites():
    return {
        0: [
            site(SPACER_A, 100),
            site(SPACER_A, 300, mismatches=2),
            site(SPACER_A, 400, pam_class="NAG", pam="CAG"),
            site(SPACER_A, 500, pam_class="none", pam="TTT"),
        ],
        1: [site(SPACER_C, 200, strand="-", pam="AGG")],
        2: [site(SPACER_G, 600), site(SPACER_G, 700, strand="-")],
        3: [
            site(SPACERS["offt"], 800, mismatches=1),
            site(SPACERS["offt"], 900, mismatches=3),
        ],
        4: [
            site(SPACERS["alt_only"], 1000, pam_class="NAG", pam="TAG"),
            site(SPACERS["alt_only"], 1100, pam_class="none", pam="TCC"),
        ],
    }


def summarized(guides, sites, orientation=VERDICTS, **kwargs):
    return summarize(
        guides,
        sites,
        orientation,
        evidence=EVIDENCE if orientation is VERDICTS else None,
        alt_pams=kwargs.pop("alt_pams", ALT_PAMS),
        **kwargs,
    )


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


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
        "offt": "mismatch_only",
        "alt_only": "no_site",
        "nothing": "no_site",
        "dup": "unique",
    }
    uniq = by_alias["uniq"]
    # The NAG site and the PAM-less site neither count as NGG nor make it multi.
    assert (uniq.n_sites, uniq.n_alt_pam) == ((1, 0, 1, 0), (1, 0, 0, 0))
    assert (uniq.guide_chr, uniq.guide_start, uniq.guide_end, uniq.strand) == (
        "chr1",
        100,
        120,
        "+",
    )
    assert (uniq.site_pam, uniq.cut_site, uniq.has_perfect_site) == ("TGG", 117, True)
    assert (uniq.mapped_sequence, uniq.input_sequence) == (SPACER_A, SPACER_A)
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
    offt = by_alias["offt"]
    assert (offt.has_perfect_site, offt.min_mismatches, offt.n_sites) == (
        False,
        1,
        (0, 1, 0, 1),
    )
    alt_only = by_alias["alt_only"]
    assert (alt_only.min_mismatches, alt_only.n_alt_pam) == (None, (1, 0, 0, 0))
    assert by_alias["nothing"].n_sites == (0, 0, 0, 0)


def test_orientation_and_its_evidence_reach_the_summaries(guides, sites):
    summaries = summarized(guides, sites)

    assert [s.orientation for s in summaries] == [
        "forward",
        "forward",
        "reverse_complement",
        "no_pam",
        "no_perfect_hit",
        "forward",
        "forward",
    ]
    # Evidence is only recorded for a reverse_complement entry.
    assert [s.orientation_evidence for s in summaries] == [
        None,
        None,
        "ambiguous",
        None,
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


def test_entries_sharing_a_sequence_have_their_sites_written_once(guides, sites):
    fanned = fan_out_sites(guides, sites)

    # 'uniq' and 'dup' are one guide: one guide_id, one set of site rows.
    assert guides[6].guide_id == guides[0].guide_id
    assert [s.start for s in fanned if s.guide_id == SPACER_A + "NGG"] == [
        100,
        300,
        400,
        500,
    ]
    assert len(fanned) == 4 + 1 + 2 + 2 + 2


def test_without_an_orientation_check_every_guide_is_not_checked(guides, sites):
    summaries = summarize(guides, sites, alt_pams=ALT_PAMS)

    assert {summary.orientation for summary in summaries} == {"not_checked"}
    assert {summary.orientation_evidence for summary in summaries} == {None}


def test_max_mismatches_sets_the_count_columns(guides):
    summaries = summarize(guides, {0: [site(SPACER_A, 100)]}, max_mismatches=0)

    assert summaries[0].n_sites == (1,)
    assert summaries[0].n_alt_pam == (0,)
    # Without alternative PAMs there are no alternative-PAM columns at all.
    assert guide_columns(0)[16:] == ["n_sites_0mm", "input_sequence"]
    assert guide_columns(0, ALT_PAMS)[16:] == [
        "n_sites_0mm",
        "n_alt_pam_0mm",
        "input_sequence",
    ]


@pytest.mark.parametrize(
    ("sites", "kwargs", "message"),
    [
        ({0: [site(SPACER_A, 1, mismatches=4)]}, {}, "more than max_mismatches=3"),
        ({0: [site(SPACER_A, 1, pam_class="NAG")]}, {"alt_pams": ()}, "not one of"),
        ({}, {"orientation": {0: "forward"}}, "No orientation verdict"),
        ({}, {"max_mismatches": -1}, "non-negative integer"),
    ],
)
def test_summarize_rejects_inconsistent_inputs(guides, sites, kwargs, message):
    with pytest.raises(ValueError, match=message):
        summarize(guides, sites, **kwargs)


def test_guides_tsv(tmp_path, guides, sites):
    summaries = summarized(guides, sites)

    write_guides_tsv(
        tmp_path / "guides.tsv", summaries, max_mismatches=3, alt_pams=ALT_PAMS
    )

    lines = (tmp_path / "guides.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == [
        "guide_id",
        "alias",
        "spacer",
        "mapped_sequence",
        "site_pam",
        "has_perfect_site",
        "guide_chr",
        "guide_start",
        "guide_end",
        "strand",
        "cut_site",
        "site_class",
        "orientation",
        "orientation_evidence",
        "g_added",
        "min_mismatches",
        "n_sites_0mm",
        "n_sites_1mm",
        "n_sites_2mm",
        "n_sites_3mm",
        "n_alt_pam_0mm",
        "n_alt_pam_1mm",
        "n_alt_pam_2mm",
        "n_alt_pam_3mm",
        "input_sequence",
    ]
    assert lines[1].split("\t") == (
        [SPACER_A + "NGG", "uniq", SPACER_A, SPACER_A, "TGG", "TRUE"]
        + ["chr1", "100", "120", "+", "117"]
        + ["unique", "forward", "", "FALSE", "0"]
        + ["1", "0", "1", "0", "1", "0", "0", "0"]
        + [SPACER_A]
    )
    assert lines[3].split("\t")[:16] == (
        [SPACER_G + "NGG", "multi", SPACER_G, "", "", "TRUE", "", "", "", "", ""]
        + ["multi", "reverse_complement", "ambiguous", "FALSE", "0"]
    )
    dup = lines[7].split("\t")
    assert (dup[0], dup[1], dup[-1]) == (SPACER_A + "NGG", "dup", SPACER_A)
    assert len(lines) == 1 + len(SPACERS)


def test_guides_tsv_rejects_counts_of_another_length(tmp_path, guides, sites):
    summaries = summarize(guides, sites, max_mismatches=4, alt_pams=ALT_PAMS)

    with pytest.raises(ValueError, match="max_mismatches=3"):
        write_guides_tsv(tmp_path / "guides.tsv", summaries, max_mismatches=3)


def test_sites_tsv(tmp_path, guides, sites):
    write_sites_tsv(tmp_path / "sites.tsv", guides, sites)

    table = read_table(tmp_path / "sites.tsv")
    assert list(table.columns) == [
        "guide_id",
        "chr",
        "start",
        "end",
        "strand",
        "pam",
        "pam_class",
        "n_mismatches",
        "mismatch_positions",
        "mapped_sequence",
        "genomic_protospacer",
        "added_g_base",
    ]
    # 'uniq' and 'dup' are one guide, so their four sites appear once.
    assert len(table) == 11
    second = table.iloc[1].to_dict()
    assert second["start"] == "300"
    assert second["mismatch_positions"] == "1,2"
    assert second["added_g_base"] == ""
    assert second["mapped_sequence"] == SPACER_A
    assert table[table.guide_id == SPACER_A + "NGG"]["start"].tolist() == [
        "100",
        "300",
        "400",
        "500",
    ]


def test_the_bed_holds_the_unique_guides(tmp_path, guides, sites):
    summaries = summarized(guides, sites)

    write_guides_bed(tmp_path / "guides.bed", summaries)

    assert (tmp_path / "guides.bed").read_text(encoding="utf-8").splitlines() == [
        "chr1\t100\t120\tdup\t0\t+",
        "chr1\t100\t120\tuniq\t0\t+",
        "chr1\t200\t220\tuniq_minus\t0\t-",
    ]


def test_the_cut_site_column_is_three_bases_from_the_pam(tmp_path, guides, sites):
    # Still worth pinning now that cut_sites.bed is gone: the boundary is
    # end - 3 on + and start + 3 on -, counted from the PAM either way.
    by_alias = {s.alias: s for s in summarized(guides, sites)}

    assert by_alias["uniq"].strand == "+"
    assert by_alias["uniq"].guide_end - by_alias["uniq"].cut_site == 3
    assert by_alias["uniq_minus"].strand == "-"
    assert by_alias["uniq_minus"].cut_site - by_alias["uniq_minus"].guide_start == 3
    assert (by_alias["uniq"].cut_site, by_alias["uniq_minus"].cut_site) == (117, 203)


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
    assert classes == {"unique": 3, "multi": 1, "mismatch_only": 1, "no_site": 2}
    verdicts = {name: int(v) for (s, name), v in rows.items() if s == "orientation"}
    assert verdicts == {
        "forward": 4,
        "reverse_complement": 1,
        "no_pam": 1,
        "no_perfect_hit": 1,
        "not_checked": 0,
    }
    assert sum(classes.values()) == sum(verdicts.values()) == len(guides)
    assert rows[("scan", "n_no_pam")] == "2"


def test_qc_tsv_without_an_orientation_check(tmp_path, guides, sites):
    write_qc_tsv(tmp_path / "qc.tsv", summarize(guides, sites, alt_pams=ALT_PAMS))

    table = read_table(tmp_path / "qc.tsv")
    verdicts = table[table.section == "orientation"].set_index("name")["value"]
    assert verdicts.to_dict() == {
        "forward": "0",
        "reverse_complement": "0",
        "no_pam": "0",
        "no_perfect_hit": "0",
        "not_checked": "7",
    }


def test_qc_tsv_rejects_a_clashing_section(tmp_path, guides, sites):
    with pytest.raises(ValueError, match="built-in section names: class"):
        write_qc_tsv(
            tmp_path / "qc.tsv",
            summarize(guides, sites, alt_pams=ALT_PAMS),
            counts={"class": {"n": 1}},
        )


def test_coordinates_and_cut_sites_from_a_hand_built_reference(tmp_path):
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
    assert [guide.g_added for guide in guides] == [True] * 4
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
    assert [scan.sites[row][0].added_g_base for row in range(4)] == ["A", "T", "C", "G"]
    # The added G is kept in mapped_sequence only where the genome also has a G.
    assert [s.mapped_sequence for s in summaries] == [
        SPACER_G,
        SPACER_C,
        SPACER_A,
        "G" + SPACER_T,
    ]
    # The 3 PAM-proximal spacer bases lie on the PAM side of each cut.
    assert contig.sequence[117:120] == SPACER_G[-3:]
    assert reverse_complement(contig.sequence[153:156]) == SPACER_C[-3:]
    assert contig.sequence[218:221] == SPACER_A[-3:]
    assert contig.sequence[252:255] == SPACER_T[-3:]

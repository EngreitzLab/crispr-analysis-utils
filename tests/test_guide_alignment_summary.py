import pandas as pd
import pytest
from guide_alignment_fixtures import (
    SPACER_A,
    SPACER_C,
    SPACER_G,
    Contig,
    hit,
    write_fasta,
    write_sam,
)

from crispr_analysis_utils.guide_alignment import (
    OrientationResult,
    PerfectHit,
    Site,
    fan_out_sites,
    guide_columns,
    read_guides,
    read_sites,
    summarize,
    write_cut_sites_bed,
    write_guides_bed,
    write_guides_tsv,
    write_orientation_tsv,
    write_sites_tsv,
    write_summary_tsv,
)
from crispr_analysis_utils.guide_alignment.iupac import reverse_complement

SPACERS = {
    "uniq": SPACER_A,
    "uniq_minus": SPACER_C,
    "multi": SPACER_G,
    "offt": "TTGCCAGTACGATCAGGTCA",
    "alt_only": "CAGTTCGAACTGGTACCATG",
    "nothing": "GGTACCTTAGCAATCGCTAG",
    "dup": SPACER_A,
}


def site(guide_id, start, strand="+", mismatches=0, pam_class="NGG", pam="TGG"):
    return Site(
        guide_id=guide_id,
        chr="chr1",
        start=start,
        end=start + 20,
        strand=strand,
        pam=pam,
        pam_class=pam_class,
        n_mismatches=mismatches,
        mismatch_positions=tuple(range(1, mismatches + 1)),
        genomic_protospacer="A" * 20,
        leading_g_base=None,
    )


@pytest.fixture
def guides():
    return read_guides(pd.DataFrame(list(SPACERS.items())))


@pytest.fixture
def sites():
    return {
        0: [
            site("uniq", 100),
            site("uniq", 300, mismatches=2),
            site("uniq", 400, pam_class="NAG", pam="CAG"),
            site("uniq", 500, pam_class="none", pam="TTT"),
        ],
        1: [site("uniq_minus", 200, strand="-", pam="AGG")],
        2: [site("multi", 600), site("multi", 700, strand="-")],
        3: [site("offt", 800, mismatches=1), site("offt", 900, mismatches=3)],
        4: [
            site("alt_only", 1000, pam_class="NAG", pam="TAG"),
            site("alt_only", 1100, pam_class="none", pam="TCC"),
        ],
    }


VERDICTS = {0: "ok", 1: "ok", 2: "reversed", 3: "no_pam", 4: "no_perfect_hit", 5: "ok"}


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def test_classes_counts_and_coordinates(guides, sites):
    summaries = summarize(guides, sites, VERDICTS)

    by_id = {summary.guide_id: summary for summary in summaries}
    assert [summary.guide_id for summary in summaries] == list(SPACERS)
    assert {name: s.alignment_class for name, s in by_id.items()} == {
        "uniq": "unique",
        "uniq_minus": "unique",
        "multi": "multi",
        "offt": "off_target_only",
        "alt_only": "no_site",
        "nothing": "no_site",
        "dup": "unique",
    }
    uniq = by_id["uniq"]
    # The NAG site and the PAM-less site neither count as NGG nor make it multi.
    assert (uniq.n_ngg, uniq.n_alt_pam) == ((1, 0, 1, 0), (1, 0, 0, 0))
    assert (uniq.guide_chr, uniq.guide_start, uniq.guide_end, uniq.strand) == (
        "chr1",
        100,
        120,
        "+",
    )
    assert (uniq.pam, uniq.cut_site, uniq.targeting) == ("TGG", 117, True)
    assert by_id["uniq_minus"].cut_site == 203
    assert by_id["uniq_minus"].pam == "AGG"
    multi = by_id["multi"]
    assert (multi.targeting, multi.guide_chr, multi.cut_site, multi.pam) == (
        True,
        None,
        None,
        None,
    )
    offt = by_id["offt"]
    assert (offt.targeting, offt.min_ngg_mismatches, offt.n_ngg) == (
        False,
        1,
        (0, 1, 0, 1),
    )
    alt_only = by_id["alt_only"]
    assert (alt_only.min_ngg_mismatches, alt_only.n_alt_pam) == (None, (1, 0, 0, 0))
    assert by_id["nothing"].n_ngg == (0, 0, 0, 0)


def test_a_duplicate_gets_the_first_guide_s_sites_and_verdict(guides, sites):
    summaries = summarize(guides, sites, VERDICTS)

    uniq, dup = summaries[0], summaries[6]
    assert dup.duplicate_of == "uniq"
    assert uniq.duplicate_of is None
    assert dup.orientation == uniq.orientation == "ok"
    for field in ("alignment_class", "guide_start", "cut_site", "n_ngg", "pam"):
        assert getattr(dup, field) == getattr(uniq, field)

    fanned = fan_out_sites(guides, sites)
    assert [s.start for s in fanned if s.guide_id == "dup"] == [100, 300, 400, 500]
    assert [s.start for s in fanned if s.guide_id == "uniq"] == [100, 300, 400, 500]
    assert len(fanned) == 4 + 1 + 2 + 2 + 2 + 4


def test_without_an_orientation_check_every_guide_is_not_checked(guides, sites):
    summaries = summarize(guides, sites)

    assert {summary.orientation for summary in summaries} == {"not_checked"}


def test_max_mismatches_sets_the_count_columns(guides):
    summaries = summarize(
        guides, {0: [site("uniq", 100)]}, max_mismatches=0, alt_pams=()
    )

    assert summaries[0].n_ngg == (1,)
    assert summaries[0].n_alt_pam == (0,)
    assert guide_columns(0)[13:15] == ["n_ngg_0mm", "n_alt_pam_0mm"]


@pytest.mark.parametrize(
    ("sites", "kwargs", "message"),
    [
        ({0: [site("uniq", 1, mismatches=4)]}, {}, "more than max_mismatches=3"),
        ({0: [site("uniq", 1, pam_class="NAG")]}, {"alt_pams": ()}, "not one of"),
        ({}, {"orientation": {0: "ok"}}, "No orientation verdict"),
    ],
)
def test_summarize_rejects_inconsistent_inputs(guides, sites, kwargs, message):
    with pytest.raises(ValueError, match=message):
        summarize(guides, sites, **kwargs)


def test_guides_tsv(tmp_path, guides, sites):
    summaries = summarize(guides, sites, VERDICTS)

    write_guides_tsv(tmp_path / "guides.tsv", summaries, max_mismatches=3)

    lines = (tmp_path / "guides.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == [
        "guide_id",
        "spacer",
        "pam",
        "targeting",
        "guide_chr",
        "guide_start",
        "guide_end",
        "strand",
        "cut_site",
        "alignment_class",
        "orientation",
        "leading_g",
        "min_ngg_mismatches",
        "n_ngg_0mm",
        "n_ngg_1mm",
        "n_ngg_2mm",
        "n_ngg_3mm",
        "n_alt_pam_0mm",
        "n_alt_pam_1mm",
        "n_alt_pam_2mm",
        "n_alt_pam_3mm",
        "duplicate_of",
        "input_sequence",
    ]
    assert lines[1].split("\t") == (
        ["uniq", SPACER_A, "TGG", "TRUE", "chr1", "100", "120", "+", "117"]
        + ["unique", "ok", "FALSE", "0", "1", "0", "1", "0", "1", "0", "0", "0"]
        + ["", SPACER_A]
    )
    assert lines[3].split("\t")[:13] == (
        ["multi", SPACER_G, "", "TRUE", "", "", "", "", ""]
        + ["multi", "reversed", "FALSE", "0"]
    )
    assert lines[7].split("\t")[-2:] == ["uniq", SPACER_A]
    assert len(lines) == 1 + len(SPACERS)


def test_guides_tsv_rejects_counts_of_another_length(tmp_path, guides, sites):
    summaries = summarize(guides, sites, max_mismatches=4)

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
        "genomic_protospacer",
        "leading_g_base",
    ]
    assert len(table) == 15
    second = table.iloc[1].to_dict()
    assert second["start"] == "300"
    assert second["mismatch_positions"] == "1,2"
    assert second["leading_g_base"] == ""
    assert table[table.guide_id == "dup"]["start"].tolist() == [
        "100",
        "300",
        "400",
        "500",
    ]


def test_orientation_tsv(tmp_path, guides):
    orientation = OrientationResult(
        hits={
            0: [
                PerfectHit("chr1", 100, 120, "+", "AAA", "TGG"),
                PerfectHit("chr2", 5, 25, "-", "CC", "AGG"),
            ]
        },
        verdicts=VERDICTS,
    )

    write_orientation_tsv(tmp_path / "orientation.tsv", guides, orientation)

    table = read_table(tmp_path / "orientation.tsv")
    assert list(table.columns) == [
        "guide_id",
        "chr",
        "start",
        "end",
        "strand",
        "flank_5p",
        "flank_3p",
        "verdict",
    ]
    assert table.iloc[1].tolist() == ["uniq", "chr2", "5", "25", "-", "CC", "AGG", "ok"]
    assert table.iloc[2].tolist() == ["uniq_minus", "", "", "", "", "", "", "ok"]
    assert table.guide_id.tolist().count("dup") == 2
    # One row per hit, or one row for a guide without hits.
    assert len(table) == 2 + 5 + 2


def test_orientation_tsv_needs_every_verdict(tmp_path, guides):
    with pytest.raises(ValueError, match="No orientation verdict for guide 'multi'"):
        write_orientation_tsv(
            tmp_path / "orientation.tsv",
            guides,
            OrientationResult(verdicts={0: "ok", 1: "ok"}),
        )


def test_bed_files_hold_the_unique_guides(tmp_path, guides, sites):
    summaries = summarize(guides, sites, VERDICTS)

    write_guides_bed(tmp_path / "guides.bed", summaries)
    write_cut_sites_bed(tmp_path / "cut_sites.bed", summaries)

    assert (tmp_path / "guides.bed").read_text(encoding="utf-8").splitlines() == [
        "chr1\t100\t120\tdup\t0\t+",
        "chr1\t100\t120\tuniq\t0\t+",
        "chr1\t200\t220\tuniq_minus\t0\t-",
    ]
    assert (tmp_path / "cut_sites.bed").read_text(encoding="utf-8").splitlines() == [
        "chr1\t117\t118\tdup\t0\t+",
        "chr1\t117\t118\tuniq\t0\t+",
        "chr1\t203\t204\tuniq_minus\t0\t-",
    ]


def test_summary_tsv_counts_sum_to_the_guides(tmp_path, guides, sites):
    summaries = summarize(guides, sites, VERDICTS)

    write_summary_tsv(
        tmp_path / "summary.tsv",
        summaries,
        counts={"site_search": {"n_records": 12, "n_over_max_mismatches": 2}},
        parameters={
            "pam": "NGG",
            "alt_pams": ["NAG", "NGA"],
            "add_leading_g": False,
            "contigs": None,
        },
    )

    table = read_table(tmp_path / "summary.tsv")
    assert list(table.columns) == ["section", "name", "value"]
    rows = {(r.section, r.name): r.value for r in table.itertuples()}
    assert rows[("guides", "n_guides")] == "7"
    assert rows[("guides", "n_distinct_spacers")] == "6"
    assert rows[("guides", "n_duplicates")] == "1"
    classes = {name: int(v) for (s, name), v in rows.items() if s == "class"}
    assert classes == {"unique": 3, "multi": 1, "off_target_only": 1, "no_site": 2}
    verdicts = {name: int(v) for (s, name), v in rows.items() if s == "orientation"}
    assert verdicts == {
        "ok": 4,
        "reversed": 1,
        "no_pam": 1,
        "no_perfect_hit": 1,
        "not_checked": 0,
    }
    assert sum(classes.values()) == sum(verdicts.values()) == len(guides)
    assert rows[("site_search", "n_over_max_mismatches")] == "2"
    assert rows[("parameter", "alt_pams")] == "NAG,NGA"
    assert rows[("parameter", "add_leading_g")] == "FALSE"
    assert rows[("parameter", "contigs")] == ""


def test_summary_tsv_without_an_orientation_check(tmp_path, guides, sites):
    write_summary_tsv(tmp_path / "summary.tsv", summarize(guides, sites))

    table = read_table(tmp_path / "summary.tsv")
    verdicts = table[table.section == "orientation"].set_index("name")["value"]
    assert verdicts.to_dict() == {
        "ok": "0",
        "reversed": "0",
        "no_pam": "0",
        "no_perfect_hit": "0",
        "not_checked": "7",
    }


def test_summary_tsv_rejects_a_clashing_section(tmp_path, guides, sites):
    with pytest.raises(ValueError, match="built-in section names: class"):
        write_summary_tsv(
            tmp_path / "summary.tsv",
            summarize(guides, sites),
            counts={"class": {"n": 1}},
        )


def test_cut_sites_from_hand_built_alignments(tmp_path):
    guides = read_guides(
        pd.DataFrame([("plus", SPACER_G), ("minus", SPACER_C), ("plus_g", SPACER_A)]),
        add_leading_g=True,
    )
    assert [guide.leading_g for guide in guides] == [False, True, True]
    contig = Contig()
    contig.add("A" * 100)
    assert contig.add(SPACER_G + "TGG") == 100
    contig.add("A" * 27)
    assert contig.add_minus("T" + SPACER_C + "AGG") == 150
    contig.add("A" * 26)
    assert contig.add("C" + SPACER_A + "GGG") == 200
    contig.add("A" * 10)
    fasta = write_fasta(tmp_path / "ref.fa", {"chr1": contig.sequence})
    sam = write_sam(
        tmp_path / "sites.sam",
        {"chr1": len(contig.sequence)},
        [
            hit("0:NGG", "chr1", 100, "+", 23),
            hit("1:NGG", "chr1", 150, "-", 24),
            hit("2:NGG", "chr1", 200, "+", 24),
        ],
    )

    result = read_sites(sam, guides, fasta)
    summaries = summarize(guides, result.sites)
    write_cut_sites_bed(tmp_path / "cut_sites.bed", summaries)

    # + : protospacer [100, 120), cut between 116 and 117 (b = 117).
    # - : PAM [150, 153), protospacer [153, 173), b = 156.
    # + with G: G at 200, protospacer [201, 221), b = 218.
    assert [(s.guide_start, s.guide_end, s.cut_site) for s in summaries] == [
        (100, 120, 117),
        (153, 173, 156),
        (201, 221, 218),
    ]
    # The 3 PAM-proximal spacer bases lie on the PAM side of each cut.
    assert contig.sequence[117:120] == SPACER_G[-3:]
    assert reverse_complement(contig.sequence[153:156]) == SPACER_C[-3:]
    assert contig.sequence[218:221] == SPACER_A[-3:]
    assert (tmp_path / "cut_sites.bed").read_text(encoding="utf-8").splitlines() == [
        "chr1\t117\t118\tplus\t0\t+",
        "chr1\t156\t157\tminus\t0\t-",
        "chr1\t218\t219\tplus_g\t0\t+",
    ]

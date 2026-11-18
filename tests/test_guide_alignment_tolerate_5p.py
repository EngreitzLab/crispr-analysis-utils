"""`tolerate_5p`: a PAM-distal protospacer base the genome need not hold.

A designer with no G-starting protospacer at a target writes the G anyway, so
the intended site is a 20-base window pairing at 19 of its positions. Every
test plants such a site and reads it back.
"""

import csv
import types

import pandas as pd
import pytest
from guide_alignment_fixtures import SPACER_A, SPACER_C, Contig, write_fasta

from crispr_analysis_utils.guide_alignment.library import read_guides
from crispr_analysis_utils.guide_alignment.pipeline import run
from crispr_analysis_utils.guide_alignment.scan import scan_guides
from crispr_analysis_utils.guide_alignment.summary import summarize

SUBSTITUTED_A = "G" + SPACER_A[1:]
SUBSTITUTED_C = "G" + SPACER_C[1:]


@pytest.fixture
def reference(tmp_path):
    """SPACER_A on the plus strand and SPACER_C on the minus, both with a PAM."""
    chr1 = Contig()
    chr1.add("TTTT")
    plus = chr1.add(SPACER_A)
    chr1.add("TGG")
    chr1.add("TTTT")
    # add_minus lays down the reverse complement, so the PAM comes first.
    minus = chr1.add_minus(SPACER_C + "TGG") + 3
    chr1.add("TTTT")
    return types.SimpleNamespace(
        path=write_fasta(tmp_path / "ref.fa", {"chr1": chr1.sequence}),
        plus=plus,
        minus=minus,
    )


def test_a_substituted_5p_base_hides_the_site_until_it_is_tolerated(reference):
    guides = read_guides(pd.DataFrame([("sub", SUBSTITUTED_A)]))

    assert scan_guides(guides, reference.path).sites == {}

    scan = scan_guides(guides, reference.path, tolerate_5p=1)
    (site,) = scan.sites[0]
    assert (site.chr, site.start, site.end, site.strand) == (
        "chr1",
        reference.plus,
        reference.plus + 20,
        "+",
    )
    assert site.genomic_spacer_5p == SPACER_A[0]
    assert site.mapped_sequence == SPACER_A + "TGG"
    assert site.n_mismatches == 1
    assert site.mismatch_positions == (1,)
    summary = summarize(
        guides, scan.sites, scan.verdicts, max_mismatches=1, tolerate_5p=1
    )[0]
    assert (summary.site_class, summary.has_perfect_site) == ("unique", False)
    # A scan told to tolerate nothing would not count this site at all.
    assert summarize(guides, scan.sites, scan.verdicts)[0].site_class == "no_site"


def test_the_tolerated_base_moves_end_not_start_on_the_minus_strand(reference):
    guides = read_guides(pd.DataFrame([("sub", SUBSTITUTED_C)]))

    (site,) = scan_guides(guides, reference.path, tolerate_5p=1).sites[0]
    assert (site.start, site.end, site.strand) == (
        reference.minus,
        reference.minus + 20,
        "-",
    )
    assert site.genomic_spacer_5p == SPACER_C[0]


def test_tolerating_a_base_leaves_a_fully_matching_guide_where_it_was(reference):
    guides = read_guides(pd.DataFrame([("exact", SPACER_A)]))

    (plain,) = scan_guides(guides, reference.path).sites[0]
    (tolerant,) = scan_guides(guides, reference.path, tolerate_5p=1).sites[0]

    assert (tolerant.start, tolerant.end, tolerant.strand) == (
        plain.start,
        plain.end,
        plain.strand,
    )
    assert tolerant.mapped_sequence == plain.mapped_sequence
    # The column is the genome's base either way; here the design paired it.
    assert tolerant.genomic_spacer_5p == SPACER_A[0]
    assert plain.genomic_spacer_5p == ""


def test_the_added_g_is_reported_apart_from_the_tolerated_base(reference):
    # 21 bases: an added G, then a protospacer whose position 1 was written G.
    guides = read_guides(pd.DataFrame([("both", "G" + SUBSTITUTED_A)]))
    assert guides[0].n_extra == 1

    (site,) = scan_guides(guides, reference.path, tolerate_5p=1).sites[0]
    assert site.genomic_spacer_5p == SPACER_A[0]
    assert site.genomic_5p_base == "T"
    assert site.paired_5p is False
    assert (site.start, site.end) == (reference.plus, reference.plus + 20)


def test_a_site_that_pairs_outranks_one_found_only_by_tolerating(tmp_path):
    chr1 = Contig()
    chr1.add("TTTT")
    designed = chr1.add(SPACER_A)
    chr1.add("TGG" + "TTTT")
    # The same 19-base core, under a different base at position 1: a place
    # the shorter search reaches and the guide does not pair with throughout.
    chr1.add("C" + SPACER_A[1:] + "TGG" + "TTTT")
    reference = write_fasta(tmp_path / "ref.fa", {"chr1": chr1.sequence})
    guides = read_guides(pd.DataFrame([("g", SPACER_A)]))

    scan = scan_guides(guides, reference, tolerate_5p=1)
    assert len(scan.sites[0]) == 2

    summary = summarize(
        guides, scan.sites, scan.verdicts, max_mismatches=1, tolerate_5p=1
    )[0]
    assert (summary.site_class, summary.guide_start) == ("unique", designed)
    # Both sites bear a PAM; the class counts the better one.
    assert summary.n_sites == (1, 1)


def test_a_match_with_no_room_for_the_tolerated_base_is_not_a_site(tmp_path):
    # The core matches at base 0, so the protospacer would start at -1.
    reference = write_fasta(
        tmp_path / "edge.fa", {"chr1": SPACER_A[1:] + "TGG" + "TTTT"}
    )
    guides = read_guides(pd.DataFrame([("sub", SUBSTITUTED_A)]))

    scan = scan_guides(guides, reference, tolerate_5p=1)
    assert scan.sites == {}
    assert scan.n_matches == 1


@pytest.mark.parametrize("value", [-1, 20, 25])
def test_tolerate_5p_is_validated(reference, value):
    guides = read_guides(pd.DataFrame([("g", SPACER_A)]))

    with pytest.raises(ValueError):
        scan_guides(guides, reference.path, tolerate_5p=value)


def test_the_run_writes_the_tolerated_base_and_records_the_setting(reference, tmp_path):
    guides = pd.DataFrame([("sub", SUBSTITUTED_A)])

    result = run(guides, reference.path, tmp_path / "out", tolerate_5p=1)

    with open(
        tmp_path / "out" / "alignments" / "guide_genome_annotation.tsv",
        encoding="utf-8",
    ) as handle:
        (row,) = list(csv.DictReader(handle, delimiter="\t"))
    assert row["site_class"] == "unique"
    assert row["genomic_spacer_5p"] == SPACER_A[0]
    assert (row["start"], row["end"]) == (str(reference.plus), str(reference.plus + 20))
    assert result.metadata["parameters"]["tolerate_5p"] == 1

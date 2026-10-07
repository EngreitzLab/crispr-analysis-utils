import pandas as pd
import pytest
from guide_alignment_fixtures import (
    SPACER_A,
    SPACER_C,
    SPACER_G,
    Contig,
    hit,
    unmapped,
    write_fasta,
    write_sam,
)

from crispr_analysis_utils.guide_alignment import (
    PerfectHit,
    check_orientation,
    library_reversed,
    read_guides,
)
from crispr_analysis_utils.guide_alignment.orientation import (
    NO_PAM,
    NO_PERFECT_HIT,
    OK,
    REVERSED,
)

SPACER_4 = "TTGCCAGTACGATCAGGTCA"
SPACER_5 = "CAGTTCGAACTGGTACCATG"
SPACER_6 = "GGTACCTTAGCAATCGCTAG"


def library(*spacers):
    rows = [(f"g{i + 1}", spacer) for i, spacer in enumerate(spacers)]
    return read_guides(pd.DataFrame(rows))


def run(tmp_path, reference, records, guides, **kwargs):
    """Write the reference and one SAM, then check the orientation."""
    fasta = write_fasta(tmp_path / "ref.fa", reference)
    lengths = {name: len(sequence) for name, sequence in reference.items()}
    sam = write_sam(tmp_path / "spacers.sam", lengths, records)
    return check_orientation(sam, guides, fasta, **kwargs)


def test_verdicts_on_both_strands(tmp_path):
    contig = Contig()
    # Each piece is 5' flank + spacer + 3' flank, in guide orientation.
    plus_ok = contig.add("AAA" + SPACER_A + "TGG") + 3
    minus_ok = contig.add_minus("TTT" + SPACER_C + "AGG") + 3
    plus_reversed = contig.add("CCA" + SPACER_G + "AAA") + 3
    minus_reversed = contig.add_minus("CCT" + SPACER_4 + "TTT") + 3
    no_pam = contig.add("AAA" + SPACER_5 + "AAA") + 3
    assert contig.sequence[minus_ok - 3 : minus_ok] == "CCT"
    assert contig.sequence[minus_reversed + 20 : minus_reversed + 23] == "AGG"

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [
            hit("0", "chr1", plus_ok, "+", 20),
            hit("1", "chr1", minus_ok, "-", 20),
            hit("2", "chr1", plus_reversed, "+", 20),
            hit("3", "chr1", minus_reversed, "-", 20),
            hit("4", "chr1", no_pam, "+", 20),
        ],
        library(SPACER_A, SPACER_C, SPACER_G, SPACER_4, SPACER_5, SPACER_6),
    )

    assert result.verdicts == {
        0: OK,
        1: OK,
        2: REVERSED,
        3: REVERSED,
        4: NO_PAM,
        5: NO_PERFECT_HIT,
    }
    assert result.hits[0] == [PerfectHit("chr1", 3, 23, "+", "AAA", "TGG")]
    assert result.hits[1] == [PerfectHit("chr1", 29, 49, "-", "TTT", "AGG")]
    assert result.hits[2] == [PerfectHit("chr1", 55, 75, "+", "CCA", "AAA")]
    assert result.hits[3] == [PerfectHit("chr1", 81, 101, "-", "CCT", "TTT")]
    assert 5 not in result.hits


def test_ok_beats_reversed_and_reversed_beats_no_pam(tmp_path):
    contig = Contig()
    a_reversed = contig.add("CCA" + SPACER_A + "AAA") + 3
    a_ok = contig.add("AAA" + SPACER_A + "TGG") + 3
    c_no_pam = contig.add("AAA" + SPACER_C + "AAA") + 3
    c_reversed = contig.add("CCC" + SPACER_C + "AAA") + 3

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [
            hit("0", "chr1", a_reversed, "+", 20),
            hit("0", "chr1", a_ok, "+", 20, secondary=True),
            hit("1", "chr1", c_no_pam, "+", 20),
            hit("1", "chr1", c_reversed, "+", 20, secondary=True),
        ],
        library(SPACER_A, SPACER_C),
    )

    assert result.verdicts == {0: OK, 1: REVERSED}
    assert len(result.hits[0]) == len(result.hits[1]) == 2


@pytest.mark.parametrize(
    ("alt_pams", "verdicts"),
    [(("NAG", "NGA"), {0: OK, 1: REVERSED}), ((), {0: NO_PAM, 1: NO_PAM})],
)
def test_alternative_pams_count_on_either_side(tmp_path, alt_pams, verdicts):
    contig = Contig()
    nag_3prime = contig.add("AAA" + SPACER_A + "TAG") + 3
    nga_reversed = contig.add("TCA" + SPACER_C + "AAA") + 3  # TCN reverses NGA

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [
            hit("0", "chr1", nag_3prime, "+", 20),
            hit("1", "chr1", nga_reversed, "+", 20),
        ],
        library(SPACER_A, SPACER_C),
        alt_pams=alt_pams,
    )

    assert result.verdicts == verdicts


def test_a_genomic_n_in_a_flank(tmp_path):
    contig = Contig()
    n_at_n = contig.add("AAA" + SPACER_A + "NGG") + 3
    n_at_g = contig.add("AAA" + SPACER_C + "TNG") + 3

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [hit("0", "chr1", n_at_n, "+", 20), hit("1", "chr1", n_at_g, "+", 20)],
        library(SPACER_A, SPACER_C),
    )

    assert result.verdicts == {0: OK, 1: NO_PAM}


def test_a_flank_cut_by_a_contig_end_matches_no_pam(tmp_path):
    minus = Contig()
    minus.add_minus("CCA" + SPACER_G)  # its 3' flank would lie before base 0
    reference = {
        "cut_3prime": "AAA" + SPACER_A + "TG",  # NGG cut to TG
        "cut_5prime": "CC" + SPACER_C + "AAA",  # CCN cut to CC
        "minus_start": minus.sequence,
    }

    result = run(
        tmp_path,
        reference,
        [
            hit("0", "cut_3prime", 3, "+", 20),
            hit("1", "cut_5prime", 2, "+", 20),
            hit("2", "minus_start", 0, "-", 20),
        ],
        library(SPACER_A, SPACER_C, SPACER_G),
    )

    # The cut flank matches nothing; the other flank of the same hit still counts.
    assert result.verdicts == {0: NO_PAM, 1: NO_PAM, 2: REVERSED}
    assert result.hits[0] == [PerfectHit("cut_3prime", 3, 23, "+", "AAA", "TG")]
    assert result.hits[1] == [PerfectHit("cut_5prime", 2, 22, "+", "CC", "AAA")]
    assert result.hits[2] == [PerfectHit("minus_start", 0, 20, "-", "CCA", "")]


def test_imperfect_and_unmapped_reads_give_no_perfect_hit(tmp_path):
    assert SPACER_A[0] == "A"
    contig = "AAA" + "C" + SPACER_A[1:] + "TGG"

    result = run(
        tmp_path,
        {"chr1": contig},
        [hit("0", "chr1", 3, "+", 20), unmapped("1")],
        library(SPACER_A, SPACER_C, SPACER_G),
    )

    assert result.verdicts == {0: NO_PERFECT_HIT, 1: NO_PERFECT_HIT, 2: NO_PERFECT_HIT}
    assert result.hits == {}
    assert result.counts() == {
        "n_records": 2,
        "n_unmapped": 1,
        "n_imperfect": 1,
        "n_outside_contigs": 0,
        "n_off_contig_end": 0,
        "n_repeated": 0,
        "n_hits": 0,
    }


def test_duplicates_get_no_verdict_of_their_own(tmp_path):
    contig = "AAA" + SPACER_A + "TGG"

    result = run(
        tmp_path,
        {"chr1": contig},
        [hit("0", "chr1", 3, "+", 20), hit("0", "chr1", 3, "+", 20, secondary=True)],
        library(SPACER_A, SPACER_A, SPACER_C),
    )

    assert result.verdicts == {0: OK, 2: NO_PERFECT_HIT}
    assert len(result.hits[0]) == 1
    assert result.n_repeated == 1


def test_the_contig_filter_applies(tmp_path):
    reference = {"chr1": "A" * 30, "chrUn": "AAA" + SPACER_A + "TGG"}

    result = run(
        tmp_path,
        reference,
        [hit("0", "chrUn", 3, "+", 20)],
        library(SPACER_A),
        contigs=["chr1"],
    )

    assert result.verdicts == {0: NO_PERFECT_HIT}
    assert result.n_outside_contigs == 1


@pytest.mark.parametrize(
    ("record", "message"),
    [
        (hit("0", "chr1", 3, "+", 23), "all 20 bases"),
        (hit("0:NGG", "chr1", 3, "+", 20), "does not carry a guide row"),
    ],
)
def test_records_must_match_the_library(tmp_path, record, message):
    with pytest.raises(ValueError, match=message):
        run(
            tmp_path,
            {"chr1": "AAA" + SPACER_A + "TGG"},
            [record],
            library(SPACER_A),
        )


@pytest.mark.parametrize(
    ("verdicts", "expected"),
    [
        ([REVERSED, REVERSED, OK, NO_PAM, NO_PERFECT_HIT], True),
        ([REVERSED, OK], False),
        ([REVERSED, OK, OK], False),
        ([NO_PAM, NO_PERFECT_HIT], False),
        ([], False),
    ],
)
def test_library_reversed(verdicts, expected):
    assert library_reversed(verdicts) is expected

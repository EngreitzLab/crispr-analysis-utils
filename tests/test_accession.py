import itertools
import random
import re

import numpy as np
import pandas as pd
import pytest

import crispr_analysis_utils as cau
from crispr_analysis_utils import accession

# Pinned: a change here renames every guide the lab has accessioned.
GOLDEN = {
    "GAAAAGCCAACATGAATGCAGNGG": "GU602A170WJ",
    "GTTTGAGCACATAGCTCTAAGNGG": "GU6ZRJ4S7E2",
    "GTGAGCACATAGCTCTAAGTGNGG": "GU6W92CKQ1E",
    # One short enough to check by hand: GACGT packs to 0b1_10_00_01_10_11,
    # which is 1563, which is 1, 16, 27 in base 32, so "1", "G", "V".
    "GACGTNGG": "GU1GV",
}


@pytest.mark.parametrize(("sequence", "expected"), GOLDEN.items())
def test_golden_accessions(sequence, expected):
    assert cau.accession.guide_accession(sequence) == expected


@pytest.mark.parametrize(("sequence", "expected"), GOLDEN.items())
def test_golden_accessions_decode(sequence, expected):
    assert accession.guide_sequence(expected) == sequence


def test_accession_shape():
    assert re.fullmatch(r"GU[0-9A-HJKMNP-TV-Z]+", accession.guide_accession("GACGTNGG"))


@pytest.mark.parametrize("length", [19, 20, 21, 22])
def test_accession_width(length):
    sequence = "G" + "ACGT" * 8
    expected = {19: 8, 20: 9, 21: 9, 22: 9}[length]
    got = accession.guide_accession(sequence[:length] + accession.PAM)
    assert len(got) - 2 == expected


@pytest.mark.parametrize(
    "variant",
    [
        "gaaaagccaacatgaatgcagngg",
        "GaAaAgCcAaCaTgAaTgCaGnGg",
        "  GAAAAGCCAACATGAATGCAGNGG",
        "GAAAAGCCAACATGAATGCAGNGG\n",
        "\tGAAAAGCCAACATGAATGCAGNGG\r\n",
        np.str_("GAAAAGCCAACATGAATGCAGNGG"),
    ],
)
def test_case_and_surrounding_whitespace_are_ignored(variant):
    assert accession.guide_accession(variant) == "GU602A170WJ"


def test_normalize_guide_sequence():
    assert accession.normalize_guide_sequence(" gacgTngg\n") == "GACGTNGG"


# --- No two guides share an accession -------------------------------------


def test_every_guide_of_a_length_round_trips():
    # 4^7 guides: every G + 7 bases + NGG. Exhaustive, so nothing can collide.
    sequences = [
        "G" + "".join(bases) + accession.PAM
        for bases in itertools.product("ACGT", repeat=7)
    ]
    accessions = accession.guide_accessions(sequences)
    assert len(set(accessions)) == len(sequences)
    assert [accession.guide_sequence(a) for a in accessions] == sequences


def test_guides_of_different_lengths_never_collide():
    # The leading 1 bit is what keeps G + AAAA apart from G + AA.
    sequences = [
        "G" + "".join(bases) + accession.PAM
        for length in range(0, 9)
        for bases in itertools.product("ACGT", repeat=length)
    ]
    accessions = accession.guide_accessions(sequences)
    assert len(set(accessions)) == len(sequences)


def test_long_guides_round_trip():
    rng = random.Random(0)
    for _ in range(2000):
        length = rng.randint(18, 24)
        sequence = "G" + "".join(rng.choices("ACGT", k=length)) + accession.PAM
        assert accession.guide_sequence(accession.guide_accession(sequence)) == sequence


# --- The PAM --------------------------------------------------------------


@pytest.mark.parametrize(
    "sequence",
    [
        "GAAAAGCCAACATGAATGCAG",  # no PAM
        "GAAAAGCCAACATGAATGCAGAGG",  # a genomic PAM
        "GAAAAGCCAACATGAATGCAGTGG",
        "GAAAAGCCAACATGAATGCAGNAG",  # an alternative PAM
        "GAAAAGCCAACATGAATGCAGNGA",
        "GAAAAGCCAACATGAATGCAGNG",
        "GNGGAAAAGCCAACATGAATGCAG",  # a 5' PAM
    ],
)
def test_a_sequence_without_a_trailing_ngg_is_rejected(sequence):
    with pytest.raises(ValueError, match="does not end with the PAM 'NGG'"):
        accession.guide_accession(sequence)


def test_a_spacer_ending_in_gg_is_not_read_as_a_pam():
    # One spacer in 16 ends with GG; reading those bases as a PAM would give
    # the same accession as a genuinely shorter guide.
    with pytest.raises(ValueError, match="does not end with the PAM"):
        accession.guide_accession("GAAAAGCCAACATGAATGCGG")
    assert accession.guide_accession(
        "GAAAAGCCAACATGAATGCGGNGG"
    ) != accession.guide_accession("GAAAAGCCAACATGAATGCNGG")


def test_an_n_before_the_pam_is_rejected():
    with pytest.raises(ValueError, match="holds an N before the PAM"):
        accession.guide_accession("GAAAAGCCAACNTGAATGCAGNGG")


# --- The leading G --------------------------------------------------------


@pytest.mark.parametrize("sequence", ["AAAAGCCAACATGAATGCAGNGG", "TACGTNGG", "NGG"])
def test_a_sequence_without_a_leading_g_is_rejected(sequence):
    with pytest.raises(ValueError, match="does not start with G"):
        accession.guide_accession(sequence)


def test_the_g_is_part_of_the_accession():
    assert accession.guide_accession("GGACGTNGG") != accession.guide_accession(
        "GACGTNGG"
    )


# --- Other input errors ---------------------------------------------------


@pytest.mark.parametrize("sequence", ["", "   ", "\n"])
def test_empty_sequence_is_rejected(sequence):
    with pytest.raises(ValueError, match="must not be empty"):
        accession.guide_accession(sequence)


@pytest.mark.parametrize(
    ("sequence", "unknown"),
    [
        ("GACGUNGG", "U"),
        ("GACGRNGG", "R"),
        ("GAC GTNGG", " "),
        ("GAC-GTNGG", "-"),
    ],
)
def test_non_acgt_is_rejected(sequence, unknown):
    with pytest.raises(ValueError, match=f"holds {unknown!r}"):
        accession.guide_accession(sequence)


def test_non_ascii_is_rejected():
    # No character outside ASCII uppercases into A, C, G, T or N, so they all
    # land in the unknown-character message.
    with pytest.raises(ValueError, match="holds 'É'"):
        accession.guide_accession("GACGTéNGG")


@pytest.mark.parametrize("value", [None, float("nan"), 7, b"GACGTNGG", ["GACGTNGG"]])
def test_non_string_is_rejected(value):
    with pytest.raises(TypeError, match="must be a string"):
        accession.guide_accession(value)


# --- Decoding -------------------------------------------------------------


def test_guide_sequence_reads_crockford_aliases():
    assert accession.guide_sequence("GU1GV") == "GACGTNGG"
    assert accession.guide_sequence("gu1gv") == "GACGTNGG"
    assert accession.guide_sequence("GUIGV") == "GACGTNGG"  # I is 1
    assert accession.guide_sequence("GULGV") == "GACGTNGG"  # L is 1


@pytest.mark.parametrize("value", [None, 7, b"GU1GV"])
def test_guide_sequence_rejects_a_non_string(value):
    with pytest.raises(TypeError, match="must be a string"):
        accession.guide_sequence(value)


@pytest.mark.parametrize("text", ["", "GU", "1GV", "EL1GV"])
def test_guide_sequence_needs_the_prefix_and_some_bases(text):
    with pytest.raises(ValueError, match="does not begin with 'GU'"):
        accession.guide_sequence(text)


def test_guide_sequence_rejects_a_bad_character():
    with pytest.raises(ValueError, match="holds 'U', which is not a base 32"):
        accession.guide_sequence("GU1GU")


@pytest.mark.parametrize("text", ["GU0", "GU1", "GU2"])
def test_guide_sequence_rejects_what_encodes_no_guide(text):
    # 0 has no leading 1 bit; 1 is the empty guide; 2 is a guide starting with C.
    with pytest.raises(ValueError, match="encodes no guide"):
        accession.guide_sequence(text)


# --- A whole library ------------------------------------------------------


def test_guide_accessions_keeps_order_and_repeats():
    sequences = [*GOLDEN, "gaaaagccaacatgaatgcagngg"]
    assert accession.guide_accessions(sequences) == [*GOLDEN.values(), "GU602A170WJ"]


def test_guide_accessions_takes_a_series_and_a_generator():
    series = pd.Series(list(GOLDEN), index=[10, 20, 30, 40])
    assert accession.guide_accessions(series) == list(GOLDEN.values())
    assert accession.guide_accessions(s for s in GOLDEN) == list(GOLDEN.values())
    assert accession.guide_accessions([]) == []


def test_guide_accessions_names_the_bad_position():
    series = pd.Series(["GACGTNGG", "GACGTNGG", None], index=["a", "b", "c"])
    with pytest.raises(TypeError, match="position 2: A guide sequence must be"):
        accession.guide_accessions(series)
    with pytest.raises(ValueError, match="position 1: Guide sequence 'GACGT'"):
        accession.guide_accessions(["GACGTNGG", "GACGT"])

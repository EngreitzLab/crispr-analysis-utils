"""Accessions encode the protospacer, so the same target always gets the same id."""

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
    "AAAAGCCAACATGAATGCAG": "GU102A170WJ",
    "TTTGAGCACATAGCTCTAAG": "GU1ZRJ4S7E2",
    # Short enough to check by hand: ACGT packs to 0b1_00_01_10_11, which is
    # 283, which is 8, 27 in base 32, so "8" then "V".
    "ACGT": "GU8V",
}


@pytest.mark.parametrize(("spacer", "expected"), GOLDEN.items())
def test_golden_accessions(spacer, expected):
    assert cau.accession.guide_accession(spacer) == expected


@pytest.mark.parametrize(("spacer", "expected"), GOLDEN.items())
def test_golden_accessions_decode(spacer, expected):
    assert accession.spacer_of(expected) == spacer


def test_accession_shape():
    assert re.fullmatch(r"GU[0-9A-HJKMNP-TV-Z]+", accession.guide_accession("ACGT"))


@pytest.mark.parametrize(("length", "width"), [(19, 8), (20, 9), (21, 9), (22, 9)])
def test_accession_width(length, width):
    spacer = ("ACGT" * 8)[:length]
    assert len(accession.guide_accession(spacer)) - 2 == width


@pytest.mark.parametrize(
    "variant",
    [
        "aaaagccaacatgaatgcag",
        "AaAaGcCaAcAtGaAtGcAg",
        "  AAAAGCCAACATGAATGCAG",
        "AAAAGCCAACATGAATGCAG\n",
        "\tAAAAGCCAACATGAATGCAG\r\n",
        np.str_("AAAAGCCAACATGAATGCAG"),
    ],
)
def test_case_and_surrounding_whitespace_are_ignored(variant):
    assert accession.guide_accession(variant) == "GU102A170WJ"


def test_normalize_spacer():
    assert accession.normalize_spacer(" acgT\n") == "ACGT"


# --- the point of keying on the protospacer ------------------------------


def test_the_same_target_gets_one_accession_whatever_the_5prime_base():
    # One sub-library orders the bare protospacer, another prepends a G for
    # the U6 promoter. Same target, so the same accession.
    protospacer = "CGGGGCCTCCGGACAGGGCA"
    bare = protospacer
    with_g = "G" + protospacer
    assert accession.guide_accession(bare) == accession.guide_accession(with_g[1:])
    # The ordered sequences differ, and guide_id is what keeps them apart.
    assert accession.guide_accession(bare) != accession.guide_accession(with_g)


def test_a_protospacer_and_its_reverse_complement_differ():
    spacer = "AAAAGCCAACATGAATGCAG"
    reverse = spacer[::-1].translate(str.maketrans("ACGT", "TGCA"))
    assert accession.guide_accession(spacer) != accession.guide_accession(reverse)


# --- no two protospacers share an accession ------------------------------


def test_every_spacer_of_a_length_round_trips():
    # 4^7 protospacers, exhaustive, so nothing can collide unnoticed.
    spacers = ["".join(bases) for bases in itertools.product("ACGT", repeat=7)]
    accessions = accession.guide_accessions(spacers)
    assert len(set(accessions)) == len(spacers)
    assert [accession.spacer_of(a) for a in accessions] == spacers


def test_spacers_of_different_lengths_never_collide():
    # The leading 1 bit is what keeps AAAA apart from AA.
    spacers = [
        "".join(bases)
        for length in range(1, 9)
        for bases in itertools.product("ACGT", repeat=length)
    ]
    assert len(set(accession.guide_accessions(spacers))) == len(spacers)


def test_long_spacers_round_trip():
    rng = random.Random(0)
    for _ in range(2000):
        spacer = "".join(rng.choices("ACGT", k=rng.randint(18, 24)))
        assert accession.spacer_of(accession.guide_accession(spacer)) == spacer


# --- input errors ---------------------------------------------------------


@pytest.mark.parametrize("spacer", ["", "   ", "\n"])
def test_empty_spacer_is_rejected(spacer):
    with pytest.raises(ValueError, match="must not be empty"):
        accession.guide_accession(spacer)


@pytest.mark.parametrize(
    ("spacer", "unknown"),
    [
        ("ACGU", "U"),
        ("ACGN", "N"),
        ("ACGR", "R"),
        ("AC GT", " "),
        ("AC-GT", "-"),
    ],
)
def test_non_acgt_is_rejected(spacer, unknown):
    with pytest.raises(ValueError, match=f"holds {unknown!r}"):
        accession.guide_accession(spacer)


def test_a_pam_left_on_the_sequence_is_not_detected():
    # Nothing in the sequence says where the protospacer ends, so leaving the
    # PAM on gives a different accession and no error. The caller decides.
    assert accession.guide_accession(
        "ACGTACGTACGTACGTACGT"
    ) != accession.guide_accession("ACGTACGTACGTACGTACGTTGG")


def test_non_ascii_is_rejected():
    # No character outside ASCII uppercases into A, C, G or T, so they all
    # land in the unknown-character message.
    with pytest.raises(ValueError, match="holds 'É'"):
        accession.guide_accession("ACGTé")


@pytest.mark.parametrize("value", [None, float("nan"), 7, b"ACGT", ["ACGT"]])
def test_non_string_is_rejected(value):
    with pytest.raises(TypeError, match="must be a string"):
        accession.guide_accession(value)


# --- decoding -------------------------------------------------------------


def test_spacer_of_reads_crockford_aliases():
    assert accession.spacer_of("GU8V") == "ACGT"
    assert accession.spacer_of("gu8v") == "ACGT"


@pytest.mark.parametrize("value", [None, 7, b"GU8V"])
def test_spacer_of_rejects_a_non_string(value):
    with pytest.raises(TypeError, match="must be a string"):
        accession.spacer_of(value)


@pytest.mark.parametrize("text", ["", "GU", "8V", "EL8V"])
def test_spacer_of_needs_the_prefix_and_some_bases(text):
    with pytest.raises(ValueError, match="does not begin with 'GU'"):
        accession.spacer_of(text)


def test_spacer_of_rejects_a_bad_character():
    with pytest.raises(ValueError, match="holds 'U', which is not a base 32"):
        accession.spacer_of("GU8U")


@pytest.mark.parametrize("text", ["GU0", "GU1"])
def test_spacer_of_rejects_what_encodes_nothing(text):
    # 0 has no leading 1 bit; 1 is the empty protospacer.
    with pytest.raises(ValueError, match="encodes no protospacer"):
        accession.spacer_of(text)


# --- a whole library ------------------------------------------------------


def test_guide_accessions_keeps_order_and_repeats():
    spacers = [*GOLDEN, "aaaagccaacatgaatgcag"]
    assert accession.guide_accessions(spacers) == [*GOLDEN.values(), "GU102A170WJ"]


def test_guide_accessions_takes_a_series_and_a_generator():
    series = pd.Series(list(GOLDEN), index=[10, 20, 30])
    assert accession.guide_accessions(series) == list(GOLDEN.values())
    assert accession.guide_accessions(s for s in GOLDEN) == list(GOLDEN.values())
    assert accession.guide_accessions([]) == []


def test_guide_accessions_names_the_bad_position():
    series = pd.Series(["ACGT", "ACGT", None], index=["a", "b", "c"])
    with pytest.raises(TypeError, match="position 2: A protospacer must be"):
        accession.guide_accessions(series)
    with pytest.raises(ValueError, match="position 1: Protospacer 'AXGT'"):
        accession.guide_accessions(["ACGT", "AXGT"])

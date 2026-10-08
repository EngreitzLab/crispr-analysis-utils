import itertools
import re
import uuid

import numpy as np
import pandas as pd
import pytest

import crispr_analysis_utils as cau
from crispr_analysis_utils import accession

# Pinned: a change to any of these changes every accession the lab has made.
GOLDEN = {
    "GAAAAGCCAACATGAATGCAG": "GU5C664EE08629",
    "GTTTGAGCACATAGCTCTAAG": "GU0A344F9F7147",
    "GAAAAGCCAACATGAATGCA": "GU7B0BBA9859A2",
}


def test_namespace_is_the_documented_derivation():
    assert accession.GUIDE_NAMESPACE == uuid.uuid5(
        uuid.NAMESPACE_URL,
        "https://github.com/EngreitzLab/crispr-analysis-utils/accession/guide",
    )


@pytest.mark.parametrize(("sequence", "expected"), GOLDEN.items())
def test_golden_accessions(sequence, expected):
    assert cau.accession.guide_accession(sequence) == expected


@pytest.mark.parametrize("sequence", GOLDEN)
def test_accession_is_the_first_digits_of_the_uuid5(sequence):
    digits = uuid.uuid5(accession.GUIDE_NAMESPACE, sequence).hex[:12].upper()
    assert accession.guide_accession(sequence) == "GU" + digits


def test_accession_shape():
    assert re.fullmatch(r"GU[0-9A-F]{12}", accession.guide_accession("ACGT"))


@pytest.mark.parametrize(
    "variant",
    [
        "gaaaagccaacatgaatgcag",
        "GaAaAgCcAaCaTgAaTgCaG",
        "  GAAAAGCCAACATGAATGCAG",
        "GAAAAGCCAACATGAATGCAG\n",
        "\tGAAAAGCCAACATGAATGCAG\r\n",
        np.str_("GAAAAGCCAACATGAATGCAG"),
    ],
)
def test_case_and_surrounding_whitespace_are_ignored(variant):
    assert accession.guide_accession(variant) == "GU5C664EE08629"


def test_normalize_guide_sequence():
    assert accession.normalize_guide_sequence(" acgT\n") == "ACGT"


def test_the_g_and_the_pam_change_the_accession():
    spacer = "AAAAGCCAACATGAATGCAG"
    ids = {
        accession.guide_accession(spacer),
        accession.guide_accession("G" + spacer),
        accession.guide_accession(spacer + "TGG"),
    }
    assert len(ids) == 3


@pytest.mark.parametrize("sequence", ["", "   ", "\n"])
def test_empty_sequence_is_rejected(sequence):
    with pytest.raises(ValueError, match="must not be empty"):
        accession.guide_accession(sequence)


@pytest.mark.parametrize(
    ("sequence", "others"),
    [
        ("ACGU", "U"),
        ("ACGN", "N"),
        ("ACGR", "R"),
        ("AC GT", " "),
        ("ACG-T", "-"),
        ("ACGTNGG", "N"),
    ],
)
def test_non_acgt_is_rejected(sequence, others):
    with pytest.raises(ValueError, match=f"holds {others!r}"):
        accession.guide_accession(sequence)


def test_non_ascii_is_rejected():
    # No non-ASCII letter uppercases to A, C, G or T today; this keeps it so.
    with pytest.raises(ValueError, match="only A, C, G and T"):
        accession.guide_accession("ACGTé")


@pytest.mark.parametrize("value", [None, float("nan"), 7, b"ACGT", ["ACGT"]])
def test_non_string_is_rejected(value):
    with pytest.raises(TypeError, match="must be a string"):
        accession.guide_accession(value)


def test_guide_accessions_keeps_order_and_repeats():
    sequences = list(GOLDEN) + ["gaaaagccaacatgaatgcag"]
    assert accession.guide_accessions(sequences) == [*GOLDEN.values(), "GU5C664EE08629"]


def test_guide_accessions_takes_a_series_and_a_generator():
    series = pd.Series(list(GOLDEN), index=[10, 20, 30])
    assert accession.guide_accessions(series) == list(GOLDEN.values())
    assert accession.guide_accessions(s for s in GOLDEN) == list(GOLDEN.values())
    assert accession.guide_accessions([]) == []


def test_guide_accessions_names_the_bad_position():
    series = pd.Series(["ACGT", "ACGT", None], index=["a", "b", "c"])
    with pytest.raises(TypeError, match="position 2: A guide sequence must be"):
        accession.guide_accessions(series)
    with pytest.raises(ValueError, match="position 1: Guide sequence 'AXGT'"):
        accession.guide_accessions(["ACGT", "AXGT"])


def test_guide_accessions_refuses_a_collision(monkeypatch):
    # One hexadecimal digit leaves 16 accessions, so 17 sequences must collide.
    monkeypatch.setattr(accession, "_DIGITS", 1)
    sequences = ["".join(bases) for bases in itertools.product("ACGT", repeat=3)]
    with pytest.raises(
        ValueError, match=r"Two different sequences have the accession GU[0-9A-F]:"
    ):
        accession.guide_accessions(sequences[:17])


def test_guide_accessions_has_no_collision_on_all_8mers():
    sequences = ["".join(bases) for bases in itertools.product("ACGT", repeat=8)]
    assert len(set(accession.guide_accessions(sequences))) == 4**8

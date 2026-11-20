import pytest

from crispr_analysis_utils.guide_alignment import iupac

TABLE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def test_the_table_is_the_iupac_one():
    # Which genomic base satisfies a code is the Rust core's to decide, and
    # tests/test_core_scan.py checks it; this table only validates and
    # complements patterns.
    assert iupac.IUPAC_BASES == {
        code: frozenset(bases) for code, bases in TABLE.items()
    }


def test_complement_and_reverse_complement():
    assert iupac.complement("ACGTRYSWKMBDHVN") == "TGCAYRSWMKVHDBN"
    assert iupac.complement("acgtn") == "tgcan"
    assert iupac.reverse_complement("NGG") == "CCN"
    assert iupac.reverse_complement("NRG") == "CYN"
    assert iupac.reverse_complement("AACX") == "XGTT"


def test_complement_agrees_with_the_base_sets():
    for code, bases in iupac.IUPAC_BASES.items():
        assert iupac.IUPAC_BASES[iupac.complement(code)] == {
            iupac.complement(base) for base in bases
        }
        assert iupac.complement(iupac.complement(code)) == code


def test_count_ambiguous():
    assert iupac.count_ambiguous("NGG") == 1
    assert iupac.count_ambiguous("ngg") == 1
    assert iupac.count_ambiguous("NNGRRT") == 4
    assert iupac.count_ambiguous("TGG") == 0


def test_validate_pams():
    assert iupac.validate_pams(["ngg", "NAG", "nga"]) == ("NGG", "NAG", "NGA")
    assert iupac.validate_pams(("NGG",)) == ("NGG",)


@pytest.mark.parametrize(
    ("pams", "message"),
    [
        ([], "At least one"),
        (["NGG", ""], "must not be empty"),
        (["NGZ"], "not IUPAC"),
        (["NGG", "ngg"], "more than once"),
        (["NGG", "NG"], "same length"),
    ],
)
def test_validate_pams_rejects(pams, message):
    with pytest.raises(ValueError, match=message):
        iupac.validate_pams(pams)


def test_validate_pams_rejects_a_bare_string():
    with pytest.raises(TypeError, match="not a single string"):
        iupac.validate_pams("NGG")

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


@pytest.mark.parametrize(("code", "bases"), sorted(TABLE.items()))
def test_base_matches_follows_the_iupac_table(code, bases):
    for base in "ACGT":
        assert iupac.base_matches(base, code) is (base in bases)
        assert iupac.base_matches(base.lower(), code.lower()) is (base in bases)


def test_a_genomic_n_matches_only_n():
    assert iupac.base_matches("N", "N")
    for code in "ACGTRYSWKMBDHV":
        assert not iupac.base_matches("N", code)
        assert not iupac.base_matches("n", code)


def test_base_matches_rejects_what_is_not_one_iupac_code():
    with pytest.raises(ValueError, match="not an IUPAC"):
        iupac.base_matches("A", "X")
    with pytest.raises(ValueError, match="not an IUPAC"):
        iupac.base_matches("A", "NG")


def test_pattern_matches():
    assert iupac.pattern_matches("TGG", "NGG")
    assert iupac.pattern_matches("tgg", "ngg")
    assert iupac.pattern_matches("NGG", "NGG")  # a genomic N at the N position
    assert not iupac.pattern_matches("TNG", "NGG")  # and at a G position
    assert not iupac.pattern_matches("TAG", "NGG")
    assert iupac.pattern_matches("TAG", "NRG")
    assert not iupac.pattern_matches("GG", "NGG")  # a flank cut by a contig end
    assert not iupac.pattern_matches("", "NGG")


def test_pattern_matches_rejects_a_non_iupac_pattern():
    with pytest.raises(ValueError, match="IUPAC"):
        iupac.pattern_matches("TTT", "NGX")
    with pytest.raises(ValueError, match="IUPAC"):
        iupac.pattern_matches("TT", "NGX")


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

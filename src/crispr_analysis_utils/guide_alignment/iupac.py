"""IUPAC nucleotide codes: matching, complements, and PAM patterns.

A PAM pattern such as ``NGG`` is a string of IUPAC codes. A genomic base
matches a code when it is one of the bases the code stands for; ``N`` matches
anything, a genomic ``N`` included, and a genomic ``N`` matches nothing else.
"""

from __future__ import annotations

from collections.abc import Iterable

IUPAC_BASES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
"""The bases each IUPAC nucleotide code stands for."""

_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def base_matches(base: str, code: str) -> bool:
    """Whether a genomic base matches an IUPAC code.

    Parameters
    ----------
    base
        One genomic base. Case is ignored.
    code
        One IUPAC nucleotide code. Case is ignored.

    Returns
    -------
    bool
        True when `code` is ``N``, or when `base` is one of the bases `code`
        stands for. A genomic ``N``, or any character other than A, C, G and
        T, matches only ``N``.

    Raises
    ------
    ValueError
        If `code` is not a single IUPAC nucleotide code.
    """
    code = code.upper()
    if len(code) != 1 or code not in IUPAC_BASES:
        raise ValueError(f"{code!r} is not an IUPAC nucleotide code.")
    return code == "N" or base.upper() in IUPAC_BASES[code]


def pattern_matches(sequence: str, pattern: str) -> bool:
    """Whether a genomic sequence matches an IUPAC pattern, base by base.

    Parameters
    ----------
    sequence
        Genomic bases. Case is ignored.
    pattern
        IUPAC codes, such as ``NGG``. Case is ignored.

    Returns
    -------
    bool
        True when the two have the same length and every base matches the
        code at its position (see `base_matches`). A sequence shorter than the
        pattern, such as a flank cut short by the end of a contig, never
        matches.

    Raises
    ------
    ValueError
        If `pattern` holds a character that is not an IUPAC nucleotide code.
    """
    pattern = validate_pattern(pattern)
    if len(sequence) != len(pattern):
        return False
    return all(
        code == "N" or base in IUPAC_BASES[code]
        for base, code in zip(sequence.upper(), pattern, strict=True)
    )


def complement(sequence: str) -> str:
    """Complement an IUPAC sequence, keeping its case.

    Parameters
    ----------
    sequence
        IUPAC nucleotide codes, in either case.

    Returns
    -------
    str
        Each code replaced by its complement (A-T, C-G, R-Y, K-M, B-V, D-H;
        S, W and N are their own complements). Characters outside the IUPAC
        alphabet are left unchanged.
    """
    return sequence.translate(_COMPLEMENT)


def reverse_complement(sequence: str) -> str:
    """Reverse-complement an IUPAC sequence, keeping its case.

    Parameters
    ----------
    sequence
        IUPAC nucleotide codes, in either case.

    Returns
    -------
    str
        The complement (see `complement`) read in reverse.
    """
    return complement(sequence)[::-1]


def count_ambiguous(pattern: str) -> int:
    """Count the characters of a pattern that are not A, C, G or T.

    Parameters
    ----------
    pattern
        IUPAC codes, such as ``NGG``. Case is ignored.

    Returns
    -------
    int
        The number of positions holding a code other than A, C, G and T,
        such as the ``N`` of ``NGG``.
    """
    return sum(code not in "ACGT" for code in pattern.upper())


def validate_pattern(pattern: str) -> str:
    """Uppercase an IUPAC pattern and check its alphabet.

    Parameters
    ----------
    pattern
        IUPAC codes, such as ``ngg``.

    Returns
    -------
    str
        The pattern, uppercased.

    Raises
    ------
    ValueError
        If the pattern is empty or holds a character that is not an IUPAC
        nucleotide code.
    """
    upper = pattern.upper()
    if not upper:
        raise ValueError("A PAM pattern must not be empty.")
    invalid = sorted(set(upper) - IUPAC_BASES.keys())
    if invalid:
        raise ValueError(
            f"PAM pattern {pattern!r} holds characters that are not IUPAC "
            f"nucleotide codes: {', '.join(map(repr, invalid))}."
        )
    return upper


def validate_pams(pams: Iterable[str]) -> tuple[str, ...]:
    """Check a list of PAM patterns: the primary PAM first, then the alternatives.

    Parameters
    ----------
    pams
        IUPAC patterns, such as ``["NGG", "NAG", "NGA"]``. The first is the
        primary PAM.

    Returns
    -------
    tuple of str
        The patterns, uppercased, in the order given.

    Raises
    ------
    TypeError
        If `pams` is a single string rather than a list of patterns.
    ValueError
        If the list is empty, a pattern is empty or holds a character that is
        not an IUPAC nucleotide code, a pattern appears twice, or the patterns
        differ in length.
    """
    if isinstance(pams, str):
        raise TypeError("pams must be a list of PAM patterns, not a single string.")
    patterns = tuple(validate_pattern(pam) for pam in pams)
    if not patterns:
        raise ValueError("At least one PAM pattern is required.")
    repeated = sorted({p for p in patterns if patterns.count(p) > 1})
    if repeated:
        raise ValueError(f"PAM patterns appear more than once: {', '.join(repeated)}.")
    if len({len(p) for p in patterns}) > 1:
        raise ValueError(
            f"PAM patterns must all have the same length: {', '.join(patterns)}."
        )
    return patterns

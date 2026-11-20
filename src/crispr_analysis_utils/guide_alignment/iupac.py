"""IUPAC nucleotide codes: complements and PAM patterns.

A PAM pattern such as ``NGG`` is a string of IUPAC codes. A genomic base
matches a code when it is one of the bases the code stands for; ``N`` matches
anything, a genomic ``N`` included, and a genomic ``N`` matches nothing else.
"""

from __future__ import annotations

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

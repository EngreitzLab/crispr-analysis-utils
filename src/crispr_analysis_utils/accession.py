"""Accession ids for SpCas9 guides, encoded from the guide's own sequence.

A guide is written ``G...NGG``. Its accession is ``GU`` followed by the bases
before the PAM, packed two bits each behind a leading 1 bit and written in
Crockford's base 32. The packing is reversible, so two guides never share an
accession.
"""

from __future__ import annotations

from collections.abc import Iterable

PAM = "NGG"
"""The PAM every accessioned guide ends with."""

_PREFIX = "GU"
_BASES = "ACGT"
_ALPHABET = "0123456789ABCDEFGHJKMNPQRSTVWXYZ"
_VALUES = {character: value for value, character in enumerate(_ALPHABET)}
_VALUES.update({"I": 1, "L": 1, "O": 0})


def normalize_guide_sequence(sequence: str) -> str:
    """Strip and uppercase a guide sequence, and check it is a whole guide.

    Parameters
    ----------
    sequence
        A guide, 5' to 3': a G, the rest of the protospacer, then ``NGG``.

    Returns
    -------
    str
        `sequence` stripped and in uppercase.

    Raises
    ------
    TypeError
        If `sequence` is not a string.
    ValueError
        If `sequence` is empty, holds a base other than A, C, G and T outside
        the PAM's N, does not begin with G, or does not end with ``NGG``.
    """
    if not isinstance(sequence, str):
        raise TypeError(
            f"A guide sequence must be a string, not {type(sequence).__name__} "
            f"({sequence!r})."
        )
    stripped = sequence.strip()
    if not stripped:
        raise ValueError("A guide sequence must not be empty.")
    normalized = stripped.upper()

    unknown = "".join(sorted(set(normalized) - set(_BASES + "N")))
    if unknown:
        raise ValueError(
            f"Guide sequence {sequence!r} holds {unknown!r}: a guide is written "
            f"with A, C, G and T, and the N of its {PAM} PAM."
        )
    if not normalized.endswith(PAM):
        raise ValueError(
            f"Guide sequence {sequence!r} does not end with the PAM {PAM!r}. "
            "Write the PAM as 'NGG' itself, not as the bases the genome holds "
            "there."
        )
    spacer = normalized[: -len(PAM)]
    if "N" in spacer:
        raise ValueError(
            f"Guide sequence {sequence!r} holds an N before the PAM: only the "
            "PAM's first base may be N."
        )
    if not spacer.startswith("G"):
        raise ValueError(
            f"Guide sequence {sequence!r} does not start with G. Write the whole "
            "guide, including the 5' G it carries or is cloned with."
        )
    return normalized


def guide_accession(sequence: str) -> str:
    """Return a guide's accession.

    Parameters
    ----------
    sequence
        A guide, ``G...NGG``, as `normalize_guide_sequence` takes it.

    Returns
    -------
    str
        ``GU`` and the encoded bases: 9 characters for a guide of 20 to 22
        bases before the PAM, 8 for 19.

    Raises
    ------
    TypeError
        If `sequence` is not a string.
    ValueError
        If `sequence` is not a whole guide written ``G...NGG``.

    See Also
    --------
    guide_sequence : The inverse.

    Examples
    --------
    >>> guide_accession("GAAAAGCCAACATGAATGCAGNGG")
    'GU602A170WJ'
    """
    spacer = normalize_guide_sequence(sequence)[: -len(PAM)]
    # The leading 1 distinguishes a guide from the same bases after an A.
    packed = 1
    for base in spacer:
        packed = packed * 4 + _BASES.index(base)

    characters = ""
    while packed:
        packed, remainder = divmod(packed, 32)
        characters = _ALPHABET[remainder] + characters
    return _PREFIX + characters


def guide_sequence(accession: str) -> str:
    """Return the guide an accession encodes.

    Parameters
    ----------
    accession
        An accession from `guide_accession`. Case is ignored, and I, L and O
        read as 1, 1 and 0.

    Returns
    -------
    str
        The guide, ``G...NGG``.

    Raises
    ------
    TypeError
        If `accession` is not a string.
    ValueError
        If `accession` does not begin with ``GU`` and base 32 characters, or
        encodes no guide.

    Examples
    --------
    >>> guide_sequence("GU602A170WJ")
    'GAAAAGCCAACATGAATGCAGNGG'
    """
    if not isinstance(accession, str):
        raise TypeError(
            f"An accession must be a string, not {type(accession).__name__} "
            f"({accession!r})."
        )
    text = accession.strip().upper()
    if not text.startswith(_PREFIX) or len(text) == len(_PREFIX):
        raise ValueError(
            f"Accession {accession!r} does not begin with {_PREFIX!r} and its "
            "encoded bases."
        )

    packed = 0
    for character in text[len(_PREFIX) :]:
        if character not in _VALUES:
            raise ValueError(
                f"Accession {accession!r} holds {character!r}, which is not a "
                "base 32 character."
            )
        packed = packed * 32 + _VALUES[character]

    spacer = ""
    while packed > 1:
        packed, remainder = divmod(packed, 4)
        spacer = _BASES[remainder] + spacer
    if packed != 1 or not spacer.startswith("G"):
        raise ValueError(f"Accession {accession!r} encodes no guide.")
    return spacer + PAM


def guide_accessions(sequences: Iterable[str]) -> list[str]:
    """Return every guide's accession, in order.

    Parameters
    ----------
    sequences
        Guides, as `guide_accession` takes them.

    Returns
    -------
    list of str
        One accession per sequence.

    Raises
    ------
    TypeError
        If a sequence is not a string. The message names its position.
    ValueError
        If a sequence is not a whole guide written ``G...NGG``. The message
        names its position.
    """
    accessions = []
    for position, sequence in enumerate(sequences):
        try:
            accessions.append(guide_accession(sequence))
        except (TypeError, ValueError) as error:
            raise type(error)(f"Sequence at position {position}: {error}") from None
    return accessions

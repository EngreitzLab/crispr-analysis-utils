"""Accession ids for guides, encoded from the protospacer.

An accession is ``GU`` followed by the protospacer packed two bits per base
behind a leading 1 bit, in Crockford's base 32. The packing is reversible, so
two different protospacers can never share an accession -- not in a library,
and not among every guide the lab will ever order.

The protospacer is the key, not the ordered sequence, so the same target gets
the same accession whatever a lab adds at the 5' end for its promoter and
whatever it calls the guide. The exact sequence that was ordered is
`guide_alignment`'s ``guide_id``; this answers the other question, which is
whether two guides hit the same place.
"""

from __future__ import annotations

from collections.abc import Iterable

_PREFIX = "GU"
_BASES = "ACGT"
# Crockford's base 32: the digits and the capitals, without I, L, O and U.
_ALPHABET = "0123456789ABCDEFGHJKMNPQRSTVWXYZ"
_VALUES = {character: value for value, character in enumerate(_ALPHABET)}
# Crockford reads these as the digits they look like.
_VALUES.update({"I": 1, "L": 1, "O": 0})


def normalize_spacer(spacer: str) -> str:
    """Strip and uppercase a protospacer, and check it is one.

    Parameters
    ----------
    spacer
        The protospacer, 5' to 3', without the PAM and without any base added
        5' of it for the promoter.

    Returns
    -------
    str
        `spacer` stripped and in uppercase.

    Raises
    ------
    TypeError
        If `spacer` is not a string, such as the NaN pandas reads for an empty
        cell.
    ValueError
        If `spacer` is empty or holds anything but A, C, G and T.

    Notes
    -----
    The PAM is genomic, so it is not part of a protospacer and must not be
    passed here. A sequence ending in its PAM is accepted and gives a
    different accession, silently: one spacer in sixteen ends in ``GG``, so
    there is no way to tell from the sequence alone. Decide which bases are
    the protospacer before calling this.
    """
    if not isinstance(spacer, str):
        raise TypeError(
            f"A protospacer must be a string, not {type(spacer).__name__} ({spacer!r})."
        )
    stripped = spacer.strip()
    if not stripped:
        raise ValueError("A protospacer must not be empty.")
    normalized = stripped.upper()
    unknown = "".join(sorted(set(normalized) - set(_BASES)))
    if unknown:
        raise ValueError(
            f"Protospacer {spacer!r} holds {unknown!r}: a protospacer is written "
            "with A, C, G and T."
        )
    return normalized


def guide_accession(spacer: str) -> str:
    """Return the accession of a protospacer.

    Parameters
    ----------
    spacer
        The protospacer, as `normalize_spacer` takes it.

    Returns
    -------
    str
        ``GU`` and the encoded bases: 9 characters for a 20-nt protospacer,
        8 for 19.

    Raises
    ------
    TypeError
        If `spacer` is not a string.
    ValueError
        If `spacer` is empty or holds anything but A, C, G and T.

    See Also
    --------
    spacer_of : The inverse.
    guide_accessions : Accession a whole library.

    Examples
    --------
    >>> guide_accession("AAAAGCCAACATGAATGCAG")
    'GU102A170WJ'
    >>> guide_accession(" aaaagccaacatgaatgcag ")
    'GU102A170WJ'
    """
    packed = 1
    for base in normalize_spacer(spacer):
        packed = packed * 4 + _BASES.index(base)

    characters = ""
    while packed:
        packed, remainder = divmod(packed, 32)
        characters = _ALPHABET[remainder] + characters
    return _PREFIX + characters


def spacer_of(accession: str) -> str:
    """Return the protospacer an accession encodes.

    Parameters
    ----------
    accession
        An accession from `guide_accession`. Case is ignored, and I, L and O
        read as 1, 1 and 0.

    Returns
    -------
    str
        The protospacer, 5' to 3'.

    Raises
    ------
    TypeError
        If `accession` is not a string.
    ValueError
        If `accession` does not begin with ``GU`` and base 32 characters, or
        encodes nothing.

    Examples
    --------
    >>> spacer_of("GU102A170WJ")
    'AAAAGCCAACATGAATGCAG'
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
    if packed != 1 or not spacer:
        raise ValueError(f"Accession {accession!r} encodes no protospacer.")
    return spacer


def guide_accessions(spacers: Iterable[str]) -> list[str]:
    """Return every protospacer's accession, in order.

    Repeats of one protospacer share an accession; different ones never do.

    Parameters
    ----------
    spacers
        Protospacers, as `guide_accession` takes them: a list, or a pandas
        Series.

    Returns
    -------
    list of str
        One accession per protospacer.

    Raises
    ------
    TypeError
        If a protospacer is not a string. The message names its position.
    ValueError
        If a protospacer is empty or holds anything but A, C, G and T. The
        message names its position.
    """
    accessions = []
    for position, spacer in enumerate(spacers):
        try:
            accessions.append(guide_accession(spacer))
        except (TypeError, ValueError) as error:
            raise type(error)(f"Spacer at position {position}: {error}") from None
    return accessions

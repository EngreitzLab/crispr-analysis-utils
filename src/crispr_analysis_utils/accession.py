"""Stable accession ids for guides, computed from their sequence.

A guide's accession is ``GU`` followed by 12 uppercase hexadecimal digits: the
first 12 of the UUID version 5 of its sequence in `GUIDE_NAMESPACE`. The same
sequence gets the same accession on any machine, in any library, whatever the
lab calls the guide.
"""

from __future__ import annotations

import re
import uuid
from collections.abc import Iterable

GUIDE_NAMESPACE = uuid.UUID("ee28d900-f51b-5e2b-9ca7-4c5808da45cc")
"""The UUID namespace of guide accessions.

It is ``uuid.uuid5(uuid.NAMESPACE_URL,
"https://github.com/EngreitzLab/crispr-analysis-utils/accession/guide")``,
written out so that no edit to that string can change every accession.
"""

_PREFIX = "GU"
# The first 12 hexadecimal digits of a UUID are all hash: its version digit is
# the 13th.
_DIGITS = 12
_ACGT = re.compile(r"[ACGT]+")


def normalize_guide_sequence(sequence: str) -> str:
    """Return the form of a guide sequence that its accession is computed from.

    Parameters
    ----------
    sequence
        A guide sequence, 5' to 3'. Case and surrounding whitespace are
        ignored.

    Returns
    -------
    str
        `sequence` without surrounding whitespace, in uppercase.

    Raises
    ------
    TypeError
        If `sequence` is not a string, such as the NaN pandas reads for an
        empty cell.
    ValueError
        If `sequence` is empty, or holds anything but A, C, G and T.
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
    if not (stripped.isascii() and _ACGT.fullmatch(normalized)):
        others = "".join(sorted(set(normalized) - set("ACGT")))
        raise ValueError(
            f"Guide sequence {sequence!r} holds {others!r}: only A, C, G and T "
            "are allowed."
        )
    return normalized


def guide_accession(sequence: str) -> str:
    """Return the accession of a guide.

    Pass the sequence that is aligned: 5' to 3', with the 5' G when the guide
    carries one, and without the PAM. A guide written with and without its G
    gets two accessions.

    Parameters
    ----------
    sequence
        The guide sequence. Case and surrounding whitespace are ignored.

    Returns
    -------
    str
        ``GU`` followed by 12 uppercase hexadecimal digits.

    Raises
    ------
    TypeError
        If `sequence` is not a string.
    ValueError
        If `sequence` is empty, or holds anything but A, C, G and T.

    See Also
    --------
    guide_accessions : The accessions of a whole library, checked for
        collisions.

    Examples
    --------
    >>> guide_accession("GAAAAGCCAACATGAATGCAG")
    'GU5C664EE08629'
    >>> guide_accession(" gaaaagccaacatgaatgcag ")
    'GU5C664EE08629'
    """
    return _accession(normalize_guide_sequence(sequence))


def guide_accessions(sequences: Iterable[str]) -> list[str]:
    """Return the accession of each guide, and check that none collide.

    Repeats of one sequence share an accession. Two different sequences with
    the same accession are a collision, which raises an error rather than let
    one accession name two guides. To check a library against accessions
    already in use, pass all their sequences in one call.

    Parameters
    ----------
    sequences
        Guide sequences, as `guide_accession` takes them: a list, or a pandas
        Series.

    Returns
    -------
    list of str
        One accession per sequence, in order.

    Raises
    ------
    TypeError
        If a sequence is not a string. The message gives its position.
    ValueError
        If a sequence is empty or holds anything but A, C, G and T (the
        message gives its position), or two different sequences have the same
        accession.
    """
    accessions = []
    first_sequence: dict[str, str] = {}
    for position, sequence in enumerate(sequences):
        try:
            normalized = normalize_guide_sequence(sequence)
        except (TypeError, ValueError) as error:
            raise type(error)(f"Sequence at position {position}: {error}") from None
        accession = _accession(normalized)
        first = first_sequence.setdefault(accession, normalized)
        if first != normalized:
            raise ValueError(
                f"Two different sequences have the accession {accession}: "
                f"{first} and {normalized}."
            )
        accessions.append(accession)
    return accessions


def _accession(normalized: str) -> str:
    digits = uuid.uuid5(GUIDE_NAMESPACE, normalized).hex[:_DIGITS]
    return _PREFIX + digits.upper()

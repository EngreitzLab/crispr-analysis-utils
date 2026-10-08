"""The guide library: reading the guide table and writing the reads to align.

A guide is known by its row: its 0-based position in the list `read_guides`
returns. Read names carry that row, so the alignments can be matched back to
the guide without parsing ids.
"""

from __future__ import annotations

import gzip
import io
import logging
import os
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import pandas as pd

from .iupac import validate_pams

logger = logging.getLogger(__name__)

_ACGT = frozenset("ACGT")
_TABLE_BREAKS = frozenset("\t\r\n")
_QUALITY = "I"
_MAX_LISTED = 5


@dataclass(frozen=True)
class Guide:
    """One guide of the library.

    Attributes
    ----------
    guide_id
        The id from the guide table, as a string.
    input_sequence
        The sequence from the guide table, stripped of surrounding whitespace
        and uppercased.
    spacer
        The 3'-most `spacer_length` bases of `input_sequence`, 5' to 3'.
    leading_g
        Whether this guide's site reads carry an extra 5' G.
    duplicate_of
        The id of the first guide in the table with the same spacer, or None
        for that first guide itself.
    """

    guide_id: str
    input_sequence: str
    spacer: str
    leading_g: bool
    duplicate_of: str | None


def read_guides(
    source: str | os.PathLike[str] | pd.DataFrame,
    *,
    id_col: int | str = 0,
    spacer_col: int | str = 1,
    sep: str = "\t",
    header: bool = False,
    spacer_length: int = 20,
    add_leading_g: bool = False,
) -> list[Guide]:
    """Read and validate a guide table.

    Parameters
    ----------
    source
        A delimited text file, or a pandas DataFrame. A file is read with every
        column as a string, so an id such as ``0001`` stays ``0001`` and an id
        such as ``NA`` is not taken for a missing value.
    id_col
        The guide id column: a 0-based position, or a column name.
    spacer_col
        The sequence column: a 0-based position, or a column name. Sequences
        are given 5' to 3'.
    sep
        The field separator of a file. Ignored for a DataFrame.
    header
        Whether the file's first line holds column names. Ignored for a
        DataFrame.
    spacer_length
        The spacer is the 3'-most `spacer_length` bases of each sequence.
    add_leading_g
        Mark a guide whose spacer does not start with G as carrying an extra
        5' G in its site reads (`Guide.leading_g`).

    Returns
    -------
    list of Guide
        One guide per table row, in table order. A guide whose spacer equals
        an earlier guide's has `duplicate_of` set to that earlier guide's id;
        a warning is logged when there are any.

    Raises
    ------
    ValueError
        If the table has no rows; a column is missing; an id or a sequence is
        missing; an id is empty or holds a tab or a line break; ids repeat; a
        sequence, once stripped and uppercased, holds a character other than
        A, C, G and T; or a sequence is shorter than `spacer_length`.
    """
    if (
        isinstance(spacer_length, bool)
        or not isinstance(spacer_length, int)
        or spacer_length < 1
    ):
        raise ValueError("spacer_length must be a positive integer.")

    if isinstance(source, pd.DataFrame):
        table = source
        origin = "the guide DataFrame"
    else:
        origin = os.fspath(source)
        try:
            table = pd.read_csv(
                source,
                sep=sep,
                header=0 if header else None,
                dtype=str,
                keep_default_na=False,
            )
        except pd.errors.EmptyDataError:
            raise ValueError(f"No guides in {origin}.") from None
    if table.empty:
        raise ValueError(f"No guides in {origin}.")

    raw_ids = _column(table, id_col, "id_col", origin)
    raw_sequences = _column(table, spacer_col, "spacer_col", origin)
    _check_present(raw_ids, "guide id", origin)
    _check_present(raw_sequences, "guide sequence", origin)

    ids = [str(value) for value in raw_ids]
    sequences = [str(value).strip().upper() for value in raw_sequences]

    empty = [str(row + 1) for row, guide_id in enumerate(ids) if guide_id == ""]
    if empty:
        raise ValueError(
            f"Guide ids must not be empty in {origin}: data rows "
            f"{_listing(empty)} (counted from 1)."
        )

    unsafe = [repr(guide_id) for guide_id in ids if _TABLE_BREAKS & set(guide_id)]
    if unsafe:
        raise ValueError(
            f"Guide ids must not hold tabs or line breaks in {origin}: "
            f"{_listing(unsafe)}."
        )

    repeated = [
        f"{guide_id!r} ({count} times)"
        for guide_id, count in Counter(ids).items()
        if count > 1
    ]
    if repeated:
        raise ValueError(f"Guide ids must be unique in {origin}: {_listing(repeated)}.")

    invalid = [
        f"{guide_id!r} ({''.join(sorted(set(sequence) - _ACGT))!r})"
        for guide_id, sequence in zip(ids, sequences, strict=True)
        if not set(sequence) <= _ACGT
    ]
    if invalid:
        raise ValueError(
            f"Guide sequences may hold only A, C, G and T in {origin}: "
            f"{_listing(invalid)}."
        )

    short = [
        f"{guide_id!r} ({len(sequence)} bases)"
        for guide_id, sequence in zip(ids, sequences, strict=True)
        if len(sequence) < spacer_length
    ]
    if short:
        raise ValueError(
            f"Guide sequences must be at least {spacer_length} bases long in "
            f"{origin}: {_listing(short)}."
        )

    first_by_spacer: dict[str, str] = {}
    guides = []
    for guide_id, sequence in zip(ids, sequences, strict=True):
        spacer = sequence[-spacer_length:]
        first = first_by_spacer.setdefault(spacer, guide_id)
        guides.append(
            Guide(
                guide_id=guide_id,
                input_sequence=sequence,
                spacer=spacer,
                leading_g=add_leading_g and not spacer.startswith("G"),
                duplicate_of=None if first == guide_id else first,
            )
        )

    duplicates = [
        f"{guide.guide_id!r} (same spacer as {guide.duplicate_of!r})"
        for guide in guides
        if guide.duplicate_of is not None
    ]
    if duplicates:
        logger.warning(
            "%d guides in %s share their spacer with an earlier guide; each "
            "spacer is aligned once and its results are copied to every id "
            "(see duplicate_of): %s",
            len(duplicates),
            origin,
            _listing(duplicates),
        )
    return guides


def write_spacer_reads(guides: Sequence[Guide], path: str | os.PathLike[str]) -> int:
    """Write the orientation-check reads: one bare spacer per distinct spacer.

    Parameters
    ----------
    guides
        The library, as returned by `read_guides`.
    path
        The FASTQ to write. A path ending in ``.gz`` is gzip-compressed.

    Returns
    -------
    int
        The number of reads written.

    Notes
    -----
    Each guide without `duplicate_of` gives one read: its name is the guide's
    row, its sequence the spacer (no PAM, no added G), and every quality is
    ``I``.
    """
    n_reads = 0
    with _open_fastq(path) as handle:
        for row, guide in enumerate(guides):
            if guide.duplicate_of is None:
                _write_read(handle, str(row), guide.spacer)
                n_reads += 1
    return n_reads


def write_site_reads(
    guides: Sequence[Guide],
    path: str | os.PathLike[str],
    *,
    pams: Iterable[str],
    leading_g: bool,
) -> int:
    """Write the site-search reads: one read per guide and PAM.

    Parameters
    ----------
    guides
        The library, as returned by `read_guides`.
    path
        The FASTQ to write. A path ending in ``.gz`` is gzip-compressed.
    pams
        The PAM patterns, primary first, as IUPAC codes (for example
        ``["NGG", "NAG", "NGA"]``). They must share one length.
    leading_g
        Write the guides whose `Guide.leading_g` is this value: True for the
        reads that carry an extra 5' G, False for those that do not.

    Returns
    -------
    int
        The number of reads written.

    Raises
    ------
    ValueError
        If the PAM patterns fail `validate_pams`.

    Notes
    -----
    Only guides without `duplicate_of` are written. A guide's reads are
    consecutive, one per PAM in the order given. Each read is named
    ``<row>:<pam>`` and its sequence is ``G`` (when `leading_g` is True), the
    spacer, then the PAM pattern; the PAM is uppercased in both, and
    otherwise written as given, IUPAC codes included. Every quality is ``I``.
    """
    patterns = validate_pams(pams)
    prefix = "G" if leading_g else ""
    n_reads = 0
    with _open_fastq(path) as handle:
        for row, guide in enumerate(guides):
            if guide.duplicate_of is not None or guide.leading_g != leading_g:
                continue
            for pam in patterns:
                _write_read(handle, f"{row}:{pam}", f"{prefix}{guide.spacer}{pam}")
                n_reads += 1
    return n_reads


def _column(
    table: pd.DataFrame, column: int | str, argument: str, origin: str
) -> pd.Series:
    """One column of the guide table, by 0-based position or by name."""
    if isinstance(column, bool) or not isinstance(column, (int, str)):
        raise ValueError(f"{argument} must be a column position or a column name.")
    if isinstance(column, int):
        if not 0 <= column < table.shape[1]:
            raise ValueError(
                f"{argument}={column} is out of range: {origin} has "
                f"{table.shape[1]} columns."
            )
        return table.iloc[:, column]
    if column not in table.columns:
        raise ValueError(f"{argument}={column!r} is not a column of {origin}.")
    return table[column]


def _check_present(values: pd.Series, what: str, origin: str) -> None:
    """Raise when any value is missing, naming the 1-based data rows."""
    missing = [str(row + 1) for row, absent in enumerate(values.isna()) if absent]
    if missing:
        raise ValueError(
            f"Missing {what} in {origin}: data rows {_listing(missing)} "
            "(counted from 1)."
        )


def _listing(items: Sequence[str]) -> str:
    """A comma-separated list, cut after the first few items."""
    shown = ", ".join(items[:_MAX_LISTED])
    if len(items) > _MAX_LISTED:
        shown += f", and {len(items) - _MAX_LISTED} more"
    return shown


def _open_fastq(path: str | os.PathLike[str]) -> TextIO:
    """Open a FASTQ for writing; gzip-compressed when the path ends in .gz."""
    path = Path(path)
    if path.suffix == ".gz":
        return io.TextIOWrapper(
            gzip.GzipFile(path, mode="wb", mtime=0), encoding="ascii", newline="\n"
        )
    return path.open("w", encoding="ascii", newline="\n")


def _write_read(handle: TextIO, name: str, sequence: str) -> None:
    """Write one FASTQ record with a constant quality."""
    handle.write(f"@{name}\n{sequence}\n+\n{_QUALITY * len(sequence)}\n")

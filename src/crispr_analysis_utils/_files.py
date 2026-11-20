"""File helpers for the pipeline."""

from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Sequence
from pathlib import Path

_GZIP_MAGIC = b"\x1f\x8b"
_SEPARATORS = {".csv": ",", ".tsv": "\t", ".txt": "\t"}


def file_md5(path: str | os.PathLike[str]) -> str:
    """The MD5 checksum of a file, as 32 hexadecimal digits."""
    with open(path, "rb") as handle:
        digest = hashlib.file_digest(handle, lambda: hashlib.md5(usedforsecurity=False))
    return digest.hexdigest()


def is_gzip(path: str | os.PathLike[str]) -> bool:
    """Whether a file starts with the gzip magic bytes."""
    with open(path, "rb") as handle:
        return handle.read(2) == _GZIP_MAGIC


def check_plain_fasta(path: str | os.PathLike[str], what: str) -> Path:
    """The path of an existing, uncompressed FASTA, or an error naming `what`."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"The {what} {path} does not exist.")
    if is_gzip(path):
        raise ValueError(
            f"The {what} {path} is gzip-compressed; the scan reads an "
            "uncompressed FASTA."
        )
    return path


def write_json(path: str | os.PathLike[str], data: object) -> None:
    """Write JSON through a temporary file and a rename, so no reader sees half."""
    path = Path(path)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def read_columns(path: str | os.PathLike[str], columns: Sequence[str], what: str):
    """The named columns of a CSV or TSV, as text, without the rows missing one.

    The suffix says the separator: ``.csv`` a comma, ``.tsv`` and ``.txt`` a
    tab, either optionally gzipped. Every value is read as text and stripped;
    a row with an empty cell in any named column is dropped.

    Parameters
    ----------
    path
        The table.
    columns
        The header names to keep, in this order.
    what
        What the table is, for the error messages.

    Returns
    -------
    pandas.DataFrame
        Exactly `columns`.

    Raises
    ------
    FileNotFoundError
        The file does not exist.
    ValueError
        The suffix names no separator, or a column is missing.
    """
    import pandas as pd

    path = Path(path)
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes[-1:] == [".gz"]:
        suffixes.pop()
    sep = _SEPARATORS.get(suffixes[-1] if suffixes else "")
    if sep is None:
        raise ValueError(
            f"The {what} {path} must end in .csv, .tsv or .txt (optionally .gz),"
            " which says how it is separated."
        )
    if not path.is_file():
        raise FileNotFoundError(f"The {what} {path} does not exist.")
    frame = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(
            f"The {what} {path} has no {', '.join(missing)} column; its columns"
            f" are {', '.join(map(str, frame.columns))}."
        )
    frame = frame[list(columns)].apply(lambda column: column.str.strip())
    return frame[(frame != "").all(axis=1)].reset_index(drop=True)

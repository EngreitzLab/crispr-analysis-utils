"""File helpers for the pipeline, and the run record every step keeps."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import time
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from datetime import UTC, datetime
from itertools import count
from pathlib import Path

from ._version import __version__

_GZIP_MAGIC = b"\x1f\x8b"
_SEPARATORS = {".csv": ",", ".tsv": "\t", ".txt": "\t"}

METADATA_DIR = "run_metadata"
"""Where every step writes the record of each invocation."""

_PACKAGE = __name__.rsplit(".", 1)[0]
_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"


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


def truthy(value: object) -> bool:
    """TRUE, true and True all mean the same in a written table; a missing cell does not."""
    import pandas as pd

    if isinstance(value, str):
        return value.strip().upper() == "TRUE"
    if value is None or pd.isna(value):
        return False
    return bool(value)


def yaml_text(text: str) -> str:
    """`text` on one line, as a double-quoted YAML scalar.

    A JSON string is a valid double-quoted YAML scalar, quotes and backslashes
    escaped.
    """
    return json.dumps(" ".join(text.split()), ensure_ascii=False)


def now() -> str:
    """The time now, in UTC, as ISO 8601 to the second."""
    return datetime.now(UTC).isoformat(timespec="seconds")


@contextmanager
def timed(state: dict, name: str) -> Iterator[None]:
    """Record in ``state["seconds"][name]`` how long the block took."""
    started = time.perf_counter()
    try:
        yield
    finally:
        state["seconds"][name] = round(time.perf_counter() - started, 3)


@contextmanager
def run_record(
    outdir: Path,
    step: str,
    results: str,
    parameters: dict,
    logger: logging.Logger,
    title: str,
) -> Iterator[tuple[dict, Path]]:
    """One invocation of a step: its folders, its log and its record.

    Makes ``outdir/results``, ``outdir/logs`` and ``outdir/run_metadata``,
    appends the package's log from INFO up to ``logs/<step>.log``, and keeps
    ``run_metadata/<step>_<UTC timestamp>.json`` from the start. The record is
    rewritten complete when the block ends, or with the error when it raises,
    which is logged and raised again.

    Parameters
    ----------
    outdir
        The step's output folder.
    step
        The step's name, which names its log and its records.
    results
        The folder inside `outdir` the results go in.
    parameters
        The parameters to record.
    logger
        The step's logger, which logs the failure.
    title
        The step as the failure line names it, such as ``Element mapping``.

    Yields
    ------
    tuple[dict, Path]
        The record, to fill in, and its path, for a step that saves it
        before it ends.
    """
    for folder in (results, "logs", METADATA_DIR):
        (outdir / folder).mkdir(parents=True, exist_ok=True)
    state: dict = {
        "step": step,
        "complete": False,
        "error": None,
        "started": now(),
        "finished": None,
        "cau_version": __version__,
        "inputs": {},
        "parameters": parameters,
        "counts": {},
        "seconds": {},
    }
    record = _reserve(outdir / METADATA_DIR, step, state)
    started = time.perf_counter()

    handler = logging.FileHandler(outdir / "logs" / f"{step}.log", encoding="utf-8")
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(_LOG_FORMAT))
    package = logging.getLogger(_PACKAGE)
    level, propagate = package.level, package.propagate
    if package.getEffectiveLevel() > logging.INFO:
        package.setLevel(logging.INFO)
    package.addHandler(handler)
    # The run's log belongs in the file, not on the terminal. pyranges calls
    # logging.basicConfig on import, which puts a handler on the root logger.
    package.propagate = False
    try:
        yield state, record
        state["complete"] = True
    except BaseException as error:
        state["error"] = f"{type(error).__name__}: {error}"
        logger.error("%s failed: %s", title, state["error"])
        raise
    finally:
        state["finished"] = now()
        state["seconds"]["total"] = round(time.perf_counter() - started, 3)
        write_json(record, state)
        package.removeHandler(handler)
        package.setLevel(level)
        package.propagate = propagate
        handler.close()


def _reserve(folder: Path, step: str, state: dict) -> Path:
    """A record path no other invocation holds, created with `state` in it.

    Two runs in the same second get ``<step>_<stamp>.json`` and
    ``<step>_<stamp>_2.json``, which sorts after it.
    """
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    # shortcut: _10 sorts before _2; ten runs into one folder in one second is not a use.
    for n in count(1):
        path = folder / (
            f"{step}_{stamp}.json" if n == 1 else f"{step}_{stamp}_{n}.json"
        )
        try:
            with open(path, "x", encoding="utf-8") as handle:
                handle.write(json.dumps(state, indent=2) + "\n")
            return path
        except FileExistsError:
            pass

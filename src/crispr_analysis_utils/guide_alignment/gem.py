"""GEM3: building the index and mapping the guide reads.

`build_index` runs gem-indexer and writes a sidecar with the reference's
checksum next to the index. `map_reads` runs gem-mapper with fixed options
and a search budget set by the mismatches allowed, the ambiguous PAM bases and
the added G (`mapper_options`), and writes the alignments it can place inside
their contig. The guide alignment docs explain each option.
"""

from __future__ import annotations

import gzip
import json
import logging
import os
import re
import shlex
import subprocess
import tempfile
import time
from collections import Counter
from collections.abc import Iterable
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import BinaryIO

from ._files import check_plain_fasta, file_md5, write_json

logger = logging.getLogger(__name__)

MAPPER_OPTIONS: tuple[str, ...] = (
    "--mapping-mode=customed",
    "--alignment-model=hamming",
    "--clipping=none",
    "--alignment-local=never",
    "--sam-compact=false",
    "--max-reported-matches=all",
    "--alignment-max-bandwidth=0",
)
"""The gem-mapper options every pass uses, before its search budget."""

INDEX_SUFFIX = ".gem"
"""The extension gem-indexer gives the index."""
SIDECAR_SUFFIX = ".json"
"""Appended to the index file name for its sidecar."""

_SIGNAL = "Signal raised"
_CRASH_HINT = (
    "GEM 3.6 can crash on a read that aligns across the start of its index, "
    "the first base of the reference's first contig on the + strand: a read "
    "with an added 5' G does whenever its protospacer starts at that base. A "
    "first contig that starts with N, as in hg38, cannot trigger it. Index a "
    "reference whose first contig starts with N (reorder or drop contigs), or "
    "align without the added G."
)
_CIGAR_OP = re.compile(rb"(\d+)([MIDNSHP=X])")
_REFERENCE_OPS = frozenset(b"MDN=X")
_LOG_TAIL = 5


class GemError(RuntimeError):
    """gem-indexer or gem-mapper is missing, failed, or crashed.

    Attributes
    ----------
    log
        The log of the failed run, or None.
    read_name
        The read gem-mapper was mapping when it crashed, when its log names
        one; None otherwise.
    """

    def __init__(
        self,
        message: str,
        *,
        log: Path | None = None,
        read_name: str | None = None,
    ) -> None:
        super().__init__(message)
        self.log = log
        self.read_name = read_name


@dataclass(frozen=True)
class MappingResult:
    """One gem-mapper run, and what was kept of its output.

    Attributes
    ----------
    argv
        The command that ran.
    sam
        The alignments written.
    log
        gem-mapper's log.
    seconds
        The wall time of the run.
    n_records
        Alignment records gem-mapper wrote, header lines excluded.
    n_unmapped
        Records of unmapped reads, kept.
    n_out_of_range
        Records of mapped reads whose position is not a number, or whose
        aligned span does not lie inside their contig, dropped.
    mapped_per_read
        The kept records of mapped reads, by read name.
    """

    argv: tuple[str, ...]
    sam: Path
    log: Path
    seconds: float
    n_records: int
    n_unmapped: int
    n_out_of_range: int
    mapped_per_read: dict[str, int] = field(repr=False)

    def counts(self) -> dict[str, int]:
        """The record counts by attribute name."""
        return {
            "n_records": self.n_records,
            "n_unmapped": self.n_unmapped,
            "n_out_of_range": self.n_out_of_range,
        }


def gem_version(executable: str = "gem-mapper") -> str:
    """The version a GEM executable reports, such as ``v3.6.0-bundle-release``.

    Parameters
    ----------
    executable
        gem-mapper or gem-indexer, by name on the PATH or by path.

    Returns
    -------
    str
        The first line the executable prints for ``--version``.

    Raises
    ------
    GemError
        If the executable is not found, or fails to report a version.
    """
    try:
        completed = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            text=True,
            check=False,
            errors="replace",
        )
    except FileNotFoundError:
        raise GemError(_not_found(executable)) from None
    output = (completed.stdout + completed.stderr).strip()
    if completed.returncode != 0 or not output:
        raise GemError(
            f"{executable} --version failed (exit status {completed.returncode})."
        )
    return output.splitlines()[0].strip()


def index_file(index: str | os.PathLike[str]) -> Path:
    """The index file for an index prefix or an index path.

    Parameters
    ----------
    index
        ``<prefix>`` or ``<prefix>.gem``.

    Returns
    -------
    Path
        ``<prefix>.gem``.
    """
    path = Path(index)
    if path.suffix == INDEX_SUFFIX:
        return path
    return path.with_name(path.name + INDEX_SUFFIX)


def sidecar_path(index: str | os.PathLike[str]) -> Path:
    """The sidecar of an index: ``<prefix>.gem.json``."""
    path = index_file(index)
    return path.with_name(path.name + SIDECAR_SUFFIX)


def read_sidecar(index: str | os.PathLike[str]) -> dict | None:
    """The sidecar `build_index` wrote next to an index, or None if there is none.

    Raises
    ------
    ValueError
        If the sidecar is not a JSON object.
    """
    path = sidecar_path(index)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(
            f"The index sidecar {path} is not valid JSON: {error}"
        ) from None
    if not isinstance(data, dict):
        raise ValueError(f"The index sidecar {path} is not a JSON object.")
    return data


def build_index(
    reference_fasta: str | os.PathLike[str],
    prefix: str | os.PathLike[str],
    *,
    threads: int | None = None,
    gem_indexer: str = "gem-indexer",
) -> Path:
    """Build a GEM index from a reference FASTA, with a checksum sidecar.

    Parameters
    ----------
    reference_fasta
        The reference, an uncompressed FASTA. Contig names are the first word
        of each header line.
    prefix
        The index path without its extension. A trailing ``.gem`` is dropped.
        Missing folders are created.
    threads
        Threads for gem-indexer; None leaves its default, every core.
    gem_indexer
        The gem-indexer executable, by name on the PATH or by path.

    Returns
    -------
    Path
        The index, ``<prefix>.gem``.

    Raises
    ------
    FileNotFoundError
        If the reference does not exist.
    ValueError
        If the reference is gzip-compressed, or `threads` is not a positive
        integer.
    GemError
        If gem-indexer is not found or fails. The partial index is removed.

    Notes
    -----
    Next to the index go gem-indexer's ``<prefix>.info``, its output in
    ``<prefix>.log``, and the sidecar ``<prefix>.gem.json``: the reference's
    path, size and MD5 checksum, the GEM version, the command and its wall
    time. gem-indexer's temporary files go in a folder next to the index,
    removed when it finishes.
    """
    reference = check_plain_fasta(reference_fasta, "reference").resolve()
    if threads is not None:
        _check_positive(threads, "threads")
    index = index_file(prefix).resolve()
    stem = index.with_suffix("")
    folder = index.parent
    folder.mkdir(parents=True, exist_ok=True)
    sidecar = sidecar_path(index)
    sidecar.unlink(missing_ok=True)
    log = stem.with_name(stem.name + ".log")
    version = gem_version(gem_indexer)

    logger.info("Indexing %s into %s", reference, index)
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix=".gem-indexer-", dir=folder) as scratch:
        # gem-indexer joins this folder and a file name with no separator.
        argv = [
            gem_indexer,
            "-i",
            str(reference),
            "-o",
            str(stem),
            "--tmp-folder",
            scratch + os.sep,
        ]
        if threads is not None:
            argv += ["-t", str(threads)]
        returncode = _run_logged(argv, log, cwd=folder)
    seconds = time.perf_counter() - started

    text = log.read_text(encoding="utf-8", errors="replace")
    if returncode != 0 or _SIGNAL in text:
        index.unlink(missing_ok=True)
        raise GemError(
            f"gem-indexer failed on {reference} (exit status {returncode}); see "
            f"{log}:\n{_tail(text)}",
            log=log,
        )
    write_json(
        sidecar,
        {
            "reference": str(reference),
            "reference_size": reference.stat().st_size,
            "reference_md5": file_md5(reference),
            "gem_version": version,
            "argv": argv,
            "seconds": round(seconds, 3),
            "created": datetime.now(UTC).isoformat(timespec="seconds"),
        },
    )
    logger.info("Built %s in %.1f s", index, seconds)
    return index


def mapper_options(
    *, max_mismatches: int, n_ambiguous: int = 0, leading_g: bool = False
) -> list[str]:
    """The gem-mapper options of one pass: `MAPPER_OPTIONS`, then its budget.

    Parameters
    ----------
    max_mismatches
        The spacer mismatches to find, k; 0 for exact matches.
    n_ambiguous
        The highest number of non-ACGT codes among the PAMs of the reads, n.
    leading_g
        Whether the reads carry an added 5' G, g.

    Returns
    -------
    list of str
        `MAPPER_OPTIONS`, then ``--alignment-max-error=k+g``,
        ``--complete-search-error=k+n+g`` and
        ``--complete-strata-after-best=k+n+g``.

    Raises
    ------
    ValueError
        If `max_mismatches` or `n_ambiguous` is not a non-negative integer.
    """
    _check_non_negative(max_mismatches, "max_mismatches")
    _check_non_negative(n_ambiguous, "n_ambiguous")
    g = int(bool(leading_g))
    errors = max_mismatches + g
    budget = max_mismatches + n_ambiguous + g
    return [
        *MAPPER_OPTIONS,
        f"--alignment-max-error={errors}",
        f"--complete-search-error={budget}",
        f"--complete-strata-after-best={budget}",
    ]


def mapper_argv(
    index: str | os.PathLike[str],
    fastq: str | os.PathLike[str],
    *,
    max_mismatches: int,
    n_ambiguous: int = 0,
    leading_g: bool = False,
    threads: int = 1,
    gem_mapper: str = "gem-mapper",
) -> list[str]:
    """The gem-mapper command of one pass, writing SAM to standard output.

    Parameters
    ----------
    index
        The index file, ``<prefix>.gem``.
    fastq
        The reads.
    max_mismatches
        The spacer mismatches to find, k (see `mapper_options`).
    n_ambiguous
        The highest number of non-ACGT codes among the PAMs of the reads, n.
    leading_g
        Whether the reads carry an added 5' G, g.
    threads
        gem-mapper threads.
    gem_mapper
        The gem-mapper executable, by name on the PATH or by path.

    Returns
    -------
    list of str
        ``gem-mapper -I <index> -i <fastq> -t <threads>``, then the options.

    Raises
    ------
    ValueError
        If a budget value or `threads` is out of range.
    """
    _check_positive(threads, "threads")
    return [
        gem_mapper,
        "-I",
        os.fspath(index),
        "-i",
        os.fspath(fastq),
        "-t",
        str(threads),
        *mapper_options(
            max_mismatches=max_mismatches,
            n_ambiguous=n_ambiguous,
            leading_g=leading_g,
        ),
    ]


def map_reads(
    index: str | os.PathLike[str],
    fastq: str | os.PathLike[str],
    sam: str | os.PathLike[str],
    *,
    log: str | os.PathLike[str],
    max_mismatches: int,
    n_ambiguous: int = 0,
    leading_g: bool = False,
    threads: int = 1,
    gem_mapper: str = "gem-mapper",
) -> MappingResult:
    """Map reads with gem-mapper and write the alignments it places correctly.

    gem-mapper's SAM output is streamed: header lines and records of unmapped
    reads are kept, and so is each record of a mapped read whose aligned span
    lies inside its contig. Any other record (gem-mapper can place a read that
    overhangs a contig end at a position past the contig, or at a position
    that is not a valid SAM position) is dropped and counted.

    Parameters
    ----------
    index
        The index, ``<prefix>.gem`` or ``<prefix>``.
    fastq
        The reads, FASTQ (gem-mapper also reads it gzip-compressed).
    sam
        The SAM to write; gzip-compressed if the name ends in ``.gz``. It is
        written under a temporary name and renamed when gem-mapper succeeds.
    log
        gem-mapper's log: the command, then everything it prints to standard
        error. Missing folders are created.
    max_mismatches
        The spacer mismatches to find, k (see `mapper_options`).
    n_ambiguous
        The highest number of non-ACGT codes among the PAMs of the reads, n.
    leading_g
        Whether the reads carry an added 5' G, g.
    threads
        gem-mapper threads.
    gem_mapper
        The gem-mapper executable, by name on the PATH or by path.

    Returns
    -------
    MappingResult
        The command, the run time and the record counts.

    Raises
    ------
    FileNotFoundError
        If the index or the reads do not exist.
    ValueError
        If a budget value or `threads` is out of range, or a mapped record
        names a contig the SAM header does not list.
    GemError
        If gem-mapper is not found, exits with a non-zero status, or reports
        a signal (a crash) in its log. No SAM is left behind.
    """
    index = index_file(index).resolve()
    fastq = Path(fastq).resolve()
    for path, what in ((index, "index"), (fastq, "FASTQ")):
        if not path.is_file():
            raise FileNotFoundError(f"The {what} {path} does not exist.")
    argv = mapper_argv(
        index,
        fastq,
        max_mismatches=max_mismatches,
        n_ambiguous=n_ambiguous,
        leading_g=leading_g,
        threads=threads,
        gem_mapper=gem_mapper,
    )
    sam = Path(sam)
    log = Path(log)
    sam.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    partial = sam.with_name(sam.name + ".partial")

    started = time.perf_counter()
    try:
        with log.open("w", encoding="utf-8") as log_handle:
            log_handle.write(f"$ {shlex.join(argv)}\n")
            log_handle.flush()
            try:
                process = subprocess.Popen(
                    argv,
                    stdout=subprocess.PIPE,
                    stderr=log_handle,
                    cwd=log.parent,
                )
            except FileNotFoundError:
                raise GemError(_not_found(gem_mapper)) from None
            with process, _open_output(partial) as output:
                try:
                    kept = _filter_sam(process.stdout, output)
                except BaseException:
                    process.kill()
                    raise
            returncode = process.returncode
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    seconds = time.perf_counter() - started

    text = log.read_text(encoding="utf-8", errors="replace")
    if returncode != 0 or _SIGNAL in text:
        partial.unlink(missing_ok=True)
        raise _mapping_error(fastq, log, returncode, text)
    partial.replace(sam)
    result = MappingResult(
        argv=tuple(argv),
        sam=sam,
        log=log,
        seconds=seconds,
        n_records=kept.n_records,
        n_unmapped=kept.n_unmapped,
        n_out_of_range=kept.n_out_of_range,
        mapped_per_read=kept.mapped_per_read,
    )
    logger.info(
        "gem-mapper mapped %s in %.1f s: %s", fastq.name, seconds, result.counts()
    )
    if result.n_out_of_range:
        logger.warning(
            "Dropped %d records of %s that gem-mapper placed outside their contig.",
            result.n_out_of_range,
            fastq.name,
        )
    return result


@dataclass
class _Kept:
    """Counts from `_filter_sam`."""

    n_records: int = 0
    n_unmapped: int = 0
    n_out_of_range: int = 0
    mapped_per_read: dict[str, int] = field(default_factory=dict)


def _filter_sam(lines: Iterable[bytes], output: BinaryIO) -> _Kept:
    """Copy SAM lines, dropping mapped records that do not fit inside their contig."""
    kept = _Kept()
    lengths: dict[bytes, int] = {}
    spans: dict[bytes, int] = {}
    per_read: Counter[bytes] = Counter()
    for line in lines:
        if line.startswith(b"@"):
            if line.startswith(b"@SQ\t"):
                name, length = _sq_fields(line)
                lengths[name] = length
            output.write(line)
            continue
        fields = line.split(b"\t", 6)
        if len(fields) < 6:
            raise ValueError(f"gem-mapper wrote a malformed SAM record: {line!r}")
        kept.n_records += 1
        if int(fields[1]) & 4:
            kept.n_unmapped += 1
            output.write(line)
            continue
        length = lengths.get(fields[2])
        if length is None:
            raise ValueError(
                f"gem-mapper placed read {fields[0].decode(errors='replace')!r} on "
                f"contig {fields[2].decode(errors='replace')!r}, which its SAM "
                "header does not list."
            )
        cigar = fields[5]
        span = spans.get(cigar)
        if span is None:
            span = spans[cigar] = _reference_span(cigar)
        try:
            position = int(fields[3])
        except ValueError:
            position = 0
        if position < 1 or position - 1 + span > length:
            kept.n_out_of_range += 1
            continue
        per_read[fields[0]] += 1
        output.write(line)
    kept.mapped_per_read = {
        name.decode(errors="replace"): count for name, count in per_read.items()
    }
    return kept


def _sq_fields(line: bytes) -> tuple[bytes, int]:
    """The contig name and length of an @SQ header line."""
    tags = dict(
        (field_[:2], field_[3:]) for field_ in line.rstrip(b"\r\n").split(b"\t")[1:]
    )
    try:
        return tags[b"SN"], int(tags[b"LN"])
    except (KeyError, ValueError):
        raise ValueError(f"gem-mapper wrote a malformed @SQ line: {line!r}") from None


def _reference_span(cigar: bytes) -> int:
    """The reference bases a CIGAR string covers (0 for ``*``)."""
    return sum(
        int(length)
        for length, op in _CIGAR_OP.findall(cigar)
        if op[0] in _REFERENCE_OPS
    )


@contextmanager
def _open_output(path: Path):
    """Open a SAM for binary writing; gzip-compressed if the name ends in .gz."""
    compressed = path.name.removesuffix(".partial").endswith(".gz")
    with path.open("wb") as raw:
        if compressed:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw, compresslevel=1, mtime=0
            ) as handle:
                yield handle
        else:
            yield raw


def _mapping_error(fastq: Path, log: Path, returncode: int, text: str) -> GemError:
    """The error for a failed gem-mapper run, naming the read it crashed on."""
    if _SIGNAL not in text:
        return GemError(
            f"gem-mapper failed on {fastq} (exit status {returncode}); see "
            f"{log}:\n{_tail(text)}",
            log=log,
        )
    read_name = _crashed_read(text)
    on_read = f" while mapping read {read_name!r}" if read_name else ""
    return GemError(
        f"gem-mapper crashed on {fastq}{on_read} (exit status {returncode}, "
        f"'{_SIGNAL}'); see {log}. {_CRASH_HINT}",
        log=log,
        read_name=read_name,
    )


def _crashed_read(text: str) -> str | None:
    """The read gem-mapper's crash report names, if any."""
    report = text.find("GEM::Input.State")
    if report < 0:
        return None
    match = re.search(r"^@(\S+)", text[report:], flags=re.MULTILINE)
    return match.group(1) if match else None


def _run_logged(argv: list[str], log: Path, *, cwd: Path) -> int:
    """Run a command with its output going to `log`; return its exit status."""
    with log.open("w", encoding="utf-8") as handle:
        handle.write(f"$ {shlex.join(argv)}\n")
        handle.flush()
        try:
            completed = subprocess.run(
                argv, stdout=handle, stderr=subprocess.STDOUT, cwd=cwd, check=False
            )
        except FileNotFoundError:
            raise GemError(_not_found(argv[0])) from None
    return completed.returncode


def _tail(text: str) -> str:
    """The last few non-empty lines of a log."""
    lines = [line for line in text.splitlines() if line.strip()]
    return "\n".join(lines[-_LOG_TAIL:])


def _not_found(executable: str) -> str:
    return (
        f"{executable} was not found. GEM3 comes with this project's pixi "
        "environments, or from bioconda as gem3-mapper."
    )


def _check_non_negative(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer.")


def _check_positive(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer.")

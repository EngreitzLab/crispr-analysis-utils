"""The guide alignment: read a guide table, scan a reference, write the tables.

One pass of the compiled core finds every exact site of every guide, on both
strands, and reads the PAM 3' of each and the base 5' of it. There is no
aligner, no index and no error budget: a guide is designed, not sequenced, so
its own site is an exact match, and everything the tables report follows from
that.

`outdir` receives ``guides_metadata.tsv``, ``sites.tsv``, ``orientation.tsv``,
``guides.bed``, ``cut_sites.bed`` and ``qc.tsv``, with this invocation's
record under ``run_metadata/`` and the log under ``logs/``.
"""

from __future__ import annotations

import logging
import os
import time
from collections import Counter
from collections.abc import Collection, Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .._version import __version__
from ._files import check_plain_fasta, file_md5, write_json
from .iupac import reverse_complement, validate_pams
from .library import Guide, read_guides
from .orientation import FORWARD, REVERSE_COMPLEMENT, VERDICTS
from .scan import ScanResult, scan_guides
from .summary import (
    CLASSES,
    GuideSummary,
    summarize,
    write_cut_sites_bed,
    write_guides_bed,
    write_guides_tsv,
    write_orientation_tsv,
    write_qc_tsv,
    write_sites_tsv,
)

logger = logging.getLogger(__name__)

METADATA_DIR = "run_metadata"
"""Where each invocation leaves its metadata, one file per invocation."""
STEP = "guide_alignment"
"""The step those metadata files are named for."""
LOG_FILE = "guide_alignment.log"
_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"

ORIENTATION_CHECKS = ("warn", "error", "off")
"""The values of `run`'s ``orientation_check``."""


class OrientationError(RuntimeError):
    """The library reads 3' to 5' and `orientation_check` is ``error``."""


@dataclass
class RunResult:
    """What `run` produced.

    Attributes
    ----------
    outdir
        The output folder.
    guides
        The library, as returned by `read_guides`.
    summaries
        One per guide, in library order, as returned by `summarize`.
    scan
        The genome scan: each guide's sites, its orientation, and the counts.
    metadata
        The record written under ``run_metadata/``.
    """

    outdir: Path
    guides: list[Guide]
    summaries: list[GuideSummary]
    scan: ScanResult
    metadata: dict = field(default_factory=dict)


def run(
    guides: str | os.PathLike[str] | pd.DataFrame,
    reference: str | os.PathLike[str],
    outdir: str | os.PathLike[str],
    *,
    pam: str = "NGG",
    spacer_length: int = 20,
    add_leading_g: bool = False,
    orientation_check: str = "error",
    fix_orientation: bool = False,
    contigs: Collection[str] | None = None,
    id_col: int | str = 0,
    spacer_col: int | str = 1,
    sep: str = "\t",
    header: bool = False,
) -> RunResult:
    """Find every guide's genomic site and write the tables.

    Parameters
    ----------
    guides
        The guide table, as `read_guides` takes it.
    reference
        The reference FASTA, uncompressed.
    outdir
        Where the tables go. Created if it does not exist; nothing is ever
        written outside it.
    pam
        The PAM, in IUPAC codes, 3' of the protospacer.
    spacer_length
        See `read_guides`.
    add_leading_g
        See `read_guides`.
    orientation_check
        ``error``: stop when entries reading 3' to 5' outnumber those reading
        5' to 3', which is a table given reverse-complemented throughout. A
        handful is reported either way, since one guide can look reversed for
        want of nothing more than a PAM that differs between assemblies.
        ``warn``: only report. ``off``: skip the check.
    fix_orientation
        Align the reverse complement of the entries that read 3' to 5',
        instead of leaving them as given.
    contigs
        If given, read only these contigs.
    id_col, spacer_col, sep, header
        See `read_guides`.

    Returns
    -------
    RunResult
        The library, the summaries, the scan and this invocation's metadata.

    Raises
    ------
    FileNotFoundError
        If `reference` does not exist.
    ValueError
        If `reference` is compressed, `pam` is not IUPAC codes,
        `orientation_check` is not one of `ORIENTATION_CHECKS`, or the guide
        table fails `read_guides`.
    OrientationError
        If the library reads 3' to 5' and `orientation_check` is ``error``.
    """
    if isinstance(contigs, str):
        raise TypeError("contigs must be a collection of contig names, not a string.")
    if orientation_check not in ORIENTATION_CHECKS:
        raise ValueError(
            f"orientation_check must be one of {', '.join(ORIENTATION_CHECKS)}, "
            f"not {orientation_check!r}."
        )
    patterns = validate_pams([pam])
    reference = check_plain_fasta(reference, "reference")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)
    (outdir / METADATA_DIR).mkdir(exist_ok=True)
    record = outdir / METADATA_DIR / f"{STEP}_{_stamp()}.json"

    parameters = {
        "pam": patterns[0],
        "spacer_length": spacer_length,
        "add_leading_g": add_leading_g,
        "orientation_check": orientation_check,
        "fix_orientation": fix_orientation,
        "contigs": None if contigs is None else list(contigs),
        "id_col": id_col,
        "spacer_col": spacer_col,
        "sep": sep,
        "header": header,
    }
    state: dict = {
        "step": STEP,
        "complete": False,
        "error": None,
        "started": _now(),
        "finished": None,
        "cau_version": __version__,
        "inputs": {},
        "parameters": parameters,
        "counts": {},
        "classes": {},
        "orientation": {},
        "orientation_fixed": 0,
        "seconds": {},
    }
    started = time.perf_counter()

    with _log_to(outdir / "logs" / LOG_FILE):
        logger.info("Guide alignment into %s, with %s", outdir, parameters)
        try:
            result = _run(guides, reference, outdir, parameters, state, record, started)
        except Exception as error:
            state["error"] = f"{type(error).__name__}: {error}"
            state["finished"] = _now()
            state["seconds"]["total"] = round(time.perf_counter() - started, 3)
            write_json(record, state)
            logger.error("Guide alignment failed: %s", state["error"])
            raise
    return result


def _run(
    guides: str | os.PathLike[str] | pd.DataFrame,
    reference: Path,
    outdir: Path,
    parameters: dict,
    state: dict,
    record: Path,
    started: float,
) -> RunResult:
    """The run itself, inside the log handler and the metadata record."""
    with _timed(state, "guides"):
        source = guides if isinstance(guides, pd.DataFrame) else os.fspath(guides)
        library = read_guides(
            guides,
            id_col=parameters["id_col"],
            spacer_col=parameters["spacer_col"],
            sep=parameters["sep"],
            header=parameters["header"],
            spacer_length=parameters["spacer_length"],
            add_leading_g=parameters["add_leading_g"],
            pam=parameters["pam"],
        )
    state["inputs"] = {
        "guides": "DataFrame" if isinstance(guides, pd.DataFrame) else source,
        "reference": str(reference),
    }
    logger.info(
        "%d entries, %d distinct guides, %d with a base 5' of the spacer",
        len(library),
        len({guide.guide_id for guide in library}),
        sum(guide.g_added for guide in library),
    )

    with _timed(state, "checksums"):
        state["inputs"]["reference_md5"] = file_md5(reference)
    write_json(record, state)

    with _timed(state, "scan"):
        scan = scan_guides(
            library,
            reference,
            pam=parameters["pam"],
            contigs=parameters["contigs"],
        )
    state["counts"]["scan"] = scan.counts()
    write_json(record, state)

    library, scan = _check_orientation(library, scan, reference, parameters, state)

    with _timed(state, "outputs"):
        summaries = summarize(
            library,
            scan.sites,
            scan.verdicts,
            evidence=scan.evidence,
            pam=parameters["pam"],
            max_mismatches=0,
        )
        write_guides_tsv(outdir / "guides_metadata.tsv", summaries, max_mismatches=0)
        write_sites_tsv(outdir / "sites.tsv", library, scan.sites)
        write_orientation_tsv(outdir / "orientation.tsv", library, scan)
        write_guides_bed(outdir / "guides.bed", summaries)
        write_cut_sites_bed(outdir / "cut_sites.bed", summaries)
        write_qc_tsv(outdir / "qc.tsv", summaries, counts=state["counts"])

    state["classes"] = _ordered((summary.site_class for summary in summaries), CLASSES)
    state["orientation"] = _ordered(
        (summary.orientation for summary in summaries), VERDICTS
    )
    state["complete"] = True
    state["finished"] = _now()
    state["seconds"]["total"] = round(time.perf_counter() - started, 3)
    write_json(record, state)
    logger.info(
        "Done in %.1f s: classes %s; outputs in %s",
        state["seconds"]["total"],
        state["classes"],
        outdir,
    )
    return RunResult(
        outdir=outdir,
        guides=list(library),
        summaries=summaries,
        scan=scan,
        metadata=state,
    )


def _check_orientation(
    library: list[Guide],
    scan: ScanResult,
    reference: Path,
    parameters: dict,
    state: dict,
) -> tuple[list[Guide], ScanResult]:
    """Act on the orientation verdicts, flipping entries when asked to."""
    check = parameters["orientation_check"]
    if check == "off":
        logger.info("Orientation check: off")
        return library, scan

    counts = Counter(scan.verdicts.values())
    backwards = [
        row for row, verdict in scan.verdicts.items() if verdict == REVERSE_COMPLEMENT
    ]
    if not backwards:
        return library, scan

    aliases = ", ".join(repr(library[row].alias) for row in backwards[:5])
    if len(backwards) > 5:
        aliases += ", ..."
    message = (
        f"{len(backwards)} of {len(scan.verdicts)} entries read 3' to 5': their "
        f"reverse complement has a site with a PAM and the sequence as given "
        f"does not ({aliases})."
    )
    # One guide can look reversed because its PAM differs between assemblies,
    # so only a table that is reversed throughout stops the run.
    if counts[REVERSE_COMPLEMENT] > counts[FORWARD] and check == "error":
        raise OrientationError(
            f"{message} Pass --fix-orientation to align their reverse complement "
            "instead, or --orientation-check warn to carry on."
        )
    if not parameters["fix_orientation"]:
        logger.warning(message)
        return library, scan

    logger.warning("%s Reverse-complementing those entries.", message)
    with _timed(state, "rescan"):
        library = _reverse_complemented(library, backwards, parameters)
        scan = scan_guides(
            library,
            reference,
            pam=parameters["pam"],
            contigs=parameters["contigs"],
        )
    state["orientation_fixed"] = len(backwards)
    state["counts"]["scan"] = scan.counts()
    return library, scan


def _reverse_complemented(
    library: Sequence[Guide], rows: Sequence[int], parameters: Mapping[str, object]
) -> list[Guide]:
    """The library with only `rows` read again from their reverse complement."""
    flip = set(rows)
    table = pd.DataFrame(
        {
            "alias": [guide.alias for guide in library],
            "sequence": [
                reverse_complement(guide.input_sequence)
                if row in flip
                else guide.input_sequence
                for row, guide in enumerate(library)
            ],
        }
    )
    return read_guides(
        table,
        id_col="alias",
        spacer_col="sequence",
        spacer_length=parameters["spacer_length"],
        add_leading_g=parameters["add_leading_g"],
        pam=parameters["pam"],
    )


def _ordered(values: Iterable[str], order: Sequence[str]) -> dict[str, int]:
    counts = Counter(values)
    return {name: counts[name] for name in order}


def _stamp() -> str:
    """The UTC timestamp that names this invocation's metadata file."""
    return time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


@contextmanager
def _timed(state: dict, name: str) -> Iterator[None]:
    """Record how long a step took, in seconds."""
    start = time.perf_counter()
    try:
        yield
    finally:
        state["seconds"][name] = round(time.perf_counter() - start, 3)


@contextmanager
def _log_to(path: Path) -> Iterator[None]:
    """Append this package's log records, from INFO up, to `path`."""
    package = logging.getLogger(__package__)
    handler = logging.FileHandler(path, mode="a", encoding="utf-8")
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(_LOG_FORMAT))
    level = package.level
    if package.getEffectiveLevel() > logging.INFO:
        package.setLevel(logging.INFO)
    package.addHandler(handler)
    try:
        yield
    finally:
        package.removeHandler(handler)
        handler.close()
        package.setLevel(level)

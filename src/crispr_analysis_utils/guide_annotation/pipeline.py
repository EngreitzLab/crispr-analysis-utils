"""The gene annotation: read the aligned sites, place them in genes, write the table.

`outdir` gets three folders and nothing else: ``annotations/`` with the one
results table and its ``guides_gene_annotation_columns.yaml``, ``run_metadata/`` with this
invocation's record, and ``logs/``.
"""

from __future__ import annotations

import logging
import os
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .._files import file_md5, write_json
from .._version import __version__
from .annotate import AnnotationError, annotate_sites
from .db import Store, load_store
from .promoters import promoters_from_store, read_promoters
from .summary import TABLE, counts, write_annotation_tsv, write_columns_yaml

logger = logging.getLogger(__name__)

RESULTS_DIR = "annotations"
"""Where the step's results go, so nothing of its own sits beside the inputs."""
METADATA_DIR = "run_metadata"
"""Where each invocation leaves its metadata, one file per invocation."""
STEP = "guide_annotation"
"""The step those metadata files are named for."""
LOG_FILE = "guide_annotation.log"
_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"

PROMOTER_WINDOW = 1000
"""How far either side of a TSS counts as its promoter, by default."""


@dataclass
class RunResult:
    """What `run` produced.

    Attributes
    ----------
    outdir
        The output folder.
    annotated
        The annotation table: one row per site per overlapping gene.
    store
        The annotation store the run used.
    metadata
        The record written under ``run_metadata/``.
    """

    outdir: Path
    annotated: pd.DataFrame
    store: Store
    metadata: dict = field(default_factory=dict)


def run(
    alignment: str | os.PathLike[str] | pd.DataFrame,
    store: str | os.PathLike[str] | Store,
    outdir: str | os.PathLike[str],
    *,
    promoters: str | os.PathLike[str] | None = None,
    promoter_window: int = PROMOTER_WINDOW,
) -> RunResult:
    """Annotate every aligned site against the genes of the store.

    Parameters
    ----------
    alignment
        ``guides_alignment_metadata.tsv`` from `guide_alignment.run`, or the
        same table as a DataFrame. Only the sites with a PAM are annotated,
        which is what it means for a guide to map.
    store
        The annotation store, as `build_store` wrote it.
    outdir
        Where the table goes. Created if it does not exist; nothing is ever
        written outside it.
    promoters
        A TSS BED. Without one, each gene's reference transcript start is
        used instead.
    promoter_window
        How far either side of a TSS counts as its promoter.

    Returns
    -------
    RunResult
        The table, the store and the record.

    Raises
    ------
    AnnotationError
        The sites and the store name contigs differently, or the alignment
        table holds no mapped site.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)
    (outdir / METADATA_DIR).mkdir(exist_ok=True)
    (outdir / RESULTS_DIR).mkdir(exist_ok=True)
    record = outdir / METADATA_DIR / f"{STEP}_{_stamp()}.json"

    parameters = {
        "promoters": None if promoters is None else str(promoters),
        "promoter_window": promoter_window,
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
        "store": {},
        "counts": {},
        "seconds": {},
    }
    started = time.perf_counter()

    with _log_to(outdir / "logs" / LOG_FILE):
        logger.info("Gene annotation into %s, with %s", outdir, parameters)
        try:
            return _run(
                alignment, store, outdir, promoters, parameters, state, record, started
            )
        except Exception as error:
            state["error"] = f"{type(error).__name__}: {error}"
            state["finished"] = _now()
            state["seconds"]["total"] = round(time.perf_counter() - started, 3)
            write_json(record, state)
            logger.error("Gene annotation failed: %s", state["error"])
            raise


def _run(
    alignment: str | os.PathLike[str] | pd.DataFrame,
    store: str | os.PathLike[str] | Store,
    outdir: Path,
    promoters: str | os.PathLike[str] | None,
    parameters: dict,
    state: dict,
    record: Path,
    started: float,
) -> RunResult:
    """The run itself, inside the log handler and the metadata record."""
    with _timed(state, "store"):
        loaded = store if isinstance(store, Store) else load_store(store)
    state["store"] = loaded.provenance
    state["inputs"]["store"] = "Store" if isinstance(store, Store) else str(store)
    logger.info(
        "Store: %d genes, %d contigs",
        loaded.provenance.get("genes", 0),
        len(loaded.contigs),
    )

    with _timed(state, "sites"):
        sites = _read_sites(alignment)
    state["inputs"]["alignment"] = (
        "DataFrame" if isinstance(alignment, pd.DataFrame) else str(alignment)
    )
    if not isinstance(alignment, pd.DataFrame):
        state["inputs"]["alignment_md5"] = file_md5(alignment)
    logger.info("%d mapped sites to annotate", len(sites))

    with _timed(state, "promoters"):
        if promoters is None:
            tss = promoters_from_store(loaded)
            logger.info("No promoter file: %d TSS from the store", len(tss))
        else:
            tss = read_promoters(promoters)
            state["inputs"]["promoters_md5"] = file_md5(promoters)
            logger.info("%d TSS from %s", len(tss), promoters)
    write_json(record, state)

    with _timed(state, "annotate"):
        annotated = annotate_sites(
            sites,
            loaded,
            tss,
            promoter_window=parameters["promoter_window"],
        )

    with _timed(state, "outputs"):
        tables = outdir / RESULTS_DIR
        write_annotation_tsv(tables / TABLE, annotated)
        write_columns_yaml(tables / "guides_gene_annotation_columns.yaml")

    state["counts"] = counts(annotated)
    state["complete"] = True
    state["finished"] = _now()
    state["seconds"]["total"] = round(time.perf_counter() - started, 3)
    write_json(record, state)
    logger.info("Gene annotation complete: %s", state["counts"])
    return RunResult(outdir=outdir, annotated=annotated, store=loaded, metadata=state)


def _read_sites(alignment: str | os.PathLike[str] | pd.DataFrame) -> pd.DataFrame:
    """The alignment table's mapped sites: the rows whose genomic PAM matched."""
    if isinstance(alignment, pd.DataFrame):
        frame = alignment.copy()
    else:
        path = Path(alignment)
        if not path.is_file():
            raise FileNotFoundError(f"The alignment table {path} does not exist.")
        frame = pd.read_csv(path, sep="\t", dtype={"chr": str}, low_memory=False)

    if "pam_matches" in frame:
        matched = frame["pam_matches"]
        frame = frame[matched.map(_truthy).fillna(False)]
    frame = frame[frame["chr"].notna()] if "chr" in frame else frame
    if frame.empty:
        raise AnnotationError(
            "The alignment table holds no mapped site, so there is nothing to"
            " annotate. Every guide matched nowhere, or had no PAM."
        )
    for column in ("start", "end", "cut_site"):
        frame[column] = frame[column].astype("int64")
    return frame.reset_index(drop=True)


def _truthy(value: object) -> bool:
    """TRUE, true and True all mean the same in a written table.

    An empty cell is not a PAM. It arrives as NaN, which is truthy, so it has
    to be rejected before `bool`.
    """
    if value is None or (isinstance(value, float) and value != value):
        return False
    if isinstance(value, str):
        return value.strip().upper() == "TRUE"
    return bool(value)


def _stamp() -> str:
    """A UTC timestamp for a file name."""
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def _now() -> str:
    """The time now, in UTC, as ISO 8601."""
    return datetime.now(UTC).isoformat(timespec="seconds")


@contextmanager
def _timed(state: dict, name: str) -> Iterator[None]:
    """Record how long a step took."""
    started = time.perf_counter()
    try:
        yield
    finally:
        state["seconds"][name] = round(time.perf_counter() - started, 3)


@contextmanager
def _log_to(path: Path) -> Iterator[None]:
    """Send this package's log to a file for the length of a run."""
    handler = logging.FileHandler(path, encoding="utf-8")
    handler.setFormatter(logging.Formatter(_LOG_FORMAT))
    package = logging.getLogger(__name__.rsplit(".", 2)[0])
    package.addHandler(handler)
    previous = package.level
    if previous == logging.NOTSET:
        package.setLevel(logging.INFO)
    # The run's log belongs in the file, not on the terminal. pyranges calls
    # logging.basicConfig on import, which puts a handler on the root logger.
    propagate = package.propagate
    package.propagate = False
    try:
        yield
    finally:
        package.removeHandler(handler)
        package.setLevel(previous)
        package.propagate = propagate
        handler.close()

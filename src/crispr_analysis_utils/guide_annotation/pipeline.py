"""The gene annotation: read the aligned sites, place them in genes, write the table.

`outdir` gets three folders and nothing else: ``annotations/`` with the one
results table and its ``guide_gene_annotation_columns.yaml``, ``run_metadata/`` with this
invocation's record, and ``logs/``.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .._files import file_md5, run_record, timed, truthy, write_json
from .annotate import MIDPOINT, AnnotationError, annotate_sites
from .db import Store, load_store
from .promoters import promoters_from_store, read_promoters
from .summary import TABLE, counts, write_annotation_tsv, write_columns_yaml

logger = logging.getLogger(__name__)

RESULTS_DIR = "annotations"
"""Where the step's results go, so nothing of its own sits beside the inputs."""
STEP = "guide_annotation"
"""The step its run records and log are named for."""

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
    distance_from: str = MIDPOINT,
) -> RunResult:
    """Annotate every aligned site against the genes of the store.

    Parameters
    ----------
    alignment
        ``guide_genome_annotation.tsv`` from `guide_alignment.run`, or the
        same table as a DataFrame. Only the sites with a PAM are annotated,
        which is what it means for a guide to map.
    store
        The annotation store, as `build_store` wrote it.
    outdir
        Where the table goes. Created if it does not exist; nothing is ever
        written outside it.
    promoters
        A TSS BED. Without one, a site per transcript is taken from the store.
    promoter_window
        How far either side of a TSS counts as its promoter.
    distance_from
        Where distances are measured from: the protospacer midpoint, or the
        target site 3 bases 5' of the PAM.

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
    parameters = {
        "promoters": None if promoters is None else str(promoters),
        "promoter_window": promoter_window,
        "distance_from": distance_from,
    }
    with run_record(
        outdir, STEP, RESULTS_DIR, parameters, logger, "Gene annotation"
    ) as (state, record):
        state["store"] = {}
        logger.info("Gene annotation into %s, with %s", outdir, parameters)
        return _run(alignment, store, outdir, promoters, parameters, state, record)


def _run(
    alignment: str | os.PathLike[str] | pd.DataFrame,
    store: str | os.PathLike[str] | Store,
    outdir: Path,
    promoters: str | os.PathLike[str] | None,
    parameters: dict,
    state: dict,
    record: Path,
) -> RunResult:
    """The run itself, inside the log handler and the metadata record."""
    with timed(state, "store"):
        loaded = store if isinstance(store, Store) else load_store(store)
    state["store"] = loaded.provenance
    state["inputs"]["store"] = "Store" if isinstance(store, Store) else str(store)
    logger.info(
        "Store: %d genes, %d contigs",
        loaded.provenance.get("genes", 0),
        len(loaded.contigs),
    )

    with timed(state, "sites"):
        sites = _read_sites(alignment)
    state["inputs"]["alignment"] = (
        "DataFrame" if isinstance(alignment, pd.DataFrame) else str(alignment)
    )
    if not isinstance(alignment, pd.DataFrame):
        state["inputs"]["alignment_md5"] = file_md5(alignment)
    logger.info("%d mapped sites to annotate", len(sites))

    with timed(state, "promoters"):
        if promoters is None:
            tss = promoters_from_store(loaded)
            logger.info("No promoter file: %d TSS from the store", len(tss))
        else:
            tss = read_promoters(promoters)
            state["inputs"]["promoters_md5"] = file_md5(promoters)
            logger.info("%d TSS from %s", len(tss), promoters)
    write_json(record, state)

    with timed(state, "annotate"):
        annotated = annotate_sites(
            sites,
            loaded,
            tss,
            promoter_window=parameters["promoter_window"],
            distance_from=parameters["distance_from"],
        )

    with timed(state, "outputs"):
        tables = outdir / RESULTS_DIR
        write_annotation_tsv(tables / TABLE, annotated)
        write_columns_yaml(tables / "guide_gene_annotation_columns.yaml")

    state["counts"] = counts(annotated)
    logger.info("Gene annotation complete: %s", state["counts"])
    return RunResult(outdir=outdir, annotated=annotated, store=loaded, metadata=state)


def _read_sites(alignment: str | os.PathLike[str] | pd.DataFrame) -> pd.DataFrame:
    """The alignment table's mapped sites: those the guide's class counts.

    Those are the rows `best_site` marks, or, in a table written before that
    column, every row whose genomic PAM matched.
    """
    if isinstance(alignment, pd.DataFrame):
        frame = alignment.copy()
    else:
        path = Path(alignment)
        if not path.is_file():
            raise FileNotFoundError(f"The alignment table {path} does not exist.")
        frame = pd.read_csv(path, sep="\t", dtype={"chr": str}, low_memory=False)

    flag = "best_site" if "best_site" in frame else "pam_matches"
    if flag in frame:
        frame = frame[frame[flag].map(truthy)]
    frame = frame[frame["chr"].notna()] if "chr" in frame else frame
    if frame.empty:
        raise AnnotationError(
            "The alignment table holds no mapped site, so there is nothing to"
            " annotate. Every guide matched nowhere, or had no PAM."
        )
    for column in ("start", "end", "target_site"):
        frame[column] = frame[column].astype("int64")
    return frame.reset_index(drop=True)

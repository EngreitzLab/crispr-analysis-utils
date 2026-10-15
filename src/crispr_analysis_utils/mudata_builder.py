"""Assemble a MuData from a gene and a guide count matrix.

Assembly only. Nothing is derived from the counts, no cell is dropped, no
guide is called and no pair is proposed: the result says what was measured
and what each feature is, and leaves every judgement to the steps after it.

The layout follows the IGVF CRISPR pipeline's own object, so the two can be
compared field by field: modalities ``gene`` and ``guide``, a cell index of
``{barcode}_{measurement_set}``, raw counts in ``X``, and the annotation of
each feature in ``var``.

Cells are the union of the two matrices, not their intersection. A barcode
seen in one modality and not the other is kept, and MuData's own ``obsm``
masks record which modality holds it. Filtering is a later, visible decision.
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from ._files import file_md5, write_json
from ._version import __version__
from .guide_annotation.db import _unversioned

GENE_COLUMNS = ("symbol", "gene_type", "gene_chr", "gene_start", "gene_end")
"""What `var` holds for a gene, beside the id that indexes it."""

MODALITIES = ("gene", "guide")
"""The two modalities, in the order the object carries them."""

TABLE = "mudata.h5mu"
"""The one output."""

RESULTS_DIR = "mudata"
METADATA_DIR = "run_metadata"
STEP = "mudata_builder"
LOG_FILE = "mudata_builder.log"
_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"

logger = logging.getLogger(__name__)


def _now() -> str:
    """The moment, as the run record writes it."""
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _stamp() -> str:
    """The run record's file stamp."""
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


class BuildError(RuntimeError):
    """The inputs cannot be assembled into one object."""


def read_counts(path: str | os.PathLike[str]):
    """Read one count matrix.

    Parameters
    ----------
    path
        An ``.h5ad`` written by kb, or anything anndata can read.

    Returns
    -------
    anndata.AnnData
        The matrix as it was stored, layers and all.

    Raises
    ------
    BuildError
        The file is missing or unreadable.
    """
    import anndata as ad

    path = Path(path)
    if not path.is_file():
        raise BuildError(f"The count matrix {path} does not exist.")
    try:
        return ad.read_h5ad(path)
    except (OSError, KeyError, ValueError) as error:
        raise BuildError(f"{path} is not a readable h5ad: {error}") from error


def gene_var(store, ids: pd.Index) -> pd.DataFrame:
    """The annotation of every gene in the matrix, in the matrix's own order.

    Parameters
    ----------
    store
        A loaded annotation `Store`. Every gene type it holds is used, not
        only the ones guides are annotated against.
    ids
        The matrix's gene ids, versioned or not.

    Returns
    -------
    pandas.DataFrame
        Indexed by the ids as given, so it can be attached to ``var``
        directly. A gene the store does not hold keeps its row with the
        annotation empty, rather than being dropped.
    """
    # Every gene the store holds, not only the annotated types: the counts
    # come from a transcriptome built on the whole GTF, so a lncRNA feature
    # is a real row of the matrix and deserves its name.
    genes = store.rows("gene", primary=False).drop_duplicates("ensg").set_index("ensg")
    annotation = pd.DataFrame(
        {
            "symbol": genes["gene_symbol"],
            "gene_type": genes["gene_type"].astype("string"),
            "gene_chr": genes["Chromosome"].astype("string"),
            "gene_start": genes["Start"],
            "gene_end": genes["End"],
        }
    )
    # kb indexes genes as the GTF names them, with the version; the store
    # strips it. The rule has to be the store's own, or a PAR_Y gene, whose
    # id ends in neither, comes back annotated with nothing.
    out = annotation.reindex(_unversioned(pd.Series(ids, dtype="string")))
    out.index = pd.Index(ids, name="gene_id")
    return out


def guide_var(path: str | os.PathLike[str], ids: pd.Index) -> pd.DataFrame:
    """The guide design, in the matrix's own order.

    Parameters
    ----------
    path
        ``guide_metadata.tsv`` from `cau igvf-export`, one row per guide.
    ids
        The matrix's guide ids.

    Returns
    -------
    pandas.DataFrame
        Every column of the file except the id that indexes it. A guide the
        file does not describe keeps its row with the annotation empty.

    Raises
    ------
    BuildError
        The file is missing, or has no ``guide_id`` column.
    """
    path = Path(path)
    if not path.is_file():
        raise BuildError(f"The guide metadata {path} does not exist.")
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if "guide_id" not in frame.columns:
        raise BuildError(
            f"{path} has no guide_id column, so it is not a guide metadata"
            " table. `cau igvf-export run` writes one."
        )
    frame = frame.drop_duplicates("guide_id").set_index("guide_id")
    out = frame.reindex(pd.Index(ids, dtype="string"))
    out.index = pd.Index(ids, name="guide_id")
    return out


def reference_block(
    store,
    *,
    gtf_accession: str | None = None,
    genome: str | os.PathLike[str] | None = None,
    genome_accession: str | None = None,
) -> dict:
    """What the object says about the annotation it was built against.

    Carried in ``uns['reference']`` so the assembled file answers "which
    GENCODE, which genome" without the run record beside it.

    Parameters
    ----------
    store
        A loaded `Store`; its provenance names the GTF it was built from.
    gtf_accession
        The GTF's IGVF portal id, which no local file carries. The CRISPR
        pipeline's inputs have it; without it the field falls back to the
        file's name, so it is never empty and never a guess.
    genome
        The reference FASTA, for its checksum. Nothing is read from it.
    genome_accession
        The genome's IGVF portal id, falling back to the file's name.

    Returns
    -------
    dict
        ``gtf``, ``gtf_md5``, ``gtf_accession``, ``genome``, ``genome_md5``,
        each file named rather than pathed, since the object is shared and a
        path is not,
        ``genome_accession``, ``gene_types`` and ``cau_version``. A checksum
        is empty when the file is no longer where the store recorded it: a
        path is a claim, a checksum is evidence. An accession falls back to
        the file's name, so the field always says which file was meant.
    """
    gtf = str(store.provenance.get("gtf", ""))
    # The name, never the path: this travels with the object, and a path
    # carries the filesystem it was built on.
    block = {
        "gtf": Path(gtf).name if gtf else "",
        "gtf_md5": file_md5(gtf) if gtf and Path(gtf).is_file() else "",
        "gtf_accession": gtf_accession or (Path(gtf).name if gtf else ""),
        "genome": Path(genome).name if genome else "",
        "genome_md5": (file_md5(genome) if genome and Path(genome).is_file() else ""),
        "genome_accession": (genome_accession or (Path(genome).name if genome else "")),
        "gene_types": list(store.provenance.get("gene_types", [])),
        "cau_version": str(store.provenance.get("cau_version", "")),
    }
    if genome and not Path(genome).is_file():
        raise BuildError(f"The genome {genome} does not exist.")
    return block


def check_same_store(
    guide_metadata: str | os.PathLike[str], store_path: str | os.PathLike[str]
) -> str:
    """Refuse a store other than the one the guide annotation was built on.

    The guide table carries no provenance of its own, but the run record
    beside it does. Annotating guides against one GENCODE and genes against
    another produces an object whose two halves disagree, and nothing
    downstream would notice.

    Parameters
    ----------
    guide_metadata
        ``guide_metadata.tsv``; its ``../run_metadata`` is where the record
        of the run that wrote it lives.
    store_path
        The store this build was given.

    Returns
    -------
    str
        What was compared, for the log. Empty when no record was found, which
        is not an error: the table may have been copied away from its run.

    Raises
    ------
    BuildError
        The record names a different store.
    """
    records = Path(guide_metadata).resolve().parent.parent / "run_metadata"
    if not records.is_dir():
        return ""
    written = sorted(records.glob("igvf_export_*.json"))
    if not written:
        return ""
    state = json.loads(written[-1].read_text(encoding="utf-8"))
    used = state.get("inputs", {}).get("store", "")
    if not used:
        return ""
    if Path(used).resolve() != Path(store_path).resolve():
        raise BuildError(
            f"The guide metadata was written against {used}, but this build "
            f"was given {store_path}. Genes and guides would be annotated "
            "against different GENCODE releases. Pass the same store, or "
            "re-run `cau igvf-export run` against this one."
        )
    return used


def build(
    rna,
    guide,
    *,
    gene_table: pd.DataFrame,
    guide_table: pd.DataFrame,
    measurement_set: str,
    reference: dict | None = None,
    covariates: pd.DataFrame | None = None,
    tables: dict[str, pd.DataFrame] | None = None,
):
    """Assemble the two matrices into one MuData.

    Parameters
    ----------
    rna, guide
        The count matrices, as `read_counts` returns them.
    gene_table, guide_table
        What `gene_var` and `guide_var` return, already in each matrix's
        order.
    measurement_set
        The channel this pair of matrices is, which the cell index carries so
        two channels can be concatenated without their barcodes colliding.
    reference
        What `reference_block` returns, stored in ``uns['reference']``.
    covariates
        One row per cell, indexed by barcode. Its columns join the root
        ``obs``. Nothing is computed; whatever is given is carried.
    tables
        Whole tables to carry in ``mod['guide'].uns``, by name. The gene
        annotation is one row per guide per site per gene, so it cannot go in
        ``var`` without collapsing it; keeping it whole loses nothing and
        decides nothing.

    Returns
    -------
    mudata.MuData
        Modalities ``gene`` and ``guide``, cells the union of the two.

    Raises
    ------
    BuildError
        A matrix and its annotation disagree on shape.
    """
    import mudata as md

    for name, matrix, table in (
        ("gene", rna, gene_table),
        ("guide", guide, guide_table),
    ):
        if len(table) != matrix.n_vars:
            raise BuildError(
                f"The {name} annotation has {len(table)} rows for "
                f"{matrix.n_vars} features. They must line up."
            )

    rna = rna.copy()
    guide = guide.copy()
    rna.var = gene_table
    guide.var = guide_table
    for matrix in (rna, guide):
        # {barcode}_{measurement_set}, as the pipeline indexes cells, so a
        # concatenation across channels cannot collide.
        barcodes = pd.Index(matrix.obs_names, dtype="string")
        matrix.obs["barcode"] = barcodes
        matrix.obs["measurement_set"] = measurement_set
        matrix.obs_names = pd.Index(barcodes + "_" + measurement_set)
        # kb names the index "barcode", and anndata refuses an index whose
        # name is also a column holding different values.
        matrix.obs.index.name = None
        matrix.var.index.name = None

    # MuData takes the union of the two indexes and records membership in
    # obsm, so a barcode seen in one modality only is kept and visible.
    for name, table in (tables or {}).items():
        guide.uns[name] = table

    data = md.MuData({"gene": rna, "guide": guide})
    data.uns["reference"] = dict(reference or {})
    data.uns["measurement_set"] = measurement_set

    if covariates is not None and not covariates.empty:
        joined = covariates.reindex(
            pd.Index(data.obs_names, dtype="string").str.rsplit("_", n=1).str[0]
        )
        joined.index = data.obs_names
        for column in joined.columns:
            data.obs[column] = joined[column]
    return data


@dataclass
class RunResult:
    """What `run` produced.

    Attributes
    ----------
    outdir
        The output folder.
    data
        The assembled MuData.
    metadata
        The record written under ``run_metadata/``.
    """

    outdir: Path
    data: object
    metadata: dict = field(default_factory=dict)


def run(
    rna: str | os.PathLike[str],
    guide: str | os.PathLike[str],
    guide_metadata: str | os.PathLike[str],
    store: str | os.PathLike[str],
    outdir: str | os.PathLike[str],
    *,
    measurement_set: str,
    genome: str | os.PathLike[str] | None = None,
    gtf_accession: str | None = None,
    genome_accession: str | None = None,
    covariates: str | os.PathLike[str] | None = None,
    guide_annotation: str | os.PathLike[str] | None = None,
    element_mapping: str | os.PathLike[str] | None = None,
) -> RunResult:
    """Assemble one channel's matrices into a MuData, and write it.

    Parameters
    ----------
    rna, guide
        The two count matrices, as ``.h5ad``.
    guide_metadata
        ``guide_metadata.tsv`` from `cau igvf-export`.
    store
        The annotation store, which must be the one the guide metadata was
        built against.
    outdir
        Where the object goes. Created if it does not exist; nothing is ever
        written outside it.
    measurement_set
        The channel these matrices are.
    genome, gtf_accession, genome_accession
        Recorded in ``uns['reference']``; nothing is read from the genome
        but its checksum.
    covariates
        A table of per-cell covariates, indexed by barcode in its first
        column. Carried into the root ``obs`` unchanged.
    guide_annotation, element_mapping
        ``guide_gene_annotation.tsv`` and ``guide_element_mapping.tsv``.
        Carried whole in ``mod['guide'].uns``, since neither is one row per
        guide. ``var`` keeps the IGVF columns alone.

    Returns
    -------
    RunResult
        The object, and the record of the run.
    """
    from .guide_annotation import load_store

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)
    (outdir / METADATA_DIR).mkdir(exist_ok=True)
    (outdir / RESULTS_DIR).mkdir(exist_ok=True)
    record = outdir / METADATA_DIR / f"{STEP}_{_stamp()}.json"

    state: dict = {
        "step": STEP,
        "complete": False,
        "error": None,
        "started": _now(),
        "finished": None,
        "cau_version": __version__,
        "inputs": {},
        "parameters": {"measurement_set": measurement_set},
        "counts": {},
    }
    handler = logging.FileHandler(outdir / "logs" / LOG_FILE, encoding="utf-8")
    handler.setFormatter(logging.Formatter(_LOG_FORMAT))
    package = logging.getLogger(__name__.rsplit(".", 1)[0])
    package.addHandler(handler)
    if package.level == logging.NOTSET:
        package.setLevel(logging.INFO)
    propagate = package.propagate
    package.propagate = False
    try:
        compared = check_same_store(guide_metadata, store)
        if compared:
            logger.info("Guide metadata and genes share the store %s", compared)
        else:
            logger.warning(
                "No run record beside the guide metadata, so the store it was "
                "built against could not be checked."
            )

        loaded = load_store(store)
        rna_counts = read_counts(rna)
        guide_counts = read_counts(guide)
        for name, path in (
            ("rna", rna),
            ("guide", guide),
            ("guide_metadata", guide_metadata),
            ("store", store),
        ):
            state["inputs"][name] = str(path)
            state["inputs"][f"{name}_md5"] = file_md5(path)

        covariate_table = None
        if covariates is not None:
            covariate_table = pd.read_csv(covariates, sep="\t", index_col=0)
            state["inputs"]["covariates"] = str(covariates)
            state["inputs"]["covariates_md5"] = file_md5(covariates)

        tables = {}
        for name, path in (
            ("annotation", guide_annotation),
            ("element_mapping", element_mapping),
        ):
            if path is None:
                continue
            tables[name] = pd.read_csv(path, sep="\t", low_memory=False)
            state["inputs"][name] = str(path)
            state["inputs"][f"{name}_md5"] = file_md5(path)

        reference = reference_block(
            loaded,
            gtf_accession=gtf_accession,
            genome=genome,
            genome_accession=genome_accession,
        )
        state["parameters"]["reference"] = reference

        data = build(
            rna_counts,
            guide_counts,
            gene_table=gene_var(loaded, rna_counts.var_names),
            guide_table=guide_var(guide_metadata, guide_counts.var_names),
            measurement_set=measurement_set,
            reference=reference,
            covariates=covariate_table,
            tables=tables,
        )
        data.write(outdir / RESULTS_DIR / TABLE)

        # mudata may hand the membership mask back as a column vector, so
        # flatten before counting rather than trusting its shape.
        in_gene = np.asarray(data.obsm["gene"]).ravel()
        in_guide = np.asarray(data.obsm["guide"]).ravel()
        state["counts"] = {
            "cells": int(data.n_obs),
            "cells_gene_only": int((in_gene & ~in_guide).sum()),
            "cells_guide_only": int((in_guide & ~in_gene).sum()),
            "cells_both": int((in_gene & in_guide).sum()),
            "genes": int(data["gene"].n_vars),
            "guides": int(data["guide"].n_vars),
            "genes_unannotated": int(data["gene"].var["symbol"].isna().sum()),
            "guides_unannotated": int(data["guide"].var["spacer"].isna().sum()),
            "uns_tables": {k: int(len(v)) for k, v in tables.items()},
        }
        state["complete"] = True
        state["finished"] = _now()
        write_json(record, state)
        logger.info("MuData complete: %s", state["counts"])
        return RunResult(outdir=outdir, data=data, metadata=state)
    except Exception as error:
        state["error"] = f"{type(error).__name__}: {error}"
        state["finished"] = _now()
        write_json(record, state)
        logger.error("MuData build failed: %s", state["error"])
        raise
    finally:
        package.removeHandler(handler)
        package.propagate = propagate
        handler.close()

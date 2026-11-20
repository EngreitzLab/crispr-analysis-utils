"""One MuData an inference reads, laid out like the IGVF CRISPR pipeline's.

Takes the QC'd per-channel objects of `mudata_builder` and `qc`, each already
carrying guide calls in ``mod['guide'].layers['guide_assignment']``, and the
pairs of `pairs`, and writes the object the pipeline hands to SCEPTRE:

- one cell index across both modalities, the ``{barcode}_{measurement_set}``
  ids `mudata_builder` gave, refused if two channels share one;
- unversioned gene ids, so they match the pairs;
- ``intended_target_key`` on every guide, and the controls bucketed as
  ``non-targeting|N`` by the pipeline's own rule;
- the pipeline's covariate names in ``obs``;
- ``moi`` and ``capture_method`` in ``mod['guide'].uns``;
- ``pairs_to_test`` in ``uns``, restricted to the genes and guides present.

Guide calls are not made here: how a guide is called is the analysis's
decision, and the layer records it.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ._files import file_md5, run_record, truthy

logger = logging.getLogger(__name__)

ASSIGNMENT = "guide_assignment"
"""The layer of binary guide calls, as the pipeline names it."""

NON_TARGETING = "non-targeting"
CONTROL_TYPES = ("safe-targeting", "non-targeting", "negative control")
"""Guide types the pipeline treats as controls, as `intended_target_key_utils` does."""

MOI = ("high", "low")
TARGET_COLUMNS = (
    "intended_target_name",
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end",
)
PAIR_COLUMNS = (
    "guide_id",
    "gene_id",
    *TARGET_COLUMNS,
    "intended_target_key",
    "pair_type",
)
"""What ``uns['pairs_to_test']`` holds, in the pipeline's order; extra pair columns follow."""

COVARIATES = (
    "batch",
    "total_gene_umis",
    "total_guide_umis",
    "num_expressed_genes",
    "percent_mito",
)
"""Root ``obs``: the pipeline's design frame under its own names."""

RESULTS_DIR = "inference_input"
STEP = "inference_input"


class InputError(RuntimeError):
    """The channels cannot be made into one inference object."""


@dataclass
class RunResult:
    """What `run` produced."""

    outdir: Path
    path: Path
    report: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)


def target_keys(var: pd.DataFrame) -> pd.DataFrame:
    """Add ``intended_target_key``, and bucket the controls as the pipeline does.

    A targeting guide's key is ``coords::name::chr::start::end``, or
    ``name_only::name`` without coordinates. Control guides (``targeting``
    false, or a control ``type``) are sorted by guide id and cut into groups
    the size of the median number of guides per targeting element, named
    ``non-targeting|N``, with chromosome ``non-targeting`` and start and end
    ``N``. That is the rule of the pipeline's
    ``annotate_intended_target_groups``; an id-preserving ``element_id`` it
    also honours is absent from the IGVF table.

    Parameters
    ----------
    var
        ``mod['guide'].var``, indexed by guide id, with ``targeting``,
        ``type`` and the four ``intended_target_*`` columns.

    Returns
    -------
    pandas.DataFrame
        A copy, with ``targeting`` as bool, coordinates as ``Int64`` and the
        new key column.
    """
    out = var.copy()
    out["targeting"] = out["targeting"].map(truthy).astype(bool)
    for column in ("intended_target_name", "intended_target_chr"):
        out[column] = out[column].astype("string").replace("", pd.NA)
    for column in ("intended_target_start", "intended_target_end"):
        out[column] = pd.to_numeric(
            out[column].replace("", pd.NA), errors="coerce"
        ).astype("Int64")

    kind = out["type"].astype("string").fillna("").str.strip().str.lower()
    name = out["intended_target_name"].fillna("").str.strip().str.lower()
    control = (
        ~out["targeting"]
        | kind.isin(CONTROL_TYPES)
        | (name == NON_TARGETING)
        | name.str.startswith(f"{NON_TARGETING}|")
    )
    keys = _keys(out)
    sizes = keys[~control].value_counts()
    size = max(1, int(np.median(sizes.to_numpy()))) if len(sizes) else 1

    order = pd.Index(sorted(out.index[control].astype(str)))
    group = pd.Series(np.arange(len(order)) // size + 1, index=order, dtype="Int64")
    out.loc[order, "intended_target_name"] = [f"{NON_TARGETING}|{n}" for n in group]
    out.loc[order, "intended_target_chr"] = NON_TARGETING
    out.loc[order, "intended_target_start"] = group
    out.loc[order, "intended_target_end"] = group
    out["intended_target_key"] = _keys(out)
    return out


def _keys(frame: pd.DataFrame) -> pd.Series:
    """The pipeline's target key, row by row."""
    coords = (
        frame[["intended_target_chr", "intended_target_start", "intended_target_end"]]
        .notna()
        .all(axis=1)
    )
    name = frame["intended_target_name"].astype("string").fillna("<NA>")
    full = (
        "coords::"
        + name
        + "::"
        + frame["intended_target_chr"].astype("string").fillna("")
        + "::"
        + frame["intended_target_start"].astype("string").fillna("")
        + "::"
        + frame["intended_target_end"].astype("string").fillna("")
    )
    return full.where(coords, "name_only::" + name).astype(str)


def cell_covariates(gene, guide) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-cell totals under the pipeline's names.

    Returns
    -------
    tuple of pandas.DataFrame
        For ``mod['gene'].obs``: ``total_gene_umis``, ``num_expressed_genes``,
        ``percent_mito`` and ``pct_counts_ribo``, percentages from 0 to 100,
        mitochondrial by contig and ribosomal as `qc.gene_classes` decides.
        For ``mod['guide'].obs``: ``total_guide_umis`` and
        ``num_expressed_guides``.
    """
    from .qc import gene_classes

    classes = gene_classes(gene.var)
    x = gene.X.tocsr()
    x.eliminate_zeros()
    total = np.asarray(x.sum(axis=1)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        share = {
            name: 100 * np.asarray(x[:, mask].sum(axis=1)).ravel() / total
            for name, mask in classes.items()
        }
    genes = pd.DataFrame(
        {
            "total_gene_umis": total,
            "num_expressed_genes": np.diff(x.indptr),
            "percent_mito": np.nan_to_num(share["mito"]),
            "pct_counts_ribo": np.nan_to_num(share["ribo"]),
        },
        index=gene.obs_names,
    )
    g = guide.X.tocsr()
    g.eliminate_zeros()
    guides = pd.DataFrame(
        {
            "total_guide_umis": np.asarray(g.sum(axis=1)).ravel(),
            "num_expressed_guides": np.diff(g.indptr),
        },
        index=guide.obs_names,
    )
    return genes, guides


def pairs_frame(
    pairs: pd.DataFrame, guide_var: pd.DataFrame, gene_ids: pd.Index
) -> tuple[pd.DataFrame, dict]:
    """The pairs the object can test, keyed as the pipeline keys them.

    Parameters
    ----------
    pairs
        `pairs.run` output: ``guide_id``, ``gene_id``, ``pair_type`` and the
        four ``intended_target_*`` columns, among others.
    guide_var
        `target_keys` output.
    gene_ids
        The object's gene index.

    Returns
    -------
    tuple
        The kept pairs, with ``intended_target_key`` taken from the guide and
        low-cardinality text as categoricals; and a count of what was dropped,
        by reason, since the pipeline drops such rows without a word.
    """
    missing = {"guide_id", "gene_id", "pair_type", *TARGET_COLUMNS} - set(pairs)
    if missing:
        raise InputError(f"The pairs table lacks {', '.join(sorted(missing))}.")
    out = pairs.copy()
    out["gene_id"] = out["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    no_guide = ~out["guide_id"].isin(guide_var.index)
    no_gene = ~out["gene_id"].isin(gene_ids)
    dropped = {
        "pairs_in": len(out),
        "guide_absent": int(no_guide.sum()),
        "gene_absent": int((no_gene & ~no_guide).sum()),
        "guides_absent": sorted(out.loc[no_guide, "guide_id"].unique().tolist()),
        "genes_absent": sorted(out.loc[no_gene, "gene_id"].unique().tolist()),
    }
    out = out[~no_guide & ~no_gene].copy()
    if dropped["pairs_in"] and out.empty:
        raise InputError(
            f"None of the {dropped['pairs_in']} pairs names a guide and a gene the"
            f" counts carry (guides such as {dropped['guides_absent'][:3]}, genes"
            f" such as {dropped['genes_absent'][:3]}): the pairs were built for"
            " another library or annotation."
        )
    # Controls are renamed by bucketing; a targeting guide's element must agree,
    # or the pairs were built from another guide table than the counts carry.
    stated = guide_var.loc[out["guide_id"]]
    for column in TARGET_COLUMNS:
        if column.endswith(("_start", "_end")):
            mine = pd.to_numeric(out[column], errors="coerce").astype("Int64")
            theirs = stated[column].astype("Int64")
        else:
            mine = out[column].astype("string")
            theirs = stated[column].astype("string")
        mine, theirs = mine.reset_index(drop=True), theirs.reset_index(drop=True)
        same = (mine == theirs).fillna(False) | (mine.isna() & theirs.isna())
        differ = int((~same.astype(bool)).sum())
        if differ:
            raise InputError(
                f"{differ} pairs disagree with the guides on {column}: the pairs and the"
                " count objects were built from different guide tables."
            )
        out[column] = stated[column].to_numpy()
    out["intended_target_key"] = guide_var.loc[
        out["guide_id"], "intended_target_key"
    ].to_numpy()
    out = out[[*PAIR_COLUMNS, *[c for c in out.columns if c not in PAIR_COLUMNS]]]
    for column in (
        "gene_id",
        "guide_id",
        "intended_target_name",
        "intended_target_chr",
        "intended_target_key",
        "pair_type",
    ):
        out[column] = out[column].astype("category")
    dropped["pairs_out"] = len(out)
    return out.reset_index(drop=True), dropped


def build(
    channels: list,
    pairs: pd.DataFrame,
    *,
    moi: str,
    capture_method: str,
    min_guides: int | None = None,
    max_guides: int | None = None,
):
    """Assemble channels into one inference object.

    Parameters
    ----------
    channels
        QC'd per-channel MuData objects, each with ``gene`` and ``guide`` and
        a binary ``mod['guide'].layers['guide_assignment']``.
    pairs
        The pairs to test.
    moi
        ``high`` or ``low``.
    capture_method
        How guides were captured, as the pipeline records it (e.g.
        ``CROP-seq`` or ``direct capture``). Carried, not read.
    min_guides, max_guides
        Keep cells with at least, and at most, this many assigned guides.
        None keeps every cell.

    Returns
    -------
    tuple
        The MuData, and a report of what was kept and dropped.

    Raises
    ------
    InputError
        A channel lacks a modality or the guide calls, the channels disagree
        on their features, or the gene ids collide once unversioned.
    """
    import anndata as ad
    import mudata as md
    import scipy.sparse as sp

    if moi not in MOI:
        raise InputError(f"moi must be one of {', '.join(MOI)}, not {moi!r}.")
    if not channels:
        raise InputError("No channels given.")
    first = channels[0]
    report: dict = {"channels": [], "cells": {}}
    genes, guides = [], []
    for data in channels:
        for modality in ("gene", "guide"):
            if modality not in data.mod:
                raise InputError(f"A channel has no {modality!r} modality.")
        gene, guide = data["gene"], data["guide"]
        if ASSIGNMENT not in guide.layers:
            raise InputError(
                f"A channel has no mod['guide'].layers[{ASSIGNMENT!r}]. Call guides first;"
                " this step does not."
            )
        calls = sp.csr_matrix(guide.layers[ASSIGNMENT])
        calls.eliminate_zeros()
        guide.layers[ASSIGNMENT] = calls
        if calls.nnz and not np.isin(calls.data, (0, 1)).all():
            raise InputError(f"mod['guide'].layers[{ASSIGNMENT!r}] is not binary.")
        if not gene.var_names.equals(
            first["gene"].var_names
        ) or not guide.var_names.equals(first["guide"].var_names):
            raise InputError(
                "The channels do not share their genes and guides, in order."
            )
        both = gene.obs_names.intersection(guide.obs_names, sort=False)
        name = _measurement_set(data)
        report["channels"].append(
            {
                "measurement_set": name,
                "gene_cells": gene.n_obs,
                "guide_cells": guide.n_obs,
                "both": len(both),
            }
        )
        gene_part = gene[both].copy()
        guide_part = guide[both].copy()
        descriptors = data.obs.loc[both].drop(
            columns=[c for c in data.obs.columns if c in gene_part.obs], errors="ignore"
        )
        if "batch" in descriptors:
            descriptors = descriptors.rename(columns={"batch": "channel_batch"})
        gene_part.obs = gene_part.obs.join(descriptors)
        gene_part.obs["batch"] = gene_part.obs.get("measurement_set", name)
        genes.append(gene_part)
        guides.append(guide_part)

    gene = ad.concat(genes, merge="first", uns_merge="first")
    guide = ad.concat(guides, merge="first", uns_merge="first")
    if gene.obs_names.has_duplicates:
        raise InputError(
            "Two channels share a cell id; the index needs its measurement set."
        )

    calls = sp.csr_matrix(guide.layers[ASSIGNMENT])
    per_cell = np.diff(calls.indptr)
    keep = np.ones(guide.n_obs, dtype=bool)
    if min_guides is not None:
        keep &= per_cell >= min_guides
    if max_guides is not None:
        keep &= per_cell <= max_guides
    report["cells"] = {
        "in": int(guide.n_obs),
        "below_min_guides": int((per_cell < min_guides).sum())
        if min_guides is not None
        else 0,
        "above_max_guides": int((per_cell > max_guides).sum())
        if max_guides is not None
        else 0,
        "out": int(keep.sum()),
    }
    gene, guide = gene[keep].copy(), guide[keep].copy()

    from .guide_annotation.db import _unversioned

    versioned = gene.var_names
    unversioned = pd.Index(_unversioned(pd.Series(versioned)).to_numpy())
    if unversioned.has_duplicates:
        raise InputError(
            f"{unversioned.duplicated().sum()} gene ids collide once unversioned, e.g."
            f" {', '.join(unversioned[unversioned.duplicated()][:3])}."
        )
    gene.var.insert(0, "gene_id_versioned", np.asarray(versioned))
    gene.var_names = unversioned
    gene.var.index.name = None

    guide.var = target_keys(guide.var)
    guide.var.insert(0, "guide_id", guide.var_names.astype(str))

    gene_cov, guide_cov = cell_covariates(gene, guide)
    for column in gene_cov:
        gene.obs[column] = gene_cov[column]
    for column in guide_cov:
        guide.obs[column] = guide_cov[column]

    for part in (gene, guide):
        part.X = _counts(part.X)
    guide.layers[ASSIGNMENT] = sp.csr_matrix(guide.layers[ASSIGNMENT]).astype(np.uint8)
    guide.uns["moi"] = np.array([moi], dtype=object)
    guide.uns["capture_method"] = np.array([capture_method], dtype=object)

    # SCEPTRE fits these as log terms and the pipeline refuses a zero, but
    # filtering is the analysis's: count such cells, keep them.
    report["cells"]["zero_covariate"] = {
        column: int((frame[column] <= 0).sum())
        for column, frame in (
            ("total_gene_umis", gene_cov),
            ("num_expressed_genes", gene_cov),
            ("total_guide_umis", guide_cov),
        )
    }

    data = md.MuData({"gene": gene, "guide": guide})
    data.obs = pd.DataFrame(
        {
            "batch": gene.obs["batch"].astype(str).to_numpy(),
            "total_gene_umis": gene.obs["total_gene_umis"].to_numpy(),
            "total_guide_umis": guide.obs["total_guide_umis"].to_numpy(),
            "num_expressed_genes": gene.obs["num_expressed_genes"].to_numpy(),
            "percent_mito": gene.obs["percent_mito"].to_numpy(),
        },
        index=gene.obs_names,
    )
    kept, dropped = pairs_frame(pairs, guide.var, gene.var_names)
    data.uns["pairs_to_test"] = kept
    report["pairs"] = {k: v for k, v in dropped.items() if not isinstance(v, list)}
    report["pairs"]["guides_absent"] = dropped["guides_absent"][:50]
    report["pairs"]["genes_absent"] = dropped["genes_absent"][:50]
    for key in ("reference",):
        if key in first.uns:
            data.uns[key] = first.uns[key]
    data.uns["qc"] = {
        _measurement_set(channel): dict(channel.uns.get("qc", {}))
        for channel in channels
    }
    data.uns["inference_input"] = {"cells": report["cells"], "pairs": report["pairs"]}
    return data, report


def run(
    inputs: list[str | os.PathLike[str]],
    pairs: str | os.PathLike[str],
    outdir: str | os.PathLike[str],
    *,
    name: str,
    moi: str = "high",
    capture_method: str,
    min_guides: int | None = None,
    max_guides: int | None = None,
) -> RunResult:
    """Read channels and pairs, build the object, and write it.

    Parameters
    ----------
    inputs
        QC'd per-channel ``.h5mu`` files with guide calls.
    pairs
        ``pairs_to_test.csv`` from `pairs.run`.
    outdir
        Where the object goes. Nothing is written outside it.
    name
        The object's file stem, such as ``day4``.
    moi, capture_method, min_guides, max_guides
        As in `build`.

    Returns
    -------
    RunResult
        The path written, the report and the run record.
    """
    import mudata as md

    outdir = Path(outdir)
    parameters = {
        "name": name,
        "moi": moi,
        "capture_method": capture_method,
        "min_guides": min_guides,
        "max_guides": max_guides,
    }
    with run_record(
        outdir, STEP, RESULTS_DIR, parameters, logger, "Inference input"
    ) as (state, _):
        for path in [*inputs, pairs]:
            if not Path(path).is_file():
                raise FileNotFoundError(f"{path} does not exist.")
        state["inputs"]["channels"] = [str(p) for p in inputs]
        state["inputs"]["pairs"] = str(pairs)
        state["inputs"]["pairs_md5"] = file_md5(pairs)

        channels = []
        for path in inputs:
            logger.info("Reading %s", path)
            channels.append(md.read_h5mu(path))
        # Ids and names are text: a guide called 0001 must still match its counts.
        table = pd.read_csv(
            pairs,
            dtype={
                column: str
                for column in (
                    "guide_id",
                    "gene_id",
                    "gene_symbol",
                    "intended_target_name",
                    "intended_target_chr",
                )
            },
            keep_default_na=False,
            na_values=[""],
        )
        data, report = build(
            channels,
            table,
            moi=moi,
            capture_method=capture_method,
            min_guides=min_guides,
            max_guides=max_guides,
        )
        for line in _report_lines(report):
            logger.info(line)
        path = outdir / RESULTS_DIR / f"{name}.h5mu"
        data.write_h5mu(path)

        state["counts"] = {
            "channels": report["channels"],
            "cells": report["cells"],
            "pairs": report["pairs"],
            "genes": int(data["gene"].n_vars),
            "guides": int(data["guide"].n_vars),
        }
        return RunResult(outdir=outdir, path=path, report=report, metadata=state)


def _report_lines(report: dict) -> list[str]:
    """The report, as the log and the command print it."""
    cells, pairs = report["cells"], report["pairs"]
    return [
        f"{len(report['channels'])} channel{'s' * (len(report['channels']) != 1)},"
        f" {cells['in']} cells in both modalities",
        f"{cells['below_min_guides']} below the guide minimum, {cells['above_max_guides']} above the maximum,"
        f" {cells['out']} kept",
        f"{pairs['pairs_out']} of {pairs['pairs_in']} pairs kept: {pairs['guide_absent']} name a guide"
        f" the counts lack, {pairs['gene_absent']} a gene",
    ] + [
        f"{n} cells have {column} 0, which SCEPTRE cannot take the log of: filter them"
        " before the inference"
        for column, n in cells.get("zero_covariate", {}).items()
        if n
    ]


def _measurement_set(data) -> str:
    """The channel's measurement set, as `mudata_builder` records it."""
    value = data.uns.get("measurement_set")
    if isinstance(value, np.ndarray):
        value = value.item() if value.size == 1 else str(value)
    if value is None:
        sets = data["gene"].obs.get("measurement_set")
        value = str(sets.iloc[0]) if sets is not None and len(sets) else "unknown"
    return str(value)


def _counts(x):
    """Raw counts as CSR in the narrowest unsigned integer type that holds them."""
    import scipy.sparse as sp

    x = sp.csr_matrix(x)
    if x.nnz and not np.array_equal(x.data, np.round(x.data)):
        raise InputError(
            "A count matrix holds non-integer values; give raw UMI counts."
        )
    top = x.data.max() if x.nnz else 0
    dtype = np.uint16 if top <= np.iinfo(np.uint16).max else np.uint32
    return x.astype(dtype)

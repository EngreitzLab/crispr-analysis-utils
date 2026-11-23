"""How many genes a pair's link passes over, and whether those are read out.

Ported from the intervening-gene classification in the Voineagu lab's
``EnhGenePairs.R`` (`astrocyte_crispri
<https://github.com/Voineagulab/astrocyte_crispri/blob/main/5_FunctionalAnnotation/EnhGenePairs/EnhGenePairs.R>`_),
using this repository's own annotation store and readout panel in place of
its external gene table. Its TAD analysis is not ported: it needs a TAD call
file this package has no model for.

For each row of :mod:`crispr_analysis_utils.pairs`' output, every annotated
gene on the element's chromosome is ranked by distance from the element's
midpoint. A gene is "skipped" when it is at least as close as the pair's own
gene and sits on the same side of the element, in genome coordinates (not
flipped to the gene's orientation, the way ``pairs.distance_to_tss`` is): the
pair's link passes over it on the way to its own gene. Whether a skipped gene
is in the readout panel is what tells a "known effect" pair apart from one
resting on an untested gene.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ._files import file_md5, run_record, yaml_text

logger = logging.getLogger(__name__)

NEAREST = "nearest"
NO_SKIP = "no_skip"
SKIP_IN_PANEL = "skip_in_panel"
SKIP_OUTSIDE_PANEL = "skip_outside_panel"
TRANS = "trans"
CATEGORIES = (NEAREST, NO_SKIP, SKIP_IN_PANEL, SKIP_OUTSIDE_PANEL, TRANS)

COUNT_COLUMNS = (
    "n_nearer",
    "n_skipped",
    "n_skipped_in_panel",
    "n_skipped_outside_panel",
)

COLUMNS = {
    "intended_target_name": "The element, as pairs_to_test.csv names it.",
    "intended_target_chr": "Chromosome of the element.",
    "intended_target_start": "Start of the element, 0-based.",
    "intended_target_end": "End of the element, half-open.",
    "gene_id": "ENSEMBL gene id, without its version.",
    "gene_symbol": "The gene's symbol in the annotation store.",
    "pair_type": "positive_control or discovery, as pairs_to_test.csv has it.",
    "distance_to_tss": "Bases from the element's midpoint to the gene's TSS, from pairs_to_test.csv.",
    "in_panel": "TRUE when the gene is in the readout panel.",
    "n_nearer": (
        "How many other annotated genes have a TSS at least as close to the"
        " element's midpoint as this gene's. Empty for a trans pair."
    ),
    "n_skipped": (
        "Of those nearer genes, how many sit on the same side of the element"
        " as this gene: the ones this pair's link passes over. Empty for a"
        " trans pair."
    ),
    "n_skipped_in_panel": "Of the skipped genes, how many are in the readout panel.",
    "n_skipped_outside_panel": "Of the skipped genes, how many are outside the readout panel.",
    "category": (
        "nearest: no gene is closer. no_skip: nearer genes exist but none on"
        " this gene's side. skip_in_panel: at least one skipped gene is in"
        " the readout panel. skip_outside_panel: every skipped gene is"
        " outside it. trans: the gene is not on the element's chromosome, so"
        " nothing is ranked."
    ),
}
"""Every column of the skipped-genes table, and what it holds."""

TABLE = "skipped_genes.csv"
RESULTS_DIR = "skipped_genes"
STEP = "skipped_genes"


class SkippedGenesError(RuntimeError):
    """The inputs cannot be classified."""


@dataclass
class RunResult:
    """What `run` produced."""

    outdir: Path
    table: pd.DataFrame
    metadata: dict = field(default_factory=dict)


def _candidate_genes(elements: pd.DataFrame, tss: pd.DataFrame) -> pd.DataFrame:
    """Every annotated gene on each element's chromosome, with its signed distance.

    The distance is the element's midpoint minus the gene's TSS, in genome
    coordinates: not flipped to the gene's orientation, since what matters
    here is which side of the element a gene sits on, the same for a plus-
    or minus-strand gene.
    """
    from .pairs import TARGET_COLUMNS

    key = list(TARGET_COLUMNS)
    readout = tss.rename_axis("gene_id").reset_index()[["gene_id", "chr", "tss"]]
    midpoint = (
        elements["intended_target_start"] + elements["intended_target_end"]
    ) // 2
    candidates = elements.assign(midpoint=midpoint).merge(
        readout, left_on="intended_target_chr", right_on="chr"
    )
    candidates["distance"] = candidates["midpoint"] - candidates["tss"]
    return candidates[[*key, "gene_id", "distance"]]


def classify(pairs: pd.DataFrame, tss: pd.DataFrame, panel: pd.Index) -> pd.DataFrame:
    """Classify every unique element-gene pair in `pairs`.

    Parameters
    ----------
    pairs
        A `crispr_analysis_utils.pairs.build_pairs`-shaped table, or
        ``pairs_to_test.csv`` read back: at least `pairs.TARGET_COLUMNS` and
        ``gene_id``, with numeric ``intended_target_start``/``_end``.
    tss
        `crispr_analysis_utils.pairs.gene_tss` output: every annotated
        gene's TSS, not only the readout panel.
    panel
        The readout genes, unversioned, as `pairs.read_panel` returns.

    Returns
    -------
    pandas.DataFrame
        One row per unique element and gene of `pairs`: `pairs.TARGET_COLUMNS`,
        ``gene_id``, the four `COUNT_COLUMNS` and ``category``. A gene on a
        different chromosome than its element (a trans positive control)
        gets ``category="trans"`` and empty counts: there is nothing to rank
        it against.

    Raises
    ------
    SkippedGenesError
        A required column is missing, or a pair's gene is not in `tss`.
    """
    from .pairs import TARGET_COLUMNS

    key = list(TARGET_COLUMNS)
    missing = {*key, "gene_id"} - set(pairs)
    if missing:
        raise SkippedGenesError(f"The pairs table lacks {', '.join(sorted(missing))}.")

    unique_pairs = pairs[[*key, "gene_id"]].drop_duplicates().reset_index(drop=True)
    unknown = sorted(set(unique_pairs["gene_id"]) - set(tss.index))
    if unknown:
        raise SkippedGenesError(
            f"{len(unknown)} pair genes are not in the annotation, e.g."
            f" {', '.join(unknown[:3])}. The pairs table and the annotation"
            " must be the same GENCODE release."
        )

    elements = unique_pairs[key].drop_duplicates()
    candidates = _candidate_genes(elements, tss)

    targets = unique_pairs.rename(columns={"gene_id": "target_gene_id"}).merge(
        candidates.rename(
            columns={"gene_id": "target_gene_id", "distance": "target_distance"}
        ),
        on=[*key, "target_gene_id"],
        how="left",
    )
    trans = targets[targets["target_distance"].isna()].copy()
    same = targets[targets["target_distance"].notna()].copy()

    joined = same.merge(candidates, on=key, how="left")
    joined = joined[joined["gene_id"] != joined["target_gene_id"]]
    nearer = joined[joined["distance"].abs() <= joined["target_distance"].abs()].copy()
    nearer["same_side"] = np.sign(nearer["distance"]) == np.sign(
        nearer["target_distance"]
    )
    skipped = nearer[nearer["same_side"]].copy()
    skipped["in_panel"] = skipped["gene_id"].isin(panel)

    group = [*key, "target_gene_id"]
    counted = pd.concat(
        [
            nearer.groupby(group).size().rename("n_nearer"),
            skipped.groupby(group).size().rename("n_skipped"),
            skipped[skipped["in_panel"]]
            .groupby(group)
            .size()
            .rename("n_skipped_in_panel"),
            skipped[~skipped["in_panel"]]
            .groupby(group)
            .size()
            .rename("n_skipped_outside_panel"),
        ],
        axis=1,
    )
    same = same.merge(counted, on=group, how="left")
    same[list(COUNT_COLUMNS)] = same[list(COUNT_COLUMNS)].fillna(0).astype("Int64")
    same["category"] = np.select(
        [
            same["n_nearer"] == 0,
            same["n_skipped"] == 0,
            same["n_skipped_in_panel"] > 0,
        ],
        [NEAREST, NO_SKIP, SKIP_IN_PANEL],
        default=SKIP_OUTSIDE_PANEL,
    )

    trans[list(COUNT_COLUMNS)] = pd.NA
    for column in COUNT_COLUMNS:
        trans[column] = trans[column].astype("Int64")
    trans["category"] = TRANS

    table = pd.concat([same, trans], ignore_index=True).rename(
        columns={"target_gene_id": "gene_id"}
    )
    return table[[*key, "gene_id", *COUNT_COLUMNS, "category"]]


def counts(table: pd.DataFrame) -> dict:
    """How many pairs fall in each category, for the run record and qc.tsv."""
    return {
        category: int((table["category"] == category).sum()) for category in CATEGORIES
    }


def write_qc_tsv(path: str | os.PathLike[str], found: dict) -> None:
    """Counts, one per category, in the section/name/value shape the other steps use."""
    rows = [("category", category, value) for category, value in found.items()]
    pd.DataFrame(rows, columns=["section", "name", "value"]).to_csv(
        path, sep="\t", index=False
    )


def write_columns_yaml(path: str | os.PathLike[str]) -> None:
    """Write skipped_genes_columns.yaml, generated from `COLUMNS` so it cannot drift."""
    lines = ["# What each column of this folder's output holds.", "", f"{TABLE}:"]
    lines += [f"  {name}: {yaml_text(text)}" for name, text in COLUMNS.items()]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_table(path: str | os.PathLike[str], table: pd.DataFrame) -> None:
    """Write the table, with booleans as TRUE and FALSE like the other outputs."""
    out = table.copy()
    out["in_panel"] = out["in_panel"].map({True: "TRUE", False: "FALSE"})
    out.to_csv(path, index=False)


def run(
    pairs: str | os.PathLike[str],
    store: str | os.PathLike[str],
    panel: str | os.PathLike[str],
    outdir: str | os.PathLike[str],
    *,
    gene_column: str = "gene_id",
) -> RunResult:
    """Classify every pair and write the table.

    Parameters
    ----------
    pairs
        ``pairs_to_test.csv`` from `crispr_analysis_utils.pairs.run`.
    store
        The annotation store the pairs were built with.
    panel
        The readout genes: the same file given to `pairs.run` as its panel.
    outdir
        Where the table goes. Nothing is written outside it.
    gene_column
        The panel's column of ENSEMBL gene ids.

    Returns
    -------
    RunResult
        The table and the run record.
    """
    from .guide_annotation.db import load_store
    from .pairs import TARGET_COLUMNS, gene_tss, read_panel

    outdir = Path(outdir)
    parameters = {"gene_column": gene_column}
    with run_record(outdir, STEP, RESULTS_DIR, parameters, logger, "Skipped genes") as (
        state,
        _,
    ):
        for name, path in (("pairs", pairs), ("store", store), ("panel", panel)):
            if not Path(path).is_file():
                raise FileNotFoundError(f"The {name} {path} does not exist.")
            state["inputs"][name] = str(path)
            state["inputs"][f"{name}_md5"] = file_md5(path)

        pairs_table = pd.read_csv(pairs, dtype=str, keep_default_na=False)
        missing = {*TARGET_COLUMNS, "gene_id"} - set(pairs_table)
        if missing:
            raise SkippedGenesError(
                f"The pairs table lacks {', '.join(sorted(missing))}."
            )
        for column in ("intended_target_start", "intended_target_end"):
            pairs_table[column] = pd.to_numeric(pairs_table[column], errors="coerce")

        tss = gene_tss(load_store(store))
        genes = read_panel(panel, gene_column)
        logger.info(
            "%d pairs, %d annotated genes, %d panel genes",
            len(pairs_table),
            len(tss),
            len(genes),
        )

        classified = classify(pairs_table, tss, genes)
        context = pairs_table[
            [
                *TARGET_COLUMNS,
                "gene_id",
                "gene_symbol",
                "pair_type",
                "distance_to_tss",
                "in_panel",
            ]
        ].drop_duplicates([*TARGET_COLUMNS, "gene_id"])
        context["in_panel"] = context["in_panel"].str.strip().str.upper() == "TRUE"
        table = context.merge(classified, on=[*TARGET_COLUMNS, "gene_id"], how="left")
        table = table[list(COLUMNS)]

        write_table(outdir / RESULTS_DIR / TABLE, table)
        write_columns_yaml(outdir / RESULTS_DIR / "skipped_genes_columns.yaml")
        found = counts(table)
        write_qc_tsv(outdir / RESULTS_DIR / "qc.tsv", found)

        state["counts"] = found
        logger.info(
            "Skipped genes complete: %d of %d pairs skip a gene in the panel",
            found[SKIP_IN_PANEL],
            len(table),
        )
        return RunResult(outdir=outdir, table=table, metadata=state)

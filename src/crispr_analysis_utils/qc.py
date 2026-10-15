"""Quality-control figures, each published as three files.

A figure nobody can check is an assertion. Every one here writes the numbers
it was drawn from, a caption that says what it shows, and a PDF at print
resolution, under one stem:

    <stem>.pdf      the figure, 300 dpi, vector text
    <stem>.tsv.gz   the numbers, one row per plotted point
    <stem>.txt      the caption, as it would appear under the figure

The data file is the point: a reviewer can redraw the figure, and a reader
can find the number behind a claim without rerunning anything.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd

DPI = 300
"""Print resolution, for the raster parts of an otherwise vector PDF."""

OKABE_ITO = (
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000",
)
"""Categorical colours, colourblind-safe. `cividis` for continuous scales."""


class FigureError(RuntimeError):
    """The figure cannot be drawn from what was given."""


def style():
    """Set the colours and the defaults every figure here shares.

    Returns
    -------
    module
        `matplotlib.pyplot`, so a caller need not import it separately.
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams.update(
        {
            "axes.prop_cycle": mpl.cycler(color=list(OKABE_ITO)),
            "image.cmap": "cividis",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "0.8",
            "grid.linewidth": 0.6,
            "axes.axisbelow": True,
            "figure.dpi": 140,
            "savefig.dpi": DPI,
            "savefig.facecolor": "white",
            # Text stays text in the PDF, so a figure can be relabelled
            # without redrawing it.
            "pdf.fonttype": 42,
        }
    )
    return plt


def publish(
    stem: str | os.PathLike[str],
    figure,
    data: pd.DataFrame,
    caption: str,
) -> dict[str, Path]:
    """Write a figure, its numbers and its caption under one stem.

    Parameters
    ----------
    stem
        Path without a suffix. Three files are written beside it.
    figure
        A matplotlib figure.
    data
        The numbers the figure was drawn from, one row per plotted point.
    caption
        One or two sentences saying what the figure shows. Written as given,
        so it can go straight under the figure.

    Returns
    -------
    dict
        The three paths, by suffix.

    Raises
    ------
    FigureError
        The data is empty, or the caption is blank. Either makes the figure
        unpublishable, and failing here is cheaper than finding out later.
    """
    if data.empty:
        raise FigureError("A figure must be published with the numbers behind it.")
    caption = caption.strip()
    if not caption:
        raise FigureError("A figure must be published with a caption.")

    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    written = {
        ".pdf": stem.with_suffix(".pdf"),
        ".tsv.gz": stem.with_name(stem.name + ".tsv.gz"),
        ".txt": stem.with_suffix(".txt"),
    }
    figure.savefig(written[".pdf"], dpi=DPI, bbox_inches="tight")
    with gzip.open(written[".tsv.gz"], "wt", encoding="utf-8") as handle:
        data.to_csv(handle, sep="\t", index=False)
    written[".txt"].write_text(caption + "\n", encoding="utf-8")
    return written


def umis_per_cell(data, modality: str) -> np.ndarray:
    """Total UMIs for every cell of the object, zero where it has none.

    A cell the modality does not hold has no UMIs in it, which is not the
    same as being absent from the object: the union keeps both, and a knee
    plot that dropped them would overstate how deep the channel is.
    """
    held = np.asarray(data.obsm[modality]).ravel()
    umis = np.zeros(data.n_obs)
    umis[held] = np.asarray(data[modality].X.sum(axis=1)).ravel()
    return umis


def filter_cells(data, min_umis: int, modality: str = "gene"):
    """Keep the barcodes whose `modality` UMI count reaches `min_umis`.

    Parameters
    ----------
    data : mudata.MuData
        The object to filter.
    min_umis : int
        The floor, read off the knee of this channel. There is no default:
        where the knee falls depends on the channel's depth, so a floor
        carried over from another one is a guess.
    modality : str
        The modality whose UMIs decide, `"gene"` by default. A barcode with
        real RNA and no guide is kept, since an unassigned cell is a result.

    Returns
    -------
    tuple of (mudata.MuData, numpy.ndarray)
        The filtered copy, and the boolean mask over the original barcodes.
    """
    keep = umis_per_cell(data, modality) >= min_umis
    return data[keep].copy(), keep


def knee(data, modalities: tuple[str, ...] = ("gene", "guide")):
    """Barcodes ranked by UMI count, log on both axes.

    Parameters
    ----------
    data
        A MuData as `mudata_builder` assembles it.
    modalities
        Which modalities to draw, each its own curve.

    Returns
    -------
    tuple
        The figure, the numbers, and the caption, ready for `publish`.

    Raises
    ------
    FigureError
        No modality has a barcode with a UMI.
    """
    plt = style()
    figure, axes = plt.subplots(figsize=(6, 4.5))

    rows = []
    for modality in modalities:
        ranked = np.sort(umis_per_cell(data, modality))[::-1]
        ranked = ranked[ranked > 0]
        if not len(ranked):
            continue
        rank = np.arange(1, len(ranked) + 1)
        axes.plot(rank, ranked, label=modality, linewidth=1.5)
        rows.append(
            pd.DataFrame(
                {"modality": modality, "rank": rank, "umis": ranked.astype(int)}
            )
        )
    if not rows:
        raise FigureError("No modality has a barcode with a UMI to rank.")

    axes.set(
        xscale="log",
        yscale="log",
        xlabel="barcode rank",
        ylabel="UMIs",
    )
    axes.legend(frameon=False)
    figure.tight_layout()

    numbers = pd.concat(rows, ignore_index=True)
    counted = ", ".join(
        f"{int((numbers['modality'] == m).sum()):,} {m}"
        for m in numbers["modality"].unique()
    )
    caption = (
        "Barcodes ranked by total UMIs, both axes logarithmic. Every barcode "
        f"the object holds is ranked ({counted} barcodes with at least one "
        "UMI); a barcode a modality does not hold counts as zero and so "
        "falls outside the curve. The knee is where real cells give way to "
        "ambient background, and its position is what a count threshold has "
        "to justify itself against."
    )
    return figure, numbers, caption


def mad_outlier(
    values,
    nmads: float = 5.0,
    side: str = "both",
    log: bool = False,
) -> np.ndarray:
    """Flag the values that sit more than `nmads` MADs from the median.

    The median absolute deviation is the spread of the channel measured
    against itself, so a threshold built from it moves with the channel's
    own depth instead of being carried over from another one. The MAD is
    rescaled to a standard deviation (`scale="normal"`), so `nmads` reads
    as sigmas on a normal distribution.

    Parameters
    ----------
    values : array-like
        One number per cell. `NaN` is never flagged: a cell the modality
        does not hold was not measured.
    nmads : float
        How many MADs from the median is still acceptable, 5 by default.
    side : {"both", "upper", "lower"}
        Which tail to flag. A count metric is two-sided; a contamination
        fraction only has an upper tail worth cutting.
    log : bool
        Work on `log10(x + 1)` rather than the raw value. Off by default:
        a count distribution wants it, a percentage does not.

    Returns
    -------
    numpy.ndarray
        Boolean, true where the value is an outlier.

    Raises
    ------
    FigureError
        The MAD is zero, so every value equals the median and no threshold
        can be read off it. This is what an unfiltered barcode space looks
        like: pick the cells first, then measure their spread.
    """
    from scipy.stats import median_abs_deviation

    if side not in ("both", "lower", "upper"):
        raise FigureError(f"side must be both, upper or lower, not {side!r}")
    numbers = np.asarray(values, dtype=float)
    measured = ~np.isnan(numbers)
    scale = np.log10(numbers[measured] + 1) if log else numbers[measured]
    middle = np.median(scale)
    spread = median_abs_deviation(scale, scale="normal")
    if spread == 0:
        raise FigureError(
            "The MAD is zero, so the threshold would be the median itself. "
            "Run this on cells, not on the whole barcode space."
        )

    low, high = middle - nmads * spread, middle + nmads * spread
    flagged = np.zeros(numbers.shape, dtype=bool)
    if side in ("both", "lower"):
        flagged[measured] |= scale < low
    if side in ("both", "upper"):
        flagged[measured] |= scale > high
    return flagged


KEPT = OKABE_ITO[4]
"""Blue, for what a threshold keeps: the cells the analysis goes on with."""

REMOVED = "0.7"
"""Grey, for what it removes. Neutral, so the eye reads the kept cells first."""

LOG_METRICS = ("total_counts", "n_genes_by_counts")
"""Metrics whose distribution is read on a log axis, because they are counts."""


def metric_distributions(
    frame: pd.DataFrame,
    flagged=None,
    log: tuple[str, ...] = LOG_METRICS,
    bins: int = 60,
):
    """One histogram per QC metric, over the cells the frame describes.

    Parameters
    ----------
    frame : pandas.DataFrame
        One column per metric, one row per cell, as `cell_fractions`
        returns. Rows that are entirely `NaN` are cells the modality does
        not hold, and are dropped rather than drawn as zero.
    flagged : array-like of bool, optional
        True for the cells a threshold rejects. Drawn in a second colour on
        top of the kept ones, so what a cut removes is visible against what
        it leaves.
    log : tuple of str
        Metrics binned and drawn on a log axis. A count distribution spans
        orders of magnitude; a percentage does not.
    bins : int
        Bins per panel.

    Returns
    -------
    tuple of (matplotlib.figure.Figure, pandas.DataFrame, str)
        The figure, one row per bin per metric, and the caption.

    Raises
    ------
    FigureError
        The frame has no columns, or no cell survives the `NaN` drop.
    """
    plt = style()

    if frame.shape[1] == 0:
        raise FigureError("No metric to draw.")
    measured = frame.dropna(how="all")
    if measured.empty:
        raise FigureError("No measured cell to draw.")
    cut = (
        np.zeros(len(measured), dtype=bool)
        if flagged is None
        else np.asarray(flagged, dtype=bool)[frame.index.isin(measured.index)]
    )

    metrics = list(measured.columns)
    columns = min(3, len(metrics))
    rows = -(-len(metrics) // columns)
    figure, axes = plt.subplots(
        rows, columns, figsize=(4.0 * columns, 3.0 * rows), squeeze=False
    )

    records = []
    for axis, metric in zip(axes.ravel(), metrics, strict=False):
        values = measured[metric].to_numpy(dtype=float)
        if metric in log:
            # A zero has no place on a log axis, so the range starts at the
            # smallest value actually present.
            positive = values[values > 0]
            floor = float(positive.min()) if positive.size else 1.0
            edges = np.logspace(np.log10(floor), np.log10(np.nanmax(values)), bins + 1)
            axis.set_xscale("log")
        else:
            edges = np.linspace(np.nanmin(values), np.nanmax(values), bins + 1)

        kept_counts, _ = np.histogram(values[~cut], bins=edges)
        cut_counts, _ = np.histogram(values[cut], bins=edges)
        axis.stairs(kept_counts, edges, fill=True, color=KEPT, label="kept")
        if cut.any():
            axis.stairs(
                cut_counts,
                edges,
                baseline=kept_counts,
                fill=True,
                color=REMOVED,
                label="flagged",
            )
        # Black, so the median is not read as another grey category.
        axis.axvline(np.nanmedian(values), color="black", linewidth=1.0, linestyle="--")
        axis.set(xlabel=metric, ylabel="cells")
        axis.grid(which="major", color="0.85", linewidth=0.6)
        axis.set_axisbelow(True)

        records.append(
            pd.DataFrame(
                {
                    "metric": metric,
                    "bin_left": edges[:-1],
                    "bin_right": edges[1:],
                    "n_kept": kept_counts,
                    "n_flagged": cut_counts,
                }
            )
        )

    for axis in axes.ravel()[len(metrics) :]:
        axis.set_visible(False)
    if cut.any():
        axes.ravel()[0].legend(frameon=False, fontsize="small")
    figure.tight_layout()

    numbers = pd.concat(records, ignore_index=True)
    medians = ", ".join(
        f"{metric} {measured[metric].median():,.2f}" for metric in metrics
    )
    counted = (
        f"{int((~cut).sum()):,} kept and {int(cut.sum()):,} flagged"
        if cut.any()
        else f"{len(measured):,}"
    )
    caption = (
        f"Distribution of each quality-control metric over {counted} cells. "
        "Counts are binned and drawn on a log axis, fractions on a linear "
        "one; the dashed line is the median (" + medians + "). The flagged "
        "cells are stacked on the kept ones, so a threshold's effect is "
        "read against the distribution it was taken from."
        if cut.any()
        else f"Distribution of each quality-control metric over {counted} "
        "cells. Counts are binned and drawn on a log axis, fractions on a "
        "linear one; the dashed line is the median (" + medians + ")."
    )
    return figure, numbers, caption


MITO_CONTIG = "chrM"
"""Mitochondrial genes are the ones on the mitochondrial contig."""

RIBO_PREFIXES = ("RPS", "RPL")
"""Cytoplasmic ribosomal protein genes, by symbol."""


def gene_classes(var: pd.DataFrame) -> dict[str, np.ndarray]:
    """Which genes are mitochondrial, ribosomal and read out by the assay.

    Mitochondrial genes are taken from the contig, not an ``MT`` symbol
    prefix. On GENCODE v43 the contig and the ``MT-`` prefix agree exactly,
    37 genes, but a prefix without the hyphen also matches 529 nuclear genes,
    mostly metallothioneins, and the MTRNR2L family: nuclear copies of
    MT-RNR2 that are highly expressed and would inflate the mitochondrial
    share of exactly the cells one would then discard as dying.

    Ribosomal genes are the RPS and RPL families by symbol, where no contig
    helps, split by whether the gene is protein-coding: the pseudogenes share
    the naming, are transcribed, and belong under their own label rather than
    swelling the ribosomal share.

    Parameters
    ----------
    var
        ``mod['gene'].var``, needing ``gene_chr`` and ``symbol``, and
        ``readout`` for the targeted class.

    Returns
    -------
    dict
        A boolean mask per class, over the genes in order.
    """
    symbol = (
        var.get("symbol", pd.Series(index=var.index, dtype="string"))
        .astype("string")
        .fillna("")
    )
    contig = (
        var.get("gene_chr", pd.Series(index=var.index, dtype="string"))
        .astype("string")
        .fillna("")
    )
    # An RPS or RPL name is shared by the ribosomal proteins and by their
    # many pseudogenes, which are transcribed and would otherwise inflate the
    # ribosomal share. Counted, but under their own label.
    ribosomal = symbol.str.upper().str.startswith(RIBO_PREFIXES).to_numpy()
    kind = (
        var.get("gene_type", pd.Series(index=var.index, dtype="string"))
        .astype("string")
        .fillna("")
    )
    coding = (kind == "protein_coding").to_numpy()
    classes = {
        "mito": (contig == MITO_CONTIG).to_numpy(),
        "ribo": ribosomal & coding,
        "ribo_pseudogene": ribosomal & ~coding,
    }
    if "readout" in var:
        classes["targeted"] = var["readout"].fillna(False).to_numpy().astype(bool)
    return classes


def cell_fractions(data, modality: str = "gene") -> pd.DataFrame:
    """The share of each cell's UMIs coming from each class of gene.

    The arithmetic is scanpy's `calculate_qc_metrics`, so the columns are
    named as the IGVF pipeline names them and the two can be compared
    directly. What this adds is which genes belong to which class, and the
    handling of a cell the modality does not hold.

    Parameters
    ----------
    data
        A MuData as `mudata_builder` assembles it.
    modality
        Which modality to measure; the gene counts by default.

    Returns
    -------
    pandas.DataFrame
        One row per cell of the whole object, indexed by cell:
        ``total_counts``, ``n_genes_by_counts`` and a ``pct_counts_`` column
        per class. A cell the modality does not hold has nothing measured,
        so its percentages are NaN rather than zero.
    """
    import scanpy as sc

    counts = data[modality]
    classes = gene_classes(counts.var)

    # scanpy reads the class membership off var, by the names given here.
    measured = counts.copy()
    for name, mask in classes.items():
        measured.var[name] = mask
    sc.pp.calculate_qc_metrics(
        measured, qc_vars=list(classes), percent_top=None, log1p=False, inplace=True
    )

    wanted = ["total_counts", "n_genes_by_counts"] + [
        f"pct_counts_{name}" for name in classes
    ]
    frame = measured.obs[wanted]
    # Back onto every cell of the object: one the modality does not hold has
    # no measurement, which is not a measurement of zero.
    return frame.reindex(pd.Index(data.obs_names))

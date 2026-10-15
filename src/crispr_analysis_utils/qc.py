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
    numbers = np.asarray(values, dtype=float)
    measured = ~np.isnan(numbers)
    low, high = mad_bounds(numbers[measured], nmads=nmads, side=side, log=log)

    flagged = np.zeros(numbers.shape, dtype=bool)
    flagged[measured] |= numbers[measured] < low
    flagged[measured] |= numbers[measured] > high
    return flagged


def mad_bounds(
    values,
    nmads: float = 5.0,
    side: str = "both",
    log: bool = False,
) -> tuple[float, float]:
    """The thresholds `mad_outlier` applies, in the values' own units.

    A figure draws these and colours by them, so the line on the plot and
    the cut in the table cannot disagree.

    Parameters
    ----------
    values : array-like
        One number per cell, with no `NaN`.
    nmads, side, log
        As `mad_outlier`.

    Returns
    -------
    tuple of float
        Lower and upper threshold. A tail that `side` does not cut comes
        back as `-inf` or `inf`.

    Raises
    ------
    FigureError
        `side` is not one of the three, or the MAD is zero.
    """
    from scipy.stats import median_abs_deviation

    if side not in ("both", "lower", "upper"):
        raise FigureError(f"side must be both, upper or lower, not {side!r}")
    numbers = np.asarray(values, dtype=float)
    scale = np.log10(numbers + 1) if log else numbers
    middle = np.median(scale)
    spread = median_abs_deviation(scale, scale="normal")
    if spread == 0:
        raise FigureError(
            "The MAD is zero, so the threshold would be the median itself. "
            "Run this on cells, not on the whole barcode space."
        )

    low, high = middle - nmads * spread, middle + nmads * spread
    if log:
        low, high = 10**low - 1, 10**high - 1
    return (
        float(low) if side in ("both", "lower") else -np.inf,
        float(high) if side in ("both", "upper") else np.inf,
    )


KEPT = OKABE_ITO[4]
"""Blue, for what a threshold keeps: the cells the analysis goes on with."""

REMOVED = "0.7"
"""Grey, for what it removes. Neutral, so the eye reads the kept cells first."""

LOG_METRICS = ("total_counts", "n_genes_by_counts")
"""Metrics whose distribution is read on a log axis, because they are counts."""


def metric_distributions(
    frame: pd.DataFrame,
    nmads: float = 5.0,
    log: tuple[str, ...] = LOG_METRICS,
    sides: dict[str, str] | None = None,
    bins: int = 60,
):
    """One histogram per QC metric, split at that metric's own threshold.

    Every panel computes its MAD band, draws it, and colours by it, so the
    line on the figure is the cut the table would make.

    Parameters
    ----------
    frame : pandas.DataFrame
        One column per metric, one row per cell, as `cell_fractions`
        returns. Rows that are entirely `NaN` are cells the modality does
        not hold, and are dropped rather than drawn as zero.
    nmads : float
        Width of each band, in MADs.
    log : tuple of str
        Metrics binned, thresholded and drawn on a log axis. A count
        distribution spans orders of magnitude; a percentage does not.
    sides : dict, optional
        Which tail each metric cuts. The default is `"both"` for a count
        and `"upper"` for everything else, since a cell with unusually
        little mitochondrial signal is not a problem.
    bins : int
        Bins per panel.

    Returns
    -------
    tuple of (matplotlib.figure.Figure, pandas.DataFrame, str)
        The figure, one row per bin per metric carrying both thresholds,
        and the caption.

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

    metrics = list(measured.columns)
    sides = sides or {}
    columns = min(3, len(metrics))
    rows = -(-len(metrics) // columns)
    figure, axes = plt.subplots(
        rows, columns, figsize=(4.0 * columns, 3.0 * rows), squeeze=False
    )

    records = []
    unthresholded: list[str] = []
    for axis, metric in zip(axes.ravel(), metrics, strict=False):
        values = measured[metric].to_numpy(dtype=float)
        on_log = metric in log
        side = sides.get(metric, "both" if on_log else "upper")
        try:
            low, high = mad_bounds(values, nmads=nmads, side=side, log=on_log)
        except FigureError:
            # No spread, so no threshold. The panel still shows the
            # distribution rather than taking the whole figure down.
            low, high = -np.inf, np.inf
            unthresholded.append(metric)
        within = (values >= low) & (values <= high)

        if on_log:
            # A zero has no place on a log axis, so the range starts at the
            # smallest value actually present.
            positive = values[values > 0]
            floor = float(positive.min()) if positive.size else 1.0
            edges = np.logspace(np.log10(floor), np.log10(np.nanmax(values)), bins + 1)
            axis.set_xscale("log")
        else:
            floor = float(np.nanmin(values))
            edges = np.linspace(floor, np.nanmax(values), bins + 1)
        # logspace rounds its endpoints, so a cell sitting exactly on the
        # smallest or largest value can fall outside its own range and go
        # uncounted. Pin the outer edges back onto the data.
        edges[0] = min(edges[0], floor)
        edges[-1] = max(edges[-1], float(np.nanmax(values)))

        kept_counts, _ = np.histogram(values[within], bins=edges)
        cut_counts, _ = np.histogram(values[~within], bins=edges)
        axis.stairs(kept_counts, edges, fill=True, color=KEPT, label="kept")
        if (~within).any():
            # Stacked, so the top of the grey is kept + cut. Passing the cut
            # count alone fills between it and the baseline, which paints
            # over the kept cells instead of sitting on them.
            axis.stairs(
                kept_counts + cut_counts,
                edges,
                baseline=kept_counts,
                fill=True,
                color=REMOVED,
                label="removed",
            )
        for threshold in (low, high):
            if np.isfinite(threshold):
                axis.axvline(threshold, color="black", linewidth=1.0, linestyle="--")
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
                    "n_removed": cut_counts,
                    "threshold_low": low,
                    "threshold_high": high,
                }
            )
        )

    for axis in axes.ravel()[len(metrics) :]:
        axis.set_visible(False)
    axes.ravel()[0].legend(frameon=False, fontsize="small")
    figure.tight_layout()

    numbers = pd.concat(records, ignore_index=True)
    thresholds = ", ".join(
        f"{metric} "
        + " to ".join(
            f"{bound:,.2f}"
            for bound in (
                numbers.loc[numbers["metric"] == metric, "threshold_low"].iloc[0],
                numbers.loc[numbers["metric"] == metric, "threshold_high"].iloc[0],
            )
            if np.isfinite(bound)
        )
        for metric in metrics
    )
    removed = int(numbers["n_removed"].sum())
    caption = (
        f"Distribution of each quality-control metric over {len(measured):,} "
        "cells. Counts are binned and drawn on a log axis, fractions on a "
        "linear one. The dashed lines are that metric's own "
        f"{nmads:g}-MAD thresholds ({thresholds}); cells beyond them are "
        f"grey. Each panel cuts on its own metric, so the {removed:,} grey "
        "bars summed across panels double-count a cell that several "
        "metrics reject."
    )
    if unthresholded:
        caption += (
            " No threshold is drawn for "
            + ", ".join(unthresholded)
            + ", whose values have no spread to measure one from."
        )
    return figure, numbers, caption


def metric_scatter(
    frame: pd.DataFrame,
    x: str = "total_counts",
    y: str = "pct_counts_targeted",
    x_threshold: float | None = None,
    y_threshold: float | None = None,
    log: tuple[str, ...] = LOG_METRICS,
    marginals: bool = True,
    unit: str = "cells",
    groups=None,
):
    """Two QC metrics against each other, cut into quadrants by thresholds.

    A cell can be poor on either metric alone, and the pair says which.
    The thresholds are given, not derived: a MAD band is one way to reach
    them and a library that has been looked at closely is another.

    Parameters
    ----------
    frame : pandas.DataFrame
        One row per cell, carrying at least the two columns named.
    x, y : str
        Columns to draw. Rows with `NaN` in either are dropped.
    x_threshold, y_threshold : float, optional
        Cuts, drawn as dashed lines. A cell is kept when it is at or above
        both; a threshold left out keeps that axis whole.
    log : tuple of str
        Metrics drawn on a log axis.
    marginals : bool
        Draw each metric's own histogram against its axis.
    unit : str
        What a point is, for the caption. Cells by default; the same figure
        serves a per-gene comparison.
    groups : array-like, optional
        A label per row. Given, the colour carries the label rather than
        whether the point clears the thresholds, which stay as lines. A row
        with no label is drawn grey and behind the rest. Use it where
        nothing is being filtered and the thresholds only mark where things
        were expected to fall.

    Returns
    -------
    tuple of (matplotlib.figure.Figure, pandas.DataFrame, str)
        The figure, one row per cell with its quadrant, and the caption.

    Raises
    ------
    FigureError
        A column is missing, or no cell has both metrics.
    """
    plt = style()

    for name in (x, y):
        if name not in frame:
            raise FigureError(f"{name!r} is not a column of the frame.")
    measured = frame[[x, y]].dropna()
    if measured.empty:
        raise FigureError(
            f"No {unit[:-1] if unit.endswith('s') else unit} has both metrics."
        )

    xs = measured[x].to_numpy(dtype=float)
    ys = measured[y].to_numpy(dtype=float)
    above_x = xs >= x_threshold if x_threshold is not None else np.ones(len(xs), bool)
    above_y = ys >= y_threshold if y_threshold is not None else np.ones(len(ys), bool)
    kept = above_x & above_y
    # A log axis cannot show a zero, and matplotlib drops one without
    # saying so. Keep it in the numbers and account for it in the caption.
    drawn = np.ones(len(xs), dtype=bool)
    for name, values in ((x, xs), (y, ys)):
        if name in log:
            drawn &= values > 0

    # A few hundred genes and ten thousand cells want different points: one
    # needs to be seen, the other needs not to be a solid block.
    size, alpha = (2.0, 0.4) if len(measured) > 2_000 else (14.0, 0.7)

    if marginals:
        figure = plt.figure(figsize=(6.0, 6.0))
        grid = figure.add_gridspec(
            2, 2, width_ratios=(5, 1), height_ratios=(1, 5), wspace=0.05, hspace=0.05
        )
        axis = figure.add_subplot(grid[1, 0])
        top = figure.add_subplot(grid[0, 0], sharex=axis)
        right = figure.add_subplot(grid[1, 1], sharey=axis)
    else:
        figure, axis = plt.subplots(figsize=(5.5, 4.5))
        top = right = None

    if groups is None:
        layers = [(kept & drawn, KEPT, "kept"), (~kept & drawn, REMOVED, "removed")]
    else:
        labels = pd.Series(groups, index=frame.index).reindex(measured.index)
        named = [str(name) for name in pd.unique(labels.dropna())]
        # Unlabelled first, so the classes being compared sit on top of it.
        layers = [(labels.isna().to_numpy() & drawn, REMOVED, None)]
        layers += [
            ((labels == name).to_numpy() & drawn, OKABE_ITO[i % len(OKABE_ITO)], name)
            for i, name in enumerate(named)
        ]

    # Rasterised: one point per cell stays a vector object otherwise, and a
    # channel of 10,000 cells turns into a PDF nobody can open.
    for mask, colour, label in layers:
        if mask.any():
            axis.scatter(
                xs[mask],
                ys[mask],
                s=size,
                alpha=alpha,
                linewidths=0,
                color=colour,
                label=label,
                rasterized=True,
            )
    if x in log:
        axis.set_xscale("log")
    if y in log:
        axis.set_yscale("log")
    if x_threshold is not None:
        axis.axvline(x_threshold, color="black", linewidth=1.0, linestyle="--")
    if y_threshold is not None:
        axis.axhline(y_threshold, color="black", linewidth=1.0, linestyle="--")
    axis.set(xlabel=x, ylabel=y)
    axis.grid(which="major", color="0.85", linewidth=0.6)
    axis.set_axisbelow(True)
    # Only worth a legend when there are two groups to tell apart.
    if groups is not None or (kept.any() and not kept.all()):
        axis.legend(frameon=False, fontsize="small", loc="lower right")

    for panel, values, name, orientation in (
        (top, xs, x, "vertical"),
        (right, ys, y, "horizontal"),
    ):
        if panel is None:
            continue
        shown = values[drawn]
        edges = (
            np.logspace(np.log10(shown.min()), np.log10(shown.max()), 60)
            if name in log
            else 60
        )
        panel.hist(shown, bins=edges, color="0.8", orientation=orientation)
        panel.axis("off")

    counts = {
        "kept": int(kept.sum()),
        "low_" + x: int((~above_x & above_y).sum()),
        "low_" + y: int((above_x & ~above_y).sum()),
        "low_both": int((~above_x & ~above_y).sum()),
    }
    numbers = pd.DataFrame(
        {
            unit.rstrip("s"): measured.index,
            x: xs,
            y: ys,
            "above_" + x: above_x,
            "above_" + y: above_y,
            "kept": kept,
            "drawn": drawn,
        }
    )
    if groups is not None:
        numbers["group"] = labels.to_numpy()
    stated = ", ".join(
        f"{name} at {value:,.4g}"
        for name, value in ((x, x_threshold), (y, y_threshold))
        if value is not None
    )
    caption = (
        f"{x} against {y} for {len(measured):,} {unit}, one point each. "
        + (f"Dashed lines are the thresholds ({stated}). " if stated else "")
        + (
            "Colour is the group each point belongs to, and a point in no "
            "group is grey. Counts by quadrant: "
            if groups is not None
            else "Points at or above every threshold are blue and the rest grey: "
        )
        + ", ".join(f"{label} {number:,}" for label, number in counts.items())
        + ". Histograms on the margins are each metric's own distribution "
        + f"over the same {unit}."
    )
    if not drawn.all():
        caption += (
            f" {int((~drawn).sum()):,} of the {unit} hold a zero on a log "
            "axis and so cannot be drawn; they are in the data file, marked "
            "by the drawn column, and are counted in the quadrants above."
        )
    return figure, numbers, caption


def mean_upm(counts) -> np.ndarray:
    """Mean UMIs per million per gene, over the cells of one modality.

    Each cell is scaled to a million UMIs first, so a deep cell does not
    outweigh a shallow one, and the mean is taken across cells.

    `normalization.counts_per_million` is the same arithmetic for the dense
    guide-by-sample tables of a bulk screen; a single-cell matrix is sparse
    and far too large to densify, so it is done here instead.

    Parameters
    ----------
    counts : anndata.AnnData
        One modality. Cells in rows, genes in columns.

    Returns
    -------
    numpy.ndarray
        One value per gene. A cell with no UMIs contributes nothing rather
        than dividing by zero.
    """
    import scipy.sparse as sp

    matrix = counts.X
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    # A cell with no UMIs has no composition to measure; scaling it by zero
    # keeps it in the denominator without inventing counts for it.
    scale = np.divide(
        1e6, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0
    )
    if sp.issparse(matrix):
        scaled = sp.diags(scale) @ matrix
    else:
        scaled = np.asarray(matrix) * scale[:, None]
    return np.asarray(scaled.mean(axis=0)).ravel()


def panel_expression(data, modality: str = "gene", flag: str = "readout"):
    """How much more a targeted panel's genes are seen than everything else.

    A targeted assay works when its panel dominates. This compares the two
    groups gene by gene rather than cell by cell, which is what says whether
    the enrichment comes from the panel as a whole or from a few genes.

    Parameters
    ----------
    data
        A MuData as `mudata_builder` assembles it.
    modality : str
        Which modality holds the genes.
    flag : str
        The boolean column of `var` marking the panel.

    Returns
    -------
    tuple of (matplotlib.figure.Figure, pandas.DataFrame, str)
        The figure, one row per gene, and the caption.

    Raises
    ------
    FigureError
        `var` has no such column, or it marks no gene.
    """
    plt = style()

    counts = data[modality]
    if flag not in counts.var:
        raise FigureError(f"{flag!r} is not a column of {modality}.var.")
    in_panel = counts.var[flag].fillna(False).to_numpy().astype(bool)
    if not in_panel.any():
        raise FigureError(f"{flag!r} marks no gene, so there is no panel.")

    values = mean_upm(counts)
    numbers = pd.DataFrame(
        {
            "gene": counts.var_names,
            "symbol": counts.var.get("symbol", pd.Series(index=counts.var_names)),
            "mean_upm": values,
            "in_panel": in_panel,
        }
    ).reset_index(drop=True)

    figure, axis = plt.subplots(figsize=(4.0, 4.5))
    groups = [values[in_panel], values[~in_panel]]
    # +1 so a gene seen in no cell stays on a log axis instead of vanishing.
    boxes = axis.boxplot(
        [group + 1 for group in groups],
        tick_labels=["panel", "background"],
        widths=0.5,
        patch_artist=True,
        flierprops={"markersize": 2, "markeredgewidth": 0, "alpha": 0.3},
    )
    for patch, colour in zip(boxes["boxes"], (KEPT, REMOVED), strict=False):
        patch.set_facecolor(colour)
    for median in boxes["medians"]:
        median.set_color("black")
    axis.set_yscale("log")
    axis.set_ylabel("mean UMIs per million + 1")
    axis.grid(axis="y", which="major", color="0.85", linewidth=0.6)
    axis.set_axisbelow(True)
    figure.tight_layout()

    panel_median = float(np.median(values[in_panel]))
    background_median = float(np.median(values[~in_panel]))
    seen = int((values[~in_panel] > 0).sum())
    difference = (
        f"a {panel_median / background_median:,.0f}-fold difference"
        if background_median
        # Half the reference is never seen, so its median is zero and a
        # ratio says nothing. The share that is seen does.
        else f"the background median being zero, with {seen:,} of "
        f"{int((~in_panel).sum()):,} background genes seen at all"
    )
    caption = (
        f"Mean expression per gene, in units of UMIs per million, for the "
        f"{int(in_panel.sum()):,} genes of the panel against the "
        f"{int((~in_panel).sum()):,} others. Medians are "
        f"{panel_median:,.1f} and {background_median:,.1f}, {difference}. "
        "One is added before taking the log so a gene seen in no cell still "
        "appears. Comparing genes rather than cells is what separates "
        "enrichment across the whole panel from a few genes carrying it."
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


def cell_fractions(
    data,
    modality: str = "gene",
    extra: dict[str, np.ndarray] | None = None,
) -> pd.DataFrame:
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
    extra
        Further classes, as a mapping of name to a boolean mask over the
        modality's genes. A class defined elsewhere -- an older analysis's
        own definition of ribosomal, say -- is measured alongside the
        built-in ones, so a threshold set against it still applies.

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
    for name, mask in (extra or {}).items():
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (counts.n_vars,):
            raise FigureError(
                f"The mask for {name!r} covers {mask.shape[0]:,} genes, "
                f"but the modality holds {counts.n_vars:,}."
            )
        classes[name] = mask

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

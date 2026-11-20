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


class FilterError(RuntimeError):
    """The filter cannot be applied to what was given."""


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
    # Appended, not with_suffix: that would cut "knee_v1.2" back to "knee_v1".
    written = {
        suffix: stem.with_name(stem.name + suffix)
        for suffix in (".pdf", ".tsv.gz", ".txt")
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
    held = data[modality]
    sums = pd.Series(np.asarray(held.X.sum(axis=1)).ravel(), index=held.obs_names)
    # By name: a modality need not hold its cells in the object's order.
    return sums.reindex(data.obs_names, fill_value=0).to_numpy(dtype=float)


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
    filtered = data[keep].copy()
    # The step that cuts writes down what it cut, so the record cannot be
    # left to whoever writes the file.
    record = dict(filtered.uns.get("qc", {}))
    record["floor"] = {
        "min_umis": int(min_umis),
        "modality": modality,
        "n_cells_in": int(data.n_obs),
        "n_cells_out": int(filtered.n_obs),
    }
    filtered.uns["qc"] = record
    return filtered, keep


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


KNEE_POINTS = 400
"""Ranks a stored knee curve keeps, log-spaced: enough to redraw it exactly."""

GROUP_COLOURS = (OKABE_ITO[4], OKABE_ITO[0], OKABE_ITO[2], OKABE_ITO[6])
"""Fixed slots for up to four groups. No yellow: a thin line on white loses it."""


def knee_curve(umis, n_points: int = KNEE_POINTS) -> pd.DataFrame:
    """A channel's barcode-rank curve at log-spaced ranks.

    A million barcodes are too many to keep per channel and too many to draw
    110 times; a few hundred ranks spaced evenly on the log axis draw the same
    curve. Each kept point is the exact UMI count at that rank: the curve is
    sampled, not smoothed.

    Parameters
    ----------
    umis : array-like
        Total UMIs per barcode, as `umis_per_cell` returns.
    n_points : int
        How many ranks to keep, at most.

    Returns
    -------
    pandas.DataFrame
        ``rank`` (1 is the deepest barcode) and ``umis``, over the barcodes
        with at least one UMI. Empty when there are none.
    """
    values = np.asarray(umis, dtype=float)
    ranked = np.sort(values[values > 0])[::-1]
    if not len(ranked):
        return pd.DataFrame(
            {"rank": pd.Series(dtype=int), "umis": pd.Series(dtype=int)}
        )
    ranks = np.unique(np.geomspace(1, len(ranked), n_points).round().astype(int))
    return pd.DataFrame({"rank": ranks, "umis": ranked[ranks - 1].astype(int)})


def knee_overlay(
    curves: pd.DataFrame,
    group: str | None = None,
    floor: float | None = None,
    ylabel: str = "gene UMIs per barcode",
):
    """Many channels' knees on one panel, and each group's mean beside it.

    Small multiples of a hundred knees are a hundred figures to read; one
    panel shows how they spread, and the second says where a group sits on
    average and how much its channels vary.

    Parameters
    ----------
    curves : pandas.DataFrame
        One row per point: ``channel``, ``rank``, ``umis``, and the `group`
        column, as `knee_curve` gives per channel.
    group : str, optional
        Column to colour by, such as ``day``. At most four groups.
    floor : float, optional
        A UMI floor, drawn as a dashed line on both panels.
    ylabel : str
        What the UMIs are, for the axis.

    Returns
    -------
    tuple
        The figure, the numbers (every channel's points, then each group's
        mean and spread by rank), and the caption, ready for `publish`.

    Raises
    ------
    FigureError
        No curve, a column missing, or more than four groups.

    Notes
    -----
    The mean is of log10 UMIs, so it is a geometric mean, and the band is one
    standard deviation of log10 UMIs either side. A group's mean is drawn
    only at ranks every one of its channels reaches, so the line never
    averages over a channel that has run out of barcodes.
    """
    plt = style()
    needed = {"channel", "rank", "umis"} | ({group} if group else set())
    missing = needed - set(curves.columns)
    if missing:
        raise FigureError(f"The curves lack {', '.join(sorted(missing))}.")
    curves = curves[curves["umis"] > 0]
    if curves.empty:
        raise FigureError("No curve has a barcode with a UMI.")

    labels = (
        curves[group].astype(str) if group else pd.Series("all", index=curves.index)
    )
    if group and isinstance(curves[group].dtype, pd.CategoricalDtype):
        names = [str(c) for c in curves[group].cat.categories if str(c) in set(labels)]
    else:
        names = sorted(pd.unique(labels), key=_natural)
    if len(names) > len(GROUP_COLOURS):
        raise FigureError(
            f"{len(names)} groups; at most {len(GROUP_COLOURS)} colours stay"
            " distinguishable. Facet, or fold the rest into one group."
        )
    colour = dict(zip(names, GROUP_COLOURS, strict=False))
    # One group needs no colour: grey curves and a dark mean say it all.
    line = (
        {name: ("0.55", 0.35) for name in names}
        if not group
        else {name: (colour[name], 0.45) for name in names}
    )
    if not group:
        colour = {name: "0.15" for name in names}

    figure, (left, right) = plt.subplots(
        1, 2, figsize=(7.2, 3.4), sharey=True, gridspec_kw={"wspace": 0.08}
    )
    figure.subplots_adjust(left=0.1, right=0.98, bottom=0.15, top=0.9)
    channels = curves.assign(_group=labels).groupby("channel", sort=False)
    n_channels = {name: 0 for name in names}
    for _, part in channels:
        name = part["_group"].iloc[0]
        n_channels[name] += 1
        left.plot(
            part["rank"],
            part["umis"],
            color=line[name][0],
            linewidth=0.5,
            alpha=line[name][1],
        )

    top = curves["rank"].max()
    grid = np.unique(np.geomspace(1, top, 200).round().astype(int))
    summaries = []
    for name in names:
        logs = []
        for _, part in channels:
            if part["_group"].iloc[0] != name:
                continue
            reach = grid <= part["rank"].max()
            values = np.full(len(grid), np.nan)
            values[reach] = np.interp(
                np.log10(grid[reach]),
                np.log10(part["rank"].to_numpy(dtype=float)),
                np.log10(part["umis"].to_numpy(dtype=float)),
            )
            logs.append(values)
        stack = np.vstack(logs)
        every = ~np.isnan(stack).any(axis=0)
        mean = stack[:, every].mean(axis=0)
        sd = stack[:, every].std(axis=0, ddof=1) if len(stack) > 1 else 0 * mean
        ranks = grid[every]
        right.fill_between(
            ranks,
            10 ** (mean - sd),
            10 ** (mean + sd),
            color=colour[name] if group else "0.6",
            alpha=0.25 if group else 0.4,
            linewidth=0,
        )
        right.plot(ranks, 10**mean, color=colour[name], linewidth=1.5, label=name)
        summaries.append(
            pd.DataFrame(
                {
                    "panel": "mean",
                    "group": name,
                    "rank": ranks,
                    "mean_log10_umis": mean,
                    "sd_log10_umis": sd,
                    "n_channels": len(stack),
                }
            )
        )

    for axis in (left, right):
        axis.set(xscale="log", yscale="log", xlabel="barcode rank")
        axis.grid(axis="x", visible=False)
        if floor is not None:
            axis.axhline(floor, color="0.3", linewidth=0.8, linestyle="--")
    left.set_ylabel(ylabel)
    if floor is not None:
        left.annotate(
            f"floor {floor:,g}",
            xy=(1, floor),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=8,
            color="0.3",
        )
    if group:
        right.legend(frameon=False, fontsize=8)
    for letter, axis in (("a", left), ("b", right)):
        axis.text(
            -0.02 if axis is right else -0.16,
            1.04,
            letter,
            transform=axis.transAxes,
            fontweight="bold",
            fontsize=11,
            va="bottom",
        )

    points = curves.assign(panel="channels", group=labels)[
        ["panel", "group", "channel", "rank", "umis"]
    ]
    numbers = pd.concat([points, *summaries], ignore_index=True)
    counted = ", ".join(f"{n_channels[name]} in {name}" for name in names)
    by = f" by {group}" if group else ""
    total = sum(n_channels.values())
    caption = (
        f"(a) Barcode-rank curves of {total} channel{'s' * (total != 1)}"
        f"{f' ({counted}), coloured{by}' if group else ''}, both axes logarithmic; each curve is the"
        " channel's barcodes with at least one UMI. (b) The mean of log10 UMIs"
        f" at each rank{by}, with one standard deviation either side shaded,"
        " drawn only at ranks every channel of the group reaches."
        + (
            f" The dashed line is the floor of {floor:,g} UMIs."
            if floor is not None
            else ""
        )
    )
    return figure, numbers, caption


def _natural(value: str):
    """Sort "2" before "10": a number as a number, anything else as text."""
    try:
        return (0, float(value), "")
    except ValueError:
        return (1, 0.0, value)


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


def _label(value: float) -> str:
    """A number short enough to sit on a threshold line."""
    return f"{value:,.0f}" if abs(value) >= 100 else f"{value:.3g}"


KEPT = OKABE_ITO[4]
"""Blue, for what a threshold keeps: the cells the analysis goes on with."""

REMOVED = "0.7"
"""Grey, for what it removes. Neutral, so the eye reads the kept cells first."""

EDGE = "0.2"
"""Near-black, for the outline a light fill needs to be seen on a white page."""

EDGE_WIDTH = 0.4
"""Points, thin enough that it outlines the fill rather than replacing it."""

LOG_METRICS = ("total_counts", "n_genes_by_counts")
"""Metrics whose distribution is read on a log axis, because they are counts."""


def _band(values, metric, nmads, log, sides):
    """One metric's MAD band, `(-inf, inf)` when it has none.

    The side rule and the NaN drop live here rather than in each caller,
    so the panel and the filter cut at the same place.
    """
    on_log = metric in log
    side = sides.get(metric, "both" if on_log else "upper")
    measured = values[np.isfinite(values)]
    if measured.size == 0:
        return -np.inf, np.inf
    try:
        return mad_bounds(measured, nmads=nmads, side=side, log=on_log)
    except FigureError:
        # No spread, so no threshold. That is a metric to report, not a
        # reason to take the figure or the filter down.
        return -np.inf, np.inf


def metric_distributions(
    frame: pd.DataFrame,
    nmads: float | None = 5.0,
    limits: dict[str, tuple[float | None, float | None]] | None = None,
    band: tuple[str, ...] | None = None,
    log: tuple[str, ...] = LOG_METRICS,
    sides: dict[str, str] | None = None,
    bins: int = 60,
):
    """One histogram per QC metric, split at the rules it is cut on.

    The rules come from `cell_filters`, so pass this the same `nmads`,
    `limits` and `band` as the filter: the lines on the figure are then
    the cuts the file was made with, not a second opinion.

    Parameters
    ----------
    frame : pandas.DataFrame
        One column per metric, one row per cell, as `cell_fractions`
        returns. Rows that are entirely `NaN` are cells the modality does
        not hold, and are dropped rather than drawn as zero.
    nmads : float, optional
        Width of each band, in MADs. `None` draws only `limits`.
    limits : dict, optional
        `{metric: (low, high)}`, as `cell_filters` takes them.
    band : tuple of str, optional
        The metrics to measure a band for, every column by default.
    log : tuple of str
        Metrics binned, thresholded and drawn on a log axis. A count
        distribution spans orders of magnitude; a percentage does not.
    sides : dict, optional
        Which tail each band cuts. The default is `"both"` for a count and
        `"upper"` for everything else, since a cell with unusually little
        mitochondrial signal is not a problem.
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
    from matplotlib import ticker

    plt = style()

    if frame.shape[1] == 0:
        raise FigureError("No metric to draw.")
    measured = frame.dropna(how="all")
    if measured.empty:
        raise FigureError("No measured cell to draw.")

    metrics = list(measured.columns)
    # The same call the filter makes, so a panel cannot draw a line the
    # file was not cut on.
    passes, rules = cell_filters(
        measured, nmads=nmads, limits=limits, band=band, log=log, sides=sides
    )
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
        mine = rules[rules["metric"] == metric]
        if mine.empty:
            unthresholded.append(metric)
            within = np.ones(len(values), dtype=bool)
        else:
            within = passes[list(mine["rule"])].to_numpy(dtype=bool).all(axis=1)
        # What the panel cut on, which is the tightest of its rules.
        low = float(mine["low"].max()) if not mine.empty else -np.inf
        high = float(mine["high"].min()) if not mine.empty else np.inf

        if on_log:
            # A zero has no place on a log axis, so the range starts at the
            # smallest value actually present.
            positive = values[values > 0]
            floor = float(positive.min()) if positive.size else 1.0
            edges = np.logspace(np.log10(floor), np.log10(np.nanmax(values)), bins + 1)
            axis.set_xscale("log")
            # Over a span this narrow matplotlib labels minor ticks too, and
            # four superscripts in an inch run into each other. Label 1, 2
            # and 5 per decade instead, in the digits the caption uses.
            decades = np.log10(edges[-1]) - np.log10(edges[0])
            axis.xaxis.set_major_locator(
                ticker.LogLocator(subs=(1.0,) if decades > 3 else (1.0, 2.0, 5.0))
            )
            axis.xaxis.set_major_formatter(
                ticker.FuncFormatter(lambda value, _: _label(value))
            )
            axis.xaxis.set_minor_formatter(ticker.NullFormatter())
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
        edges_of = [
            (bound, upper, rule.source)
            for rule in mine.itertuples()
            for bound, upper in ((rule.low, False), (rule.high, True))
        ]
        for threshold, upper, source in edges_of:
            if not np.isfinite(threshold):
                continue
            axis.axvline(threshold, color="black", linewidth=1.0, linestyle="--")
            # A line nobody can read is a line nobody can check. The label
            # sits inside the band, so it never leaves the panel.
            axis.annotate(
                f"{source}  {_label(threshold)}",
                xy=(threshold, 1.0),
                xycoords=axis.get_xaxis_transform(),
                xytext=(-3 if upper else 3, -3),
                textcoords="offset points",
                rotation=90,
                ha="right" if upper else "left",
                va="top",
                fontsize="x-small",
                clip_on=True,
            )
        # The panel's own result, since a cut this small leaves no grey
        # tall enough to see: with no count, one colour says nothing.
        axis.set_title(
            f"{int(cut_counts.sum()):,} of {int(kept_counts.sum() + cut_counts.sum()):,} beyond"
            if metric not in unthresholded
            else "no rule",
            loc="right",
            fontsize="small",
            color="0.35",
        )
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
    figure.tight_layout()

    numbers = pd.concat(records, ignore_index=True)
    thresholds = ", ".join(
        f"{rule.metric} {rule.source} "
        + " to ".join(
            f"{bound:,.2f}" for bound in (rule.low, rule.high) if np.isfinite(bound)
        )
        for rule in rules.itertuples()
    )
    removed = int(numbers["n_removed"].sum())
    caption = (
        f"Distribution of each quality-control metric over {len(measured):,} "
        "cells. Counts are binned and drawn on a log axis, fractions on a "
        "linear one. The dashed lines are the rules each panel is cut on, "
        "a band measured from this channel or a limit given to it "
        f"({thresholds}); cells beyond them are grey, and each panel's "
        "title counts its own. Each panel cuts on its own metric, so those "
        f"counts add to {removed:,}, which double-counts a cell that "
        "several metrics reject."
    )
    if unthresholded:
        caption += (
            " No threshold is drawn for "
            + ", ".join(unthresholded)
            + ", which no rule covers."
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
        layers = [
            (kept & drawn, KEPT, "kept", alpha, 0.0),
            (~kept & drawn, REMOVED, "removed", alpha, 0.0),
        ]
    else:
        labels = pd.Series(groups, index=frame.index).reindex(measured.index)
        # Deterministic, so a label keeps its colour across figures. A
        # Categorical says the order it wants; anything else is sorted,
        # since row order is not a meaning.
        named = (
            [str(name) for name in labels.cat.categories]
            if isinstance(labels.dtype, pd.CategoricalDtype)
            else sorted(str(name) for name in pd.unique(labels.dropna()))
        )
        # Unlabelled first and see-through, so the classes being compared
        # sit on top of it. Those are drawn solid: a handful of points at
        # the alpha that keeps a cloud of ten thousand readable is barely
        # visible.
        # The outline is what makes a pale fill findable: Okabe-Ito's yellow
        # on a white page has almost no edge of its own.
        layers = [(labels.isna().to_numpy() & drawn, REMOVED, None, alpha, 0.0)]
        layers += [
            (
                (labels == name).to_numpy() & drawn,
                OKABE_ITO[i % len(OKABE_ITO)],
                name,
                1.0,
                EDGE_WIDTH,
            )
            for i, name in enumerate(named)
        ]

    # Rasterised: one point per cell stays a vector object otherwise, and a
    # channel of 10,000 cells turns into a PDF nobody can open.
    for mask, colour, label, opacity, edge in layers:
        if mask.any():
            axis.scatter(
                xs[mask],
                ys[mask],
                s=size,
                alpha=opacity,
                linewidths=edge,
                edgecolors=EDGE,
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


def metric_density(
    frame: pd.DataFrame,
    x: str,
    y: str,
    x_threshold: float | None = None,
    y_threshold: float | None = None,
    bins: int = 120,
    xlabel: str | None = None,
    ylabel: str | None = None,
):
    """Two count metrics against each other as a density, both axes log.

    `metric_scatter` draws a point per cell, which at a hundred thousand
    cells is a solid block that hides where the cells are. This counts the
    cells in each bin and colours the bin, so the populations and the
    valleys between them show.

    Parameters
    ----------
    frame : pandas.DataFrame
        One row per cell, with columns `x` and `y`. Rows with ``NaN`` in
        either are dropped.
    x, y : str
        The two metrics, both counts.
    x_threshold, y_threshold : float, optional
        Cuts, drawn as dashed lines; the caption counts the cells in each
        quadrant they make.
    bins : int
        Bins per axis, evenly spaced on the log scale. On an axis whose
        values are all whole numbers the edges move to the nearest half
        integer, so a low count never falls in a bin no integer can reach.
    xlabel, ylabel : str, optional
        Axis labels; the column names by default.

    Returns
    -------
    tuple
        The figure, the numbers (one row per non-empty bin: its edges and
        the cells in it), and the caption, ready for `publish`.

    Raises
    ------
    FigureError
        A column is missing, or no cell has both metrics above zero.

    Notes
    -----
    A zero cannot sit on a log axis. Such cells are left out of the bins and
    the margins, counted in the caption, and still counted against the cuts.

    Moving the edges to half integers widens the bins at low counts, so each
    bin's count, and each margin bar, is divided by its size relative to an
    even bin. That keeps the colour a density rather than a bin width.
    """
    import matplotlib.colors as mcolors

    plt = style()
    for name in (x, y):
        if name not in frame:
            raise FigureError(f"{name!r} is not a column of the frame.")
    measured = frame[[x, y]].dropna()
    xs = measured[x].to_numpy(dtype=float)
    ys = measured[y].to_numpy(dtype=float)
    drawn = (xs > 0) & (ys > 0)
    if not drawn.any():
        raise FigureError("No cell has both metrics above zero.")

    lx, ly = np.log10(xs[drawn]), np.log10(ys[drawn])
    x_edges, x_size = _density_edges(xs[drawn], lx, bins)
    y_edges, y_size = _density_edges(ys[drawn], ly, bins)
    counts, _, _ = np.histogram2d(lx, ly, bins=(x_edges, y_edges))
    scaled = counts / np.outer(x_size, y_size)

    figure = plt.figure(figsize=(6.0, 6.0))
    grid = figure.add_gridspec(
        2, 3, width_ratios=(5, 1, 0.18), height_ratios=(1, 5), wspace=0.05, hspace=0.05
    )
    axis = figure.add_subplot(grid[1, 0])
    top = figure.add_subplot(grid[0, 0], sharex=axis)
    right = figure.add_subplot(grid[1, 1], sharey=axis)
    bar_axis = figure.add_subplot(grid[1, 2])

    shown = np.ma.masked_where(counts.T == 0, scaled.T)
    mesh = axis.pcolormesh(
        10**x_edges,
        10**y_edges,
        shown,
        cmap="cividis",
        norm=mcolors.LogNorm(vmin=min(1, shown.min()), vmax=max(2, shown.max())),
        rasterized=True,
    )
    colourbar = figure.colorbar(mesh, cax=bar_axis)
    colourbar.set_label("cells per bin")
    axis.set(xscale="log", yscale="log", xlabel=xlabel or x, ylabel=ylabel or y)
    axis.grid(False)
    top.stairs(counts.sum(axis=1) / x_size, 10**x_edges, fill=True, color="0.55")
    right.stairs(
        counts.sum(axis=0) / y_size,
        10**y_edges,
        fill=True,
        color="0.55",
        orientation="horizontal",
    )
    for margin in (top, right):
        margin.grid(False)
        margin.tick_params(labelbottom=False, labelleft=False)
    top.spines["left"].set_visible(False)
    top.set_yticks([])
    right.spines["bottom"].set_visible(False)
    right.set_xticks([])

    above_x = xs >= x_threshold if x_threshold is not None else np.ones(len(xs), bool)
    above_y = ys >= y_threshold if y_threshold is not None else np.ones(len(ys), bool)
    for threshold, line in ((x_threshold, axis.axvline), (y_threshold, axis.axhline)):
        if threshold is not None:
            line(threshold, color="white", linewidth=1.0, linestyle="--")
    if x_threshold is not None:
        top.axvline(x_threshold, color="0.3", linewidth=0.8, linestyle="--")
    if y_threshold is not None:
        right.axhline(y_threshold, color="0.3", linewidth=0.8, linestyle="--")

    i, j = np.nonzero(counts)
    numbers = pd.DataFrame(
        {
            f"{x}_low": 10 ** x_edges[i],
            f"{x}_high": 10 ** x_edges[i + 1],
            f"{y}_low": 10 ** y_edges[j],
            f"{y}_high": 10 ** y_edges[j + 1],
            "cells": counts[i, j].astype(int),
            "cells_per_even_bin": scaled[i, j],
        }
    )
    quadrants = {
        "both above": int((above_x & above_y).sum()),
        f"{x} below only": int((~above_x & above_y).sum()),
        f"{y} below only": int((above_x & ~above_y).sum()),
        "both below": int((~above_x & ~above_y).sum()),
    }
    cut = []
    if x_threshold is not None:
        cut.append(f"{x} {x_threshold:,g}")
    if y_threshold is not None:
        cut.append(f"{y} {y_threshold:,g}")
    caption = (
        f"{y} against {x} for {len(measured):,} cells, both axes logarithmic;"
        f" colour is the number of cells per bin ({bins} bins per axis), and the"
        " margins are each metric's own distribution."
    )
    if (x_size > 1).any() or (y_size > 1).any():
        caption += (
            " At low counts a bin spans whole numbers, so it is wider than the rest and"
            " its count is divided by its size relative to an even bin."
        )
    if len(cut) == 2:
        caption += (
            f" Dashed lines are the cuts ({', '.join(cut)}). Cells per quadrant: "
            + ", ".join(f"{name} {n:,}" for name, n in quadrants.items())
            + "."
        )
    elif cut:
        above = above_x if x_threshold is not None else above_y
        caption += (
            f" The dashed line is the cut ({cut[0]}): {int(above.sum()):,} cells at or"
            f" above it, {int((~above).sum()):,} below."
        )
    zeros = int((~drawn).sum())
    if zeros:
        caption += (
            f" {zeros:,} {'cell holds' if zeros == 1 else 'cells hold'} a zero, which a log axis cannot"
            " show; they are left out of the bins but counted against the cuts."
        )
    return figure, numbers, caption


def _density_edges(values: np.ndarray, logs: np.ndarray, bins: int):
    """Log-spaced bin edges for `metric_density`, and each bin's size.

    For whole-number values the edges move to the nearest half integer, and
    duplicates merge. The size is each bin's log width over an even bin's.
    """
    even = np.linspace(logs.min(), logs.max(), bins + 1)
    if not np.array_equal(values, np.round(values)):
        return even, np.ones(bins)
    linear = np.unique(np.append(np.round(10**even) - 0.5, values.max() + 0.5))
    edges = np.log10(linear)
    width = (even[-1] - even[0]) / bins if even[-1] > even[0] else 1.0
    return edges, np.diff(edges) / width


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


def cell_filters(
    frame: pd.DataFrame,
    nmads: float | None = 5.0,
    limits: dict[str, tuple[float | None, float | None]] | None = None,
    band: tuple[str, ...] | None = None,
    log: tuple[str, ...] = LOG_METRICS,
    sides: dict[str, str] | None = None,
):
    """Which cells pass each quality-control rule.

    Two kinds of rule: a band this channel's own spread decides, and a
    limit someone else did. A metric can carry both, so a threshold the
    screen settled on elsewhere sits next to the band measured here.

    Parameters
    ----------
    frame : pandas.DataFrame
        One column per metric, one row per cell, as `cell_fractions`
        returns, measured on cells already past the count floor. The MAD
        describes whatever population it is given, and on raw barcodes
        that population is ambient.
    nmads : float, optional
        Width of each band, in MADs. `None` measures no band, leaving only
        `limits`.
    limits : dict, optional
        `{metric: (low, high)}` in the metric's own units, `None` for an
        open side. A threshold from outside this channel goes here.
    band : tuple of str, optional
        The metrics to measure a band for. Every column by default; name a
        subset when only some of them are worth cutting on.
    log : tuple of str
        Metrics whose band is measured on a log scale.
    sides : dict, optional
        Which tail each band cuts. The default is `"both"` for a count and
        `"upper"` for everything else.

    Returns
    -------
    tuple of (pandas.DataFrame, pandas.DataFrame)
        A boolean frame, one row per cell and one column per rule, `True`
        where the cell passes it; and the rules themselves, one row each,
        carrying the metric, where the bounds came from and how many cells
        the rule rejected.

    Raises
    ------
    FilterError
        A limit or a band names a metric the frame does not hold.

    Notes
    -----
    A cell whose metric is `NaN` fails that rule: the modality did not
    measure it, and a metric a cell was never measured for is not one it
    passed.
    """
    sides = sides or {}
    limits = limits or {}
    unknown = [
        metric
        for metric in list(limits) + list(band or ())
        if metric not in frame.columns
    ]
    if unknown:
        raise FilterError(
            "No such metric: "
            + ", ".join(unknown)
            + ". The frame holds "
            + ", ".join(frame.columns)
            + "."
        )

    passes: dict[str, np.ndarray] = {}
    rules = []
    for metric in frame.columns:
        values = frame[metric].to_numpy(dtype=float)
        bounds = []
        if nmads is not None and (band is None or metric in band):
            low, high = _band(values, metric, nmads, log, sides)
            if (low, high) != (-np.inf, np.inf):
                bounds.append((f"{metric}_mad", f"{nmads:g} MAD", low, high))
        if metric in limits:
            low, high = limits[metric]
            bounds.append(
                (
                    f"{metric}_limit",
                    "limit",
                    -np.inf if low is None else float(low),
                    np.inf if high is None else float(high),
                )
            )
        for name, source, low, high in bounds:
            within = (values >= low) & (values <= high)
            passes[name] = within
            rules.append(
                {
                    "rule": name,
                    "metric": metric,
                    "source": source,
                    "low": low,
                    "high": high,
                    "n_failed": int((~within).sum()),
                }
            )

    return (
        pd.DataFrame(passes, index=frame.index),
        pd.DataFrame(
            rules, columns=["rule", "metric", "source", "low", "high", "n_failed"]
        ),
    )


def apply_filters(data, passes: pd.DataFrame, rules: pd.DataFrame, path=None):
    """Keep the cells that pass every rule, and record the rules in `uns`.

    Parameters
    ----------
    data : mudata.MuData
        The object the rules were measured on.
    passes, rules : pandas.DataFrame
        What `cell_filters` returned.
    path : str or pathlib.Path, optional
        Where to write the filtered object. Nothing is written without one.

    Returns
    -------
    mudata.MuData
        The filtered copy, with every threshold that decided it under
        `uns["qc"]`, so the file says how it was made.

    Raises
    ------
    FilterError
        There is no rule, the rules do not cover exactly this object's
        cells, or no cell passes them all.
    """
    if list(passes.index) != list(data.obs_names):
        raise FilterError(
            f"The rules cover {len(passes.index):,} cells and the object "
            f"holds {data.n_obs:,}. Measure them on the object being "
            "filtered: a cell no rule saw would pass every one of them."
        )
    if passes.shape[1] == 0:
        raise FilterError(
            "No rule to apply: no metric has the spread to measure a band, "
            "and no limit was given. Filtering on nothing would keep every "
            "cell and say it passed."
        )
    keep = passes.to_numpy(dtype=bool).all(axis=1)
    if not keep.any():
        raise FilterError("No cell passes every rule.")

    filtered = data[keep].copy()
    record = dict(filtered.uns.get("qc", {}))
    record["n_cells_in"] = int(data.n_obs)
    record["n_cells_out"] = int(filtered.n_obs)
    record["rules"] = {
        row.rule: {
            "metric": row.metric,
            "source": row.source,
            "low": float(row.low),
            "high": float(row.high),
            "n_failed": int(row.n_failed),
        }
        for row in rules.itertuples()
    }
    filtered.uns["qc"] = record

    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        filtered.write_h5mu(path)
    return filtered


def guide_counts(data, modality: str = "guide") -> pd.DataFrame:
    """Every guide's UMIs and the number of cells it was seen in.

    Parameters
    ----------
    data
        A MuData as `mudata_builder` assembles it.
    modality
        The modality holding the guide counts.

    Returns
    -------
    pandas.DataFrame
        One row per guide, indexed by guide id, with `umis` and `n_cells`.
        `attrs` carries the modality and how many cells it was summed over,
        for a caption to quote.

    Notes
    -----
    Both numbers are summed over the cells the modality holds, which under a
    union is not every cell of the object.
    """
    counts = data[modality]
    matrix = counts.X
    frame = pd.DataFrame(
        {
            "umis": np.asarray(matrix.sum(axis=0)).ravel(),
            "n_cells": np.asarray((matrix > 0).sum(axis=0)).ravel(),
        },
        index=pd.Index(counts.var_names, name="guide_id"),
    )
    frame.attrs = {"modality": modality, "n_cells": int(counts.n_obs)}
    return frame


#: What each `guide_counts` column is called on an axis.
GUIDE_LABELS = {"umis": "UMIs per guide", "n_cells": "cells per guide"}


def guide_abundance(
    counts: pd.DataFrame,
    column: str = "umis",
    groups=None,
    bins: int = 60,
):
    """How much of the library this channel saw, and how often.

    The guides never seen get a bar of their own left of the axis. A zero
    has no place on a log scale, and dropping it would answer a different
    question than the one the figure is asked.

    Parameters
    ----------
    counts : pandas.DataFrame
        One row per guide, as `guide_counts` returns.
    column : str
        Which of its columns to bin.
    groups : optional
        A label per guide, aligned to `counts.index`. Each class is drawn
        as its own outline, over the share of that class rather than the
        count: a library's classes differ in size by two orders of
        magnitude, and raw counts would show only the largest.
    bins : int
        Logarithmic bins between one and the largest count, rounded to
        integer edges so no bin is narrower than a count.

    Returns
    -------
    tuple
        The figure, the numbers, and the caption, ready for `publish`.

    Raises
    ------
    FigureError
        No guide was seen at all, so there is nothing to bin.
    """
    from matplotlib import ticker

    plt = style()
    if column not in counts.columns:
        raise FigureError(f"The frame has no {column!r} to bin.")
    values = counts[column].to_numpy(dtype=float)
    seen = values[values > 0]
    if seen.size == 0:
        raise FigureError(f"No guide has a {column} above zero to bin.")

    labels = pd.Series("every guide" if groups is None else groups, index=counts.index)
    named = (
        [str(name) for name in labels.cat.categories]
        if isinstance(labels.dtype, pd.CategoricalDtype)
        else sorted(str(name) for name in pd.unique(labels.dropna()))
    )
    labels = labels.astype(str)

    # A count is an integer, so a log bin narrower than one alternates
    # full and empty and draws a comb. Rounding the edges collapses those
    # into single-count bins and leaves the rest logarithmic.
    edges = np.unique(np.round(np.logspace(0.0, np.log10(seen.max()), bins + 1)))
    edges[-1] = max(edges[-1], float(seen.max()))

    figure, axis = plt.subplots(figsize=(6.5, 4.0))
    # The never-seen bars share the linear strip the symlog axis keeps for
    # zero, side by side so one class does not hide another.
    span = 0.8 / len(named)
    records = []
    for index, name in enumerate(named):
        mine = values[(labels == name).to_numpy()]
        if mine.size == 0:
            continue
        colour = KEPT if groups is None else OKABE_ITO[index % len(OKABE_ITO)]
        scale = float(mine.size) if groups is not None else 1.0
        # Counted as integers, then scaled: a count scaled and back is 14.999.
        n_guides = np.concatenate(
            [[(mine == 0).sum()], np.histogram(mine[mine > 0], bins=edges)[0]]
        )
        never, binned = n_guides[0] / scale, n_guides[1:] / scale
        if groups is None:
            axis.stairs(binned, edges, fill=True, color=colour)
        else:
            axis.stairs(binned, edges, color=colour, linewidth=1.4, label=name)
        axis.bar(
            0.1 + span * (index + 0.5),
            never,
            width=span * 0.9,
            color=REMOVED if groups is None else colour,
            edgecolor=EDGE,
            linewidth=EDGE_WIDTH,
        )
        records.append(
            pd.DataFrame(
                {
                    "class": name,
                    "panel": ["never seen"] + ["binned"] * (len(edges) - 1),
                    "bin_left": np.concatenate([[np.nan], edges[:-1]]),
                    "bin_right": np.concatenate([[np.nan], edges[1:]]),
                    "n_guides": n_guides,
                    "share": n_guides / mine.size,
                }
            )
        )

    axis.set_xscale("symlog", linthresh=1.0, linscale=0.35)
    axis.set_yscale("log")
    axis.set(
        xlabel=GUIDE_LABELS.get(column, column.replace("_", " ")),
        ylabel="guides" if groups is None else "share of the class",
    )
    axis.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: _label(v)))
    axis.xaxis.set_minor_formatter(ticker.NullFormatter())
    if groups is not None:
        axis.legend(frameon=False, fontsize="small")

    never = int((values == 0).sum())
    axis.set_title(
        f"{never:,} of {len(values):,} never seen",
        loc="right",
        fontsize="small",
        color="0.35",
    )
    figure.tight_layout()

    numbers = pd.concat(records, ignore_index=True)
    over = counts.attrs.get("n_cells")
    what = GUIDE_LABELS.get(column, column.replace("_", " "))
    caption = (
        f"{what[0].upper()}{what[1:]}, binned logarithmically"
        + (f" over {over:,} cells" if over else "")
        + f". The bar left of the axis is the {never:,} of {len(values):,} "
        "guides never seen, which a log scale cannot place at zero. A "
        "channel holds far fewer cells than a library holds guides, so a "
        "guide missing here is not a guide missing from the library: that "
        "question needs the channels pooled."
    )
    return figure, numbers, caption

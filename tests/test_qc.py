"""QC: what the classes mean, what a floor keeps, and what a figure writes."""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pytest

from crispr_analysis_utils.qc import (
    FigureError,
    cell_fractions,
    filter_cells,
    gene_classes,
    knee,
    publish,
    umis_per_cell,
)

BARCODES = ["AAA", "CCC", "GGG"]


@pytest.fixture(autouse=True)
def _close_figures():
    """pyplot keeps every figure until it is closed, and this file draws many."""
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


def mudata():
    """Three barcodes: the first gene-only, the last guide-only."""
    import anndata as ad
    import mudata as md

    gene_var = pd.DataFrame(
        {
            # MT-CO1 is on chrM; MT1A is a nuclear metallothionein that the
            # symbol prefix would sweep in and the contig does not.
            "gene_chr": ["chrM", "chr16", "chr19", "chr1", "chr11"],
            "symbol": ["MT-CO1", "MT1A", "RPS2", "RPL5P1", "VWF"],
            "gene_type": [
                "protein_coding",
                "protein_coding",
                "protein_coding",
                "processed_pseudogene",
                "protein_coding",
            ],
            "readout": [False, False, False, False, True],
        },
        index=pd.Index([f"ENSG{n}" for n in range(5)], name="gene_id"),
    )
    gene = ad.AnnData(
        X=np.array([[10, 10, 20, 10, 50], [0, 0, 0, 0, 100]], dtype="float32"),
        obs=pd.DataFrame(index=pd.Index(BARCODES[:2], name="barcode")),
        var=gene_var,
    )
    guide = ad.AnnData(
        X=np.array([[3, 0], [0, 7]], dtype="float32"),
        obs=pd.DataFrame(index=pd.Index(BARCODES[1:], name="barcode")),
        var=pd.DataFrame(index=pd.Index(["g1", "g2"], name="guide_id")),
    )
    return md.MuData({"gene": gene, "guide": guide})


def test_a_barcode_the_modality_does_not_hold_counts_as_zero_umis():
    data = mudata()
    assert umis_per_cell(data, "gene").tolist() == [100.0, 100.0, 0.0]
    assert umis_per_cell(data, "guide").tolist() == [0.0, 3.0, 7.0]


def test_mito_is_the_contig_not_the_symbol():
    classes = gene_classes(mudata()["gene"].var)
    # MT-CO1 only: MT1A shares the prefix and is nuclear.
    assert classes["mito"].tolist() == [True, False, False, False, False]


def test_ribosomal_pseudogenes_are_counted_under_their_own_label():
    classes = gene_classes(mudata()["gene"].var)
    assert classes["ribo"].tolist() == [False, False, True, False, False]
    assert classes["ribo_pseudogene"].tolist() == [False, False, False, True, False]


def test_the_readout_panel_becomes_a_class_only_when_var_declares_one():
    var = mudata()["gene"].var
    assert gene_classes(var)["targeted"].tolist() == [False] * 4 + [True]
    assert "targeted" not in gene_classes(var.drop(columns="readout"))


def test_fractions_are_percentages_and_a_cell_outside_the_modality_has_none():
    frame = cell_fractions(mudata())
    assert frame.index.tolist() == BARCODES
    assert frame.loc["AAA", "total_counts"] == 100.0
    assert frame.loc["AAA", "pct_counts_mito"] == pytest.approx(10.0)
    assert frame.loc["AAA", "pct_counts_targeted"] == pytest.approx(50.0)
    # Not measured, rather than measured as zero.
    assert frame.loc["GGG"].isna().all()


def test_a_floor_keeps_a_cell_with_rna_and_no_guide():
    filtered, keep = filter_cells(mudata(), min_umis=50)
    assert keep.tolist() == [True, True, False]
    assert filtered.n_obs == 2
    # AAA has no guide counts at all and is still a cell.
    assert filtered["guide"].n_obs == 1


def test_the_floor_has_no_default_because_the_knee_is_per_channel():
    with pytest.raises(TypeError):
        filter_cells(mudata())


def test_a_figure_writes_its_numbers_and_its_caption_beside_it(tmp_path):
    figure, numbers, caption = knee(mudata())
    written = publish(tmp_path / "knee", figure, numbers, caption)
    assert set(written) == {".pdf", ".tsv.gz", ".txt"}
    assert all(path.exists() for path in written.values())
    with gzip.open(written[".tsv.gz"], "rt") as handle:
        assert pd.read_csv(handle, sep="\t").equals(numbers)
    assert written[".txt"].read_text().strip() == caption.strip()


def test_a_figure_without_numbers_or_a_caption_is_refused(tmp_path):
    figure, numbers, caption = knee(mudata())
    with pytest.raises(FigureError):
        publish(tmp_path / "a", figure, numbers.iloc[:0], caption)
    with pytest.raises(FigureError):
        publish(tmp_path / "b", figure, numbers, "   ")


def test_mad_flags_both_tails_and_never_an_unmeasured_cell():
    from crispr_analysis_utils.qc import mad_outlier

    values = np.array([100.0, 110.0, 105.0, 95.0, 100.0, 1e6, np.nan])
    flagged = mad_outlier(values, nmads=3)
    assert flagged[5]  # the deep outlier
    assert not flagged[6]  # not measured, so not an outlier
    assert not flagged[:5].any()


def test_mad_can_be_held_to_one_tail():
    from crispr_analysis_utils.qc import mad_outlier

    values = np.array([100.0, 110.0, 105.0, 95.0, 100.0, 1.0, 1e6])
    assert mad_outlier(values, nmads=3, side="upper").tolist()[-2:] == [False, True]
    assert mad_outlier(values, nmads=3, side="lower").tolist()[-2:] == [True, False]


def test_a_zero_mad_is_refused_because_it_is_the_unfiltered_barcode_space():
    from crispr_analysis_utils.qc import mad_outlier

    # What the whole barcode space looks like: the median barcode is ambient.
    with pytest.raises(FigureError, match="not on the whole barcode space"):
        mad_outlier(np.array([2.0] * 1000 + [5000.0]))


def test_the_mad_is_rescaled_to_a_sigma_so_the_band_is_wider_than_the_raw_mad():
    from crispr_analysis_utils.qc import mad_outlier

    # Median 0, raw MAD 1, so scale="normal" puts one MAD at 1.4826.
    values = np.array([-1.2] + [-1.0] * 500 + [0.0] + [1.0] * 500 + [1.2])
    flagged = mad_outlier(values, nmads=1)
    # A raw MAD would flag both probes; a sigma-scaled one does not.
    assert not flagged.any()
    assert mad_outlier(values, nmads=0.8).tolist()[0] is True


def test_the_distributions_draw_one_panel_per_metric_and_bin_every_cell(tmp_path):
    from crispr_analysis_utils.qc import metric_distributions

    frame = cell_fractions(mudata())
    figure, numbers, caption = metric_distributions(frame, bins=4)

    metrics = list(frame.columns)
    assert sorted(numbers["metric"].unique()) == sorted(metrics)
    # The unmeasured cell is dropped, the other two are binned exactly once.
    per_metric = numbers.groupby("metric")[["n_kept", "n_removed"]].sum().sum(axis=1)
    assert (per_metric == 2).all()
    assert "MAD" in caption
    written = publish(tmp_path / "metrics", figure, numbers, caption)
    assert written[".pdf"].exists()


def test_the_distributions_need_a_metric_and_a_measured_cell():
    from crispr_analysis_utils.qc import metric_distributions

    frame = cell_fractions(mudata())
    with pytest.raises(FigureError, match="No metric"):
        metric_distributions(frame[[]])
    with pytest.raises(FigureError, match="No measured cell"):
        metric_distributions(frame.loc[["GGG"]])


def test_a_log_panel_bins_from_the_smallest_value_not_from_one():
    from crispr_analysis_utils.qc import metric_distributions

    frame = pd.DataFrame({"total_counts": [1000.0, 2000.0, 8000.0]})
    _, numbers, _ = metric_distributions(frame, bins=4)
    # Otherwise three of the four bins would sit under the data.
    assert numbers["bin_left"].min() == pytest.approx(1000.0)
    assert numbers["n_kept"].sum() == 3


def test_a_cell_on_the_smallest_or_largest_value_is_still_counted():
    from crispr_analysis_utils.qc import metric_distributions

    # logspace rounds these endpoints just inside the data.
    frame = pd.DataFrame({"total_counts": [271.0, 1000.0, 2321.0, 12150.0]})
    _, numbers, _ = metric_distributions(frame, bins=60)
    assert numbers["n_kept"].sum() == 4


def test_the_removed_layer_sits_on_the_kept_one_rather_than_over_it():
    from matplotlib.patches import StepPatch

    from crispr_analysis_utils.qc import metric_distributions

    # One outlier among many: the grey must be a sliver on top.
    rng = np.random.default_rng(0)
    frame = pd.DataFrame({"pct_counts_mito": np.r_[rng.normal(1.0, 0.1, 50), 9.0]})
    figure, numbers, _ = metric_distributions(frame, bins=4)

    axis = figure.axes[0]
    patches = [a for a in axis.get_children() if isinstance(a, StepPatch)]
    assert len(patches) == 2
    kept = patches[0].get_data().values
    stacked = patches[1].get_data().values
    baseline = patches[1].get_data().baseline
    # The grey starts where the blue ends and only ever goes up from there.
    assert np.allclose(baseline, kept)
    assert (stacked >= kept).all()
    assert np.allclose(stacked - kept, numbers["n_removed"].to_numpy())


def test_the_line_is_the_threshold_and_the_colour_changes_at_it():
    from crispr_analysis_utils.qc import mad_bounds, metric_distributions

    rng = np.random.default_rng(0)
    values = np.r_[rng.normal(1.0, 0.1, 50), 9.0]
    frame = pd.DataFrame({"pct_counts_mito": values})
    figure, numbers, caption = metric_distributions(frame, bins=4)

    # Upper tail only for a fraction, so one line, at the upper bound.
    _, high = mad_bounds(values, side="upper")
    lines = [x.get_xdata()[0] for x in figure.axes[0].get_lines()]
    assert lines == pytest.approx([high])
    assert numbers["threshold_high"].iloc[0] == pytest.approx(high)
    assert not np.isfinite(numbers["threshold_low"].iloc[0])
    # Everything grey lies beyond the line, everything blue short of it.
    beyond = numbers[numbers["n_removed"] > 0]
    assert (beyond["bin_right"] > high).all()
    assert f"{high:,.2f}" in caption


def test_every_line_is_labelled_with_the_rule_and_the_value_it_drew():
    from crispr_analysis_utils.qc import mad_bounds, metric_distributions

    rng = np.random.default_rng(0)
    values = np.r_[rng.normal(1.0, 0.1, 50), 9.0]
    figure, _, _ = metric_distributions(
        pd.DataFrame({"pct_counts_mito": values}), nmads=3.0, bins=4
    )

    _, high = mad_bounds(values, nmads=3.0, side="upper")
    labels = [t.get_text() for t in figure.axes[0].texts]
    assert labels == [f"3 MAD  {high:.3g}"]


def test_a_line_label_sits_inside_the_band_not_against_the_panel_edge():
    from crispr_analysis_utils.qc import metric_distributions

    rng = np.random.default_rng(0)
    figure, _, _ = metric_distributions(
        pd.DataFrame({"total_counts": rng.lognormal(7.0, 0.4, 500)}), bins=10
    )
    axis = figure.axes[0]
    figure.canvas.draw()
    lower, upper = sorted(line.get_xdata()[0] for line in axis.get_lines())
    # The loop draws the lower line first, so the labels come in that order.
    below, above = (text.get_window_extent() for text in axis.texts)
    # Both fall between the lines, where the kept cells are: outside them
    # the panel edge is a few points away and the label is cut in half.
    assert below.x0 >= axis.transData.transform((lower, 0))[0]
    assert above.x1 <= axis.transData.transform((upper, 0))[0]


def test_the_panel_says_how_many_it_put_beyond_the_band():
    from crispr_analysis_utils.qc import metric_distributions

    rng = np.random.default_rng(0)
    frame = pd.DataFrame({"pct_counts_mito": np.r_[rng.normal(1.0, 0.1, 50), 9.0, 9.1]})
    figure, numbers, _ = metric_distributions(frame, bins=4)

    assert figure.axes[0].get_title(loc="right") == "2 of 52 beyond"
    assert int(numbers["n_removed"].sum()) == 2


def test_a_panel_with_one_layer_carries_no_legend_naming_it():
    from crispr_analysis_utils.qc import metric_distributions

    # Nothing is beyond the band, so a legend would name only the blue.
    figure, _, _ = metric_distributions(
        pd.DataFrame({"pct_counts_mito": [1.0, 1.1, 1.2, 1.3]}), bins=4
    )
    assert all(axis.get_legend() is None for axis in figure.axes)


def test_a_metric_with_no_spread_loses_its_threshold_not_the_whole_figure():
    from crispr_analysis_utils.qc import metric_distributions

    frame = pd.DataFrame(
        {"pct_counts_mito": [1.0] * 50, "pct_counts_ribo": [1.0, 2.0] * 25}
    )
    figure, numbers, caption = metric_distributions(frame, bins=4)

    assert len(figure.axes) >= 2
    assert figure.axes[0].get_title(loc="right") == "no threshold"
    flat = numbers[numbers["metric"] == "pct_counts_mito"]
    assert not np.isfinite(flat["threshold_high"]).any()
    assert "pct_counts_mito" in caption.split("No threshold is drawn for")[1]


def test_the_scatter_cuts_into_quadrants_and_keeps_only_the_upper_right():
    from crispr_analysis_utils.qc import metric_scatter

    frame = pd.DataFrame(
        {
            "total_counts": [500.0, 2000.0, 500.0, 2000.0],
            "pct_counts_targeted": [10.0, 10.0, 60.0, 60.0],
        }
    )
    figure, numbers, caption = metric_scatter(frame, x_threshold=1000, y_threshold=50)
    assert numbers["kept"].tolist() == [False, False, False, True]
    assert "kept 1" in caption
    assert "low_both 1" in caption
    # One line per threshold, on the scatter panel.
    scatter_axis = figure.axes[0]
    assert len(scatter_axis.get_lines()) == 2


def test_the_scatter_drops_a_cell_missing_either_metric_and_refuses_an_empty_one():
    from crispr_analysis_utils.qc import metric_scatter

    frame = pd.DataFrame(
        {"total_counts": [1000.0, np.nan], "pct_counts_targeted": [60.0, 60.0]}
    )
    _, numbers, _ = metric_scatter(frame, marginals=False)
    assert len(numbers) == 1

    with pytest.raises(FigureError, match="not a column"):
        metric_scatter(frame, x="missing")
    with pytest.raises(FigureError, match="No cell has both"):
        metric_scatter(frame.iloc[1:2], y="total_counts", x="pct_counts_targeted")


def test_the_scatter_drops_the_legend_when_nothing_is_removed():
    from crispr_analysis_utils.qc import metric_scatter

    frame = pd.DataFrame(
        {"total_counts": [2000.0, 3000.0], "pct_counts_targeted": [60.0, 61.0]}
    )
    figure, _, _ = metric_scatter(frame, x_threshold=1000, y_threshold=25)
    assert figure.axes[0].get_legend() is None


def test_an_extra_gene_class_is_measured_alongside_the_built_in_ones():
    data = mudata()
    # Another analysis's own definition, wider than this module's.
    wider = np.array([False, True, True, True, False])
    frame = cell_fractions(data, extra={"ribo_wide": wider})

    assert "pct_counts_ribo_wide" in frame
    # AAA: 10 + 20 + 10 of 100 UMIs.
    assert frame.loc["AAA", "pct_counts_ribo_wide"] == pytest.approx(40.0)
    # The built-in classes are untouched.
    assert frame.loc["AAA", "pct_counts_ribo"] == pytest.approx(20.0)


def test_an_extra_mask_of_the_wrong_length_is_refused():
    with pytest.raises(FigureError, match="covers 2 genes"):
        cell_fractions(mudata(), extra={"wrong": np.array([True, False])})


def test_mean_upm_scales_each_cell_before_averaging():
    from crispr_analysis_utils.qc import mean_upm

    data = mudata()
    values = mean_upm(data["gene"])
    # AAA is 10/10/20/10/50 of 100; CCC is 0/0/0/0/100 of 100.
    assert values.tolist() == pytest.approx(
        [50_000.0, 50_000.0, 100_000.0, 50_000.0, 750_000.0]
    )


def test_a_cell_with_no_umis_does_not_divide_by_zero():
    import anndata as ad
    import mudata as md

    from crispr_analysis_utils.qc import mean_upm

    empty = ad.AnnData(
        X=np.array([[1.0, 1.0], [0.0, 0.0]], dtype="float32"),
        obs=pd.DataFrame(index=pd.Index(["a", "b"], name="barcode")),
        var=pd.DataFrame(index=pd.Index(["g1", "g2"], name="gene_id")),
    )
    values = mean_upm(md.MuData({"gene": empty})["gene"])
    # The empty cell contributes zero, so the mean is halved, not infinite.
    assert values.tolist() == pytest.approx([250_000.0, 250_000.0])


def test_the_panel_boxplot_splits_on_the_readout_flag(tmp_path):
    from crispr_analysis_utils.qc import panel_expression

    figure, numbers, caption = panel_expression(mudata())
    assert numbers["in_panel"].sum() == 1
    assert len(numbers) == 5
    assert "1 genes of the panel against the 4 others" in caption
    written = publish(tmp_path / "panel", figure, numbers, caption)
    assert written[".tsv.gz"].exists()


def test_the_panel_boxplot_needs_a_flag_that_marks_something():
    from crispr_analysis_utils.qc import panel_expression

    data = mudata()
    with pytest.raises(FigureError, match="not a column"):
        panel_expression(data, flag="missing")
    data["gene"].var["readout"] = False
    with pytest.raises(FigureError, match="marks no gene"):
        panel_expression(data)


def test_the_panel_caption_says_what_a_zero_background_median_means():
    from crispr_analysis_utils.qc import panel_expression

    data = mudata()
    # Most of an unfiltered reference is never seen, so its median is zero.
    data["gene"].X = np.array([[0, 0, 0, 0, 50], [0, 0, 0, 0, 100]], dtype="float32")
    _, _, caption = panel_expression(data)
    assert "inf" not in caption
    assert "background median being zero" in caption
    assert "0 of 4 background genes seen" in caption


def test_the_scatter_serves_a_per_gene_comparison_too():
    from crispr_analysis_utils.qc import metric_scatter

    frame = pd.DataFrame(
        {"reference_upm": [5.0, 100.0], "mean_umi": [0.2, 3.0]},
        index=pd.Index(["ENSG1", "ENSG2"], name="gene"),
    )
    _, numbers, caption = metric_scatter(
        frame,
        x="reference_upm",
        y="mean_umi",
        x_threshold=30,
        y_threshold=1,
        log=("reference_upm", "mean_umi"),
        unit="genes",
    )
    assert "2 genes, one point each" in caption
    assert "over the same genes" in caption
    assert numbers.columns[0] == "gene"
    assert numbers["kept"].tolist() == [False, True]


def test_a_zero_on_a_log_axis_is_reported_rather_than_silently_dropped():
    from crispr_analysis_utils.qc import metric_scatter

    frame = pd.DataFrame(
        {
            "total_counts": [0.0, 1000.0, 3000.0],
            "pct_counts_targeted": [60.0, 60.0, 60.0],
        }
    )
    _, numbers, caption = metric_scatter(frame, x_threshold=500, y_threshold=25)

    assert numbers["drawn"].tolist() == [False, True, True]
    # Still counted: it is below the threshold, and the quadrant says so.
    assert numbers["kept"].tolist() == [False, True, True]
    assert "1 of the cells hold a zero on a log axis" in caption


def test_the_point_size_follows_how_many_points_there_are():
    from crispr_analysis_utils.qc import metric_scatter

    few = pd.DataFrame(
        {"total_counts": [1000.0] * 10, "pct_counts_targeted": [60.0] * 10}
    )
    many = pd.DataFrame(
        {"total_counts": [1000.0] * 3000, "pct_counts_targeted": [60.0] * 3000}
    )
    small = metric_scatter(many, marginals=False)[0].axes[0].collections[0].get_sizes()
    large = metric_scatter(few, marginals=False)[0].axes[0].collections[0].get_sizes()
    assert large[0] > small[0]


def test_groups_colour_the_points_instead_of_the_threshold():
    from crispr_analysis_utils.qc import OKABE_ITO, REMOVED, metric_scatter

    frame = pd.DataFrame(
        {"reference_upm": [5.0, 100.0, 200.0], "mean_umi": [0.2, 3.0, 4.0]},
        index=pd.Index(["a", "b", "c"], name="gene"),
    )
    groups = pd.Series(["Day 0", "Day 4", None], index=frame.index)
    figure, numbers, caption = metric_scatter(
        frame,
        x="reference_upm",
        y="mean_umi",
        x_threshold=30,
        y_threshold=1,
        log=("reference_upm", "mean_umi"),
        unit="genes",
        groups=groups,
        marginals=False,
    )

    assert numbers["group"].tolist()[:2] == ["Day 0", "Day 4"]
    assert pd.isna(numbers["group"].iloc[2])
    labels = [t.get_text() for t in figure.axes[0].get_legend().get_texts()]
    assert labels == ["Day 0", "Day 4"]
    # Nothing is being filtered, so the caption must not say so.
    assert "removed" not in caption
    assert "Colour is the group" in caption
    # The unlabelled point is grey and drawn first, behind the classes.
    collections = figure.axes[0].collections
    assert collections[0].get_facecolor()[0][0] == pytest.approx(float(REMOVED))
    assert len(collections) == 3
    assert OKABE_ITO[0].lstrip("#")


def test_groups_still_report_the_quadrant_counts():
    from crispr_analysis_utils.qc import metric_scatter

    frame = pd.DataFrame({"reference_upm": [5.0, 100.0], "mean_umi": [0.2, 3.0]})
    _, _, caption = metric_scatter(
        frame,
        x="reference_upm",
        y="mean_umi",
        x_threshold=30,
        y_threshold=1,
        log=("reference_upm", "mean_umi"),
        unit="genes",
        groups=["low", "high"],
        marginals=False,
    )
    assert "kept 1" in caption and "low_both 1" in caption


def test_a_label_keeps_its_colour_whatever_order_the_rows_arrive_in():
    from crispr_analysis_utils.qc import metric_scatter

    def legend(order):
        frame = pd.DataFrame(
            {"reference_upm": [10.0, 20.0, 30.0], "mean_umi": [1.0, 2.0, 3.0]}
        )
        figure, _, _ = metric_scatter(
            frame,
            x="reference_upm",
            y="mean_umi",
            log=("reference_upm", "mean_umi"),
            groups=order,
            marginals=False,
        )
        return [t.get_text() for t in figure.axes[0].get_legend().get_texts()]

    assert legend(["Day 4", "Day 0", "Day 2"]) == legend(["Day 0", "Day 2", "Day 4"])


def test_a_categorical_sets_its_own_legend_order():
    from crispr_analysis_utils.qc import metric_scatter

    frame = pd.DataFrame({"reference_upm": [10.0, 20.0], "mean_umi": [1.0, 2.0]})
    wanted = ["Housekeeping", "Day 0"]
    groups = pd.Categorical(["Day 0", "Housekeeping"], categories=wanted)
    figure, _, _ = metric_scatter(
        frame,
        x="reference_upm",
        y="mean_umi",
        log=("reference_upm", "mean_umi"),
        groups=groups,
        marginals=False,
    )
    assert [t.get_text() for t in figure.axes[0].get_legend().get_texts()] == wanted


def test_a_grouped_point_is_drawn_solid_and_the_background_is_not():
    from crispr_analysis_utils.qc import metric_scatter

    frame = pd.DataFrame({"reference_upm": [10.0, 20.0], "mean_umi": [1.0, 2.0]})
    figure, _, _ = metric_scatter(
        frame,
        x="reference_upm",
        y="mean_umi",
        log=("reference_upm", "mean_umi"),
        groups=[None, "Day 4"],
        marginals=False,
    )
    background, named = figure.axes[0].collections
    assert background.get_alpha() < 1.0
    assert named.get_alpha() == 1.0


def test_a_grouped_point_is_outlined_so_a_pale_fill_still_reads():
    from crispr_analysis_utils.qc import metric_scatter

    frame = pd.DataFrame({"reference_upm": [10.0, 20.0], "mean_umi": [1.0, 2.0]})
    figure, _, _ = metric_scatter(
        frame,
        x="reference_upm",
        y="mean_umi",
        log=("reference_upm", "mean_umi"),
        groups=[None, "Day 4"],
        marginals=False,
    )
    background, named = figure.axes[0].collections
    assert background.get_linewidths()[0] == 0
    assert named.get_linewidths()[0] > 0
    assert named.get_edgecolors()[0].tolist() != named.get_facecolors()[0].tolist()

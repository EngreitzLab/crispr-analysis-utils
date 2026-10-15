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
    flagged = np.array([False, True, False])
    figure, numbers, caption = metric_distributions(frame, flagged=flagged, bins=4)

    metrics = list(frame.columns)
    assert sorted(numbers["metric"].unique()) == sorted(metrics)
    # The unmeasured cell is dropped, the other two are binned exactly once.
    per_metric = numbers.groupby("metric")[["n_kept", "n_flagged"]].sum().sum(axis=1)
    assert (per_metric == 2).all()
    assert "flagged" in caption
    written = publish(tmp_path / "metrics", figure, numbers, caption)
    assert written[".pdf"].exists()


def test_the_distributions_need_a_metric_and_a_measured_cell():
    from crispr_analysis_utils.qc import metric_distributions

    frame = cell_fractions(mudata())
    with pytest.raises(FigureError, match="No metric"):
        metric_distributions(frame[[]])
    with pytest.raises(FigureError, match="No measured cell"):
        metric_distributions(frame.loc[["GGG"]])

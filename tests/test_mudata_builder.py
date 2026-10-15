"""The MuData builder: what it assembles, and what it refuses."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from guide_annotation_fixtures import write_gtf

from crispr_analysis_utils.guide_annotation import build_store
from crispr_analysis_utils.mudata_builder import (
    BuildError,
    build,
    check_same_store,
    gene_var,
    guide_var,
    read_counts,
    reference_block,
    run,
)

MEASUREMENT_SET = "IGVFDS0000TEST"


@pytest.fixture
def store(tmp_path):
    return build_store(write_gtf(tmp_path / "test.gtf"), tmp_path / "store.parquet")


def counts(path, barcodes, features, index_name):
    """A minimal h5ad, shaped like the ones kb writes."""
    import anndata as ad

    frame = pd.DataFrame(index=pd.Index(features, name=index_name))
    obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))
    matrix = ad.AnnData(
        X=np.arange(len(barcodes) * len(features), dtype="float32").reshape(
            len(barcodes), len(features)
        ),
        obs=obs,
        var=frame,
    )
    matrix.write_h5ad(path)
    return path


def guide_table(path, guide_ids):
    pd.DataFrame(
        {
            "guide_id": guide_ids,
            "spacer": ["A" * 20] * len(guide_ids),
            "type": ["targeting", "safe-targeting", "non-targeting"][: len(guide_ids)],
        }
    ).to_csv(path, sep="\t", index=False)
    return path


def test_every_gene_type_is_annotated_not_only_the_coding_ones(store):
    # DDD is lncRNA. The counts come from the whole GTF, so it is a real row.
    ids = pd.Index(["ENSG00000000001.5", "ENSG00000000004.2"])
    out = gene_var(store, ids)
    assert list(out.index) == list(ids)
    assert out["symbol"].tolist() == ["AAA", "DDD"]
    assert out["gene_type"].tolist() == ["protein_coding", "lncRNA"]


def test_a_gene_the_store_does_not_hold_keeps_its_row(store):
    out = gene_var(store, pd.Index(["ENSG00000000001.5", "ENSG99999999999.1"]))
    assert len(out) == 2
    assert pd.isna(out["symbol"].iloc[1])


def test_the_version_rule_is_the_stores_so_par_y_keeps_its_suffix(store):
    # Splitting on the first dot would give this the chrX gene's annotation.
    out = gene_var(store, pd.Index(["ENSG00000000001.5_PAR_Y"]))
    assert pd.isna(out["symbol"].iloc[0])


def test_guide_var_keeps_every_column_and_the_matrix_order(tmp_path):
    path = guide_table(tmp_path / "guides.tsv", ["g2", "g1", "g3"])
    out = guide_var(path, pd.Index(["g1", "g2"]))
    assert list(out.index) == ["g1", "g2"]
    assert "spacer" in out.columns and "type" in out.columns
    assert "guide_id" not in out.columns


def test_the_four_guide_classes_survive_untouched(tmp_path):
    path = guide_table(tmp_path / "guides.tsv", ["g1", "g2", "g3"])
    out = guide_var(path, pd.Index(["g1", "g2", "g3"]))
    assert out["type"].tolist() == ["targeting", "safe-targeting", "non-targeting"]


def test_a_table_without_guide_id_is_refused(tmp_path):
    path = tmp_path / "bad.tsv"
    pd.DataFrame({"name": ["g1"]}).to_csv(path, sep="\t", index=False)
    with pytest.raises(BuildError, match="no guide_id column"):
        guide_var(path, pd.Index(["g1"]))


def test_cells_are_the_union_and_obsm_says_which_modality_holds_each(store, tmp_path):
    rna = read_counts(
        counts(tmp_path / "r.h5ad", ["c1", "c2"], ["ENSG00000000001.5"], "gene_id")
    )
    guide = read_counts(counts(tmp_path / "g.h5ad", ["c2", "c3"], ["g1"], "guide_id"))
    data = build(
        rna,
        guide,
        gene_table=gene_var(store, rna.var_names),
        guide_table=guide_var(
            guide_table(tmp_path / "gt.tsv", ["g1"]), guide.var_names
        ),
        measurement_set=MEASUREMENT_SET,
    )
    assert data.n_obs == 3
    assert list(data.obs_names) == [f"c{i}_{MEASUREMENT_SET}" for i in (1, 2, 3)]
    assert np.asarray(data.obsm["gene"]).ravel().tolist() == [True, True, False]
    assert np.asarray(data.obsm["guide"]).ravel().tolist() == [False, True, True]


def test_the_cell_index_carries_the_measurement_set(store, tmp_path):
    rna = read_counts(
        counts(tmp_path / "r.h5ad", ["c1"], ["ENSG00000000001.5"], "gene_id")
    )
    guide = read_counts(counts(tmp_path / "g.h5ad", ["c1"], ["g1"], "guide_id"))
    data = build(
        rna,
        guide,
        gene_table=gene_var(store, rna.var_names),
        guide_table=guide_var(
            guide_table(tmp_path / "gt.tsv", ["g1"]), guide.var_names
        ),
        measurement_set=MEASUREMENT_SET,
    )
    assert data["gene"].obs["barcode"].tolist() == ["c1"]
    assert data["gene"].obs["measurement_set"].tolist() == [MEASUREMENT_SET]


def test_an_annotation_of_the_wrong_length_is_refused(store, tmp_path):
    rna = read_counts(
        counts(tmp_path / "r.h5ad", ["c1"], ["ENSG00000000001.5"], "gene_id")
    )
    guide = read_counts(counts(tmp_path / "g.h5ad", ["c1"], ["g1"], "guide_id"))
    with pytest.raises(BuildError, match="must line up"):
        build(
            rna,
            guide,
            gene_table=pd.DataFrame(index=["a", "b"]),
            guide_table=guide_var(
                guide_table(tmp_path / "gt.tsv", ["g1"]), guide.var_names
            ),
            measurement_set=MEASUREMENT_SET,
        )


def test_the_reference_block_names_files_rather_than_pathing_them(store, tmp_path):
    block = reference_block(store, gtf_accession="IGVFFI9561BASO")
    assert "/" not in block["gtf"]
    assert block["gtf_accession"] == "IGVFFI9561BASO"
    assert len(block["gtf_md5"]) == 32


def test_an_accession_falls_back_to_the_file_name(store):
    block = reference_block(store)
    assert block["gtf_accession"] == block["gtf"] != ""


def test_a_genome_that_does_not_exist_is_refused(store, tmp_path):
    with pytest.raises(BuildError, match="does not exist"):
        reference_block(store, genome=tmp_path / "absent.fa")


def test_a_different_store_than_the_guides_were_built_on_is_refused(tmp_path):
    results = tmp_path / "results"
    (results / "igvf").mkdir(parents=True)
    (results / "run_metadata").mkdir()
    guides = results / "igvf" / "guide_metadata.tsv"
    guide_table(guides, ["g1"])
    (results / "run_metadata" / "igvf_export_20260101T000000Z.json").write_text(
        json.dumps({"inputs": {"store": str(tmp_path / "other.parquet")}}),
        encoding="utf-8",
    )
    with pytest.raises(BuildError, match="different GENCODE"):
        check_same_store(guides, tmp_path / "mine.parquet")


def test_no_run_record_is_not_an_error(tmp_path):
    guides = guide_table(tmp_path / "guide_metadata.tsv", ["g1"])
    assert check_same_store(guides, tmp_path / "mine.parquet") == ""


def test_a_run_writes_the_object_and_its_record(store, tmp_path):
    rna = counts(tmp_path / "r.h5ad", ["c1", "c2"], ["ENSG00000000001.5"], "gene_id")
    guide = counts(tmp_path / "g.h5ad", ["c2"], ["g1"], "guide_id")
    result = run(
        rna,
        guide,
        guide_table(tmp_path / "gt.tsv", ["g1"]),
        tmp_path / "store.parquet",
        tmp_path / "out",
        measurement_set=MEASUREMENT_SET,
    )
    assert (tmp_path / "out" / "mudata" / "mudata.h5mu").is_file()
    assert result.metadata["complete"] is True
    assert result.metadata["counts"]["cells"] == 2
    assert result.metadata["counts"]["cells_gene_only"] == 1
    assert result.metadata["counts"]["cells_both"] == 1
    assert result.metadata["parameters"]["reference"]["gtf"] != ""


def test_a_missing_matrix_is_refused(store, tmp_path):
    with pytest.raises(BuildError, match="does not exist"):
        read_counts(tmp_path / "absent.h5ad")

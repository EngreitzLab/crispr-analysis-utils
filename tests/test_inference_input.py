"""Joining QC'd channels into one inference object, and ``cau inference-input``.

Each channel is built here: three cells, five genes and eight guides. The
guides are two targeting elements of two guides each, three non-targeting
guides and one safe-targeting guide, so the median element has two guides
and the four controls make two buckets.
"""

import json

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from crispr_analysis_utils import cli
from crispr_analysis_utils.inference_input import (
    ASSIGNMENT,
    COVARIATES,
    PAIR_COLUMNS,
    InputError,
    build,
    run,
    target_keys,
)

GENES = [
    "ENSG00000000001.5",
    "ENSG00000000002.3",
    "ENSG00000000003.7",
    "ENSG00000000003.7_PAR_Y",
    "ENSG00000000004.1",
]


def guide_var():
    rows = [
        ("e1a", "TRUE", "targeting", "chr1:100-400", "chr1", "100", "400"),
        ("e1b", "TRUE", "targeting", "chr1:100-400", "chr1", "100", "400"),
        ("e2a", "TRUE", "positive control", "chr2:500-900", "chr2", "500", "900"),
        ("e2b", "TRUE", "positive control", "chr2:500-900", "chr2", "500", "900"),
        ("nt3", "FALSE", "non-targeting", "", "", "", ""),
        ("nt1", "FALSE", "non-targeting", "", "", "", ""),
        ("nt2", "FALSE", "non-targeting", "", "", "", ""),
        ("s1", "FALSE", "safe-targeting", "", "", "", ""),
    ]
    frame = pd.DataFrame(
        rows,
        columns=[
            "guide_id",
            "targeting",
            "type",
            "intended_target_name",
            "intended_target_chr",
            "intended_target_start",
            "intended_target_end",
        ],
    ).set_index("guide_id")
    frame.index.name = None
    return frame


def channel(ms, *, guide_cells=None, calls=None):
    """One QC'd channel, cells indexed ``{barcode}_{measurement_set}``."""
    cells = [f"{bc}_{ms}" for bc in ("AAA", "CCC", "GGG")]
    counts = np.array(
        [
            [10, 0, 5, 1, 4],
            [3, 2, 0, 0, 5],
            [0, 7, 1, 0, 2],
        ]
    )
    gene = ad.AnnData(
        X=sp.csr_matrix(counts.astype(np.float32)),
        obs=pd.DataFrame({"measurement_set": ms}, index=cells),
        var=pd.DataFrame(
            {
                "symbol": ["A", "RPL1", "C", "C", "MT-X"],
                "gene_chr": ["chr1", "chr2", "chrX", "chrY", "chrM"],
                "gene_type": ["protein_coding"] * 5,
            },
            index=GENES,
        ),
    )
    guide_counts = np.zeros((3, 8))
    guide_counts[:, 0] = [30, 0, 25]
    guide_counts[:, 4] = [0, 40, 1]
    called = (guide_counts >= 20).astype(np.uint8) if calls is None else calls
    guide = ad.AnnData(
        X=sp.csr_matrix(guide_counts),
        obs=pd.DataFrame(index=cells),
        var=guide_var(),
        layers={ASSIGNMENT: sp.csr_matrix(called)},
    )
    if guide_cells is not None:
        guide = guide[guide_cells].copy()
    data = md.MuData({"gene": gene, "guide": guide})
    data.uns["measurement_set"] = ms
    data.uns["qc"] = {"min_umis": 1}
    return data


def pairs():
    return pd.DataFrame(
        [
            ("e1a", "ENSG00000000001", "discovery"),
            ("e1b", "ENSG00000000002", "discovery"),
            ("e2a", "ENSG00000000004", "positive_control"),
            ("gone", "ENSG00000000001", "discovery"),
            ("e1a", "ENSG00000000009", "discovery"),
        ],
        columns=["guide_id", "gene_id", "pair_type"],
    ).assign(
        intended_target_name=lambda f: f.guide_id.map(
            {"e1a": "chr1:100-400", "e1b": "chr1:100-400", "e2a": "chr2:500-900"}
        ).fillna("chr9:1-2"),
        intended_target_chr=lambda f: f.intended_target_name.str.split(":").str[0],
        intended_target_start=lambda f: f.guide_id.map(
            {"e1a": 100, "e1b": 100, "e2a": 500}
        ).fillna(1),
        intended_target_end=lambda f: f.guide_id.map(
            {"e1a": 400, "e1b": 400, "e2a": 900}
        ).fillna(2),
    )


def built(*channels, **kwargs):
    kwargs.setdefault("moi", "high")
    kwargs.setdefault("capture_method", "CROP-seq")
    return build(list(channels) or [channel("MS1")], pairs(), **kwargs)


# --- target keys --------------------------------------------------------------


def test_controls_are_sorted_by_id_and_cut_into_buckets_the_size_of_the_median_element():
    var = target_keys(guide_var())
    names = var["intended_target_name"].to_dict()
    assert {g: names[g] for g in ("nt1", "nt2", "nt3", "s1")} == {
        "nt1": "non-targeting|1",
        "nt2": "non-targeting|1",
        "nt3": "non-targeting|2",
        "s1": "non-targeting|2",
    }
    assert var.loc["s1", "intended_target_key"] == (
        "coords::non-targeting|2::non-targeting::2::2"
    )


def test_a_targeting_guide_is_keyed_by_its_element():
    var = target_keys(guide_var())
    assert (
        var.loc["e1a", "intended_target_key"] == "coords::chr1:100-400::chr1::100::400"
    )
    assert var["targeting"].dtype == bool
    assert str(var["intended_target_start"].dtype) == "Int64"


# --- refusals -----------------------------------------------------------------


def test_a_channel_without_guide_calls_is_refused():
    data = channel("MS1")
    del data["guide"].layers[ASSIGNMENT]
    with pytest.raises(InputError, match="Call guides first"):
        built(data)


def test_counts_in_place_of_calls_are_refused():
    calls = np.zeros((3, 8))
    calls[0, 0] = 30
    with pytest.raises(InputError, match="not binary"):
        built(channel("MS1", calls=calls))


def test_channels_must_share_their_features():
    other = channel("MS2")
    other.mod["guide"] = other["guide"][:, ::-1].copy()
    with pytest.raises(InputError, match="do not share"):
        built(channel("MS1"), other)


def test_an_unknown_moi_is_refused():
    with pytest.raises(InputError, match="moi"):
        built(channel("MS1"), moi="medium")


def test_pairs_from_another_guide_table_are_refused():
    table = pairs()
    table.loc[0, "intended_target_end"] = 401
    with pytest.raises(InputError, match="different guide tables"):
        build([channel("MS1")], table, moi="high", capture_method="CROP-seq")


# --- the object ---------------------------------------------------------------


def test_only_cells_held_by_both_modalities_are_kept():
    data, report = built(channel("MS1", guide_cells=["AAA_MS1", "GGG_MS1"]))
    assert list(data.obs_names) == ["AAA_MS1", "GGG_MS1"]
    assert list(data["gene"].obs_names) == list(data["guide"].obs_names)
    assert report["channels"][0]["both"] == 2


def test_channels_concatenate_and_batch_is_the_measurement_set():
    data, _ = built(channel("MS1"), channel("MS2"))
    assert data.n_obs == 6
    assert list(data.obs.columns) == list(COVARIATES)
    assert list(data.obs["batch"]) == ["MS1"] * 3 + ["MS2"] * 3


def test_covariates_carry_the_pipeline_names_and_values():
    data, _ = built()
    obs = data.obs.loc["AAA_MS1"]
    assert obs["total_gene_umis"] == 20
    assert obs["num_expressed_genes"] == 4
    assert obs["percent_mito"] == pytest.approx(20.0)
    assert obs["total_guide_umis"] == 30
    assert data["gene"].obs.loc["CCC_MS1", "pct_counts_ribo"] == pytest.approx(20.0)
    assert data["guide"].obs.loc["CCC_MS1", "num_expressed_guides"] == 1


def test_gene_ids_lose_their_version_and_a_par_y_copy_stays_distinct():
    data, _ = built()
    assert list(data["gene"].var_names) == [
        "ENSG00000000001",
        "ENSG00000000002",
        "ENSG00000000003",
        "ENSG00000000003.7_PAR_Y",
        "ENSG00000000004",
    ]
    assert list(data["gene"].var["gene_id_versioned"]) == GENES


def test_gene_ids_that_collide_once_unversioned_are_refused():
    data = channel("MS1")
    data["gene"].var_names = [*GENES[:4], "ENSG00000000001.6"]
    with pytest.raises(InputError, match="collide"):
        built(data)


def test_the_guide_modality_carries_the_pipeline_fields():
    data, _ = built()
    guide = data["guide"]
    assert list(guide.var["guide_id"]) == list(guide.var_names)
    assert "intended_target_key" in guide.var
    assert list(guide.uns["moi"]) == ["high"]
    assert list(guide.uns["capture_method"]) == ["CROP-seq"]
    assert guide.layers[ASSIGNMENT].dtype == np.uint8
    assert np.issubdtype(guide.X.dtype, np.unsignedinteger)


def test_the_guide_filters_count_what_they_drop():
    # Guides called per cell: AAA 1 (e1a), CCC 1 (nt3), GGG 1 (e1a). Add one.
    calls = np.zeros((3, 8), dtype=np.uint8)
    calls[0, [0, 1]] = 1
    calls[1, 0] = 1
    data, report = built(channel("MS1", calls=calls), min_guides=1, max_guides=1)
    assert {k: v for k, v in report["cells"].items() if k != "zero_covariate"} == {
        "in": 3,
        "below_min_guides": 1,
        "above_max_guides": 1,
        "out": 1,
    }
    assert list(data.obs_names) == ["CCC_MS1"]


def test_pairs_the_object_cannot_test_are_dropped_and_counted():
    data, report = built()
    kept = data.uns["pairs_to_test"]
    assert list(kept.columns[: len(PAIR_COLUMNS)]) == list(PAIR_COLUMNS)
    assert set(zip(kept["guide_id"], kept["gene_id"], strict=True)) == {
        ("e1a", "ENSG00000000001"),
        ("e1b", "ENSG00000000002"),
        ("e2a", "ENSG00000000004"),
    }
    assert report["pairs"]["guide_absent"] == 1
    assert report["pairs"]["gene_absent"] == 1
    assert report["pairs"]["guides_absent"] == ["gone"]
    assert kept["intended_target_key"].iloc[0].startswith("coords::chr1:100-400")


# --- end to end ---------------------------------------------------------------


@pytest.fixture
def files(tmp_path):
    paths = []
    for ms in ("MS1", "MS2"):
        path = tmp_path / f"{ms}.h5mu"
        channel(ms).write_h5mu(path)
        paths.append(path)
    table = tmp_path / "pairs_to_test.csv"
    pairs().to_csv(table, index=False)
    return paths, table


def test_run_writes_the_object_and_its_record(tmp_path, files):
    paths, table = files
    result = run(
        paths,
        table,
        tmp_path / "out",
        name="day4",
        capture_method="CROP-seq",
        min_guides=1,
    )
    assert result.path == tmp_path / "out" / "inference_input" / "day4.h5mu"
    data = md.read_h5mu(result.path)
    assert data.n_obs == 6
    assert len(data.uns["pairs_to_test"]) == 3
    assert set(data.uns["qc"]) == {"MS1", "MS2"}

    written = sorted(
        p.relative_to(tmp_path / "out").as_posix()
        for p in (tmp_path / "out").rglob("*")
        if p.is_file()
    )
    assert written[:2] == ["inference_input/day4.h5mu", "logs/inference_input.log"]
    (record,) = written[2:]
    state = json.loads((tmp_path / "out" / record).read_text())
    assert state["complete"]
    assert state["parameters"]["min_guides"] == 1
    assert state["counts"]["pairs"]["pairs_out"] == 3


def test_guide_ids_that_look_like_numbers_still_match_their_pairs(tmp_path):
    numbers = {old: f"{n:04d}" for n, old in enumerate(guide_var().index, start=1)}
    data = channel("MS1")
    data.mod["guide"].var_names = data.mod["guide"].var_names.map(numbers)
    data.update()
    data.write_h5mu(tmp_path / "MS1.h5mu")
    table = pairs().assign(guide_id=lambda f: f.guide_id.map(numbers).fillna("9999"))
    table.to_csv(tmp_path / "pairs.csv", index=False)
    result = run(
        [tmp_path / "MS1.h5mu"],
        tmp_path / "pairs.csv",
        tmp_path / "out",
        name="day4",
        capture_method="CROP-seq",
    )
    assert result.report["pairs"]["pairs_out"] == 3


def test_pairs_that_match_nothing_are_refused():
    with pytest.raises(InputError, match="None of the 5 pairs"):
        build(
            [channel("MS1")],
            pairs().assign(guide_id="elsewhere"),
            moi="high",
            capture_method="CROP-seq",
        )


def test_the_command_runs(tmp_path, files, capsys):
    paths, table = files
    status = cli.main(
        [
            "inference-input",
            "run",
            "--inputs",
            *map(str, paths),
            "--pairs",
            str(table),
            "--outdir",
            str(tmp_path / "out"),
            "--name",
            "day4",
            "--capture-method",
            "CROP-seq",
        ]
    )
    assert status == 0
    out = capsys.readouterr().out
    assert "2 channels, 6 cells in both modalities" in out
    assert "3 of 5 pairs kept" in out


def test_the_command_reports_a_missing_file(tmp_path, files, capsys):
    _, table = files
    with pytest.raises(SystemExit):
        cli.main(
            [
                "inference-input",
                "run",
                "--inputs",
                str(tmp_path / "nope.h5mu"),
                "--pairs",
                str(table),
                "--outdir",
                str(tmp_path / "out"),
                "--name",
                "day4",
                "--capture-method",
                "CROP-seq",
            ]
        )
    assert "does not exist" in capsys.readouterr().err


def test_a_stored_zero_is_not_a_call():
    calls = sp.csr_matrix(np.eye(3, 8, dtype=np.uint8))
    calls.data[0] = 0  # AAA's call, kept as an explicit zero
    _, report = built(channel("MS1", calls=calls), min_guides=1)
    assert report["cells"]["below_min_guides"] == 1


def test_a_cell_without_guide_umis_is_kept_and_counted_since_filtering_comes_later():
    data = channel("MS1")
    data["guide"].X = sp.csr_matrix(data["guide"].X.multiply(np.array([[1], [1], [0]])))
    data["guide"].layers[ASSIGNMENT] = sp.csr_matrix((3, 8), dtype=np.uint8)
    out, report = built(data)
    assert out.n_obs == 3
    assert report["cells"]["zero_covariate"]["total_guide_umis"] == 1

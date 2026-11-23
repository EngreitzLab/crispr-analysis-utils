"""Classifying pairs by the genes they skip.

One element, midpoint 10,000 on chr1, with genes at fixed distances from it
(negative upstream of the midpoint, matching the element's own coordinates,
not the gene's orientation):

=========  ========  ====
Gene       Distance  Side
=========  ========  ====
NEAR40     40        +
NEAR70     70        +
FAR200     200       +   (in the readout panel)
FAR500     500       +
OPPOSITE   -300       -
AWAY1000   1000      +
FAR5000    -5000      -
TRANS      n/a (chr2)
=========  ========  ====
"""

import json
import re

import pandas as pd
import pytest
from guide_annotation_fixtures import write_gtf

from crispr_analysis_utils import cli
from crispr_analysis_utils.guide_annotation import build_store
from crispr_analysis_utils.pairs import TARGET_COLUMNS, build_pairs
from crispr_analysis_utils.skipped_genes import (
    CATEGORIES,
    COLUMNS,
    NEAREST,
    NO_SKIP,
    SKIP_IN_PANEL,
    SKIP_OUTSIDE_PANEL,
    TRANS,
    SkippedGenesError,
    classify,
)
from crispr_analysis_utils.skipped_genes import run as classify_run

NEAR40 = "ENSG00000000040"
NEAR70 = "ENSG00000000070"
FAR200 = "ENSG00000000200"
FAR500 = "ENSG00000000500"
OPPOSITE = "ENSG00000000300"
AWAY1000 = "ENSG00000001000"
FAR5000 = "ENSG00000005000"
TRANSGENE = "ENSG00000099999"

ELEMENT = ("peak1", "chr1", 9_000, 11_000)  # midpoint 10,000


def tss_frame():
    return pd.DataFrame(
        {
            "gene_symbol": [
                "NEAR40",
                "NEAR70",
                "FAR200",
                "FAR500",
                "OPPOSITE",
                "AWAY1000",
                "FAR5000",
                "TRANS",
            ],
            "chr": ["chr1"] * 7 + ["chr2"],
            "tss": [9_960, 9_930, 9_800, 9_500, 10_300, 9_000, 15_000, 9_960],
            "strand": ["+"] * 8,
        },
        index=pd.Index(
            [NEAR40, NEAR70, FAR200, FAR500, OPPOSITE, AWAY1000, FAR5000, TRANSGENE],
            name="gene_id",
        ),
    )


PANEL = pd.Index([FAR200])


def pair(gene_id):
    name, chrom, start, end = ELEMENT
    return {
        "intended_target_name": name,
        "intended_target_chr": chrom,
        "intended_target_start": start,
        "intended_target_end": end,
        "gene_id": gene_id,
    }


ALL_GENES = [NEAR40, NEAR70, FAR200, FAR500, OPPOSITE, AWAY1000, FAR5000, TRANSGENE]


def classified():
    # Ranking is always against every annotated gene on the chromosome, not
    # only the genes passed in, so one call with the whole set is enough.
    table = classify(pd.DataFrame([pair(g) for g in ALL_GENES]), tss_frame(), PANEL)
    return table.set_index("gene_id")


def test_the_nearest_gene_has_no_nearer_gene():
    row = classified().loc[NEAR40]
    assert row["n_nearer"] == 0
    assert row["category"] == NEAREST


def test_a_nearer_gene_outside_the_panel_is_skip_outside_panel():
    # NEAR40 (dist 40, same side) is nearer than NEAR70 (dist 70) and not in the panel.
    row = classified().loc[NEAR70]
    assert (row["n_nearer"], row["n_skipped"]) == (1, 1)
    assert (row["n_skipped_in_panel"], row["n_skipped_outside_panel"]) == (0, 1)
    assert row["category"] == SKIP_OUTSIDE_PANEL


def test_a_nearer_gene_on_the_opposite_side_is_not_skipped():
    # NEAR40, NEAR70 and FAR200 are all nearer to OPPOSITE than it is to
    # itself, but all three sit on the other side of the element.
    row = classified().loc[OPPOSITE]
    assert row["n_nearer"] == 3
    assert row["n_skipped"] == 0
    assert row["category"] == NO_SKIP


def test_a_nearer_gene_in_the_panel_is_skip_in_panel():
    # FAR200 is in the panel, same side as FAR500, and nearer to the element.
    row = classified().loc[FAR500]
    assert row["n_skipped_in_panel"] == 1
    assert row["n_skipped_outside_panel"] == 2  # NEAR40 and NEAR70
    assert row["category"] == SKIP_IN_PANEL


def test_a_gene_on_another_chromosome_is_trans():
    row = classified().loc[TRANSGENE]
    assert row["category"] == TRANS
    assert pd.isna(row["n_nearer"])


def test_counting_includes_a_far_opposite_gene_once_it_is_the_nearer_one():
    # FAR5000 is nearer than nothing on chr1 except itself is farthest of
    # all, so every other chr1 gene counts as nearer to it.
    row = classified().loc[FAR5000]
    assert row["n_nearer"] == 6
    assert row["n_skipped"] == 1  # only OPPOSITE shares its side
    assert row["category"] == SKIP_OUTSIDE_PANEL


def test_a_pair_gene_missing_a_required_column_is_refused():
    with pytest.raises(SkippedGenesError, match="lacks gene_id"):
        classify(pd.DataFrame([{c: 0 for c in TARGET_COLUMNS}]), tss_frame(), PANEL)


def test_a_pair_gene_absent_from_the_annotation_is_refused():
    with pytest.raises(SkippedGenesError, match="not in the annotation"):
        classify(pd.DataFrame([pair("ENSG00000000077")]), tss_frame(), PANEL)


def test_classify_returns_the_geometry_columns_in_order():
    from crispr_analysis_utils.pairs import TARGET_COLUMNS
    from crispr_analysis_utils.skipped_genes import COUNT_COLUMNS

    table = classify(pd.DataFrame([pair(NEAR40)]), tss_frame(), PANEL)
    assert list(table.columns) == [
        *TARGET_COLUMNS,
        "gene_id",
        *COUNT_COLUMNS,
        "category",
    ]


def test_every_category_is_reachable():
    assert set(classified()["category"]) == set(CATEGORIES)


# --- end to end -------------------------------------------------------------


@pytest.fixture
def inputs(tmp_path):
    store = tmp_path / "store.parquet"
    build_store(write_gtf(tmp_path / "test.gtf"), store)
    panel = tmp_path / "panel.csv"
    pd.DataFrame({"gene_id": ["ENSG00000000001.5"]}).to_csv(panel, index=False)
    guides = pd.DataFrame(
        [
            {
                "guide_id": "t1",
                "targeting": "TRUE",
                "type": "targeting",
                "intended_target_name": "peak1",
                "intended_target_chr": "chr1",
                "intended_target_start": "400",
                "intended_target_end": "600",
                "putative_target_genes": "",
            }
        ]
    )
    pairs_path = tmp_path / "pairs" / "pairs_to_test.csv"
    pairs_path.parent.mkdir()
    from crispr_analysis_utils.guide_annotation.db import load_store
    from crispr_analysis_utils.pairs import gene_tss, read_panel

    built = build_pairs(
        guides, read_panel(panel), gene_tss(load_store(store)), max_distance=5_000
    )
    built.to_csv(pairs_path, index=False)
    return pairs_path, store, panel


def test_run_writes_the_table_and_its_record(tmp_path, inputs):
    pairs_path, store, panel = inputs
    result = classify_run(pairs_path, store, panel, tmp_path / "out")

    table = pd.read_csv(tmp_path / "out" / "skipped_genes" / "skipped_genes.csv")
    assert list(table.columns) == list(COLUMNS)
    assert len(table) == len(result.table)

    written = sorted(
        p.relative_to(tmp_path / "out").as_posix()
        for p in (tmp_path / "out").rglob("*")
        if p.is_file()
    )
    assert [w for w in written if not w.startswith("run_metadata/")] == [
        "logs/skipped_genes.log",
        "skipped_genes/qc.tsv",
        "skipped_genes/skipped_genes.csv",
        "skipped_genes/skipped_genes_columns.yaml",
    ]
    (record,) = [w for w in written if w.startswith("run_metadata/")]
    assert re.fullmatch(r"run_metadata/skipped_genes_\d{8}T\d{6}Z\.json", record)
    state = json.loads((tmp_path / "out" / record).read_text())
    assert state["complete"]


def test_the_command_runs(tmp_path, inputs, capsys):
    pairs_path, store, panel = inputs
    status = cli.main(
        [
            "skipped-genes",
            "run",
            "--pairs",
            str(pairs_path),
            "--store",
            str(store),
            "--genes",
            str(panel),
            "--outdir",
            str(tmp_path / "out"),
        ]
    )
    assert status == 0
    out = capsys.readouterr().out
    assert "pairs" in out and "Outputs in" in out


def test_the_command_reports_a_missing_file(tmp_path, inputs, capsys):
    _, store, panel = inputs
    with pytest.raises(SystemExit):
        cli.main(
            [
                "skipped-genes",
                "run",
                "--pairs",
                str(tmp_path / "nope.csv"),
                "--store",
                str(store),
                "--genes",
                str(panel),
                "--outdir",
                str(tmp_path / "out"),
            ]
        )
    assert "does not exist" in capsys.readouterr().err

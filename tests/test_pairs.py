"""Building the discovery and positive-control pairs, and the ``cau pairs`` command.

Coordinates are 0-based half-open. An element's midpoint is
``(start + end) // 2``, and a TSS is a transcript's first base: its start on
``+``, its end minus one on ``-``.
"""

import json
import re

import pandas as pd
import pytest
from guide_annotation_fixtures import write_gtf

from crispr_analysis_utils import cli
from crispr_analysis_utils.guide_annotation import build_store
from crispr_analysis_utils.guide_annotation.db import Store
from crispr_analysis_utils.pairs import (
    COLUMNS,
    DISCOVERY,
    POSITIVE_CONTROL,
    TABLE,
    PairsError,
    build_pairs,
    gene_tss,
    read_panel,
    run,
)

PLUS = "ENSG00000000001"
MINUS = "ENSG00000000002"
FAR = "ENSG00000000003"
OUTSIDE = "ENSG00000000009"


def tss_frame():
    """Three readout genes and one the panel leaves out."""
    return pd.DataFrame(
        {
            "gene_symbol": ["PLUS", "MINUS", "FAR", "OUTSIDE"],
            "chr": ["chr1", "chr1", "chr2", "chr1"],
            "tss": [10_000, 20_000, 10_000, 30_000],
            "strand": ["+", "-", "+", "+"],
        },
        index=pd.Index([PLUS, MINUS, FAR, OUTSIDE], name="gene_id"),
    )


def guide(
    guide_id, start, end, *, kind="targeting", chrom="chr1", genes="", targeting="TRUE"
):
    """One row of the IGVF per-guide table."""
    return {
        "guide_id": guide_id,
        "targeting": targeting,
        "type": kind,
        "intended_target_name": f"{chrom}:{start}-{end}",
        "intended_target_chr": chrom,
        "intended_target_start": str(start),
        "intended_target_end": str(end),
        "putative_target_genes": genes,
    }


def control(guide_id, kind):
    """A control guide: no element at all."""
    row = guide(guide_id, 0, 0, kind=kind, targeting="FALSE")
    row.update(intended_target_chr="", intended_target_start="", intended_target_end="")
    return row


PANEL = pd.Index([PLUS, MINUS, FAR])


def pairs_of(rows, max_distance=1_000):
    return build_pairs(
        pd.DataFrame(rows), PANEL, tss_frame(), max_distance=max_distance
    )


def keyed(pairs):
    return set(
        zip(pairs["guide_id"], pairs["gene_id"], pairs["pair_type"], strict=True)
    )


def test_discovery_pairs_a_targeting_guide_with_each_readout_gene_in_reach():
    # Midpoint 10,500: 500 past PLUS's TSS, 9,500 short of MINUS's.
    pairs = pairs_of([guide("g1", 10_400, 10_600)])
    assert keyed(pairs) == {("g1", PLUS, DISCOVERY)}


def test_the_distance_limit_is_inclusive():
    # Midpoint 11,000, exactly 1,000 from PLUS's TSS.
    assert keyed(pairs_of([guide("g1", 10_900, 11_100)])) == {("g1", PLUS, DISCOVERY)}
    assert keyed(pairs_of([guide("g1", 10_902, 11_102)])) == set()


def test_discovery_never_crosses_chromosomes():
    # FAR's TSS has the same coordinate on chr2.
    pairs = pairs_of([guide("g1", 10_400, 10_600)])
    assert FAR not in set(pairs["gene_id"])


def test_discovery_uses_only_panel_genes():
    # Midpoint 30,000 is OUTSIDE's TSS, and OUTSIDE is not in the panel.
    assert keyed(pairs_of([guide("g1", 29_900, 30_100)])) == set()


def test_distance_is_signed_in_the_gene_orientation():
    pairs = pairs_of([guide("up", 9_400, 9_600), guide("down", 20_400, 20_600)])
    distance = dict(zip(pairs["guide_id"], pairs["distance_to_tss"], strict=True))
    # 500 before PLUS's TSS on +; 500 past MINUS's TSS in coordinates, which
    # is upstream on -.
    assert distance == {"up": -500, "down": -500}


def test_a_positive_control_pairs_with_the_genes_it_names():
    pairs = pairs_of(
        [
            guide(
                "p1", 10_400, 10_600, kind="positive control", genes=f"{PLUS},{OUTSIDE}"
            )
        ]
    )
    assert keyed(pairs) == {
        ("p1", PLUS, POSITIVE_CONTROL),
        ("p1", OUTSIDE, POSITIVE_CONTROL),
    }
    in_panel = dict(zip(pairs["gene_id"], pairs["in_panel"], strict=True))
    assert in_panel == {PLUS: True, OUTSIDE: False}


def test_a_positive_element_gene_is_positive_for_every_guide_and_other_genes_stay_discovery():
    # The element is within reach of PLUS and MINUS; PLUS is the named gene.
    pairs = pairs_of(
        [
            guide("p1", 14_900, 15_100, kind="positive control", genes=PLUS),
            guide("t1", 14_900, 15_100),
        ],
        max_distance=10_000,
    )
    assert keyed(pairs) == {
        ("p1", PLUS, POSITIVE_CONTROL),
        ("p1", MINUS, DISCOVERY),
        ("t1", PLUS, POSITIVE_CONTROL),
        ("t1", MINUS, DISCOVERY),
    }


def test_a_peak_near_a_tss_is_discovery_not_positive_control():
    # A targeting guide whose element sits on PLUS's TSS, typed targeting.
    assert keyed(pairs_of([guide("t1", 9_900, 10_100, genes=PLUS)])) == {
        ("t1", PLUS, DISCOVERY)
    }


def test_controls_get_no_pairs():
    pairs = pairs_of(
        [control("nt", "non-targeting"), control("safe", "safe-targeting")]
    )
    assert pairs.empty


def test_versions_are_stripped_from_named_genes():
    pairs = pairs_of(
        [guide("p1", 10_400, 10_600, kind="positive control", genes=f"{PLUS}.7")]
    )
    assert list(pairs["gene_id"]) == [PLUS]


def test_columns_are_in_order():
    assert list(pairs_of([guide("g1", 10_400, 10_600)]).columns) == list(COLUMNS)


def test_a_panel_gene_missing_from_the_annotation_is_refused():
    with pytest.raises(PairsError, match="panel genes are not in the annotation"):
        build_pairs(
            pd.DataFrame([guide("g1", 1, 2)]),
            pd.Index(["ENSG00000000077"]),
            tss_frame(),
        )


def test_a_positive_control_gene_missing_from_the_annotation_is_refused():
    rows = [
        guide("p1", 10_400, 10_600, kind="positive control", genes="ENSG00000000077")
    ]
    with pytest.raises(
        PairsError, match="positive-control genes are not in the annotation"
    ):
        pairs_of(rows)


def test_read_panel_takes_the_first_column_of_ensembl_ids(tmp_path):
    path = tmp_path / "panel.csv"
    pd.DataFrame(
        {"symbol": ["A", "B", "A"], "id": [f"{PLUS}.3", MINUS, f"{PLUS}.3"]}
    ).to_csv(path, index=False)
    assert list(read_panel(path)) == [PLUS, MINUS]


def test_read_panel_refuses_a_table_without_ensembl_ids(tmp_path):
    path = tmp_path / "panel.csv"
    pd.DataFrame({"symbol": ["A", "B"]}).to_csv(path, index=False)
    with pytest.raises(PairsError, match="no column of ENSEMBL gene ids"):
        read_panel(path)


def test_gene_tss_uses_the_reference_transcript_not_the_gene_span():
    features = pd.DataFrame(
        [
            # The gene spans 100-900; its reference transcript starts at 300.
            (
                "chr1",
                100,
                900,
                "+",
                "gene",
                "ENSG00000000001",
                "AAA",
                None,
                True,
                "protein_coding",
            ),
            (
                "chr1",
                100,
                900,
                "+",
                "transcript",
                "ENSG00000000001",
                "AAA",
                "ENST1",
                False,
                "protein_coding",
            ),
            (
                "chr1",
                300,
                900,
                "+",
                "transcript",
                "ENSG00000000001",
                "AAA",
                "ENST2",
                True,
                "protein_coding",
            ),
            (
                "chr1",
                100,
                800,
                "-",
                "transcript",
                "ENSG00000000002",
                "BBB",
                "ENST3",
                True,
                "protein_coding",
            ),
            # A gene of another type keeps only its gene row.
            (
                "chr2",
                50,
                60,
                "-",
                "gene",
                "ENSG00000000004",
                "DDD",
                None,
                False,
                "lncRNA",
            ),
        ],
        columns=[
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "Feature",
            "ensg",
            "gene_symbol",
            "enst",
            "is_reference",
            "gene_type",
        ],
    )
    tss = gene_tss(Store(features=features, provenance={}))
    assert tss.loc["ENSG00000000001", "tss"] == 300
    assert tss.loc["ENSG00000000002", "tss"] == 799
    assert tss.loc["ENSG00000000004", "tss"] == 59


# --- end to end -------------------------------------------------------------


@pytest.fixture
def inputs(tmp_path):
    store = tmp_path / "store.parquet"
    build_store(write_gtf(tmp_path / "test.gtf"), store)
    metadata = tmp_path / "guide_metadata.tsv"
    # AAA's TSS is 1000 on chr1 +; the element's midpoint is 1500.
    pd.DataFrame(
        [
            guide("pc", 1_400, 1_600, kind="positive control", genes="ENSG00000000001"),
            guide("t", 1_400, 1_600),
            control("nt", "non-targeting"),
        ]
    ).to_csv(metadata, sep="\t", index=False)
    panel = tmp_path / "panel.csv"
    pd.DataFrame({"gene_id": ["ENSG00000000001.5", "ENSG00000000002.5"]}).to_csv(
        panel, index=False
    )
    return metadata, panel, store


def test_run_writes_the_table_and_its_record(tmp_path, inputs):
    metadata, panel, store = inputs
    result = run(metadata, panel, store, tmp_path / "out", max_distance=5_000)

    table = pd.read_csv(tmp_path / "out" / "pairs" / TABLE)
    assert list(table.columns) == list(COLUMNS)
    assert set(table["in_panel"]) <= {True, False}
    assert keyed(table) == keyed(result.pairs)
    assert ("pc", "ENSG00000000001", POSITIVE_CONTROL) in keyed(table)
    assert ("t", "ENSG00000000001", POSITIVE_CONTROL) in keyed(table)

    written = sorted(
        p.relative_to(tmp_path / "out").as_posix()
        for p in (tmp_path / "out").rglob("*")
        if p.is_file()
    )
    assert [w for w in written if not w.startswith("run_metadata/")] == [
        "logs/pairs.log",
        "pairs/pairs_to_test.csv",
        "pairs/pairs_to_test_columns.yaml",
        "pairs/panel_genes_without_pairs.tsv",
        "pairs/qc.tsv",
    ]
    (record,) = [w for w in written if w.startswith("run_metadata/")]
    assert re.fullmatch(r"run_metadata/pairs_\d{8}T\d{6}Z\.json", record)
    state = json.loads((tmp_path / "out" / record).read_text())
    assert state["complete"] and state["parameters"] == {"max_distance": 5_000}


def test_the_command_runs(tmp_path, inputs, capsys):
    metadata, panel, store = inputs
    status = cli.main(
        [
            "pairs",
            "run",
            "--guide-metadata",
            str(metadata),
            "--genes",
            str(panel),
            "--store",
            str(store),
            "--outdir",
            str(tmp_path / "out"),
        ]
    )
    assert status == 0
    assert "positive_control: 2 pairs" in capsys.readouterr().out


def test_the_command_reports_a_missing_file(tmp_path, inputs, capsys):
    _, panel, store = inputs
    with pytest.raises(SystemExit):
        cli.main(
            [
                "pairs",
                "run",
                "--guide-metadata",
                str(tmp_path / "nope.tsv"),
                "--genes",
                str(panel),
                "--store",
                str(store),
                "--outdir",
                str(tmp_path / "out"),
            ]
        )
    assert "does not exist" in capsys.readouterr().err

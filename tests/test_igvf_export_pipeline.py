"""``igvf_export.run`` end to end, and the ``cau igvf-export`` command.

A run writes into ``outdir`` and nowhere else: ``igvf/`` with the table and
its data dictionary, ``run_metadata/`` with this invocation's record, and
``logs/``.
"""

import json
import logging
import re

import pandas as pd
import pytest
import yaml
from element_mapping_fixtures import write_bed
from guide_annotation_fixtures import write_gtf
from igvf_export_fixtures import (
    ELEMENTS,
    GENE_PAIRS,
    GUIDES,
    WINDOW,
    mapping_frame,
    write_alignment,
    write_genes,
)
from pandas.testing import assert_frame_equal

from crispr_analysis_utils import cli
from crispr_analysis_utils.element_mapping import read_elements, write_mapping_tsv
from crispr_analysis_utils.guide_annotation import build_store
from crispr_analysis_utils.igvf_export import (
    COLUMNS,
    NON,
    POSITIVE,
    SAFE,
    TABLE,
    TARGETING,
    ExportError,
    run,
)

AAA = "ENSG00000000001"
BBB = "ENSG00000000002"

EXPECTED_FILES = [
    "igvf/columns.yaml",
    f"igvf/{TABLE}",
    "logs/igvf_export.log",
    re.compile(r"run_metadata/igvf_export_\d{8}T\d{6}Z\.json"),
]
"""Everything a finished run leaves behind, in the order `written` lists it."""


@pytest.fixture
def store(tmp_path):
    build_store(write_gtf(tmp_path / "test.gtf"), tmp_path / "store.parquet")
    return tmp_path / "store.parquet"


@pytest.fixture
def bed(tmp_path):
    return write_bed(tmp_path / "elements.bed", ELEMENTS)


@pytest.fixture
def alignment(tmp_path):
    return write_alignment(tmp_path / "guides_alignment_metadata.tsv")


@pytest.fixture
def mapping(tmp_path, bed):
    path = tmp_path / "guide_element_mapping.tsv"
    write_mapping_tsv(path, mapping_frame(read_elements(bed)))
    return path


@pytest.fixture
def genes(tmp_path):
    return write_genes(tmp_path / "pairs.tsv")


@pytest.fixture
def inputs(alignment, mapping, bed, store):
    """The four inputs every run takes, in the order `run` takes them."""
    return (alignment, mapping, bed, store)


def written(outdir):
    """Every file under `outdir`, as a relative path."""
    return sorted(
        str(path.relative_to(outdir)) for path in outdir.rglob("*") if path.is_file()
    )


def matches(paths, expected):
    """Each path against its name or its pattern, in order."""
    if len(paths) != len(expected):
        return False
    return all(
        want.fullmatch(got) if isinstance(want, re.Pattern) else want == got
        for got, want in zip(paths, expected, strict=True)
    )


def record_of(outdir):
    """This invocation's metadata."""
    records = sorted((outdir / "run_metadata").glob("igvf_export_*.json"))
    assert len(records) == 1
    return json.loads(records[0].read_text())


def table_of(outdir):
    """The written table, read back as text columns."""
    return pd.read_csv(
        outdir / "igvf" / TABLE, sep="\t", dtype=str, keep_default_na=False
    )


def row_for(table, alias):
    """The one row of a library entry."""
    rows = table[table["guide_id"] == alias]
    assert len(rows) == 1
    return rows.iloc[0]


# --- what a run writes ----------------------------------------------------


def test_a_run_writes_four_files_and_nothing_else(inputs, tmp_path):
    run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    assert matches(written(tmp_path / "out"), EXPECTED_FILES)


def test_nothing_is_written_beside_the_inputs(inputs, genes, tmp_path):
    run(*inputs, tmp_path / "out", genes=genes, promoter_window=WINDOW)
    outside = sorted(
        str(path.relative_to(tmp_path))
        for path in tmp_path.rglob("*")
        if path.is_file() and not path.is_relative_to(tmp_path / "out")
    )
    assert outside == [
        "elements.bed",
        "guide_element_mapping.tsv",
        "guides_alignment_metadata.tsv",
        "pairs.tsv",
        "store.parquet",
        "test.gtf",
    ]


def test_the_outdir_is_made_however_deep(inputs, tmp_path):
    result = run(*inputs, tmp_path / "a" / "b", promoter_window=WINDOW)
    assert result.outdir == tmp_path / "a" / "b"
    assert (tmp_path / "a" / "b" / "igvf" / TABLE).is_file()


def test_the_written_table_is_the_one_the_result_carries(inputs, tmp_path):
    result = run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    written_back = table_of(tmp_path / "out")
    assert list(written_back.columns) == list(COLUMNS) == list(result.table.columns)
    assert written_back["guide_id"].tolist() == result.table["guide_id"].tolist()
    assert len(written_back) == len(GUIDES)


def test_the_table_has_every_column_of_the_data_dictionary(inputs, tmp_path):
    run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    folder = tmp_path / "out" / "igvf"
    parsed = yaml.safe_load((folder / "columns.yaml").read_text())
    assert list(table_of(tmp_path / "out").columns) == list(parsed[TABLE])


def test_the_written_table_holds_no_decimal_point(inputs, tmp_path):
    run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    assert ".0" not in (tmp_path / "out" / "igvf" / TABLE).read_text()


def test_the_run_is_logged_to_a_file_in_the_outdir(inputs, tmp_path):
    run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    log = (tmp_path / "out" / "logs" / "igvf_export.log").read_text()
    assert "IGVF export complete" in log
    assert POSITIVE in log


def test_the_run_log_does_not_reach_the_root_logger(inputs, tmp_path):
    """pyranges calls logging.basicConfig on import, which would print the run's
    log to stderr. The log belongs in the file."""
    root = logging.getLogger()
    seen: list[logging.LogRecord] = []

    class Catch(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            seen.append(record)

    catch = Catch()
    root.addHandler(catch)
    level = root.level
    root.setLevel(logging.INFO)
    try:
        run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    finally:
        root.removeHandler(catch)
        root.setLevel(level)

    assert [record for record in seen if record.name.startswith("crispr_")] == []
    assert (
        "IGVF export complete"
        in (tmp_path / "out" / "logs" / "igvf_export.log").read_text()
    )


# --- the table a run produces ---------------------------------------------


def test_every_entry_of_the_library_is_exported_once(inputs, tmp_path):
    table = table_of(run(*inputs, tmp_path / "out", promoter_window=WINDOW).outdir)
    assert table["guide_id"].tolist() == [alias for alias, *_ in GUIDES]


def test_the_types_are_read_off_the_genome(inputs, tmp_path):
    table = table_of(run(*inputs, tmp_path / "out", promoter_window=WINDOW).outdir)
    kinds = dict(zip(table["guide_id"], table["type"], strict=True))
    assert kinds["promoter"] == POSITIVE
    assert kinds["exon"] == TARGETING
    assert (kinds["safe"], kinds["non_targeting"]) == (SAFE, NON)


def test_a_non_targeting_entry_is_written_with_empty_cells(inputs, tmp_path):
    run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    row = row_for(table_of(tmp_path / "out"), "non_targeting")
    assert row["type"] == NON
    assert row["targeting"] == "FALSE"
    for column in ("guide_chr", "guide_start", "guide_end", "strand", "pam"):
        assert row[column] == ""


def test_a_pam_is_the_genomes_bases(inputs, tmp_path):
    run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    assert row_for(table_of(tmp_path / "out"), "promoter")["pam"] == "TGG"


# --- the genes a run was told about ---------------------------------------


def test_without_a_genes_table_a_promoter_takes_its_nearest_gene(inputs, tmp_path):
    run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    row = row_for(table_of(tmp_path / "out"), "promoter")
    assert row["putative_target_genes"] == AAA
    assert row["intended_target_name"] == AAA


def test_a_genes_table_keeps_every_measured_gene(inputs, genes, tmp_path):
    run(*inputs, tmp_path / "out", genes=genes, promoter_window=WINDOW)
    row = row_for(table_of(tmp_path / "out"), "promoter")
    assert row["putative_target_genes"] == f"{AAA},{BBB}"
    assert row["intended_target_name"] == AAA


def test_the_first_column_of_ensembl_ids_is_the_gene(inputs, tmp_path):
    # A pairs table written for another tool is used as it is: the first
    # column names the element, and the first ENSEMBL column the gene.
    path = write_genes(
        tmp_path / "pairs.tsv",
        [(element, "measured", ensg) for element, ensg in GENE_PAIRS],
        columns=("element_id", "how", "ensembl_id"),
    )
    run(*inputs, tmp_path / "out", genes=path, promoter_window=WINDOW)
    row = row_for(table_of(tmp_path / "out"), "promoter")
    assert row["putative_target_genes"] == f"{AAA},{BBB}"


def test_a_comma_separated_pairs_file_is_read_as_one(inputs, tmp_path):
    path = tmp_path / "pairs.csv"
    pd.DataFrame(list(GENE_PAIRS), columns=["element_id", "ensembl_id"]).to_csv(
        path, index=False
    )
    run(*inputs, tmp_path / "out", genes=path, promoter_window=WINDOW)
    row = row_for(table_of(tmp_path / "out"), "promoter")
    assert row["putative_target_genes"] == f"{AAA},{BBB}"


def test_a_genes_table_without_ensembl_ids_is_refused(inputs, tmp_path):
    path = write_genes(
        tmp_path / "pairs.tsv",
        [(element, "a gene") for element, _ in GENE_PAIRS],
        columns=("element_id", "symbol"),
    )
    with pytest.raises(ExportError, match="no column of ENSEMBL gene ids"):
        run(*inputs, tmp_path / "out", genes=path, promoter_window=WINDOW)


# --- the promoters a run measures against ---------------------------------


def test_without_a_bed_the_stores_reference_starts_are_the_tss(inputs, tmp_path):
    result = run(*inputs, tmp_path / "out", promoter_window=WINDOW)
    kinds = dict(
        zip(
            result.elements["element_id"],
            result.elements["genomic_element"],
            strict=True,
        )
    )
    assert kinds["P1"] == "promoter"


def test_a_tss_bed_stands_in_for_the_stores_starts(inputs, tmp_path):
    # One TSS, 39 kb from P1, so nothing the store would call a promoter is
    # one any more. P1 falls in AAA's first exon instead.
    path = tmp_path / "tss.bed"
    path.write_text(
        "chr1\t40000\t40001\tFAR\t0\t+\tENSG00000000009\n", encoding="utf-8"
    )
    result = run(*inputs, tmp_path / "out", promoters=path, promoter_window=WINDOW)
    kinds = dict(
        zip(
            result.elements["element_id"],
            result.elements["genomic_element"],
            strict=True,
        )
    )
    assert kinds["P1"] == "exon"
    assert row_for(table_of(tmp_path / "out"), "promoter")["type"] == TARGETING


def test_the_window_widens_what_counts_as_a_promoter(inputs, tmp_path):
    narrow = run(*inputs, tmp_path / "narrow", promoter_window=WINDOW)
    wide = run(*inputs, tmp_path / "wide", promoter_window=2000)
    assert narrow.table.loc[narrow.table["guide_id"] == "exon", "type"].item() == (
        TARGETING
    )
    assert wide.table.loc[wide.table["guide_id"] == "exon", "type"].item() == POSITIVE
    assert wide.metadata["parameters"] == {"promoter_window": 2000}


# --- the record of the run ------------------------------------------------


def test_the_record_says_the_run_finished_and_what_it_found(inputs, genes, tmp_path):
    result = run(*inputs, tmp_path / "out", genes=genes, promoter_window=WINDOW)
    state = record_of(tmp_path / "out")
    assert state == result.metadata
    assert state["step"] == "igvf_export"
    assert state["complete"] is True
    assert state["error"] is None
    assert state["parameters"] == {"promoter_window": WINDOW}
    assert state["counts"]["guides"] == len(GUIDES) == len(result.table)
    assert state["counts"]["type"] == {
        TARGETING: 7,
        SAFE: 1,
        NON: 1,
        POSITIVE: 2,
    }
    assert state["counts"]["elements"] == len(ELEMENTS)
    assert state["counts"]["with_a_gene"] == 3
    assert state["counts"]["missing_pam"] == 0


def test_the_record_carries_the_inputs_and_their_checksums(inputs, tmp_path):
    run(*(str(path) for path in inputs), tmp_path / "out", promoter_window=WINDOW)
    state = record_of(tmp_path / "out")
    alignment, mapping, bed, store = inputs
    assert state["inputs"]["alignment"] == str(alignment)
    assert state["inputs"]["mapping"] == str(mapping)
    assert state["inputs"]["elements"] == str(bed)
    assert state["inputs"]["store"] == str(store)
    assert state["inputs"]["genes"] is None
    assert len(state["inputs"]["alignment_md5"]) == 32
    assert len(state["inputs"]["elements_md5"]) == 32
    assert state["cau_version"]


def test_a_failed_run_still_leaves_its_record(mapping, bed, store, tmp_path):
    with pytest.raises(FileNotFoundError):
        run(tmp_path / "nowhere.tsv", mapping, bed, store, tmp_path / "out")
    state = record_of(tmp_path / "out")
    assert state["complete"] is False
    assert state["error"].startswith("FileNotFoundError: ")
    assert state["counts"] == {}
    # The table is never written, so no half-finished results are left behind.
    assert matches(
        written(tmp_path / "out"), ["logs/igvf_export.log", EXPECTED_FILES[-1]]
    )


def test_a_second_run_is_the_same_table(inputs, genes, tmp_path):
    first = run(*inputs, tmp_path / "one", genes=genes, promoter_window=WINDOW)
    second = run(*inputs, tmp_path / "two", genes=genes, promoter_window=WINDOW)
    assert_frame_equal(first.table, second.table)


# --- the cau command ------------------------------------------------------


def test_igvf_export_is_registered():
    assert cli.COMMANDS["igvf-export"][0] == "igvf_export"


def test_parse_run_defaults():
    args = cli.build_parser().parse_args(
        [
            "igvf-export",
            "run",
            "--alignment",
            "aln.tsv",
            "--mapping",
            "map.tsv",
            "--elements",
            "e.bed",
            "--store",
            "store.parquet",
            "--outdir",
            "out",
        ]
    )
    assert (args.command, args.action) == ("igvf-export", "run")
    assert (args.alignment, args.mapping) == ("aln.tsv", "map.tsv")
    assert (args.elements, args.store, args.outdir) == ("e.bed", "store.parquet", "out")
    assert (args.promoters, args.genes) == (None, None)
    assert args.promoter_window == 2000


def test_the_optional_arguments_are_parsed_when_given():
    args = cli.build_parser().parse_args(
        [
            "igvf-export",
            "run",
            "--alignment",
            "aln.tsv",
            "--mapping",
            "map.tsv",
            "--elements",
            "e.bed",
            "--store",
            "store.parquet",
            "--outdir",
            "out",
            "--promoters",
            "tss.bed",
            "--genes",
            "pairs.tsv",
            "--promoter-window",
            "500",
        ]
    )
    assert (args.promoters, args.genes) == ("tss.bed", "pairs.tsv")
    assert args.promoter_window == 500


@pytest.mark.parametrize(
    "argv",
    [
        ["igvf-export"],
        ["igvf-export", "export"],
        ["igvf-export", "run"],
        ["igvf-export", "run", "--alignment", "aln.tsv"],
        [
            "igvf-export",
            "run",
            "--alignment",
            "aln.tsv",
            "--mapping",
            "map.tsv",
            "--elements",
            "e.bed",
            "--store",
            "store.parquet",
        ],
    ],
)
def test_bad_arguments_are_usage_errors(argv):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(argv)
    assert exit_info.value.code == 2


def command(inputs, **extra):
    """``cau igvf-export run`` over the fixture inputs."""
    alignment, mapping, bed, store = inputs
    argv = [
        "igvf-export",
        "run",
        "--alignment",
        str(alignment),
        "--mapping",
        str(mapping),
        "--elements",
        str(bed),
        "--store",
        str(store),
        "--outdir",
        "out",
        "--promoter-window",
        str(WINDOW),
    ]
    for name, value in extra.items():
        argv += [f"--{name}", str(value)]
    return argv


def test_a_run_from_the_command_line(inputs, tmp_path, capsys):
    status = cli.main(command(inputs))

    assert status == 0
    assert capsys.readouterr().out.splitlines() == [
        "11 guides: 7 targeting, 1 safe-targeting, 1 non-targeting, 2 positive control",
        "elements: 1 promoter, 1 enhancer, 4 distal element, 2 exon, 1 intron",
        "Outputs in out",
    ]
    assert matches(written(tmp_path / "out"), EXPECTED_FILES)


@pytest.mark.parametrize(
    ("option", "message"),
    [
        pytest.param("alignment", "nowhere.tsv", id="the_alignment"),
        pytest.param("mapping", "nowhere.tsv", id="the_mapping"),
        pytest.param("genes", "nowhere.tsv", id="the_measured_genes"),
    ],
)
def test_a_missing_input_is_reported_without_a_traceback(
    inputs, tmp_path, capsys, option, message
):
    argv = command(inputs)
    if option == "genes":
        argv += ["--genes", "nowhere.tsv"]
    else:
        argv[argv.index(f"--{option}") + 1] = "nowhere.tsv"

    with pytest.raises(SystemExit) as exit_info:
        cli.main(argv)

    assert exit_info.value.code == 1
    error = capsys.readouterr().err
    assert error.startswith("cau: error: ")
    assert message in error
    assert "Traceback" not in error

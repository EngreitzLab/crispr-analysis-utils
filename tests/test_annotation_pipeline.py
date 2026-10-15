"""``guide_annotation.run`` end to end, the table it writes, and the ``cau`` command.

A run writes into ``outdir`` and nowhere else: ``annotations/`` with the one
table and its data dictionary, ``run_metadata/`` with this invocation's
record, and ``logs/``.
"""

import json
import re

import pandas as pd
import pytest
import yaml
from guide_annotation_fixtures import (
    site,
    sites_frame,
    write_genes_only_gtf,
    write_gtf,
)

from crispr_analysis_utils import cli
from crispr_analysis_utils.guide_annotation import (
    AnnotationError,
    build_store,
    run,
)
from crispr_analysis_utils.guide_annotation.summary import (
    COLUMNS,
    GROUPS,
    TABLE,
    annotation_columns,
    counts,
    write_annotation_tsv,
    write_columns_yaml,
)

EXPECTED_FILES = [
    f"annotations/{TABLE}",
    "annotations/guides_gene_annotation_columns.yaml",
    "logs/guide_annotation.log",
    re.compile(r"run_metadata/guide_annotation_\d{8}T\d{6}Z\.json"),
]

SITES = (
    site("exon", "chr1", 2030, "+"),
    site("two_genes", "chr1", 2600, "+"),
    site("intergenic", "chr1", 20000, "+"),
    site("no_pam", "chr1", 1700, "+", pam_matches=False),
)


@pytest.fixture
def store(tmp_path):
    return build_store(write_gtf(tmp_path / "test.gtf"), tmp_path / "store.parquet")


def write_alignment(path, rows=SITES, *, as_text=True):
    """``guides_alignment_metadata.tsv``, with pam_matches written as a table has it."""
    frame = pd.DataFrame(list(rows))
    if as_text:
        frame["pam_matches"] = frame["pam_matches"].map({True: "TRUE", False: "FALSE"})
    frame.to_csv(path, sep="\t", index=False)
    return path


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
    records = sorted((outdir / "run_metadata").glob("guide_annotation_*.json"))
    assert len(records) == 1
    return json.loads(records[0].read_text())


# --- what a run writes ----------------------------------------------------


def test_a_run_writes_four_files_and_nothing_else(store, tmp_path):
    run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "out")
    assert matches(written(tmp_path / "out"), EXPECTED_FILES)


def test_nothing_is_written_beside_the_inputs(store, tmp_path):
    run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "out")
    outside = [
        path.relative_to(tmp_path)
        for path in tmp_path.rglob("*")
        if path.is_file() and not path.is_relative_to(tmp_path / "out")
    ]
    assert sorted(str(path) for path in outside) == [
        "aln.tsv",
        "store.parquet",
        "test.gtf",
    ]


def test_the_outdir_is_made_however_deep(store, tmp_path):
    result = run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "a" / "b")
    assert result.outdir == tmp_path / "a" / "b"
    assert (tmp_path / "a" / "b" / "annotations" / TABLE).is_file()


def test_the_table_has_every_column_of_the_data_dictionary(store, tmp_path):
    run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "out")
    table = pd.read_csv(tmp_path / "out" / "annotations" / TABLE, sep="\t")
    assert list(table.columns) == annotation_columns()


def test_the_run_is_logged_to_a_file_in_the_outdir(store, tmp_path, capsys):
    run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "out")
    log = (tmp_path / "out" / "logs" / "guide_annotation.log").read_text()
    assert "Gene annotation complete" in log
    assert "mapped sites to annotate" in log
    assert "Gene annotation complete" not in capsys.readouterr().err


# --- only the sites that mapped -------------------------------------------


@pytest.mark.parametrize("as_text", [True, False])
def test_only_the_sites_with_a_pam_are_annotated(store, tmp_path, as_text):
    # A written table holds TRUE and FALSE, which pandas reads back as bools,
    # and a caller may hand over the booleans directly.
    result = run(
        write_alignment(tmp_path / "aln.tsv", as_text=as_text), store, tmp_path / "out"
    )
    assert "no_pam" not in result.annotated["guide_id"].tolist()
    assert result.metadata["counts"]["sites"] == 3


@pytest.mark.parametrize(
    "value", [True, "TRUE", "true", "True", " TRUE ", pytest.param(1, id="one")]
)
def test_every_spelling_of_true_maps(store, tmp_path, value):
    frame = pd.DataFrame([site("g", "chr1", 2030, "+")])
    frame["pam_matches"] = [value]
    assert len(run(frame, store, tmp_path / "out").annotated) == 1


@pytest.mark.parametrize("value", [False, "FALSE", "false", pytest.param(0, id="zero")])
def test_every_spelling_of_false_does_not(store, tmp_path, value):
    frame = pd.DataFrame([site("g", "chr1", 2030, "+")])
    frame["pam_matches"] = [value]
    with pytest.raises(AnnotationError, match="no mapped site"):
        run(frame, store, tmp_path / "out")


def test_a_table_without_a_pam_column_annotates_every_row(store, tmp_path):
    frame = sites_frame(*SITES).drop(columns=["pam_matches"])
    assert run(frame, store, tmp_path / "out").metadata["counts"]["sites"] == 4


def test_nothing_mapped_is_an_annotation_error(store, tmp_path):
    rows = [site("no_pam", "chr1", 1700, "+", pam_matches=False)]
    with pytest.raises(AnnotationError, match="holds no mapped site"):
        run(write_alignment(tmp_path / "aln.tsv", rows), store, tmp_path / "out")


def test_a_failed_run_still_leaves_its_record(store, tmp_path):
    rows = [site("no_pam", "chr1", 1700, "+", pam_matches=False)]
    with pytest.raises(AnnotationError):
        run(write_alignment(tmp_path / "aln.tsv", rows), store, tmp_path / "out")
    state = record_of(tmp_path / "out")
    assert state["complete"] is False
    assert state["error"].startswith("AnnotationError: ")
    assert state["counts"] == {}
    # The table is never written, so no half-finished results are left behind.
    assert matches(
        written(tmp_path / "out"),
        ["logs/guide_annotation.log", EXPECTED_FILES[-1]],
    )


def test_a_missing_alignment_table_is_a_file_not_found_error(store, tmp_path):
    with pytest.raises(FileNotFoundError, match="alignment table"):
        run(tmp_path / "nowhere.tsv", store, tmp_path / "out")


# --- the record of the run ------------------------------------------------


def test_the_record_says_the_run_finished_and_what_it_found(store, tmp_path):
    result = run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "out")
    state = record_of(tmp_path / "out")
    assert state == result.metadata
    assert state["step"] == "guide_annotation"
    assert state["complete"] is True
    assert state["error"] is None
    assert state["counts"]["rows"] == len(result.annotated)
    assert state["counts"]["location"] == {"exon": 2, "intron": 1, "intergenic": 1}
    assert state["counts"]["sites"] == 3
    assert state["counts"]["genes"] == 2


def test_the_record_carries_the_inputs_and_the_stores_provenance(store, tmp_path):
    alignment = write_alignment(tmp_path / "aln.tsv")
    run(alignment, str(tmp_path / "store.parquet"), tmp_path / "out")
    state = record_of(tmp_path / "out")
    assert state["inputs"]["alignment"] == str(alignment)
    assert len(state["inputs"]["alignment_md5"]) == 32
    assert state["store"]["gene_types"] == ["protein_coding"]
    assert state["parameters"] == {"promoters": None, "promoter_window": 1000}
    assert set(state["seconds"]) >= {"store", "sites", "promoters", "annotate", "total"}


def test_a_store_handed_over_in_memory_is_recorded_as_such(store, tmp_path):
    run(sites_frame(*SITES), store, tmp_path / "out")
    state = record_of(tmp_path / "out")
    assert state["inputs"] == {"store": "Store", "alignment": "DataFrame"}


# --- the promoters a run measures from ------------------------------------


def test_without_a_promoter_file_the_store_supplies_the_tss(store, tmp_path):
    run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "out")
    assert (
        "4 TSS from the store"
        in (tmp_path / "out" / "logs" / "guide_annotation.log").read_text()
    )


def test_a_promoter_file_is_read_and_recorded(store, tmp_path):
    bed = tmp_path / "tss.bed"
    bed.write_text("chr1\t2500\t2501\tXXX\t0\t+\n", encoding="utf-8")
    result = run(
        write_alignment(tmp_path / "aln.tsv"),
        store,
        tmp_path / "out",
        promoters=bed,
        promoter_window=50,
    )
    state = record_of(tmp_path / "out")
    assert state["parameters"] == {"promoters": str(bed), "promoter_window": 50}
    assert len(state["inputs"]["promoters_md5"]) == 32
    # Every distance is now measured from the one TSS the file gives.
    assert set(result.annotated["nearest_tss_gene"]) == {"XXX"}


# --- the table and its data dictionary ------------------------------------


def test_the_groups_account_for_every_column():
    assert sum(count for _, count in GROUPS) == len(COLUMNS)


def test_columns_yaml_parses_and_documents_the_table(store, tmp_path):
    run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "out")
    folder = tmp_path / "out" / "annotations"
    parsed = yaml.safe_load(
        (folder / "guides_gene_annotation_columns.yaml").read_text()
    )
    header = (folder / TABLE).read_text().splitlines()[0].split("\t")
    assert list(parsed) == [TABLE]
    assert list(parsed[TABLE]) == annotation_columns() == header
    assert all(parsed[TABLE].values())


def test_columns_yaml_describes_no_column_the_table_lacks(tmp_path):
    write_columns_yaml(tmp_path / "columns.yaml")
    parsed = yaml.safe_load((tmp_path / "columns.yaml").read_text())
    assert set(parsed[TABLE]) == set(COLUMNS)


def test_columns_yaml_keeps_the_group_headings(tmp_path):
    write_columns_yaml(tmp_path / "columns.yaml")
    text = (tmp_path / "columns.yaml").read_text()
    for heading, _ in GROUPS:
        assert f"  # {heading}" in text


def test_the_table_writes_flags_as_words_and_counts_as_whole_numbers(store, tmp_path):
    result = run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "out")
    lines = (tmp_path / "out" / "annotations" / TABLE).read_text().splitlines()
    header = lines[0].split("\t")
    cells = [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]
    # 147.0 in the frame, 147 in the file, and an empty cell where there is none.
    assert result.annotated.loc[0, "cds_position"] == 147
    assert [row["cds_position"] for row in cells] == ["147", "", "132", ""]
    assert [row["exonic_any_transcript"] for row in cells] == [
        "TRUE",
        "FALSE",
        "TRUE",
        "FALSE",
    ]


def test_a_column_annotate_sites_never_made_is_written_empty(tmp_path):
    bare = pd.DataFrame(
        {
            "chr": ["chr1"],
            "start": [10],
            "end": [30],
            "strand": ["+"],
            "location": ["intergenic"],
            "region": [""],
            "sense": [""],
            "ensg": [None],
            "promoter": [False],
            "exonic_any_transcript": [False],
        }
    )
    write_annotation_tsv(tmp_path / "bare.tsv", bare)
    lines = (tmp_path / "bare.tsv").read_text().splitlines()
    assert lines[0].split("\t") == annotation_columns()
    assert lines[1].split("\t")[0] == ""
    assert counts(bare)["guides"] == 0


def test_counts_sum_the_rows_not_the_sites(store, tmp_path):
    result = run(write_alignment(tmp_path / "aln.tsv"), store, tmp_path / "out")
    found = counts(result.annotated)
    assert sum(found["location"].values()) == found["rows"] == 4
    assert found["sites"] == 3
    assert found["guides"] == 3
    assert found["intergenic"] == 1


# --- the cau command ------------------------------------------------------


def test_guide_annotation_is_registered():
    assert cli.COMMANDS["guide-annotation"][0] == "guide_annotation"


def test_parse_build_db_defaults():
    args = cli.build_parser().parse_args(
        ["guide-annotation", "build-db", "gencode.gtf", "--out", "store.parquet"]
    )
    assert (args.command, args.action) == ("guide-annotation", "build-db")
    assert (args.gtf, args.out) == ("gencode.gtf", "store.parquet")
    assert args.gene_types == "protein_coding"


def test_parse_run_defaults():
    args = cli.build_parser().parse_args(
        [
            "guide-annotation",
            "run",
            "--alignment",
            "aln.tsv",
            "--store",
            "store.parquet",
            "--outdir",
            "out",
        ]
    )
    assert args.action == "run"
    assert (args.alignment, args.store, args.outdir) == (
        "aln.tsv",
        "store.parquet",
        "out",
    )
    assert args.promoters is None
    assert args.promoter_window == 1000


@pytest.mark.parametrize(
    "argv",
    [
        ["guide-annotation"],
        ["guide-annotation", "annotate"],
        ["guide-annotation", "build-db", "gencode.gtf"],
        ["guide-annotation", "run", "--alignment", "aln.tsv", "--store", "s.parquet"],
        ["guide-annotation", "run", "--outdir", "out"],
        [
            "guide-annotation",
            "run",
            "--alignment",
            "aln.tsv",
            "--store",
            "s.parquet",
            "--outdir",
            "out",
            "--promoter-window",
            "wide",
        ],
    ],
)
def test_bad_arguments_are_usage_errors(argv):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(argv)
    assert exit_info.value.code == 2


def test_build_db_then_run_from_the_command_line(tmp_path, capsys):
    gtf = write_gtf(tmp_path / "test.gtf")
    alignment = write_alignment(tmp_path / "aln.tsv")

    assert (
        cli.main(["guide-annotation", "build-db", str(gtf), "--out", "store.parquet"])
        == 0
    )
    assert capsys.readouterr().out.splitlines() == [
        "4 genes, 9 transcript, 21 exon, 13 CDS",
        "Store in store.parquet",
    ]

    status = cli.main(
        [
            "guide-annotation",
            "run",
            "--alignment",
            str(alignment),
            "--store",
            "store.parquet",
            "--outdir",
            "out",
        ]
    )

    assert status == 0
    assert capsys.readouterr().out.splitlines() == [
        "4 rows for 3 sites in 2 genes: 2 exon, 1 intron, 1 intergenic, 0 in a promoter",
        "Outputs in out",
    ]
    assert matches(written(tmp_path / "out"), EXPECTED_FILES)


def test_gene_types_is_split_on_commas(tmp_path, capsys):
    gtf = write_gtf(tmp_path / "test.gtf")
    status = cli.main(
        [
            "guide-annotation",
            "build-db",
            str(gtf),
            "--out",
            "store.parquet",
            "--gene-types",
            "protein_coding, lncRNA",
        ]
    )
    assert status == 0
    assert capsys.readouterr().out.startswith("5 genes")


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (
            ["guide-annotation", "build-db", "nowhere.gtf", "--out", "store.parquet"],
            "The GTF nowhere.gtf does not exist",
        ),
        (
            [
                "guide-annotation",
                "build-db",
                "test.gtf",
                "--out",
                "store.parquet",
                "--gene-types",
                " , ",
            ],
            "--gene-types was given no gene type",
        ),
        (
            [
                "guide-annotation",
                "run",
                "--alignment",
                "nowhere.tsv",
                "--store",
                "store.parquet",
                "--outdir",
                "out",
            ],
            "The alignment table nowhere.tsv does not exist",
        ),
        (
            [
                "guide-annotation",
                "run",
                "--alignment",
                "aln.tsv",
                "--store",
                "nowhere.parquet",
                "--outdir",
                "out",
            ],
            "The annotation store nowhere.parquet does not exist",
        ),
    ],
)
def test_a_missing_file_is_reported_without_a_traceback(
    tmp_path, capsys, argv, message
):
    write_gtf(tmp_path / "test.gtf")
    write_alignment(tmp_path / "aln.tsv")
    build_store(tmp_path / "test.gtf", tmp_path / "store.parquet")

    with pytest.raises(SystemExit) as exit_info:
        cli.main(argv)

    assert exit_info.value.code == 1
    error = capsys.readouterr().err
    assert error.startswith("cau: error: ")
    assert message in error
    assert "Traceback" not in error


def test_a_genes_only_gtf_is_reported_without_a_traceback(tmp_path, capsys):
    write_genes_only_gtf(tmp_path / "genes.gtf")

    with pytest.raises(SystemExit) as exit_info:
        cli.main(
            ["guide-annotation", "build-db", "genes.gtf", "--out", "store.parquet"]
        )

    assert exit_info.value.code == 1
    error = capsys.readouterr().err
    assert "cau: error: " in error
    assert "holds no exon rows" in error
    assert "Traceback" not in error

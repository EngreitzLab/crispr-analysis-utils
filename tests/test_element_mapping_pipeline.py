"""``element_mapping.run`` end to end, the table it writes, and the ``cau`` command.

A run writes into ``outdir`` and nowhere else: ``elements/`` with the one
table and its data dictionary, ``run_metadata/`` with this invocation's
record, and ``logs/``.
"""

import json
import logging
import re

import pandas as pd
import pytest
import yaml
from element_mapping_fixtures import (
    annotation_frame,
    site,
    sites_frame,
    unmapped,
    write_alignment,
    write_bed,
)

from crispr_analysis_utils import cli
from crispr_analysis_utils.element_mapping import (
    COLUMNS,
    NON_TARGETING,
    SAFE_TARGETING,
    TABLE,
    TARGETING,
    ElementError,
    counts,
    run,
)

EXPECTED_FILES = [
    "elements/columns.yaml",
    f"elements/{TABLE}",
    "elements/qc.tsv",
    "logs/element_mapping.log",
    re.compile(r"run_metadata/element_mapping_\d{8}T\d{6}Z\.json"),
]

ROWS = (
    # Inside E1.
    site("in_one", "chr1", 1000),
    # One guide, two sites: inside E3, and in no element.
    site("split", "chr2", 2000),
    site("split", "chr1", 5000),
    # The same bases as "in_one", but without a PAM, so it never maps.
    site("no_pam", "chr1", 1000, pam_matches=False),
    # In the library and nowhere in the genome.
    unmapped("never"),
)
"""Four guides, three mapped sites: two targeting, one safe, two non-targeting."""


@pytest.fixture
def alignment(tmp_path):
    return write_alignment(tmp_path / "aln.tsv", ROWS)


@pytest.fixture
def bed(tmp_path):
    return write_bed(tmp_path / "elements.bed")


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
    records = sorted((outdir / "run_metadata").glob("element_mapping_*.json"))
    assert len(records) == 1
    return json.loads(records[0].read_text())


def classes_of(mapping):
    """Each guide's classes, so a guide with two sites shows both."""
    return {
        guide_id: sorted(set(rows["target_class"]))
        for guide_id, rows in mapping.groupby("guide_id")
    }


# --- what a run writes ----------------------------------------------------


def test_a_run_writes_five_files_and_nothing_else(alignment, bed, tmp_path):
    run(alignment, bed, tmp_path / "out")
    assert matches(written(tmp_path / "out"), EXPECTED_FILES)


def test_nothing_is_written_beside_the_inputs(alignment, bed, tmp_path):
    run(alignment, bed, tmp_path / "out")
    outside = [
        path.relative_to(tmp_path)
        for path in tmp_path.rglob("*")
        if path.is_file() and not path.is_relative_to(tmp_path / "out")
    ]
    assert sorted(str(path) for path in outside) == ["aln.tsv", "elements.bed"]


def test_the_outdir_is_made_however_deep(alignment, bed, tmp_path):
    result = run(alignment, bed, tmp_path / "a" / "b")
    assert result.outdir == tmp_path / "a" / "b"
    assert (tmp_path / "a" / "b" / "elements" / TABLE).is_file()


def test_the_table_has_every_column_of_the_data_dictionary(alignment, bed, tmp_path):
    run(alignment, bed, tmp_path / "out")
    folder = tmp_path / "out" / "elements"
    table = pd.read_csv(folder / TABLE, sep="\t")
    parsed = yaml.safe_load((folder / "columns.yaml").read_text())
    assert list(table.columns) == list(COLUMNS) == list(parsed[TABLE])


def test_the_run_is_logged_to_a_file_in_the_outdir(alignment, bed, tmp_path):
    run(alignment, bed, tmp_path / "out")
    log = (tmp_path / "out" / "logs" / "element_mapping.log").read_text()
    assert "3 elements, 4 guides, 3 mapped sites" in log
    assert "Element mapping complete" in log


def test_the_written_table_holds_no_decimal_point(alignment, bed, tmp_path):
    path = tmp_path / "ann.tsv"
    annotation_frame((site("in_one", "chr1", 1000), 12)).to_csv(
        path, sep="\t", index=False
    )
    run(alignment, bed, tmp_path / "out", annotation=path)
    text = (tmp_path / "out" / "elements" / TABLE).read_text()
    assert ".0" not in text
    assert "in_one\tchr1\t1000\t1020\t+\tE1\tenhancer\ttargeting\t12" in text
    # The guide that matches nowhere keeps every site column empty.
    assert "never\t\t\t\t\t\t\tnon_targeting\t" in text


# --- the three classes, end to end ----------------------------------------


def test_each_guide_gets_the_class_of_each_of_its_sites(alignment, bed, tmp_path):
    result = run(alignment, bed, tmp_path / "out")
    assert classes_of(result.mapping) == {
        "in_one": [TARGETING],
        "never": [NON_TARGETING],
        "no_pam": [NON_TARGETING],
        "split": [SAFE_TARGETING, TARGETING],
    }


def test_a_library_that_maps_nowhere_is_all_non_targeting(bed, tmp_path):
    result = run(sites_frame(unmapped("a"), unmapped("b")), bed, tmp_path / "out")
    assert result.mapping["target_class"].tolist() == [NON_TARGETING] * 2
    assert result.metadata["complete"] is True
    assert result.metadata["counts"]["elements_with_a_guide"] == 0


def test_elements_on_other_contigs_are_refused(alignment, tmp_path):
    bed = write_bed(tmp_path / "other.bed", [("1", 1000, 1100, "E1")], columns=4)
    with pytest.raises(ElementError, match="No contig of the sites appears"):
        run(alignment, bed, tmp_path / "out")


# --- only the sites with a PAM --------------------------------------------


@pytest.mark.parametrize("as_text", [True, False])
def test_a_site_without_a_pam_is_dropped_but_its_guide_is_kept(bed, tmp_path, as_text):
    # A written table holds TRUE and FALSE, which pandas reads back as bools,
    # and a caller may hand the booleans over directly.
    path = write_alignment(tmp_path / "aln.tsv", ROWS, as_text=as_text)
    result = run(path, bed, tmp_path / "out")
    found = result.metadata["counts"]
    assert found["guides"] == 4
    assert found["classes"] == {TARGETING: 2, SAFE_TARGETING: 1, NON_TARGETING: 2}
    # "no_pam" sits on the same bases as "in_one" and is still not targeting.
    assert classes_of(result.mapping)["no_pam"] == [NON_TARGETING]


@pytest.mark.parametrize(
    "value", [True, "TRUE", "true", "True", " TRUE ", pytest.param(1, id="one")]
)
def test_every_spelling_of_true_maps(bed, tmp_path, value):
    frame = sites_frame(site("g", "chr1", 1000))
    frame["pam_matches"] = [value]
    assert run(frame, bed, tmp_path / "out").mapping["target_class"].tolist() == [
        TARGETING
    ]


@pytest.mark.parametrize(
    "value",
    [
        False,
        "FALSE",
        "false",
        pytest.param(0, id="zero"),
        pytest.param("", id="empty"),
    ],
)
def test_every_spelling_of_false_does_not(bed, tmp_path, value):
    frame = sites_frame(site("g", "chr1", 1000))
    frame["pam_matches"] = [value]
    assert run(frame, bed, tmp_path / "out").mapping["target_class"].tolist() == [
        NON_TARGETING
    ]


def test_an_empty_pam_cell_does_not_count_as_mapped(bed, tmp_path):
    path = tmp_path / "aln.tsv"
    path.write_text(
        "guide_id\tchr\tstart\tend\tstrand\tpam_matches\n"
        "yes\tchr1\t1000\t1020\t+\tTRUE\n"
        "no\tchr1\t1000\t1020\t+\tFALSE\n"
        "blank\tchr1\t1000\t1020\t+\t\n",
        encoding="utf-8",
    )
    assert classes_of(run(path, bed, tmp_path / "out").mapping) == {
        "yes": [TARGETING],
        "no": [NON_TARGETING],
        "blank": [NON_TARGETING],
    }


def test_a_table_without_a_pam_column_maps_every_row(bed, tmp_path):
    frame = sites_frame(*ROWS).drop(columns=["pam_matches"])
    result = run(frame, bed, tmp_path / "out")
    assert classes_of(result.mapping)["no_pam"] == [TARGETING]
    assert result.metadata["counts"]["classes"][TARGETING] == 3


# --- the distance the annotation carries ----------------------------------


def test_without_an_annotation_the_distance_is_empty_and_said_to_be(
    alignment, bed, tmp_path
):
    result = run(alignment, bed, tmp_path / "out")
    assert result.mapping["distance_to_tss"].isna().all()
    log = (tmp_path / "out" / "logs" / "element_mapping.log").read_text()
    assert "No annotation given" in log


def test_an_annotation_file_fills_the_distance_of_every_class(alignment, bed, tmp_path):
    path = tmp_path / "ann.tsv"
    annotation_frame(
        (site("in_one", "chr1", 1000), 12),
        (site("split", "chr2", 2000), -4000),
    ).to_csv(path, sep="\t", index=False)

    result = run(alignment, bed, tmp_path / "out", annotation=path)

    found = result.mapping.set_index(["guide_id", "chr"])["distance_to_tss"]
    # A targeting site gets its distance as readily as a safe one.
    assert found[("in_one", "chr1")] == 12
    assert found[("split", "chr2")] == -4000
    # The site the annotation never saw keeps an empty cell.
    assert pd.isna(found[("split", "chr1")])
    assert result.metadata["inputs"]["annotation"] == str(path)


def test_the_annotation_leaves_every_class_alone(alignment, bed, tmp_path):
    path = tmp_path / "ann.tsv"
    annotation_frame((site("in_one", "chr1", 1000), 12)).to_csv(
        path, sep="\t", index=False
    )
    without = run(alignment, bed, tmp_path / "a").metadata["counts"]
    with_it = run(alignment, bed, tmp_path / "b", annotation=path).metadata["counts"]
    assert with_it == without


# --- the record of the run ------------------------------------------------


def test_the_record_says_the_run_finished_and_what_it_found(alignment, bed, tmp_path):
    result = run(alignment, bed, tmp_path / "out")
    state = record_of(tmp_path / "out")
    assert state == result.metadata
    assert state["step"] == "element_mapping"
    assert state["complete"] is True
    assert state["error"] is None
    assert state["parameters"] == {}
    assert state["counts"] == counts(result.mapping)
    assert state["counts"]["rows"] == len(result.mapping) == 5
    # Three mapped sites; the two guides that match nowhere have none.
    assert state["counts"]["sites"] == 3
    assert state["counts"]["classes"] == {
        TARGETING: 2,
        SAFE_TARGETING: 1,
        NON_TARGETING: 2,
    }
    assert state["counts"]["elements_with_a_guide"] == 2


def test_the_record_carries_the_inputs_and_their_checksums(alignment, bed, tmp_path):
    run(str(alignment), str(bed), tmp_path / "out")
    state = record_of(tmp_path / "out")
    assert state["inputs"]["alignment"] == str(alignment)
    assert state["inputs"]["elements"] == str(bed)
    assert len(state["inputs"]["alignment_md5"]) == 32
    assert len(state["inputs"]["elements_md5"]) == 32


def test_a_table_handed_over_in_memory_is_recorded_as_such(bed, tmp_path):
    run(sites_frame(*ROWS), bed, tmp_path / "out")
    state = record_of(tmp_path / "out")
    assert state["inputs"]["alignment"] == "DataFrame"
    assert "alignment_md5" not in state["inputs"]


def test_a_failed_run_still_leaves_its_record(alignment, tmp_path):
    with pytest.raises(ElementError, match="does not exist"):
        run(alignment, tmp_path / "nowhere.bed", tmp_path / "out")
    state = record_of(tmp_path / "out")
    assert state["complete"] is False
    assert state["error"].startswith("ElementError: ")
    assert state["counts"] == {}
    # The table is never written, so no half-finished results are left behind.
    assert matches(
        written(tmp_path / "out"), ["logs/element_mapping.log", EXPECTED_FILES[-1]]
    )


def test_a_missing_alignment_table_is_a_file_not_found_error(bed, tmp_path):
    with pytest.raises(FileNotFoundError, match="alignment table"):
        run(tmp_path / "nowhere.tsv", bed, tmp_path / "out")


def test_an_alignment_table_without_guide_id_is_refused(bed, tmp_path):
    frame = sites_frame(site("g", "chr1", 1000)).drop(columns=["guide_id"])
    with pytest.raises(ElementError, match="no guide_id column"):
        run(frame, bed, tmp_path / "out")


# --- the cau command ------------------------------------------------------


def test_element_mapping_is_registered():
    assert cli.COMMANDS["element-mapping"][0] == "element_mapping"


def test_parse_run_defaults():
    args = cli.build_parser().parse_args(
        [
            "element-mapping",
            "run",
            "--alignment",
            "aln.tsv",
            "--elements",
            "e.bed",
            "--outdir",
            "out",
        ]
    )
    assert (args.command, args.action) == ("element-mapping", "run")
    assert (args.alignment, args.elements, args.outdir) == ("aln.tsv", "e.bed", "out")
    assert args.annotation is None


def test_the_annotation_is_optional_and_parsed_when_given():
    args = cli.build_parser().parse_args(
        [
            "element-mapping",
            "run",
            "--alignment",
            "aln.tsv",
            "--elements",
            "e.bed",
            "--outdir",
            "out",
            "--annotation",
            "ann.tsv",
        ]
    )
    assert args.annotation == "ann.tsv"


@pytest.mark.parametrize(
    "argv",
    [
        ["element-mapping"],
        ["element-mapping", "map"],
        ["element-mapping", "run"],
        ["element-mapping", "run", "--alignment", "aln.tsv", "--elements", "e.bed"],
        ["element-mapping", "run", "--elements", "e.bed", "--outdir", "out"],
        ["element-mapping", "run", "--alignment", "aln.tsv", "--outdir", "out"],
    ],
)
def test_bad_arguments_are_usage_errors(argv):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(argv)
    assert exit_info.value.code == 2


def test_a_run_from_the_command_line(alignment, bed, tmp_path, capsys):
    status = cli.main(
        [
            "element-mapping",
            "run",
            "--alignment",
            str(alignment),
            "--elements",
            str(bed),
            "--outdir",
            "out",
        ]
    )

    assert status == 0
    assert capsys.readouterr().out.splitlines() == [
        "4 guides, 3 sites: 2 targeting, 1 safe_targeting, 2 non_targeting",
        "2 elements have at least one guide",
        "Outputs in out",
    ]
    assert matches(written(tmp_path / "out"), EXPECTED_FILES)


@pytest.mark.parametrize(
    ("alignment_arg", "elements_arg", "message"),
    [
        ("nowhere.tsv", "elements.bed", "The alignment table nowhere.tsv does not"),
        ("aln.tsv", "nowhere.bed", "The elements file nowhere.bed does not exist"),
        ("aln.tsv", "two_columns.bed", "a BED needs at least three"),
    ],
)
def test_a_bad_input_is_reported_without_a_traceback(
    tmp_path, capsys, alignment_arg, elements_arg, message
):
    write_alignment(tmp_path / "aln.tsv", ROWS)
    write_bed(tmp_path / "elements.bed")
    (tmp_path / "two_columns.bed").write_text("chr1\t1000\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exit_info:
        cli.main(
            [
                "element-mapping",
                "run",
                "--alignment",
                alignment_arg,
                "--elements",
                elements_arg,
                "--outdir",
                "out",
            ]
        )

    assert exit_info.value.code == 1
    error = capsys.readouterr().err
    assert error.startswith("cau: error: ")
    assert message in error
    assert "Traceback" not in error


def test_the_run_log_does_not_reach_the_root_logger(alignment, bed, tmp_path):
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
        run(alignment, bed, tmp_path / "out")
    finally:
        root.removeHandler(catch)
        root.setLevel(level)

    assert [r for r in seen if r.name.startswith("crispr_analysis_utils")] == []
    assert (
        "Element mapping complete"
        in (tmp_path / "out" / "logs" / "element_mapping.log").read_text()
    )

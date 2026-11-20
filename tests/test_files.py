"""The run record every step keeps, and the helpers the steps share."""

import json
import logging
from datetime import UTC, datetime

import pytest
import yaml

from crispr_analysis_utils import _files, element_mapping, igvf_export, pairs
from crispr_analysis_utils._files import METADATA_DIR, run_record
from crispr_analysis_utils.guide_alignment import summary as alignment_summary
from crispr_analysis_utils.guide_annotation import summary as annotation_summary

PACKAGE = logging.getLogger("crispr_analysis_utils")
LOGGER = logging.getLogger("crispr_analysis_utils.step")


@pytest.fixture
def package_level():
    """Put the package logger back however a test leaves it."""
    level, propagate = PACKAGE.level, PACKAGE.propagate
    yield
    PACKAGE.setLevel(level)
    PACKAGE.propagate = propagate


def records(outdir):
    return sorted((outdir / METADATA_DIR).glob("step_*.json"))


def test_a_run_leaves_the_package_logger_as_it_found_it(tmp_path, package_level):
    PACKAGE.setLevel(logging.NOTSET)
    handlers = list(PACKAGE.handlers)
    with pytest.raises(FileNotFoundError):
        pairs.run(tmp_path / "absent.tsv", "panel.csv", "store.parquet", tmp_path)
    assert PACKAGE.level == logging.NOTSET
    assert PACKAGE.propagate is True
    assert PACKAGE.handlers == handlers


def test_info_reaches_the_log_however_quiet_the_package_was(tmp_path, package_level):
    PACKAGE.setLevel(logging.WARNING)
    with run_record(tmp_path, "step", "results", {}, LOGGER, "Step"):
        LOGGER.info("inside the run")
    assert "inside the run" in (tmp_path / "logs" / "step.log").read_text()
    assert PACKAGE.level == logging.WARNING


def test_the_record_exists_from_the_start(tmp_path):
    with run_record(tmp_path, "step", "results", {"a": 1}, LOGGER, "Step"):
        (record,) = records(tmp_path)
        state = json.loads(record.read_text())
        assert (state["complete"], state["finished"]) == (False, None)
    state = json.loads(record.read_text())
    assert state["complete"] is True
    assert state["parameters"] == {"a": 1}
    assert "total" in state["seconds"]


def test_a_failure_is_recorded_logged_and_raised(tmp_path):
    with (
        pytest.raises(ValueError, match="boom"),
        run_record(tmp_path, "step", "results", {}, LOGGER, "Step"),
    ):
        raise ValueError("boom")
    (record,) = records(tmp_path)
    state = json.loads(record.read_text())
    assert state["complete"] is False
    assert state["error"] == "ValueError: boom"
    assert state["finished"] is not None
    assert (
        "Step failed: ValueError: boom" in (tmp_path / "logs" / "step.log").read_text()
    )


def test_two_runs_in_one_second_keep_both_records(tmp_path, monkeypatch):
    class Frozen(datetime):
        @classmethod
        def now(cls, tz=None):
            return datetime(2026, 1, 1, tzinfo=UTC)

    monkeypatch.setattr(_files, "datetime", Frozen)
    for n in (1, 2):
        with run_record(tmp_path, "step", "results", {"n": n}, LOGGER, "Step"):
            pass
    names = [record.name for record in records(tmp_path)]
    assert names == ["step_20260101T000000Z.json", "step_20260101T000000Z_2.json"]
    # The later run sorts last, which is how its record is found.
    assert json.loads(records(tmp_path)[-1].read_text())["parameters"] == {"n": 2}


@pytest.mark.parametrize(
    ("module", "table"),
    [
        (element_mapping, element_mapping.TABLE),
        (igvf_export, igvf_export.TABLE),
        (pairs, pairs.TABLE),
        (annotation_summary, annotation_summary.TABLE),
        (alignment_summary, "guide_genome_annotation.tsv"),
    ],
    ids=lambda value: getattr(value, "__name__", value).rsplit(".", 1)[-1],
)
def test_a_quote_or_backslash_in_a_description_stays_valid_yaml(
    module, table, monkeypatch, tmp_path
):
    column = next(iter(module.COLUMNS))
    monkeypatch.setitem(module.COLUMNS, column, 'A "quoted"\n    word, and C:\\path')
    module.write_columns_yaml(tmp_path / "columns.yaml")
    parsed = yaml.safe_load((tmp_path / "columns.yaml").read_text())
    assert parsed[table][column] == 'A "quoted" word, and C:\\path'


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("TRUE", True),
        (" true ", True),
        (True, True),
        (1, True),
        ("FALSE", False),
        ("", False),
        (None, False),
        (float("nan"), False),
        (False, False),
    ],
)
def test_truthy(value, expected):
    assert _files.truthy(value) is expected


def test_truthy_takes_a_nullable_missing_value():
    import pandas as pd

    assert _files.truthy(pd.NA) is False

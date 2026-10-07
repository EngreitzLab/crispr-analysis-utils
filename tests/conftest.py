import os

import fake_gem as fake_gem_tools
import pytest


@pytest.fixture(autouse=True)
def _run_in_tmp_path(tmp_path, monkeypatch):
    """Run every test from its own temporary directory.

    Some functions write "auto" outputs relative to the working directory
    (``filter_guide_alignments`` without SAM outputs), so a test run from the
    repository root would otherwise leave files there.
    """
    monkeypatch.chdir(tmp_path)


@pytest.fixture
def fake_gem(tmp_path, monkeypatch):
    """Put the fake gem-indexer and gem-mapper first on the PATH.

    Returns the file where the fakes record their command lines (see
    `fake_gem.calls`). The ``FAKE_GEM_*`` switches start unset.
    """
    folder = fake_gem_tools.install(tmp_path / "fake-gem-bin")
    monkeypatch.setenv("PATH", f"{folder}{os.pathsep}{os.environ.get('PATH', '')}")
    record = tmp_path / "fake-gem-calls.jsonl"
    monkeypatch.setenv("FAKE_GEM_CALLS", str(record))
    for switch in (
        "FAKE_GEM_CRASH",
        "FAKE_GEM_SIGNAL_EXIT",
        "FAKE_GEM_OUT_OF_RANGE",
        "FAKE_GEM_INDEX_FAIL",
    ):
        monkeypatch.delenv(switch, raising=False)
    return record

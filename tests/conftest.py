import pytest


@pytest.fixture(autouse=True)
def _run_in_tmp_path(tmp_path, monkeypatch):
    """Run every test from its own temporary directory.

    Some functions write "auto" outputs relative to the working directory
    (``filter_guide_alignments`` without SAM outputs), so a test run from the
    repository root would otherwise leave files there.
    """
    monkeypatch.chdir(tmp_path)

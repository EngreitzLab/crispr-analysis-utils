import pytest


@pytest.fixture(autouse=True)
def _run_in_tmp_path(tmp_path, monkeypatch):
    """Run every test from its own temporary directory.

    A test that writes to a relative path, such as a ``cau`` command given
    ``--outdir out``, then writes inside ``tmp_path`` instead of the
    repository.
    """
    monkeypatch.chdir(tmp_path)

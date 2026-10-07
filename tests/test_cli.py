import argparse
import types
from importlib.metadata import entry_points

import pytest

import crispr_analysis_utils as cau
from crispr_analysis_utils import cli


def test_console_script_points_at_main():
    (script,) = entry_points(group="console_scripts", name="cau")
    assert script.value == "crispr_analysis_utils.cli:main"


def test_package_version_comes_from_installed_metadata():
    assert cau.__version__ != "0+unknown"


def test_version_flag_prints_the_package_version(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--version"])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"cau {cau.__version__}"


def test_help_lists_every_command(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--help"])
    assert exit_info.value.code == 0
    out = capsys.readouterr().out
    for name in cli.COMMANDS:
        assert name in out


def test_missing_command_is_a_usage_error():
    with pytest.raises(SystemExit) as exit_info:
        cli.main([])
    assert exit_info.value.code == 2


@pytest.fixture
def fake_command(monkeypatch):
    """Register a `fake` subcommand whose behaviour each test chooses."""
    module = types.SimpleNamespace(
        add_arguments=lambda parser: parser.add_argument("--fail", action="store_true"),
    )

    def run(args: argparse.Namespace) -> int:
        if args.fail:
            raise cli.CommandError("something expected went wrong")
        return 0

    module.run = run
    monkeypatch.setitem(cli.COMMANDS, "fake", ("fake", "A test command."))
    monkeypatch.setattr(cli, "_command_module", lambda name: module)
    return module


def test_command_exit_status_is_returned(fake_command):
    assert cli.main(["fake"]) == 0


def test_command_error_is_reported_without_traceback(fake_command, capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["fake", "--fail"])
    assert exit_info.value.code == 1
    assert capsys.readouterr().err == "cau: error: something expected went wrong\n"

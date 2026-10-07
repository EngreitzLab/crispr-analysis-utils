"""The ``cau`` command-line interface.

Each subcommand lives in its own module under ``crispr_analysis_utils.cli``
and exposes two functions:

- ``add_arguments(parser)`` declares the subcommand's arguments;
- ``run(args) -> int`` does the work and returns the exit status.

A module imports its optional dependencies inside ``run``, so one command's
dependencies never slow down or break another command, ``--help`` or
``--version``.
"""

from __future__ import annotations

import argparse
import importlib
from collections.abc import Sequence
from types import ModuleType

from .. import __version__

# Subcommand name -> (module under crispr_analysis_utils.cli, one-line help).
COMMANDS: dict[str, tuple[str, str]] = {
    "guide-alignment": (
        "guide_alignment",
        "Align a CRISPR guide library to a reference with GEM3.",
    ),
    "install-skills": (
        "install_skills",
        "Install the bundled Claude Code skills and agents.",
    ),
}


class CommandError(Exception):
    """An expected failure, reported as ``cau: error: ...`` without a traceback."""


def build_parser() -> argparse.ArgumentParser:
    """Build the ``cau`` argument parser with every registered subcommand."""
    parser = argparse.ArgumentParser(
        prog="cau",
        description="Command-line tools from crispr-analysis-utils.",
    )
    parser.add_argument("--version", action="version", version=f"cau {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)
    for name, (module_name, help_text) in COMMANDS.items():
        subparser = subparsers.add_parser(name, help=help_text, description=help_text)
        _command_module(module_name).add_arguments(subparser)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run ``cau`` with `argv` (default: the process arguments)."""
    parser = build_parser()
    args = parser.parse_args(argv)
    module_name, _ = COMMANDS[args.command]
    try:
        return _command_module(module_name).run(args)
    except CommandError as error:
        parser.exit(1, f"cau: error: {error}\n")


def _command_module(module_name: str) -> ModuleType:
    return importlib.import_module(f"{__name__}.{module_name}")

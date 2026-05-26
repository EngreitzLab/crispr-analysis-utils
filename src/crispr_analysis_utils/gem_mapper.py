"""Helpers for building GEM indices and running GEM mapping."""

from __future__ import annotations

from pathlib import Path
import shlex

from .utils import run_shell_cmd


def build_gem_index(
    reference_fasta: str | Path,
    index_prefix: str | Path,
    *,
    threads: int | None = None,
    gem_indexer_bin: str = "gem-indexer",
    log_path: str | Path | None = None,
) -> str:
    """Build a GEM index from a reference FASTA."""
    cmd = [
        gem_indexer_bin,
        "-i",
        str(reference_fasta),
        "-o",
        str(index_prefix),
    ]
    if threads is not None:
        cmd.extend(["-t", str(threads)])
    cmd_str = " ".join(shlex.quote(token) for token in cmd)
    if log_path is not None:
        cmd_str = f"{cmd_str} > {shlex.quote(str(log_path))} 2>&1"
    return run_shell_cmd(cmd_str)


def map_guides_with_gem(
    gem_index: str | Path,
    input_fastq: str | Path,
    output_sam: str | Path,
    *,
    threads: int = 8,
    mapping_mode: str = "sensitive",
    sam_compact: bool = False,
    gem_mapper_bin: str = "gem-mapper",
    log_path: str | Path | None = None,
) -> str:
    """Run GEM mapper for single-end FASTQ input.

    If `log_path` is not provided, command output is redirected to a `.log` file
    next to `output_sam`.
    """
    output_sam_path = Path(output_sam)
    if log_path is None:
        log_path = output_sam_path.with_suffix(".log")

    cmd = [
        gem_mapper_bin,
        "-I",
        str(gem_index),
        "-i",
        str(input_fastq),
        "-o",
        str(output_sam_path),
        f"--sam-compact={'true' if sam_compact else 'false'}",
        "--mapping-mode",
        mapping_mode,
        "--threads",
        str(threads),
    ]
    cmd_str = " ".join(shlex.quote(token) for token in cmd)
    cmd_str = f"{cmd_str} > {shlex.quote(str(log_path))} 2>&1"
    return run_shell_cmd(cmd_str)

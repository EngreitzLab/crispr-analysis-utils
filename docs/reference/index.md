# API reference

The Python API of `crispr_analysis_utils`, rendered from its docstrings.

## Guide alignment

`crispr_analysis_utils.guide_alignment` finds every guide's genomic site; the
[guide](../guide-alignment.md) explains how. The package re-exports the main
names, such as `guide_alignment.run`. Module by module:

- [`pipeline`](guide-alignment-run.md): `run`, which chains every step, and
  its `RunResult`.
- [`library`](guide-alignment-library.md): reading the guide table, and
  writing the reads to align.
- [`scan`](guide-alignment-scan.md): the genome pass, over the Rust core.
- [`sites`](guide-alignment-sites.md): the `Site` record.
- [`orientation`](guide-alignment-orientation.md): the orientation verdicts.
- [`summary`](guide-alignment-summary.md): the classes, and the output tables.
- [`iupac`](guide-alignment-iupac.md): IUPAC codes and PAM patterns.

## Normalization

- [`counts_per_million`](counts-per-million.md)

## Utilities

- [`run_shell_cmd`](run-shell-cmd.md)

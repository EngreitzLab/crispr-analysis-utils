# API reference

The Python API of `crispr_analysis_utils`, rendered from its docstrings.

## Guide alignment

`crispr_analysis_utils.guide_alignment` maps a guide library to its genomic
sites; the [guide](../guide-alignment.md) explains how. Its package re-exports
the main names, such as `guide_alignment.run` and
`guide_alignment.build_index`. Module by module:

- [`pipeline`](guide-alignment-run.md): `run`, which chains every step, and
  its `RunResult`.
- [`library`](guide-alignment-library.md): reading the guide table, and
  writing the reads to align.
- [`gem`](guide-alignment-gem.md): building the GEM index and running
  gem-mapper.
- [`sites`](guide-alignment-sites.md): the site search's alignments, checked
  against the reference.
- [`orientation`](guide-alignment-orientation.md): the orientation check.
- [`summary`](guide-alignment-summary.md): the classes, and the output tables.
- [`iupac`](guide-alignment-iupac.md): IUPAC codes and PAM patterns.

## Accession

- [`accession`](accession.md): accession ids for SpCas9 guides, encoded
  from the guide's own sequence.

## Normalization

- [`counts_per_million`](counts-per-million.md)

## Utilities

- [`run_shell_cmd`](run-shell-cmd.md)

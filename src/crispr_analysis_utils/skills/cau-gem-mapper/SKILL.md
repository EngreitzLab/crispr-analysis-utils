---
name: cau-gem-mapper
description: >-
  Build a GEM3 genome index and map guide-RNA reads with crispr-analysis-utils
  (`cau.gem_mapper`: build_gem_index and map_guides_with_gem). Use when a user
  wants to align CRISPR guides (a FASTQ written by
  cau.guide_qc.guides_to_fastq) to a reference genome with GEM3, index a
  reference FASTA such as hg38, choose GEM threads or mapping mode, or debug a
  gem-indexer or gem-mapper failure. Pairs with cau-guide-qc, which filters the
  SAM this module writes.
---

# GEM mapping (`cau.gem_mapper`)

Thin wrappers that run the GEM3 binaries through bash with `pipefail`. A
failing command raises `RuntimeError` carrying its exit status, stdout and
stderr.

## Getting GEM3

`gem-indexer` and `gem-mapper` come from the bioconda package `gem3-mapper`,
which is built for linux-64 and osx-64 only, not for Apple silicon. In this
repository the pixi environments include it on linux-64 (`pixi shell`, or
prefix commands with `pixi run`). Elsewhere:
`conda install -c conda-forge -c bioconda gem3-mapper`. Running `gem-mapper`
with no arguments prints its usage, a quick check that it is installed.

## Ask before running

- **Which reference FASTA**, and its genome build: guide coordinates come out
  in it.
- **Whether an index already exists.** Indexing a whole genome is slow; reuse
  an existing `<prefix>.gem`.
- **How many threads** the machine or the cluster job allows.

## Usage

```python
import crispr_analysis_utils as cau

cau.gem_mapper.build_gem_index(
    "hg38.fa",
    "gem_index/hg38",  # GEM writes gem_index/hg38.gem
    threads=8,
    log_path="gem_index/hg38.log",
)
cau.gem_mapper.map_guides_with_gem(
    "gem_index/hg38.gem",
    "guides.fastq",
    "guides_mapped.sam",  # GEM's own output goes to guides_mapped.log
    threads=8,
)
```

## Defaults and traps

- `build_gem_index` takes an index *prefix*; GEM writes `<prefix>.gem`, the
  file `map_guides_with_gem` wants. Its `threads` defaults to GEM's own choice.
- `map_guides_with_gem` defaults: `threads=8`, `mapping_mode="sensitive"`,
  `sam_compact=False`.
- Keep `sam_compact=False`. GEM's own default is compact, which folds
  subdominant matches into the `XA` tag of one record; `cau-guide-qc` counts
  alignments as separate SAM records and would miss them.
- Logs: without `log_path`, `map_guides_with_gem` sends GEM's output to the SAM
  path with a `.log` suffix. `build_gem_index` writes a log only when
  `log_path` is given, and otherwise returns GEM's stdout.
- The binaries are looked up on PATH; `gem_indexer_bin` and `gem_mapper_bin`
  override them.

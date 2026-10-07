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
failing command raises `RuntimeError` with its exit status (`RC=`). GEM's own
output goes to the log file, always for `map_guides_with_gem` and for
`build_gem_index` when `log_path` is given, so in those cases the exception's
`STDOUT` and `STDERR` are empty: **read the log**. A missing binary shows there
as `command not found` with `RC=127`.

## Getting GEM3

`gem-indexer` and `gem-mapper` come from the bioconda package `gem3-mapper`,
which is built for linux-64 and osx-64 only, not for Apple silicon. In this
repository the pixi `default` and `dev` environments include it on both. On
Apple silicon pixi installs those environments as osx-64, which runs under
Rosetta 2, Python and pysam included (`pixi shell -e dev`, or prefix commands
with `pixi run -e dev`). Elsewhere:
`conda install -c conda-forge -c bioconda gem3-mapper`, which on Apple silicon
works only in an osx-64 environment. Running `gem-mapper` with no arguments
prints its usage, a quick check that it is installed.

## Ask before running

- **Which reference FASTA**, and its genome build: guide coordinates come out
  in it.
- **Whether an index already exists.** Indexing a whole genome is slow; reuse
  an existing `<prefix>.gem`.
- **How many threads** the machine or the cluster job allows.

## Usage

```python
from pathlib import Path

import crispr_analysis_utils as cau

Path("gem_index").mkdir(parents=True, exist_ok=True)  # neither wrapper creates folders
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
- Neither wrapper creates folders. The folder of the log (by default, the SAM's
  folder) must exist, or bash fails the redirect before GEM starts: a
  `RuntimeError` whose `STDERR` says `No such file or directory`. Create the
  index folder up front as well.
- The binaries are looked up on PATH; `gem_indexer_bin` and `gem_mapper_bin`
  override them.

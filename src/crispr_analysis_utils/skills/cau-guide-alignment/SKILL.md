---
name: cau-guide-alignment
description: >-
  Align a CRISPR guide library to a reference genome with crispr-analysis-utils
  (the `cau guide-alignment` command, or
  `crispr_analysis_utils.guide_alignment.run`; GEM3 underneath) and read its
  outputs. Use when a user wants to know where their guides land in hg38 or
  another genome: which guides are unique, multi-mapping, off-target-only or
  without a site; their protospacer coordinates and cut sites, for sceptre,
  element annotation or IGVF-style guide metadata; their sites with up to 3
  mismatches and an NGG, NAG or NGA PAM; whether a guide table is
  reverse-complemented; or to build a GEM index. Also for reading guides.tsv,
  sites.tsv, orientation.tsv, guides.bed, cut_sites.bed, summary.tsv or
  run.json, and for users of the earlier guide_qc and gem_mapper workflow
  (valid_alignments.bed, alignment_summary.tsv).
---

# Guide alignment (`cau guide-alignment`)

`cau guide-alignment run` maps every guide of a library to its genomic sites
with GEM3, checks every hit against the reference, and classifies each guide:

1. It reads the guide table, ids kept as text, and validates it.
2. **Pass 1, orientation**: it aligns each bare spacer with no mismatch
   allowed and reads the flanks of every exact hit, to check that the spacers
   are given 5' to 3'.
3. **Pass 2, sites**: it aligns one read per guide and PAM (`NGG`, `NAG` and
   `NGA` by default), the spacer followed by the PAM, with a complete search
   up to `--max-mismatches`.
4. It reads every hit back from the reference, recounts its spacer
   mismatches, keeps those within the limit, and classifies its genomic PAM.
5. It classifies each guide on its primary-PAM sites and writes the tables.

## Getting GEM3

`gem-indexer` and `gem-mapper` are not Python packages: pip installs the
package without them. They come with the repository's pixi environments
(`default` and `dev`), on linux-64 and osx-64; on Apple silicon pixi installs
those environments as osx-64, which runs under Rosetta 2. From a clone of the
repository, prefix commands with `pixi run` (from elsewhere,
`pixi run -m /path/to/crispr-analysis-utils ...`). `pixi run gem-mapper --version`
checks that GEM is there. Outside pixi, GEM3 is bioconda's `gem3-mapper`.

## Usage

```bash
# Once per reference: the GEM index, plus gem_index/hg38.gem.json with the
# reference's MD5 checksum.
pixi run cau guide-alignment index hg38.fa gem_index/hg38 --threads 8

pixi run cau guide-alignment run \
    --guides guides.tsv \
    --reference hg38.fa \
    --index gem_index/hg38.gem \
    --outdir results/guide_alignment \
    --threads 8
```

The guide table defaults to tab-separated, without a header: the guide id,
then its sequence. The run prints the class counts and logs its progress to
standard error; `cau guide-alignment run --help` lists every option. In
Python:

```python
from crispr_analysis_utils import guide_alignment

guide_alignment.build_index("hg38.fa", "gem_index/hg38", threads=8)
result = guide_alignment.run(
    "guides.tsv",  # or a pandas DataFrame
    "hg38.fa",
    "gem_index/hg38.gem",
    "results/guide_alignment",
    add_leading_g=True,
)
for guide in result.summaries:
    print(guide.guide_id, guide.alignment_class, guide.cut_site)
```

## Ask before running

Ask, don't guess:

- **The reference FASTA, and its genome build.** Coordinates come out in it.
  It must be uncompressed (a gzip-compressed FASTA is rejected) and the very
  file the index was built from: `run` compares its MD5 checksum with the one
  `PREFIX.gem.json` recorded, and stops on a mismatch.
- **Whether an index already exists.** Indexing a whole genome is slow; reuse
  a `PREFIX.gem` whose `PREFIX.gem.json` names the same reference.
- **The guide table's layout**: which column holds the id and which the
  sequence (`--id-col`, `--spacer-col`: 0-based positions, or names with
  `--header`), the separator (`--sep`: `tab`, `comma` or one character), and
  whether the first line is a header. Look at the first lines before running.
- **What the sequences hold.** They must be 5' to 3' and must not include the
  PAM: the spacer is the 3'-most `--spacer-length` bases (default 20), so a
  23-mer ending in its PAM gives a wrong spacer. Strip the PAM first.
- **The PAM**: `--pam` (default `NGG`, in IUPAC codes, 3' of the protospacer)
  and `--alt-pams` (default `NAG,NGA`; `""` for none). Every PAM must have
  the same length. Only 3' PAMs are supported.
- **Whether the guides are expressed with an added 5' G** (as from a U6
  promoter): `--add-leading-g`, off by default. See the trap below.
- **The mismatch limit**: `--max-mismatches`, default 3.
- **Which contigs count**: `--contigs chr1,chr2,...` keeps only the hits on
  those contigs (default: every contig). Each name must be in the reference.
- **Where it runs, with how many threads** (`--threads`, default 8). A whole
  library against a human genome belongs on a cluster node: the
  `cau-guide-alignment-runner` agent runs it there and watches it.
- **The output folder**, `--outdir`, which is required. Nothing is written
  anywhere else, and a rerun into the same folder reuses the GEM passes whose
  inputs have not changed.

## Defaults

| CLI option | Python argument of `run` | Default |
| --- | --- | --- |
| `--max-mismatches` | `max_mismatches` | `3` |
| `--pam` | `pam` | `NGG` |
| `--alt-pams` | `alt_pams` | `NAG,NGA` (Python: `("NAG", "NGA")`) |
| `--spacer-length` | `spacer_length` | `20` |
| `--add-leading-g` | `add_leading_g` | off |
| `--orientation-check` | `orientation_check` | `warn` |
| `--contigs` | `contigs` | every contig (Python: `None`) |
| `--threads` | `threads` | `8` |
| `--id-col`, `--spacer-col` | `id_col`, `spacer_col` | `0`, `1` |
| `--sep` | `sep` | tab |
| `--header` | `header` | off |

`cau guide-alignment index --threads` defaults to 8 as well, but the Python
`build_index` defaults to `threads=None`, which leaves gem-indexer's own
default: every core.

## How a guide is classified

Only the primary PAM's sites decide the class:

| Class | Primary-PAM sites |
| --- | --- |
| `unique` | Exactly one with 0 spacer mismatches |
| `multi` | More than one with 0 mismatches |
| `off_target_only` | None with 0 mismatches, at least one with 1 to `--max-mismatches` |
| `no_site` | None within `--max-mismatches` |

`targeting` is `TRUE` for `unique` and `multi` guides. The coordinates,
`pam`, `strand` and `cut_site` are filled only for `unique` guides.
Alternative-PAM sites are counted in `n_alt_pam_*` and never decide the
class; a hit whose genomic PAM matches no pattern is listed in `sites.tsv`
with `pam_class` `none` and counts nowhere. `n_ngg_*` and
`min_ngg_mismatches` count the primary PAM's sites, whatever `--pam` is.

## Traps

- **Orientation.** Spacers must be given 5' to 3'. A reverse-complemented
  spacer still has exact hits, with the reverse complement of its PAM on the
  5' side: its verdict is `reversed`. `--orientation-check warn` (the default)
  logs a warning when `reversed` guides outnumber `ok` ones, counted per
  distinct spacer; `error` stops the run there, before the site search;
  `off` skips pass 1, and every guide is `not_checked`. The hits and their
  flanks are in `orientation.tsv`.
- **The added G.** With `--add-leading-g`, a spacer that does not start with
  G is aligned with an extra 5' G. Its genomic base is reported
  (`leading_g_base` in `sites.tsv`) but never counted as a mismatch. A spacer
  that already starts with G gets nothing added, and its own G counts like
  any spacer base. The G is part of the read, so the site needs one more
  reference base: a site whose protospacer starts at a contig's first base
  (on `+`) or ends at its last base (on `-`), or borders a run of 50 or more
  N (gem-indexer leaves those runs out of the index), is not found with it.
- **GEM 3.6 crashes at the start of its index.** A read that aligns across
  the first base of the reference's first contig, on `+`, makes gem-mapper
  crash ("Signal raised") or, on some genomes, report a bogus hit instead;
  with the added G, a protospacer starting at that base is enough. A crash
  fails the run with a `GemError` naming the FASTQ, the read and its guide. hg38's first contig
  starts with N, so it cannot happen there. The error's remedies: index a
  reference whose first contig starts with N (reorder or drop contigs), or
  align without the added G. `--contigs` does not help: it filters hits after
  GEM has run.
- **Coordinates are 0-based and half-open** on the forward strand, as in BED,
  and the protospacer excludes the PAM and the added G. Add 1 to
  `guide_start` for a 1-based start.
- **The cut site is a boundary.** `cut_site` is `b` when the cut falls
  between bases `b - 1` and `b`, 3 bases 5' of the PAM: `guide_end - 3` on
  `+`, `guide_start + 3` on `-`. `cut_sites.bed` writes it as the one base
  `[b, b + 1)`.
- **Non-targeting controls** are expected to come out `off_target_only` at
  3 mismatches in hg38, with `targeting` `FALSE`. One that comes out `unique`
  or `multi` has a perfect site: report it.
- **Duplicate spacers.** Guides sharing a spacer are aligned once; the later
  ones get `duplicate_of` and the same sites, and a warning lists them.
- **The guide table is strict.** Ids are read as text (`0001` stays `0001`),
  must be unique, and must hold no tab or line break. Sequences may hold only
  A, C, G and T (case is ignored), at least `--spacer-length` of them.
- **An index without `PREFIX.gem.json`** (built outside this package) cannot
  be checked against the reference: `run` warns, and checks only that the
  contig names and lengths match.
- **Hits GEM places badly are dropped and counted.** `n_out_of_range` in
  `summary.tsv` counts records gem-mapper placed outside their contig (a read
  overhanging a contig end); `n_over_max_mismatches` counts the hits the
  recount against the reference drops; `n_nm_disagreements` counts hits whose
  `NM` tag differs from the recount, as GEM's hits inside runs of N do.

## Reading the results

Start with `summary.tsv`: the class and orientation counts (each sums to the
number of guides), the alignment counts of both passes, and the parameters.
Then `guides.tsv`, one row per guide. `run.json` says whether the run
completed (`complete`, `error`) and holds each GEM command and its run time.
After a failure, read `logs/guide_alignment.log` and the
`logs/gem-mapper.<pass>.log` it names.

## Outputs

Every file, column by column, with the Python return value:
`references/outputs.md`.

## Coming from `guide_qc` and `gem_mapper`

Why the counts differ from the earlier workflow, and where each column of
`valid_alignments.bed`, `guide_alignment_log.tsv` and
`alignment_summary.tsv` went: `references/migration.md`.

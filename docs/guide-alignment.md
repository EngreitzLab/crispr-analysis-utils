# Guide alignment

`cau guide-alignment` maps a CRISPR guide library to its genomic sites with
GEM3, checks every hit against the reference, and classifies each guide as
unique, multi-mapping, off-target only, or without a site. It writes a guide
table with IGVF-style column names, protospacer coordinates and cut sites,
every site with up to 3 mismatches, an orientation check of the library, and
BED files of the unique guides.

The [tutorial](tutorial-guide-alignment.md) runs it on hg38 on a SLURM
cluster. In Claude Code, the `cau-guide-alignment` skill and the
`cau-guide-alignment-runner` agent run it for you
([Agents and skills](agents-and-skills.md)).

## Install

The command comes with the package, but its aligner, GEM3, is not a Python
package. The repository's pixi environments carry it, on Linux (x86-64) and
macOS, Apple silicon included through Rosetta 2 (see
[Getting started](getting-started.md)):

```bash
git clone https://github.com/EngreitzLab/crispr-analysis-utils.git
cd crispr-analysis-utils
pixi install
pixi run gem-mapper --version
```

From inside the clone, `pixi run cau ...` runs the command; from elsewhere,
`pixi run -m /path/to/crispr-analysis-utils cau ...`.

## Quick start

Index the reference once, then align the library:

```bash
pixi run cau guide-alignment index hg38.fa gem_index/hg38 --threads 8

pixi run cau guide-alignment run \
    --guides guides.tsv \
    --reference hg38.fa \
    --index gem_index/hg38.gem \
    --outdir results/guide_alignment \
    --threads 8
```

The same from Python:

```python
from crispr_analysis_utils import guide_alignment

guide_alignment.build_index("hg38.fa", "gem_index/hg38", threads=8)
result = guide_alignment.run(
    "guides.tsv",
    "hg38.fa",
    "gem_index/hg38.gem",
    "results/guide_alignment",
)
```

`run` returns a `RunResult`, whose `summaries` hold one row of `guides.tsv`
per guide.

## The guide table

One guide per line: an id and a sequence. By default the file is
tab-separated, without a header, with the id in the first column and the
sequence in the second; `--sep`, `--header`, `--id-col` and `--spacer-col`
read other layouts, such as a CSV with column names. In Python, `run` also
takes a pandas DataFrame.

- **Ids are read as text**, so `0001` stays `0001` and `NA` is an id. They
  must be unique, not empty, and free of tabs and line breaks.
- **Sequences are given 5' to 3'**, without the PAM. Once stripped of
  surrounding whitespace and uppercased, they may hold only A, C, G and T,
  and at least `--spacer-length` of them.
- **The spacer is the 3'-most `--spacer-length` bases** (20 by default). The
  whole sequence is kept as `input_sequence`. A sequence that ends in its
  PAM gives a wrong spacer: strip the PAM first.
- **Duplicate spacers** are aligned once. Every later guide with the same
  spacer gets `duplicate_of`, the id of the first, and the same results; a
  warning lists them.

## How it works

1. **Pass 1, orientation.** Each distinct spacer is aligned bare (no PAM, no
   added G) with no mismatch allowed. For every exact hit, the reference's
   PAM-length flanks are read in guide orientation. A guide is `ok` when a
   hit has a PAM on its 3' side; else `reversed` when a hit has the reverse
   complement of a PAM on its 5' side; else `no_pam`; `no_perfect_hit`
   without any hit.
2. **The library-level check.** With `--orientation-check warn`, the default,
   a warning is logged when `reversed` guides outnumber `ok` ones, counted per
   distinct spacer: the table probably gives the spacers reverse-complemented.
   `error` stops the run there, and `off` skips pass 1 (every guide is then
   `not_checked`).
3. **Pass 2, sites.** One read per distinct spacer and PAM pattern: the
   spacer, then the pattern as written (`NGG`, `NAG`, `NGA` by default),
   after an added 5' G when `--add-leading-g` applies. The reads with and
   without the added G go to separate FASTQs, each mapped with its own error
   budget (see [the GEM options](#the-gem-options)).
4. **Every hit is checked against the reference.** The tool reads the
   protospacer, the PAM and the base under the added G from the reference,
   uppercased (soft-masked bases count like any other) and in guide
   orientation. It recounts the spacer mismatches, a genomic N counting as
   one and the added G never, and keeps the site when there are at most
   `--max-mismatches`; GEM's own `NM` is only compared. The genomic PAM gets
   the first pattern it matches, `--pam` first and then `--alt-pams` in their
   order, or `none`; a genomic N matches only an N of the pattern. A guide's
   hits from its several PAM reads are merged on contig, start and strand,
   and `--contigs`, when given, keeps only the hits on those contigs.
5. **Each guide gets a class** from its primary-PAM sites: `unique` with
   exactly one perfect site, `multi` with more, `off_target_only` with none
   but some with 1 to `--max-mismatches` mismatches, `no_site` otherwise.
   Alternative-PAM sites are counted but never decide the class, and `none`
   sites count nowhere.

`run` also checks the reference against the index. Before mapping, it
compares the reference's MD5 checksum with the one `PREFIX.gem.json`
recorded when `cau guide-alignment index` built the index, and stops on a
mismatch; an index without that file gets a warning instead. Reading the
alignments, it compares the contig names and lengths of every SAM header with
the reference's FASTA index: the reference's own `.fai` when it has one, else
one built in the output folder.

## Coordinates and cut sites

Coordinates are 0-based and half-open on the reference's forward strand, as
in BED. The protospacer interval excludes the PAM and the added G, on both
strands; the PAM lies right of it on `+` and left of it on `-`.

The cut site is a boundary `b`, 3 bases 5' of the PAM, where SpCas9 cuts:
between bases `b - 1` and `b`. So `b = guide_end - 3` on `+` and
`b = guide_start + 3` on `-`. `cut_sites.bed` writes it as the one base
`[b, b + 1)`, right of the cut on the forward strand. For example:

| Strand | `guide_start` | `guide_end` | PAM | `cut_site` | `cut_sites.bed` |
| --- | --- | --- | --- | --- | --- |
| `+` | 26900 | 26920 | `[26920, 26923)` | 26917 | `26917 26918` |
| `-` | 78932 | 78952 | `[78929, 78932)` | 78935 | `78935 78936` |

## Options

| `cau guide-alignment run` | `run` argument | Default | Meaning |
| --- | --- | --- | --- |
| `--guides` | `guides` | required | The guide table (Python: a path or a DataFrame) |
| `--reference` | `reference` | required | The reference FASTA the index was built from, uncompressed |
| `--index` | `index` | required | The GEM index, `PREFIX.gem` or `PREFIX` |
| `--outdir` | `outdir` | required | The output folder, created if missing |
| `--max-mismatches` | `max_mismatches` | `3` | Spacer mismatches a site may have |
| `--pam` | `pam` | `NGG` | The primary PAM, in IUPAC codes, 3' of the protospacer |
| `--alt-pams` | `alt_pams` | `NAG,NGA` | Alternative PAMs, comma-separated, `""` for none (Python: a list); counted, never deciding the class. Every PAM must have the same length |
| `--spacer-length` | `spacer_length` | `20` | The spacer is the 3'-most this many bases |
| `--add-leading-g` | `add_leading_g` | off | Align a spacer that does not start with G with an extra 5' G, never counted as a mismatch |
| `--orientation-check` | `orientation_check` | `warn` | `warn`, `error` or `off` (above) |
| `--contigs` | `contigs` | every contig | Keep only the hits on these contigs, comma-separated (Python: a list); each must be in the reference |
| `--threads` | `threads` | `8` | gem-mapper threads |
| `--id-col` | `id_col` | `0` | The id column: a 0-based position, or a name with `--header` |
| `--spacer-col` | `spacer_col` | `1` | The sequence column, the same way |
| `--sep` | `sep` | tab | The separator: `tab`, `comma` or one character (Python: the character) |
| `--header` | `header` | off | The first line holds column names |

`cau guide-alignment index REFERENCE PREFIX` takes `--threads`, 8 by default.
Its Python counterpart, `build_index`, defaults to `threads=None`, which
leaves gem-indexer's own default: every core. Missing folders of `PREFIX`
are created.

`--contigs` filters the hits after GEM has run: it does not shorten the
search.

## The GEM options

Every gem-mapper run uses these options, which a unit test pins:

| Option | Why |
| --- | --- |
| `--mapping-mode=customed` | GEM's own modes do not search completely: `fast`, its default, never runs a complete search, and `sensitive` does not either. `customed` lets the last two options set one |
| `--alignment-model=hamming` | Mismatches only: no insertion or deletion, so no bulges |
| `--clipping=none` | No clipped ends. Without the `=`, GEM silently ignores the option |
| `--alignment-local=never` | No local alignments: every alignment covers the whole read |
| `--sam-compact=false` | Every match in its own SAM record, not folded into an `XA` tag. It needs the `=` too |
| `--max-reported-matches=all` | GEM's default, 5, caps the search itself, not only the report |
| `--alignment-max-bandwidth=0` | Without it, GEM merges nearby hits on the same strand |
| `--alignment-max-error=k+g` | The errors an alignment may have; GEM adds the read's own non-ACGT codes, the PAM's `N`, to this limit itself |
| `--complete-search-error=k+n+g` | The errors up to which the search is complete |
| `--complete-strata-after-best=k+n+g` | How many more errors than the best match's the search goes on to cover |

Here k is `--max-mismatches`; n is the highest number of non-ACGT codes among
the PAMs of the FASTQ (1 for `NGG`, `NAG` and `NGA`), because GEM counts the
`N` of a read as an error; and g is 1 for the reads with the added G, 0 for
the others. With the defaults, the reads without the added G run with
`--alignment-max-error=3 --complete-search-error=4
--complete-strata-after-best=4`, those with it with 4, 5 and 5. Pass 1 runs
with k, n and g all 0.

These options were chosen against a brute-force scan of a synthetic genome
with planted sites: repeat families, contig edges, runs of N, soft-masked
stretches and a guide with 60 perfect copies. The test that compares the two,
`tests/test_guide_alignment_planted.py`, runs with GEM in CI, and it fails
under GEM's `fast` and `sensitive` modes, which shows it can tell them apart.

## Outputs

--8<-- "src/crispr_analysis_utils/skills/cau-guide-alignment/references/outputs.md:tables"

## Known limitations

GEM 3.6 has edge cases at contig ends and in runs of N. A crash fails the
run with a clear error; the other cases drop or miss hits, and `summary.tsv`
counts the drops.

- **A crash at the start of the index.** A read that aligns across the first
  base of the reference's first contig, on `+`, makes gem-mapper crash
  ("Signal raised") or, on some genomes, report a bogus hit. With
  `--add-leading-g`, a guide whose protospacer starts at that base is
  enough. A crash fails the run with a `GemError` naming the FASTQ, the read
  and its guide, and that pass leaves no SAM. hg38's first contig starts
  with N, so hg38 cannot trigger it. The error suggests indexing a reference
  whose first contig starts with N (reorder or drop contigs), or aligning
  without the added G. `--contigs` cannot help: it filters after GEM has
  run.
- **Hits past a contig end.** gem-mapper places a read that overhangs the end
  of a contig past that end, or at a position that is not valid SAM. Those
  records are dropped as gem-mapper writes them, and counted as
  `n_out_of_range` in `summary.tsv`.
- **Runs of N.** gem-indexer replaces every run of at least 50 N with a
  single separator, so a protospacer that starts inside such a run is not in
  the index. A hit GEM reports inside a run reads as N in the reference: the
  recount drops it (`n_over_max_mismatches`) and counts its `NM` as a
  disagreement (`n_nm_disagreements`).
- **A genomic N inside a protospacer** counts as a mismatch, and GEM's search
  does not reach a site that has one plus other mismatches.
- **The added G needs its base.** With `--add-leading-g` the G is part of the
  read, so a site whose G would fall off the contig (a protospacer starting
  at a contig's first base on `+`, or ending at its last base on `-`) or
  inside a run of 50 or more N is not found.

## IGVF column names

The first eight columns of `guides.tsv`, `guide_id` to `strand`, take their
names from the lab's IGVF-style guide metadata: the tables
`validations/metadata/guide-metadata-*-intermediate.tsv` in the
EngreitzLab/DC_TAP_Paper repository, which its
`validations/analyses/README.md` describes as built on the IGVF E2G group's
column definitions. The other columns are this tool's own. Columns the tool
cannot know, such as the intended target, are left to the user.

Three shared names do not mean quite the same thing in those tables, so
check before merging the two:

| Column | DC_TAP intermediate tables | `guides.tsv` |
| --- | --- | --- |
| `guide_start`, `guide_end` | The BLAT hit of the spacer and its PAM: 23 bases, PAM included | The protospacer only, 20 bases by default, PAM excluded |
| `targeting` | `True` when the guide aligned | `TRUE` when the guide has a perfect primary-PAM site |
| `pam` | CRISPick's PAM sequence | The genomic PAM of the unique site |

Their `predicted_cutposition`, CRISPick's 1-based cut position, is this
tool's `cut_site + 1`: both name the base right of the cut on the forward
strand, the one `cut_sites.bed` holds.

## Coming from `guide_qc` and `gem_mapper`

--8<-- "src/crispr_analysis_utils/skills/cau-guide-alignment/references/migration.md:migration"

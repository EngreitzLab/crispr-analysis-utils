# Coming from the GEM version

<!-- --8<-- [start:migration] -->
The tool used to align guides with GEM3. It now scans the reference directly:
a guide is designed, not sequenced, so its site is an **exact** match of the
protospacer, and the only other questions are what sits 3' of it (the PAM) and
what sits 5' of it. There is no aligner, no index and no error budget.

## What went away

| Gone | Why |
| --- | --- |
| `cau guide-alignment index`, `--index` | Nothing is indexed. Point `--reference` at the FASTA. |
| `--max-mismatches` | Matching is exact. Off-target search is separate work. |
| `--alt-pams` | Only `--pam` is searched. |
| `--threads` | One pass, no aligner to parallelise. |
| `alignments/*.bam`, `*.fastq` | No reads are written and nothing is aligned. |
| GEM3 itself | No external tool, so every environment is native on Apple silicon. |

A whole library now takes seconds rather than an hour, and the four GEM
defects the old version worked around — the `-M 5` search cap, the `N` of
`NGG` eating the error budget, non-exhaustive search modes and garbage
coordinates at contig ends — cannot arise.

## Renamed files

| Was | Now |
| --- | --- |
| `guides.tsv` and `sites.tsv` | one `guide_genome_annotation.tsv` |
| `guides.bed` | `guide_mapping.bed` |
| `cut_sites.bed` | gone; `cut_site` is a column |
| `orientation.tsv` | gone; its verdict and site are columns |
| `summary.tsv` | `qc.tsv`, counts only |
| `run.json` | `run_metadata/guide_alignment_<timestamp>.json` |

Everything now lives under `alignments/`, `run_metadata/` and `logs/`.

## Renamed columns

| Was | Now |
| --- | --- |
| `guide_id` (the lab's name) | `alias` |
| — | `guide_id` is now the aligned sequence |
| `targeting` | gone; use `n_cut_sites` |
| `alignment_class` | `site_class`, and `off_target_only` no longer exists |
| `n_ngg_<m>mm` | `n_cut_sites`, exact matches only |
| `min_ngg_mismatches` | gone |
| `leading_g` | gone; compare the lengths of `guide_id` and `spacer` |
| `pam_class` | `pam_matches`, `TRUE`/`FALSE` |
| `site_pam`, `guide_chr/start/end` | `pam` folded into `mapped_sequence`; `chr`, `start`, `end` |
| `no_perfect_hit` (orientation) | `no_match` |

`input_sequence` and `genomic_protospacer` are gone: the first is in your
input file and in `guide_id`, and the second is always identical to `spacer`
when matching is exact.

## What this changes in a result

- **Guides are keyed by sequence.** Two entries differing only by a 5' G are
  two guides, not duplicates. The old version truncated them to the same
  spacer and collapsed them.
- **`n_matches` is new and often surprising.** The old `multi` counted only
  PAM-bearing sites and still does; `n_matches` counts every occurrence of the
  protospacer, which is usually larger.
- **Orientation is decided from both readings**, not from a `CCN` upstream.
  A library reversed throughout stops the run; a handful of entries warns.
<!-- --8<-- [end:migration] -->

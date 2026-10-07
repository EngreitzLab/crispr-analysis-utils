# Coming from `guide_qc` and `gem_mapper`

<!-- --8<-- [start:migration] -->
`cau guide-alignment` replaces the earlier three-step workflow:
`guide_qc.guides_to_fastq`, `gem_mapper.build_gem_index` and
`gem_mapper.map_guides_with_gem`, then `guide_qc.filter_guide_alignments`,
chained by `scripts/run_guide_alignment_qc.py`.

## Why the numbers change

- **More multi-mapping guides.** gem-mapper's default reports at most 5
  matches per read, and its default `fast` mapping mode (like the
  `sensitive` mode the old pipeline used) does not search completely. The new
  options report every match of a complete search, so a guide with more than
  5 sites, or with sites the incomplete search missed, now shows them all.
- **Mismatches are counted, not ignored.** The old filter kept alignments
  with any number of spacer mismatches, up to GEM's error budget, and called
  a guide unique when it had one valid alignment. Now a site has at most
  `--max-mismatches` spacer mismatches, recounted against the reference, and
  a guide is `unique` when it has exactly one perfect site with the primary
  PAM; its other sites are counted by mismatches.
- **The PAM is read from the genome**, after the alignment, instead of being
  judged from the CIGAR. The old rules differed between strands for a
  deletion at either end of the spacer; the new alignments have no gaps.
- **Old coordinates could include the added G.** The old interval started at
  the leading G whenever GEM aligned it instead of soft-clipping it. The new
  ones never include it.

## The old files

`valid_alignments.bed` (no header, one row per valid alignment) maps to
`sites.tsv`, which has one row per site and a header:

| Old column | Content | Now |
| --- | --- | --- |
| 1 | Contig | `chr` |
| 2 | Protospacer start, 0-based | `start` |
| 3 | Protospacer end, exclusive | `end` |
| 4 | Protospacer sequence, from the read | `genomic_protospacer` holds the reference's bases; the guide's own spacer is `spacer` in `guides.tsv` |
| 5 | Mapping quality | None: every site is checked against the reference instead |
| 6 | Strand | `strand` |
| 7 | `NM` tag | `n_mismatches`, the spacer mismatches only, with their `mismatch_positions`; `NM` also counted the PAM's `N` |
| 8 | `AS` tag | None |
| 9 | Alias, or the read name | `guide_id` |

The unique guides' sites are also in `guides.tsv` (`guide_chr`,
`guide_start`, `guide_end`, `strand`) and in `guides.bed`.

`guide_alignment_log.tsv` (one row per read name) maps to `guides.tsv`:

| Old column | Now |
| --- | --- |
| `guide_id` | `guide_id` |
| `n_aligned` | `alignments_per_guide` in `run.json`: the alignment records GEM reported for the guide's reads |
| `n_valid` | `n_ngg_0mm` to `n_ngg_3mm`, by spacer mismatches; the alternative-PAM counts, `n_alt_pam_*`, are new |
| `n_discarded` | No per-guide count; the `site_search` rows of `summary.tsv` count the dropped alignments by reason |
| `n_not_mapped` | `alignment_class` `no_site` |

`alignment_summary.tsv` (`metric count`) maps to the `class` rows of
`summary.tsv`, with the meanings above:

| Old metric | Now |
| --- | --- |
| `guides_unique_valid` | `unique` |
| `guides_multi_valid` | `multi` |
| `guides_aligned_none_valid`, `guides_unmapped` | `off_target_only` or `no_site` |
| `guides_one_valid_plus_invalid` | None |

`discarded_alignments.tsv`, `invalid_alignments.tsv` and `unmapped.tsv` have
no successor: the dropped alignments are counted in `summary.tsv`, and the
raw alignments stay in `alignments/`. Nor do the valid-unique and valid-multi
SAM files.

## The old options

| `run_guide_alignment_qc.py` | `cau guide-alignment` |
| --- | --- |
| `--guides-tsv` | `run --guides` (tab-separated without a header by default; see `--sep`, `--header`) |
| `--reference-fasta` | `index REFERENCE PREFIX` once, then `run --reference` and `--index` |
| `--chromsizes` | `run --contigs`, a comma-separated list of contig names (below) |
| `--outdir`, `--threads`, `--pam`, `--add-leading-g` | The same names on `run` |
| `--mapping-mode` | None: the GEM options are fixed |
| `--allow-leading-g-softclip` | None: the added G is aligned, and a mismatch there never counts |

To keep the contigs a chrom.sizes file lists, join its first column and pass
the result to `run` as `--contigs "$contigs"`:

```bash
contigs="$(cut -f1 hg38.chrom.sizes | paste -sd, -)"  # chr1,chr2,...
```
<!-- --8<-- [end:migration] -->

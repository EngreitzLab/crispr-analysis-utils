# `filter_guide_alignments` outputs

<!-- --8<-- [start:tables] -->
All tables are tab-separated. "Read name" is the SAM QNAME, which is the guide
id when the FASTQ came from `guides_to_fastq`.

Each table has its own argument (`output_valid_bed`, `output_discarded_tsv`,
...). Left at `"auto"`, a table goes to the folder of `output_unique_sam`, else
of `output_multi_sam`, else the current working directory.

## `valid_alignments.bed`

No header. One row per valid alignment whose protospacer maps to the reference.

| Column | Content |
| --- | --- |
| 1 | Contig |
| 2 | Protospacer start, 0-based |
| 3 | Protospacer end, exclusive |
| 4 | Protospacer sequence, 5' to 3': the read without the PAM, and without the leading base when it is soft-clipped |
| 5 | Mapping quality |
| 6 | Strand of the alignment, `+` or `-` |
| 7 | `NM` tag (edit distance), or `-1` when absent |
| 8 | `AS` tag (alignment score), or `-1` when absent |
| 9 | `alias_by_guide_id[read sequence]` when that key exists, else the read name |

Columns 2 and 3 span the reference positions aligned to protospacer bases, so
the PAM is outside the interval.

## `discarded_alignments.tsv`

Header: `read_name chromosome pos1 flag mapq strand cigar NM AS MD reason`.
One row per mapped alignment that failed. `pos1` is the 1-based start; `NM` and
`AS` are `-1`, and `MD` is empty, when the tag is absent.

| `reason` | Meaning |
| --- | --- |
| `non_primary_contig` | Aligned to a contig outside `chromsizes` or `primary_contigs` |
| `discarded_tail_unaligned` | PAM not fully aligned, read no longer than the PAM, or a soft clip other than the leading base |
| `discarded_tail_mismatch` | Mismatch at a non-`N` position of the PAM |
| `discarded_protospacer_indel` | Insertion or deletion inside the spacer |

## `invalid_alignments.tsv`

Header: `guide_id contig reason`; the first column holds the read name. One row
per unmapped record (contig `.`, reason `unmapped`) and per discarded
alignment, with the reason before it is grouped as above: `non_primary_contig`,
`pam_not_fully_aligned`, `query_too_short_for_pam`, `softclip_not_allowed`,
`pam_gg_mismatch`, `spacer_insertion` or `spacer_deletion`.

## `unmapped.tsv`

Header: `read_name flag`. One row per unmapped record.

## `guide_alignment_log.tsv`

Header: `guide_id n_aligned n_valid n_discarded n_not_mapped`. One row per read
name, sorted.

| Column | Content |
| --- | --- |
| `guide_id` | Read name |
| `n_aligned` | Mapped alignments seen |
| `n_valid` | Mapped alignments that passed |
| `n_discarded` | Mapped alignments that failed |
| `n_not_mapped` | `1` if an unmapped record was seen, else `0` |

## `alignment_summary.tsv`

Header: `metric count`.

| Metric | Guides (read names) with |
| --- | --- |
| `guides_unique_valid` | exactly one valid alignment |
| `guides_multi_valid` | more than one valid alignment |
| `guides_aligned_none_valid` | mapped alignments, none of them valid |
| `guides_unmapped` | an unmapped record and no mapped alignment |
| `guides_one_valid_plus_invalid` | exactly one valid alignment and at least one discarded |

## The SAM outputs

Written only when both `output_unique_sam` and `output_multi_sam` are given.
The valid alignments of guides with exactly one go to the first; those of
guides with several go to the second.

## Return value

A dict with the five metrics of `alignment_summary.tsv`, plus:

| Key | Count |
| --- | --- |
| `valid_guides` | guides with at least one valid alignment |
| `unique_guides` | guides with exactly one valid alignment |
| `multi_guides` | guides with more than one valid alignment |
| `valid_alignments` | valid alignments |
| `invalid_alignments` | rows of `invalid_alignments.tsv`, unmapped records included |

<!-- --8<-- [end:tables] -->

# Guide QC

Example usage:

```python
import pandas as pd
import crispr_analysis_utils as cau

cau.guide_qc.guides_to_fastq(
  "guides.tsv",
  "guides.fastq",
  pam="NGG",
  add_leading_g=True
)

# DataFrame input with column overrides
guides = pd.DataFrame({"id": ["g1"], "seq": ["ACGT"]})
cau.guide_qc.guides_to_fastq(
  guides,
  "guides.fastq",
  id_col="id",
  sequence_col="seq"
)

# Filter GEM SAM/BAM alignments into valid unique vs valid multi-mapping
summary = cau.guide_qc.filter_guide_alignments(
  "guides_mapped_hg38.sam",
  "guides_valid_unique.sam",
  "guides_valid_multi.sam",
  pam="NGG",
  allow_leading_g_softclip=True  # default
  # Choose one for contig filtering:
  # chromsizes="hg38.chrom.sizes"
  # chromsizes=chromsizes_df
  # primary_contigs=["chr1", "chr2", "chrX", "chrY", "chrM"]
)
print(summary)
```

Filtering logic:

- Primary assembly contigs only
- PAM must be fully aligned
- Mismatches in PAM are not allowed except the first PAM base when using `N..` PAM (for `NGG`, `N` mismatch is accepted)
- Spacer cannot contain insertions/deletions
- Leading `G` soft-clipping is allowed by default

Output notes:

- Invalid alignments are always written to TSV
- Default invalid TSV path: `invalid_alignments.tsv` in the same folder as `guides_valid_unique.sam`
- You can override with `output_invalid_tsv="your_path.tsv"`
- Additional default outputs are written in the same folder as `guides_valid_unique.sam`

FASTQ generation note:

- `add_leading_g=True` prepends `G` only when the guide sequence does not already start with `G`

### `valid_alignments.bed`

Contains valid alignments in BED format with additional columns:

| Column | Description |
|--------|-------------|
| 1. chrom | Chromosome name |
| 2. chromStart | Start position (0-based) |
| 3. chromEnd | End position (0-based, exclusive) |
| 4. name | Guide read name |
| 5. score | Mapping quality (0-255) |
| 6. strand | Strand (`+` or `-`) |
| 7. NM | Number of mismatches (`NM` tag, `-1` if missing) |
| 8. AS | Alignment score (`AS` tag, `-1` if missing) |

### `discarded_alignments.tsv`

Contains mapped alignments that failed filtering:

| Column | Description |
|--------|-------------|
| read_name | Guide read name |
| chromosome | Chromosome name |
| pos1 | Start position (1-based) |
| flag | SAM flag |
| mapq | Mapping quality |
| strand | Strand (`+` or `-`) |
| cigar | CIGAR string |
| NM | Number of mismatches |
| AS | Alignment score |
| MD | MD tag (mismatch/deletion string) |
| reason | Discard reason |

**Discard reasons:**

- `discarded_tail_unaligned`: PAM region not properly aligned
- `discarded_tail_mismatch`: PAM region contains mismatches
- `discarded_protospacer_indel`: Protospacer region contains insertions or deletions

### `unmapped.tsv`

Contains unmapped records:

| Column | Description |
|--------|-------------|
| read_name | Guide read name |
| flag | SAM flag |

### `guide_alignment_log.tsv`

Per-guide summary counts:

| Column | Description |
|--------|-------------|
| guide_id | Guide/read name |
| n_aligned | Number of mapped alignments observed |
| n_valid | Number of mapped alignments that passed filters |
| n_discarded | Number of mapped alignments that failed filters |
| n_not_mapped | Unmapped indicator count (0 or 1) |

### `alignment_summary.tsv`

Run-level guide summary:

| Metric | Description |
|--------|-------------|
| guides_unique_valid | Number of guides with exactly one valid alignment |
| guides_multi_valid | Number of guides with multiple valid alignments |
| guides_aligned_none_valid | Number of guides with mapped alignments but zero valid alignments |
| guides_unmapped | Number of guides with no mapped alignments (`n_not_mapped = 1`) |
| guides_one_valid_plus_invalid | Number of guides with exactly one valid alignment and at least one discarded mapped alignment |

::: crispr_analysis_utils.guide_qc

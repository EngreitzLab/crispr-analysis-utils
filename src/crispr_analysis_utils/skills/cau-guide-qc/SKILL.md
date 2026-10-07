---
name: cau-guide-qc
description: >-
  Build guide FASTQs and filter guide-RNA alignments with crispr-analysis-utils
  (`cau.guide_qc`: guides_to_fastq and filter_guide_alignments). Use when a
  user wants to check where CRISPR guides (spacer plus PAM) align in a genome,
  keep only alignments with an intact PAM and no indel in the spacer, separate
  uniquely from multiply aligned guides, produce a protospacer BED, or
  interpret valid_alignments.bed, discarded_alignments.tsv,
  guide_alignment_log.tsv or alignment_summary.tsv. Pairs with cau-gem-mapper,
  which produces the SAM this module filters.
---

# Guide alignment QC (`cau.guide_qc`)

Two steps of the guide alignment QC workflow:

1. `guides_to_fastq` writes each guide as a synthetic FASTQ read: the spacer,
   an optional leading G, then the PAM. The `cau-gem-mapper` skill covers
   mapping those reads.
2. `filter_guide_alignments` reads the aligner's SAM or BAM, keeps the
   alignments that look like real target sites, and writes tables describing
   every guide.

The filter needs pysam: `pip install "crispr-analysis-utils[alignment]"`. The
repository's pixi environments already include it.

## Ask before running

- **The PAM.** The filter's default is `"NGG"` (SpCas9). Pass the *same* `pam`
  to both functions: `guides_to_fastq` defaults to `pam=""` and then appends
  nothing, while the filter treats the last `len(pam)` bases of every read as
  the PAM.
- **Whether a leading G was added** for a U6 promoter. `add_leading_g=True`
  prepends a G only to spacers that do not already start with one. The filter
  lets that single base be soft-clipped (`allow_leading_g_softclip=True`, the
  default).
- **Which contigs count.** Pass `chromsizes` (a path or a DataFrame whose first
  column lists the allowed contigs) or `primary_contigs`. With neither, every
  contig is accepted.
- **The genome build** the index was built from: the BED coordinates are in it.
- **Where the outputs go**, because of the first trap below.

## What makes an alignment valid

- The PAM is fully aligned (no soft clip or insertion there), with no mismatch
  except at `N` positions of the PAM.
- No insertion or deletion inside the spacer.
- No soft clip, except the single leading-G base when allowed.
- The same rules hold mirrored on the reverse strand, where the PAM is at the
  start of the read as SAM stores it.
- Mismatches inside the spacer are **not** filtered: read `NM` and `AS` in the
  BED. PAM mismatches are found through the `MD` tag (and `X` CIGAR
  operations), so an aligner that omits `MD` makes every PAM look clean.

## Usage

```python
import crispr_analysis_utils as cau

cau.guide_qc.guides_to_fastq(
    "guides.tsv",  # tab-separated, no header: guide id, spacer sequence
    "guides.fastq",
    pam="NGG",
    add_leading_g=True,
)
# Map guides.fastq to guides_mapped.sam (see the cau-gem-mapper skill), then:
summary = cau.guide_qc.filter_guide_alignments(
    "guides_mapped.sam",
    "qc/guides_valid_unique.sam",
    "qc/guides_valid_multi.sam",
    pam="NGG",
    chromsizes="hg38.chrom.sizes",
)
print(summary)  # counts; the tables are written next to the SAM outputs, in qc/
```

## Traps

- **Pass both SAM outputs, or an explicit path for every table.** The six
  tables go to the folder of `output_unique_sam`, else of `output_multi_sam`,
  and to the **current working directory** when neither is given.
- **The SAM files are written only when both** `output_unique_sam` and
  `output_multi_sam` are given. Passing just one writes no SAM, silently.
- **`alias_by_guide_id` is keyed by the full read sequence** as written to the
  FASTQ: uppercase spacer, any added G, and the PAM. When no key matches, BED
  column 9 is the read name, which is already the guide id when the FASTQ came
  from `guides_to_fastq`.
- **Per-guide counts are keyed by read name** (the SAM QNAME), so guide ids
  must be unique.

## Outputs

Every table, column by column, with the return value: `references/outputs.md`.

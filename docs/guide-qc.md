# Guide QC

`cau.guide_qc` prepares guides for alignment and filters the alignments:

1. `guides_to_fastq` writes each guide as a synthetic FASTQ read: the spacer,
   an optional leading G, then the PAM.
2. After mapping (see [GEM mapper](gem-mapper.md)), `filter_guide_alignments`
   keeps the alignments that look like real target sites and writes tables
   describing every guide.

The filter reads SAM and BAM with pysam, which installs with the package.

## Usage

```python
import pandas as pd

import crispr_analysis_utils as cau

cau.guide_qc.guides_to_fastq(
    "guides.tsv", "guides.fastq", pam="NGG", add_leading_g=True
)

# DataFrame input with column overrides
guides = pd.DataFrame({"id": ["g1"], "seq": ["ACGT"]})
cau.guide_qc.guides_to_fastq(guides, "guides.fastq", id_col="id", sequence_col="seq")

# Filter GEM SAM/BAM alignments into valid unique vs valid multi-mapping
summary = cau.guide_qc.filter_guide_alignments(
    "guides_mapped_hg38.sam",
    "guides_valid_unique.sam",
    "guides_valid_multi.sam",
    pam="NGG",
    allow_leading_g_softclip=True,  # default
    # Choose one for contig filtering:
    # chromsizes="hg38.chrom.sizes"
    # chromsizes=chromsizes_df
    # primary_contigs=["chr1", "chr2", "chrX", "chrY", "chrM"]
)
print(summary)
```

Pass the same `pam` to both functions: `guides_to_fastq` defaults to `pam=""`
and then appends no PAM, while the filter treats the last `len(pam)` bases of
each read as the PAM. `add_leading_g=True` prepends `G` only when the guide
does not already start with `G`.

## What makes an alignment valid

- It is on an allowed contig: one listed in `chromsizes` or `primary_contigs`.
  With neither, every contig is allowed.
- The PAM is fully aligned, with no mismatch except at its `N` positions (for
  `NGG`, the first base).
- No insertion or deletion inside the spacer.
- No soft clip, except the leading `G` base, which may be soft-clipped by
  default (`allow_leading_g_softclip=True`).
- On the reverse strand the rules are mirrored, with the PAM at the start of
  the read as SAM stores it, except for a deletion at either end of the
  spacer: one between spacer and PAM is accepted on `+` but rejected on `-`,
  and one right after the read's first base is rejected on `+` but accepted on
  `-` (with the leading-G soft clip allowed).

Mismatches inside the spacer are not filtered: read the `NM` and `AS` columns
of the BED.

## Where the outputs go

The SAM files are written only when both `output_unique_sam` and
`output_multi_sam` are given. Each table can be given its own path; one left
at `"auto"` goes to the folder of `output_unique_sam`, else of
`output_multi_sam`, and to the current working directory when neither is
given. The tables, one by one:

--8<-- "src/crispr_analysis_utils/skills/cau-guide-qc/references/outputs.md:tables"

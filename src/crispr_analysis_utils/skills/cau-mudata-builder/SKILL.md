---
name: cau-mudata-builder
description: >-
  Assemble a gene and a guide count matrix into one MuData with
  crispr-analysis-utils (the `cau mudata-builder` command, or
  `crispr_analysis_utils.mudata_builder`). Use when a user has kb or Cell
  Ranger output for a channel and wants the single object the IGVF CRISPR
  pipeline works in, when they ask how to attach gene and guide annotation
  to their counts, or when they want an h5mu to hand to a QC or inference
  step. Covers what goes in var, obs and uns, why cells are the union rather
  than the intersection, and what the step deliberately does not do.
---

# Build a MuData (`cau mudata-builder`)

```bash
cau mudata-builder run \
    --rna counts_unfiltered/adata.h5ad \
    --guide counts_unfiltered_modified/adata.h5ad \
    --guide-metadata results/igvf/guide_metadata.tsv \
    --store gencode.v43.annotation.parquet \
    --measurement-set IGVFDS8842PXCC \
    --gtf-accession IGVFFI9561BASO \
    --outdir results/
```

## Ask before running

- **Which GENCODE?** It must be the same release the guides were annotated
  against, or genes and guides describe different genomes. Ask rather than
  reach for whatever store is nearest. The run refuses a store other than the
  one the guide metadata's run record names, and warns when there is no
  record to check against.
- **Which channel?** `--measurement-set` goes into the cell index as
  `{barcode}_{measurement_set}`, so two channels concatenate without their
  barcodes colliding. There is no sensible default.
- **The portal ids.** `--gtf-accession` and `--genome-accession` take IGVF
  accessions; the CRISPR pipeline's inputs carry them. Without one the object
  records the file's name instead, so the field is never empty and never a
  guess.

## What it writes

`mudata/mudata.h5mu`, plus `run_metadata/` and `logs/`.

| Where | What |
| --- | --- |
| `mod['gene'].X` | the gene counts, untouched, with the kb layers beside them |
| `mod['gene'].var` | `symbol`, `gene_type`, `gene_chr`, `gene_start`, `gene_end` |
| `mod['guide'].X` | the guide counts, untouched |
| `mod['guide'].var` | every column of `guide_metadata.tsv` |
| `obs` | `barcode` and `measurement_set`, plus any `--covariates` given |
| `obsm['gene']`, `obsm['guide']` | which modality holds each cell |
| `uns['reference']` | the GTF and genome: name, md5 and accession |

## What it does not do

**Nothing is computed from the counts.** No per-cell or per-feature
statistics, no QC, no filtering, no doublet calls, no guide assignment, no
pairs to test. Those are later steps, and each is a judgement this one has no
business making.

## Four things worth knowing

- **Cells are the union, not the intersection.** A barcode seen in one
  modality and not the other is kept, and `obsm` records which modality holds
  it. On a real channel that is a million barcodes against 200,000 seen by
  both, so an intersection would discard most of the object before anyone
  looked at it. Filtering comes later, where it is visible.
- **Every gene type is annotated**, not only protein-coding. The counts come
  from a transcriptome built on the whole GTF, so two thirds of the features
  are lncRNA, pseudogenes and small RNAs; annotating only the coding genes
  would leave them nameless.
- **PAR genes appear twice.** 47 features carry a `_PAR_Y` suffix and have a
  chrX twin with identical sequence, so kallisto cannot tell their reads
  apart. Both rows are carried through untouched. Anything downstream that
  sums or de-duplicates them is making a choice; this step does not.
- **The ids carry GENCODE versions.** kb indexes genes as `ENSG00000290825.1`
  while the store drops the version, and the two are joined on the store's own
  rule, which strips a trailing `.N` and leaves `_PAR_Y` alone. A different
  rule silently gives the chrY copy the chrX gene's coordinates.

## Then what

A QC step filters this object; after that it is comparable with the IGVF
pipeline's own. See `cau-workflow` for where this sits, and the guide design
steps that produce `guide_metadata.tsv`.

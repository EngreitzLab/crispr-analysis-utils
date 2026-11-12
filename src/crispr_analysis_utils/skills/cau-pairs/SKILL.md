---
name: cau-pairs
description: >-
  Build the guide-gene pairs a CRISPR screen tests, discovery and positive
  control, with crispr-analysis-utils (the `cau pairs` command, or
  `crispr_analysis_utils.pairs`). Use when a user wants pairs_to_test for
  sceptre, pysceptre or the IGVF CRISPR pipeline, asks which genes an element
  is tested against, which pairs are positive controls, why a readout gene has
  no pair, or wants to change the cis window. Reads the IGVF per-guide table
  from `cau igvf-export`, so the pairs name the same elements the submission
  does.
---

# Pairs to test (`cau pairs`)

```bash
cau pairs run \
    --guide-metadata results/igvf/guide_metadata.tsv \
    --genes readout_genes.csv \
    --store gencode.v43.annotation.parquet \
    --outdir results/
```

## Ask before running

- **Which genes are read out?** `--genes` is the panel: a CSV or TSV whose
  first column of ENSEMBL ids is used, versions dropped. For a
  whole-transcriptome screen, ask which genes the user means to test rather
  than passing all of them.
- **Which annotation?** The store must be the GENCODE release the guide
  table was built with. A panel gene the store lacks stops the run.
- **How far is cis?** `--max-distance` defaults to 1,000,000 bases from the
  element's midpoint to the gene's TSS. Ask if the user has a different
  window in mind.
- **Are the positive controls right?** They come from the guide table's
  `type` and `putative_target_genes`, which `cau igvf-export` filled from its
  `--genes` pairs file, falling back to the nearest TSS's gene for a
  promoter. A validated enhancer missing from that file names no gene and
  gives no pair; the log warns. Fix the file and re-export: this step will
  not invent a gene.

## The two pair types

| `pair_type` | From |
| --- | --- |
| `positive_control` | a guide typed `positive control`, with every gene in its `putative_target_genes` |
| `discovery` | a targeting guide, with every panel gene whose TSS is within `--max-distance` of its element's midpoint |

- **A pair is decided per element.** Where an element and gene are a
  positive control's, every guide of that element pairs with that gene as
  `positive_control`. The same element still pairs with its other
  neighbours as `discovery`.
- **A peak near a TSS is discovery.** The IGVF table calls an element within
  2 kb of a TSS a `promoter`, but only a guide typed `positive control` gives a
  positive pair. Proximity is not a known effect.
- **Positive-control genes outside the panel are kept**, with `in_panel`
  FALSE. Whether they can be measured is the inference's question.
- **Controls get no pairs.** The inference calibrates on them.
- **The TSS is the reference transcript's**: MANE Select, else Ensembl
  canonical, else the longest CDS, the same transcript `cau guide-annotation`
  reports. `distance_to_tss` is signed in the gene's orientation.

## What it writes

`pairs/pairs_to_test.csv`, comma-separated because the CRISPR pipeline's
`INFERENCE_predefined_pairs_to_test` reads it with `read_csv`. One row per
guide and gene: `guide_id`, `gene_id`, `gene_name`, the four
`intended_target_*` columns, `pair_type`, `distance_to_tss`, `in_panel`.
Beside it: `pairs_to_test_columns.yaml`, `qc.tsv` (counts per pair type)
and `panel_genes_without_pairs.tsv`. Plus `run_metadata/` and `logs/`.

## Then what

The pairs key on `guide_id` and an unversioned `gene_id`, so the count
matrix must use the same ids. Rows whose guide or gene the matrix lacks are
the inference step's to report, not this one's.

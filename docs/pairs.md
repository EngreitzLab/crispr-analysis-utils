# Pairs to test

`cau pairs` turns the [IGVF per-guide table](igvf-export.md) into the
guide-gene pairs an inference tests: positive controls, whose effect is
known, and discovery pairs, whose effect is the question.

In Claude Code the `cau-pairs` skill teaches Claude to run it and read the
output ([Agents and skills](agents-and-skills.md)).

## Quick start

```bash
cau pairs run \
    --guide-metadata results/igvf/guide_metadata.tsv \
    --genes readout_genes.csv \
    --store gencode.v43.annotation.parquet \
    --outdir results/
```

`--genes` is the readout panel, a CSV or TSV with a column of ENSEMBL ids.
`--max-distance` sets the cis window, 1,000,000 bases by default. Nothing is
written outside `--outdir`.

## How a pair is decided

- **Positive control**: a guide typed `positive control`, with every gene in
  its `putative_target_genes`. `cau igvf-export` types a guide that way only
  at a design promoter the genome confirms or at a validated enhancer. The
  genes are the ones its `--genes` pairs file measured the element against;
  a promoter without such a pair gets its nearest TSS's gene.
- **Discovery**: a targeting guide, with every panel gene whose TSS lies
  within the window of its element's midpoint.
- Pairs are decided per element. An element and gene that are a positive
  control's are `positive_control` for every guide of that element, and the
  element's other neighbours stay `discovery`.
- A peak within 2 kb of a TSS is `promoter` in the IGVF table but gives
  discovery pairs, not positive ones.
- The TSS is the reference transcript's (MANE Select, then Ensembl
  canonical, then the longest CDS), as in [gene annotation](guide-annotation.md).

## Outputs

`pairs/pairs_to_test.csv` is comma-separated, the format the CRISPR
pipeline's `INFERENCE_predefined_pairs_to_test` reads. One row per guide and
gene:

| Column | Holds |
| --- | --- |
| `guide_id` | the guide, as `guide_metadata.tsv` names it |
| `gene_id` | ENSEMBL gene id, without its version |
| `gene_name` | the gene's symbol in the store |
| `intended_target_name`, `_chr`, `_start`, `_end` | the element, 0-based half-open |
| `pair_type` | `positive_control` or `discovery` |
| `distance_to_tss` | element midpoint to TSS, signed in the gene's orientation |
| `in_panel` | whether the gene is in `--genes` |

Beside it are `pairs_to_test_columns.yaml`, `qc.tsv` (pairs, guides,
elements and genes per pair type) and `panel_genes_without_pairs.tsv`.

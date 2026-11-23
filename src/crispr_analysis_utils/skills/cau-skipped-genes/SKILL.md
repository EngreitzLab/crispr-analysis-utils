---
name: cau-skipped-genes
description: >-
  Classify every discovery and positive-control pair by the genes its link
  passes over, with crispr-analysis-utils (the `cau skipped-genes` command,
  or `crispr_analysis_utils.skipped_genes`). Use after `cau pairs` when a
  user asks whether a pair's gene is really the nearest one to the element,
  how many other genes sit between the element and the paired gene, or
  whether a skipped gene is in the readout panel. Reads pairs_to_test.csv and
  the same annotation store and readout panel `cau pairs` used.
---

# Skipped genes (`cau skipped-genes`)

```bash
cau skipped-genes run \
    --pairs results/pairs/pairs_to_test.csv \
    --store gencode.v43.annotation.parquet \
    --genes readout_genes.csv \
    --outdir results/
```

Ported from the intervening-gene classification in the Voineagu lab's
`EnhGenePairs.R` (astrocyte_crispri), using this repository's own annotation
store and readout panel instead of its external gene table. Its TAD analysis
is not ported: this package has no model for a TAD call file.

## Ask before running

- **Which pairs?** `--pairs` is `pairs_to_test.csv` from `cau pairs run`.
- **Which annotation and panel?** `--store` and `--genes` must be the same
  ones `cau pairs` used to build those pairs; a different release or panel
  gives a classification that does not match the pairs it labels.

## How a pair is classified

For each pair's element, every annotated gene on its chromosome is ranked by
distance from the element's midpoint, in genome coordinates (not flipped to
the gene's orientation, the way `distance_to_tss` is in `pairs_to_test.csv`).
A gene is "skipped" when it is at least as close as the pair's own gene and
on the same side of the element: the pair's link passes over it.

| `category` | Meaning |
| --- | --- |
| `nearest` | No annotated gene has a closer TSS. |
| `no_skip` | Nearer genes exist, but none on this gene's side. |
| `skip_in_panel` | At least one skipped gene is in the readout panel. |
| `skip_outside_panel` | Every skipped gene is outside the readout panel. |
| `trans` | The gene is not on the element's chromosome (a trans positive control); nothing is ranked. |

- **Ranking uses every annotated gene**, not only the readout panel: a
  skipped gene outside the panel is still a gene the design could have
  measured against instead.
- **"In panel" means in the file given to `--genes`**, the same readout
  panel `cau pairs` used for `in_panel` in `pairs_to_test.csv`.

## What it writes

`skipped_genes/skipped_genes.csv`, comma-separated. One row per element and
gene pair: the four `intended_target_*` columns, `gene_id`, `gene_symbol`,
`pair_type`, `distance_to_tss`, `in_panel` (all from `pairs_to_test.csv`),
then `n_nearer`, `n_skipped`, `n_skipped_in_panel`, `n_skipped_outside_panel`
and `category`. The four counts are empty for a `trans` pair.
Beside it: `skipped_genes_columns.yaml` and `qc.tsv` (pairs per category).
Plus `run_metadata/` and `logs/`.

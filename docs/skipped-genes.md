# Skipped genes

`cau skipped-genes` classifies every pair from [pairs to test](pairs.md) by
the genes its link passes over: for each element, every annotated gene on its
chromosome is ranked by distance, and a nearer gene on the same side as the
pair's own gene is "skipped".

Ported from the intervening-gene classification in the Voineagu lab's
[`EnhGenePairs.R`](https://github.com/Voineagulab/astrocyte_crispri/blob/main/5_FunctionalAnnotation/EnhGenePairs/EnhGenePairs.R),
using this repository's own annotation store and readout panel instead of its
external gene table. The script's TAD analysis is not ported: it needs a TAD
call file this package has no model for.

In Claude Code the `cau-skipped-genes` skill teaches Claude to run it and read
the output ([Agents and skills](agents-and-skills.md)).

## Quick start

```bash
cau skipped-genes run \
    --pairs results/pairs/pairs_to_test.csv \
    --store gencode.v43.annotation.parquet \
    --genes readout_genes.csv \
    --outdir results/
```

`--store` and `--genes` must be the same annotation store and readout panel
`cau pairs` used to build `--pairs`. Nothing is written outside `--outdir`.

## How a pair is classified

For each pair's element, every annotated gene on its chromosome is ranked by
distance from the element's midpoint, in genome coordinates: not flipped to
the gene's orientation, the way `distance_to_tss` is. A gene is "skipped"
when it is at least as close as the pair's own gene and sits on the same
side of the element.

| `category` | Meaning |
| --- | --- |
| `nearest` | No annotated gene has a closer TSS. |
| `no_skip` | Nearer genes exist, but none on this gene's side. |
| `skip_in_panel` | At least one skipped gene is in the readout panel. |
| `skip_outside_panel` | Every skipped gene is outside the readout panel. |
| `trans` | The gene is not on the element's chromosome; nothing is ranked. |

Ranking uses every annotated gene, not only the readout panel: a skipped
gene outside the panel is still a gene the design could have measured
against instead.

## Outputs

`skipped_genes/skipped_genes.csv` is comma-separated. One row per unique
element and gene pair:

| Column | Holds |
| --- | --- |
| `intended_target_name`, `_chr`, `_start`, `_end` | the element, from `pairs_to_test.csv` |
| `gene_id`, `gene_symbol` | the gene |
| `pair_type`, `distance_to_tss`, `in_panel` | carried through from `pairs_to_test.csv` |
| `n_nearer` | annotated genes at least as close to the element as this one |
| `n_skipped` | of those, how many are on this gene's side |
| `n_skipped_in_panel`, `n_skipped_outside_panel` | the skipped genes, split by readout panel membership |
| `category` | `nearest`, `no_skip`, `skip_in_panel`, `skip_outside_panel` or `trans` |

Beside it are `skipped_genes_columns.yaml` and `qc.tsv` (pairs per category).

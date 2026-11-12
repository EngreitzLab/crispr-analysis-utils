# Inference input

`cau inference-input` joins QC'd channels that already carry guide calls into
one MuData laid out like the IGVF CRISPR pipeline's inference input, with the
[pairs to test](pairs.md) attached. One object per call: a day, a condition
or the whole screen, as the analysis groups it.

In Claude Code the `cau-inference-input` skill teaches Claude to run it and
read the output ([Agents and skills](agents-and-skills.md)).

## Quick start

```bash
cau inference-input run \
    --inputs qc/R1-D4-A-b1.h5mu qc/R2-D4-A-b1.h5mu \
    --pairs results/pairs/pairs_to_test.csv \
    --outdir results/ \
    --name day4 \
    --capture-method CROP-seq
```

Each input is a channel from [`cau mudata-builder`](mudata-builder.md) after
the QC pass, with binary guide calls in
`mod['guide'].layers['guide_assignment']`. This step does not call guides.
`--moi` defaults to `high`. `--min-guides` and `--max-guides` filter cells by
the number of guides called in them and have no default; the pipeline uses 1
and 15. Nothing is written outside `--outdir`.

## The object

| Where | What |
| --- | --- |
| cell index | `{barcode}_{measurement_set}`, the cells both modalities hold |
| `obs` | `batch` (the measurement set), `total_gene_umis`, `total_guide_umis`, `num_expressed_genes`, `percent_mito` |
| `mod['gene'].var_names` | ENSEMBL ids without their version, except a `_PAR_Y` id, left whole; `gene_id_versioned` keeps the original |
| `mod['gene'].obs` | the covariates above, plus `pct_counts_ribo` and the channel's own columns |
| `mod['guide'].var` | `guide_id`, the IGVF guide columns, and `intended_target_key` |
| `mod['guide'].obs` | `total_guide_umis`, `num_expressed_guides` |
| `mod['guide'].layers['guide_assignment']` | the guide calls, as given |
| `mod['guide'].uns` | `moi`, `capture_method` |
| `uns['pairs_to_test']` | the pairs whose guide and gene the object holds, keyed by `intended_target_key` |
| `uns['qc']` | each channel's QC record, by measurement set |
| `uns['inference_input']` | the cells and pairs this step dropped |

Control guides (`targeting` false, or typed non-targeting, safe-targeting or
negative control)
are grouped as `non-targeting|N` by the pipeline's own rule: sorted by guide
id and cut into groups the size of the median targeting element. On one
channel of a 43,766-guide screen, the keys and groups match the pipeline's
`annotate_intended_target_groups` guide for guide, and the pipeline's prepare
step keeps the same pairs and covariates.

A pair whose element disagrees with its guide's stops the run: the pairs and
the counts were built from different guide tables. So does a cell with no gene
UMIs, no expressed genes or no guide UMIs, since SCEPTRE fits the log of each
and the pipeline refuses a zero.

## Then what

pysceptre reads the object with its pairs. The pipeline's
`--INFERENCE_input_mudata` entry also takes it, but that entry requires the
`default` pairing strategy, which builds its own pairs from distance, so the
attached ones are not used.

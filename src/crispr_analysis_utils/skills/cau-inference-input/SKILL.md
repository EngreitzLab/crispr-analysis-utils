---
name: cau-inference-input
description: >-
  Join QC'd CRISPR screen channels that already carry guide calls into one
  MuData laid out like the IGVF CRISPR pipeline's inference input, with
  crispr-analysis-utils (the `cau inference-input` command, or
  `crispr_analysis_utils.inference_input`). Use when a user wants the object
  sceptre or pysceptre runs on, wants a day's or a condition's channels in
  one h5mu, asks how controls are grouped as non-targeting|N, asks why a
  pair was dropped, or wants to filter cells by how many guides they carry.
  Takes the pairs from `cau pairs` and the channels from
  `cau mudata-builder` and the `cau-qc` pass.
---

# Inference input (`cau inference-input`)

```bash
cau inference-input run \
    --inputs qc/R1-D4-A-b1.h5mu qc/R2-D4-A-b1.h5mu \
    --pairs results/pairs/pairs_to_test.csv \
    --outdir results/ \
    --name day4 \
    --capture-method CROP-seq
```

## Ask before running

- **Which channels make one object?** Every file in `--inputs` goes into one
  h5mu, and they must share their genes and guides, in order. Ask how the
  user groups them (per day, per condition, all together).
- **Are the guides called?** Each channel needs a binary
  `mod['guide'].layers['guide_assignment']`. This step does not call guides,
  and refuses a channel without the layer. How to call them (a UMI threshold,
  a mixture model) is the user's decision; ask.
- **How were guides captured?** `--capture-method` has no default. It is
  carried into `mod['guide'].uns`, not read.
- **High or low MOI?** `--moi` defaults to `high`.
- **Filter cells by guide count?** `--min-guides` and `--max-guides` have no
  default: without them every cell is kept. The CRISPR pipeline uses 1 and
  15. A minimum of 1 drops every cell with no call, which at a strict
  threshold can be a large share; ask, and report how many each would drop.
  A cell with no guide UMIs is kept but counted in the log: the inference fits
  `log(total_guide_umis)`, so such cells must be filtered before it runs.
- **Same guide table?** The pairs must come from the guide table the
  channels carry. A pair whose element disagrees with its guide's stops the
  run.

## What it does

- Keeps the cells both modalities hold, then concatenates the channels.
  Cell ids are `{barcode}_{measurement_set}` already, so they cannot clash.
- Drops the version from gene ids (`gene_id_versioned` keeps it). A
  `_PAR_Y` id is left whole, version included, so the index stays unique.
- Adds `intended_target_key` to every guide and buckets the controls as
  `non-targeting|N`, by the pipeline's rule: control guides (`targeting`
  false, or typed non-targeting, safe-targeting or negative control) sorted by guide id and
  cut into groups the size of the median targeting element.
- Writes the covariates under the pipeline's names. Root `obs`: `batch`
  (the measurement set), `total_gene_umis`, `total_guide_umis`,
  `num_expressed_genes` and `percent_mito`. `pct_counts_ribo` goes in
  `mod['gene'].obs`, and `num_expressed_guides` in `mod['guide'].obs`.
- Sets `mod['guide'].uns['moi']` and `['capture_method']`.
- Stores `uns['pairs_to_test']`, restricted to the guides and genes present.
  The rows it drops are counted in the log, the run record and
  `uns['inference_input']`, since the pipeline drops them silently.

## What it writes

`inference_input/<name>.h5mu`, plus `run_metadata/` and `logs/`.

## Then what

pysceptre reads the object with its pairs. The pipeline's
`--INFERENCE_input_mudata` entry also takes it, but that entry requires the
`default` pairing strategy and builds its own pairs, so ours are not used.

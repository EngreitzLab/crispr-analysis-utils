---
name: cau-element-mapping
description: >-
  Map a CRISPR guide library to the elements a screen targets with
  crispr-analysis-utils (the `cau element-mapping` command, or
  `crispr_analysis_utils.element_mapping`). Use after `cau guide-alignment`
  and `cau guide-annotation` when a user wants a guide-to-element mapping for
  sceptre or another test: which guides target which element, which are
  non-targeting because they match nowhere, which are safe-targeting because
  they sit far from every element and every promoter, and which map but hit
  none of the elements given. Also for reading guide_element_mapping.tsv, and
  for checking a library's own control labels against where its guides
  actually land.
---

# Element mapping (`cau element-mapping`)

The third step. Alignment says where a guide is, annotation says what gene it
is in, and this says which element of the screen it tests.

```bash
cau element-mapping run \
    --alignment results/alignments/guide_genome_annotation.tsv \
    --annotation results/annotations/guide_gene_annotation.tsv \
    --elements elements.bed \
    --outdir results/
```

## Ask before running

- **Which element list, and is it the whole one?** A screen usually has
  several design files. Pass the combined list: with a partial one, every
  guide for a missing element comes back `safe_targeting`, which looks like a
  result rather than a missing input. Check that count against what the user
  expects before believing the table.
- **Is `--annotation` available?** Without it `distance_to_tss` is empty.
  The classes do not depend on it, so the command warns but still runs.

## The three classes

Per **site**, not per guide, so a guide whose sites disagree reports each one.

| Class | Means |
| --- | --- |
| `targeting` | The protospacer overlaps this element |
| `safe_targeting` | The guide maps, but overlaps no element |
| `non_targeting` | The guide matches nowhere in the genome |

**`safe_targeting` is not a promise.** It says only that the site is outside
every element given. How far from a gene it has to be before it is safe
depends on what the screen reads out, so that is left as `distance_to_tss`
rather than fixed by a threshold here. Threshold it yourself, and say which
threshold you used.

A large `safe_targeting` count can also mean the element list was incomplete,
so check it against what the user expects before believing the table.

## Reading the result

`--outdir` gets `elements/`, `run_metadata/` and `logs/`. The results are one
table, `guide_element_mapping.tsv`, the `guide_element_mapping_columns.yaml`
that describes it, and
`qc.tsv`. See [`references/outputs.md`](references/outputs.md) for the columns
in full.

- **One row per site per element.** Elements may overlap each other, so a
  single site can appear against two of them.
- **Any overlap counts, down to one base.** For CRISPRi that is right: the
  dCas9-KRAB footprint covers the whole protospacer, so a guide hanging off
  an element's edge still represses it. For a cutting screen the target site
  would be the better test, and roughly 4% of rows are partial overlaps, so
  say which kind of screen it is before reading the counts.
- **`--genes` fills `putative_target_genes`**, from the same pairs table
  `cau igvf-export run` takes and through the same lookup, so the two tables
  name the same genes. Without it the column is empty; with it, only elements
  the pairs table covers get one.
- **`element_id` comes from the BED's fourth column** when it has one, so an
  id used elsewhere keeps its spelling. Otherwise it is `chr:start-end`,
  1-based, as a locus is read aloud.
- **That id is 1-based over a 0-based interval, and the two cannot be
  concatenated.** `chr:start-end` is what a genome browser takes, so its
  start is the BED's `Start` plus one. A table that pastes a BED's own
  columns back together is 0-based and will miss every element by a single
  base: all 26 validated enhancers of one library matched nothing that way,
  silently, until someone asked where the target genes had gone.
- **It keys back** to `guide_genome_annotation.tsv` on `guide_id` and the
  four site columns.

## Checking a library's own labels

This is the step that tests a design's control sets, so expect disagreement
and report it rather than smoothing it over:

- a guide the design calls non-targeting that lands inside an element is a
  real finding, not a bug;
- a guide the design calls safe-targeting whose `distance_to_tss` is small
  usually means the design measured safety against its own readout genes, not
  against every promoter. Report the distance rather than the label.

Say which is which, with counts, instead of reporting only the totals.

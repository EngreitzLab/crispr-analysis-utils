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
    --alignment results/alignments/guides_alignment_metadata.tsv \
    --annotation results/annotations/guides_gene_annotation.tsv \
    --elements elements.bed \
    --outdir results/
```

## Ask before running

- **Which element list, and is it the whole one?** A screen usually has
  several design files. Pass the combined list: with a partial one, every
  guide for a missing element comes back `unassigned`, which looks like a
  result rather than a missing input. Check the `unassigned` count against
  what the user expects before believing the table.
- **Is `--annotation` available?** Without it nothing can be called
  `safe_targeting`, and those guides land in `unassigned` instead. The
  command warns but still runs.
- **What safe distance?** `--safe-distance` is 1 Mb by default. Ask: a
  library's own safe-targeting controls are often much closer than that, so
  the default can disagree with the design on purpose.

## The four classes

Per **site**, not per guide, so a guide whose sites disagree reports each one.

| Class | Means |
| --- | --- |
| `targeting` | The protospacer overlaps this element |
| `safe_targeting` | Overlaps no element, and further than `--safe-distance` from any protein-coding TSS |
| `non_targeting` | The guide matches nowhere in the genome |
| `unassigned` | Maps, overlaps no element, but is nearer than the safe distance |

`unassigned` is the honest answer for "maps, but to none of the elements you
gave me". A large count means the element list was incomplete, the safe
distance is too strict, or the library targets more than this screen.

## Reading the result

`--outdir` gets `elements/`, `run_metadata/` and `logs/`. The results are one
table, `guide_element_mapping.tsv`, and the `columns.yaml` that describes it.

- **One row per site per element.** Elements may overlap each other, so a
  single site can appear against two of them.
- **`element_id` comes from the BED's fourth column** when it has one, so an
  id used elsewhere keeps its spelling. Otherwise it is `chr:start-end`,
  1-based, as a locus is read aloud.
- **It keys back** to `guides_alignment_metadata.tsv` on `guide_id` and the
  four site columns.

## Checking a library's own labels

This is the step that tests a design's control sets, so expect disagreement
and report it rather than smoothing it over:

- a guide the design calls non-targeting that lands inside an element is a
  real finding, not a bug;
- a guide the design calls safe-targeting that is `unassigned` here usually
  means the design's idea of safe was a shorter distance than
  `--safe-distance`, not that either is wrong.

Say which is which, with counts, instead of reporting only the totals.

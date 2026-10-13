# Element mapping

`cau element-mapping` says which element of a screen each guide tests. It is
the third step: [alignment](guide-alignment.md) finds where a guide is,
[annotation](guide-annotation.md) says what gene it falls in, and this maps it
to the regions the experiment targets.

In Claude Code the `cau-element-mapping` skill teaches Claude to run it and
read the output ([Agents and skills](agents-and-skills.md)).

## Quick start

```bash
cau element-mapping run \
    --alignment results/alignments/guides_alignment_metadata.tsv \
    --annotation results/annotations/guides_gene_annotation.tsv \
    --elements elements.bed \
    --outdir results/
```

`--annotation` is optional and fills `distance_to_tss`; the classes do not
depend on it. Nothing is written outside `--outdir`.

## The three classes

Decided per **site**, not per guide, so a guide whose sites disagree reports
each of them.

| Class | Means |
| --- | --- |
| `targeting` | The protospacer overlaps this element |
| `safe_targeting` | The guide maps, but the site overlaps no element |
| `non_targeting` | The guide matches nowhere in the genome |

There is deliberately no fourth class and no distance threshold. An earlier
draft called a guide safe only beyond 1 Mb from any protein-coding TSS, which
put most of a real library's safe-targeting controls somewhere else: that
design had measured safety against the screen's own readout genes, not every
promoter. Both definitions are defensible, which is the point — the tool
reports `distance_to_tss` and lets the analysis choose.

## Outputs

--8<-- "src/crispr_analysis_utils/skills/cau-element-mapping/references/outputs.md:tables"

## Use it to check the design

A library's own labels are a hypothesis. Mapping a real 43,778-guide library
against its element list disagreed with the shipped annotation in two ways,
in opposite directions:

- **18 guides labelled non-targeting cut inside a targeted element.** All 18
  are exact 20-nt protospacers with a real NGG, read straight back out of the
  reference. The aligner the old pipeline used had missed them, so genuine
  cutters were used as controls.
- **45 guides assigned to an element do not reach it.** That pipeline stripped
  a 5' G before aligning and so searched a 19-mer, which matches places the
  real 20-mer does not.

Both contaminate a downstream test, in opposite directions, and neither is
visible without redoing the mapping from the sequences.

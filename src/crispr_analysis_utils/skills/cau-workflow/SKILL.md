---
name: cau-workflow
description: >-
  Route a CRISPR guide library through crispr-analysis-utils from designed
  spacers to an IGVF submission, in the right order. Use when a user names a
  goal rather than a step, such as submitting a library to the IGVF portal,
  preparing the guide_metadata.tsv the CRISPR Perturb-seq pipeline reads,
  finding out what their guides target, or checking a design's own labels.
  Also when a step's inputs are not on disk yet and the step that produces
  them has to run first, or when a user asks which of the cau commands they
  need. Covers the order, what each step consumes and produces, the keys the
  tables join on, and what to settle with the user before anything runs.
---

# The crispr-analysis-utils workflow

Four commands, run in order. Each reads the previous one's output, so the
order is a data dependency, not a convention. The per-step skills hold the
detail; this one decides which step the user is actually asking for.

| Step | Command | Skill |
| --- | --- | --- |
| 1. Where does each guide match the genome? | `cau guide-alignment run` | `cau-guide-alignment` |
| 2. What gene is each site in? | `cau guide-annotation run` | `cau-guide-annotation` |
| 3. Which element of the screen does it test? | `cau element-mapping run` | `cau-element-mapping` |
| 4. Write the IGVF table | `cau igvf-export run` | `cau-igvf-export` |

Two modules sit outside the chain and need no predecessor:
`cau-accession` (a lab-independent id encoded from the protospacer) and
`cau-normalization` (counts per million, for a screen's count matrix).

## The chain

```bash
cau guide-alignment run \
    --guides guides.tsv --reference hg38.fa --outdir results/

cau guide-annotation build-db gencode.v43.annotation.gtf.gz \
    --out gencode.v43.protein_coding.parquet

cau guide-annotation run \
    --alignment results/alignments/guide_genome_annotation.tsv \
    --store gencode.v43.protein_coding.parquet \
    --promoters gencode.v43.protein_coding.TSS1bp.bed6 \
    --outdir results/

cau element-mapping run \
    --alignment results/alignments/guide_genome_annotation.tsv \
    --annotation results/annotations/guide_gene_annotation.tsv \
    --elements elements.bed \
    --outdir results/

cau igvf-export run \
    --alignment results/alignments/guide_genome_annotation.tsv \
    --mapping results/elements/guide_element_mapping.tsv \
    --elements elements.bed \
    --store gencode.v43.protein_coding.parquet \
    --promoters gencode.v43.protein_coding.TSS1bp.bed6 \
    --genes measured_element_gene_pairs.tsv \
    --outdir results/
```

`build-db` is a one-off per GENCODE release, not per library. Reuse the
parquet store and the TSS file across runs and across steps 2 and 4.

## What flows between the steps

`guide_genome_annotation.tsv` is the spine: one row per library entry per
site, and the only table with the whole picture. Everything downstream reads
it, and steps 3 and 4 read it directly rather than through step 2's output.

A guide with two sites has two rows, so joining on `guide_id` alone is
ambiguous. The annotation table keys back on `guide_id` plus the four site
columns; use all five.

## Settle these before running anything

Ask, do not guess:

- **The reference and the annotation must be the same assembly**, and the
  GENCODE release decides the ENSEMBL ids the IGVF table carries.
- **Does the library have an added 5' G?** It changes what is aligned and
  how a duplicate is judged. `cau-guide-alignment` has the rule.
- **Is the element BED the screen's real element set?** An incomplete one
  inflates `safe_targeting` silently rather than failing, so the count is
  checked against the design before it is believed.
- **Are there validated enhancers?** They are the one case where a design's
  own label for an element outranks what the genome says.

## Picking the step from what the user has

Work back from the output they want, then forward from the files that exist.
A step whose inputs are missing is a request for the step before it.

- Only a guide table: start at 1.
- Alignment output but no gene calls: 2.
- Alignment plus an element BED, and the question is what each guide tests
  or which guides are controls: 3.
- All three, and the goal is a submission or the pipeline's input: 4.
- A question about one guide's identity across libraries: `cau-accession`,
  no chain needed.
- A count matrix to normalize: `cau-normalization`, no chain needed.

---
name: cau-igvf-export
description: >-
  Write a CRISPR guide design in the IGVF per-guide metadata format with
  crispr-analysis-utils (the `cau igvf-export` command, or
  `crispr_analysis_utils.igvf_export`). Use when a user wants to submit a
  guide library to the IGVF portal as a tabular file with content_type
  "guide RNA sequences", or wants the guide_metadata.tsv the IGVF CRISPR
  Perturb-seq pipeline takes as input. Covers the columns the format
  requires, what each one holds, how a guide's genomic_element and
  intended_target_name are decided, and the traps that make a submission
  wrong but valid-looking.
---

# IGVF export (`cau igvf-export`)

The fourth step. Alignment, annotation and element mapping between them know
everything the format asks for; this assembles it.

```bash
cau igvf-export run \
    --alignment results/alignments/guide_genome_annotation.tsv \
    --mapping results/elements/guide_element_mapping.tsv \
    --elements elements.bed \
    --store gencode.v43.protein_coding.parquet \
    --promoters gencode.v43.protein_coding.TSS1bp.bed6 \
    --genes measured_element_gene_pairs.tsv \
    --outdir results/
```

## The format

`Format_for_Per-Guide_Metadata.pdf`, IGVF document
`f6a915a9-a6b8-4015-8aeb-236764976979`. Tabular files with `content_type`
"guide RNA sequences" must follow it. Its conventions are the ones this
package already uses: 0-based, half-open, ACTG, and **ENSEMBL ids at GENCODE
release 43**, so build the store from v43 unless the user says otherwise.

## Ask before running

- **Which GENCODE?** The format pins ENSEMBL ids to release 43. A store built
  from another release produces ids that will not validate.
- **Is there a file of measured element-to-gene pairs?** `--genes` takes one.
  Without it, a promoter's gene is inferred from the nearest TSS, which is
  reasonable but is not the same as a measured pair.
- **Is the screen CRISPRi, CRISPRa or cutting?** It does not change the file,
  but it changes whether "a guide overlapping an element" is the right
  definition of targeting, which the element mapping decided earlier.

## What is decided, and how

**`genomic_element` comes from the genome, not from the element's label.** A
design's own category for an element is a hypothesis; the annotation is
evidence. In order: a validated enhancer stays an enhancer, then within 2 kb
of a protein-coding TSS is a `promoter`, then overlapping an exon is `exon`,
then inside a gene is `intron`, and the rest is `distal element`. The 2 kb
is the format's own definition of a distal element.

The one exception is `e2g_validated`: that an element was shown to regulate a
gene is experimental knowledge no annotation can supply, so that label is
trusted.

**`type` must agree with `genomic_element`.** `positive control` means a
guide expected to produce a known measurable effect, so a row claiming it
while `genomic_element` says `exon` contradicts itself. A guide is a positive
control only when its element is a validated enhancer, or the design calls it
a promoter and the genome agrees.

**`intended_target_name` depends on the element:**

| `genomic_element` | holds |
| --- | --- |
| `promoter` | the regulated gene's ENSG |
| `exon` | the ENSEMBL exon id |
| `intron` | the intron's coordinates |
| `enhancer`, `distal element` | the element's coordinates |

## genomic_element does not say genic

`genomic_element` is decided by precedence: promoter, then exon, then intron,
then distal element. Promoter and enhancer win before the gene body is ever
tested, so the class cannot say whether the element sits in a gene.
`element_location` does, from the same overlap. On one library 499 promoter
and 110 enhancer guides are intergenic, while every exon and intron guide is
genic.

## Reading the result

`--outdir` gets `igvf/`, `run_metadata/` and `logs/`. The results are one
table, `guide_metadata.tsv`, and the `guide_metadata_columns.yaml` that
describes it. See
[`references/outputs.md`](references/outputs.md) for the columns in full and
the three checks worth running before submitting.

## Traps that make a submission wrong but valid-looking

- **`pam` is the genome's bases, not the pattern.** A file with `NGG` in
  every row has not reported a PAM. The format asks for the motif "as it
  appears on the targeted strand".
- **A guide that keeps its 5' G.** If a pipeline strips a leading G before
  aligning, the spacer is 19 nt while the coordinates span 20, and the two no
  longer describe the same thing. Check `len(spacer) == guide_end -
  guide_start` before submitting.
- **One row per `guide_id`.** A guide inside two overlapping elements still
  gets one row: the element sharing the most bases with the protospacer. The
  others are not lost, because every associated gene goes in
  `putative_target_genes`.
- **`targeting` is FALSE for both control kinds.** The pipeline uses those
  cells as the comparison set for differential expression, so mislabelling a
  real cutter as a control contaminates the null.

## Feeding the CRISPR pipeline

The IGVF CRISPR Perturb-seq pipeline (`pinellolab/CRISPR_Pipeline`,
Nextflow) takes a `guide_metadata.tsv` with a **subset** of these columns:
everything except `genomic_element`, `putative_target_genes`, `reporter` and
`imperfect`. This table is a superset, so it can be passed as it is.

`intended_target_name` is what joins the guide table to the pipeline's
`pairs_to_test.csv`, so the two files must name elements the same way.

# IGVF export

`cau igvf-export` writes a guide design in the IGVF per-guide metadata
format. It is the fourth step: [alignment](guide-alignment.md) finds each
guide's site, [annotation](guide-annotation.md) says what gene it falls in,
[element mapping](element-mapping.md) says which element it tests, and this
assembles the result into the table IGVF asks for.

In Claude Code the `cau-igvf-export` skill teaches Claude to run it and read
the output ([Agents and skills](agents-and-skills.md)).

## Quick start

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

`--genes` is optional and takes a table of measured element-to-gene pairs:
the first column names the element, and the first column of ENSEMBL ids names
the gene. Without it, a promoter's gene is inferred from the nearest TSS.

## The format

The table follows
`Format_for_Per-Guide_Metadata.pdf`, the IGVF document that files with
`content_type` "guide RNA sequences" must abide by. Its conventions are the
ones the rest of this package already uses — 0-based, half-open, ACTG — and
it pins ENSEMBL ids to **GENCODE release 43**, so build the annotation store
from v43.

## What a guide targets is decided from the genome

A design's category for an element is a hypothesis; the annotation is
evidence. In order:

1. a validated enhancer stays `enhancer`;
2. within 2 kb of a protein-coding TSS, `promoter`;
3. overlapping an exon, `exon`;
4. inside a gene, `intron`;
5. otherwise, `distal element`.

The 2 kb is the format's own definition of a distal element ("outside of a
gene body, >2 kb from any annotated TSS").

The exception is a validated enhancer: that an element was shown to regulate
a particular gene is experimental knowledge no annotation can supply, so that
label is trusted and its gene comes from the measured pairs.

This matters in practice. In a real 43,778-guide screen the design labelled
344 elements `promoter`; 300 of those have a measured gene and the genome
agreed with every one, but of the other 44, **39 contained no transcription
start site of any biotype** — checked against all 62,757 genes of v43 — and
one sat 109 kb from the nearest. Derived rather than believed, those 44 came
out as 25 `exon`, 18 `promoter` and 1 `intron`.

## `type` has to agree

`positive control` means a guide expected to produce a known measurable
effect. A row claiming that while its `genomic_element` is `exon`
contradicts itself, so a guide is a positive control only when its element is
a validated enhancer, or the design calls it a promoter **and** the genome
agrees.

## One row per guide

The format requires `guide_id` to be unique, so a guide inside two
overlapping elements still gets one row: the element sharing the most bases
with the protospacer, ties broken by element id. Nothing is lost, because
`putative_target_genes` keeps every gene of every element the guide hit —
which is also how a bidirectional promoter keeps both of its genes.

## Outputs

--8<-- "src/crispr_analysis_utils/skills/cau-igvf-export/references/outputs.md:tables"

## Feeding the CRISPR pipeline

The IGVF CRISPR Perturb-seq pipeline takes a `guide_metadata.tsv` with a
subset of these columns — everything except `genomic_element`,
`putative_target_genes`, `reporter` and `imperfect` — so this table can be
passed to it as it is.

`intended_target_name` is what joins the guide table to the pipeline's
pairs-to-test file, so the two have to name elements the same way.

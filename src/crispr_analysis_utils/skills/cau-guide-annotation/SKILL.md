---
name: cau-guide-annotation
description: >-
  Annotate aligned CRISPR guides against genes with crispr-analysis-utils (the
  `cau guide-annotation` command, or `crispr_analysis_utils.guide_annotation`).
  Use after `cau guide-alignment` when a user wants to know what each guide
  hits: exon, intron or intergenic, whether it falls in a promoter, which gene
  and transcript, how far it is from the nearest protein-coding TSS, whether it
  reads sense or antisense, where in the CDS it cuts, its GC content and its
  poly-T run. Also for building the GENCODE annotation store, for choosing a
  promoter file, and for reading guides_gene_annotation.tsv or columns.yaml.
---

# Gene annotation (`cau guide-annotation`)

Alignment answers where in the genome a guide is. This answers where in a
gene. It reads the alignment table and a store built once from GENCODE.

```bash
cau guide-annotation build-db gencode.v43.annotation.gtf.gz \
    --out gencode.v43.protein_coding.parquet

cau guide-annotation run \
    --alignment results/alignments/guides_alignment_metadata.tsv \
    --store gencode.v43.protein_coding.parquet \
    --promoters gencode.v43.protein_coding.TSS1bp.bed6 \
    --outdir results/
```

Building the store reads the whole GTF and takes about a minute. A run reads
the store in about a second, and a 43,000-site library annotates in another.

## Ask before running

- **Which GENCODE?** The default is protein-coding v43. Ask which release the
  rest of the analysis uses: a guide's gene and exon call changes between
  releases, so mixing them silently disagrees.
- **A full annotation GTF, not a genes-only one.** The store needs
  `transcript` and `exon` rows. A genes-only GTF gives gene and intergenic
  calls and nothing else, and `build-db` refuses it.
- **Which promoter file?** Ask for the lab's TSS BED. Without `--promoters`,
  each gene's reference transcript start is used instead, which is a
  defensible default but not the same file the rest of the lab uses.
- **Which gene types?** `--gene-types` defaults to `protein_coding`, so a
  guide in a lncRNA comes back `intergenic`. Say so rather than let it
  surprise someone.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--gene-types` | `protein_coding` | Which genes go in the store |
| `--promoters` | the store's TSS | A TSS BED, one site per gene |
| `--promoter-window` | `1000` | How far either side of a TSS is its promoter |

## Reading the result

`--outdir` gets `annotations/`, `run_metadata/` and `logs/`. The results are
one table, `guides_gene_annotation.tsv`, and the `columns.yaml` that
describes it, generated from the same definitions. See
[`references/outputs.md`](references/outputs.md) for the columns in full.

Four things trip people up:

- **One row per site per gene.** A site in two overlapping genes gets two
  rows, so the row count is not the site count. Only the sites with a PAM are
  annotated, so a multi-mapper appears once per cuttable site.
- **`location` and `promoter` are separate.** A cut 300 bp into exon 1 is
  both exonic and in a promoter; one column would have to hide one of them.
- **The table keys back to the alignment table** on `guide_id` plus the four
  site columns. The alias, accession and spacer are not repeated. A guide
  with two sites needs all five columns to join unambiguously.
- **`intergenic` means "not in a protein-coding gene"**, not "in no gene at
  all", whenever `--gene-types` is left at its default.

## The reference transcript

A gene has many transcripts, so one is chosen for `enst`, `exon_number` and
the CDS position: **MANE Select, else Ensembl canonical, else the longest
CDS**. `location` is called against it, and `exonic_any_transcript` says
whether any other transcript has an exon there, so a row called `intron` that
is someone else's exon is still visible.

## What is not here

**Off-target annotation.** Every site in the table is an exact match the
alignment found. Where else a guide might cut is a different question.

**Non-coding genes, enhancers and blacklists.** The store holds what
`--gene-types` kept.

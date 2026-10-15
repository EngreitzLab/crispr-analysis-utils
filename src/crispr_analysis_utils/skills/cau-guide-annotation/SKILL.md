---
name: cau-guide-annotation
description: >-
  Annotate aligned CRISPR guides against genes with crispr-analysis-utils (the
  `cau guide-annotation` command, or `crispr_analysis_utils.guide_annotation`).
  Use after `cau guide-alignment` when a user wants to know what each guide
  hits: genic or intergenic, exon or intron, whether it falls in a promoter
  and whose, which gene
  and transcript, how far it is from the nearest protein-coding TSS, whether it
  reads sense or antisense, where in the CDS it cuts, its GC content and its
  poly-T run. Also for building the GENCODE annotation store, for choosing a
  promoter file, and for reading guide_gene_annotation.tsv or the columns file beside it.
---

# Gene annotation (`cau guide-annotation`)

Alignment answers where in the genome a guide is. This answers where in a
gene. It reads the alignment table and a store built once from GENCODE.

```bash
cau guide-annotation build-db gencode.v43.annotation.gtf.gz \
    --out gencode.v43.protein_coding.parquet

cau guide-annotation run \
    --alignment results/alignments/guide_genome_annotation.tsv \
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
- **Which promoter file?** Ask for the lab's TSS BED. A file with one site
  per gene loses alternative promoters. Without `--promoters`,
  each gene's reference transcript start is used instead, which is a
  defensible default but not the same file the rest of the lab uses.
- **Which gene types?** `--gene-types` defaults to `protein_coding`, so a
  guide in a lncRNA comes back `intergenic`. Say so rather than let it
  surprise someone.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--gene-types` | `protein_coding` | Which genes go in the store |
| `--promoters` | the store's TSS | A TSS BED, one row per site. A seventh column is the gene id, an eighth the transcript id, which is what fills `promoter_enst` |
| `--promoter-window` | `1000` | How far either side of a TSS is its promoter |
| `--distance-from` | `midpoint` | Where distances are measured from: the protospacer midpoint, or `target_site` 3 bases 5' of the PAM |

## Reading the result

`--outdir` gets `annotations/`, `run_metadata/` and `logs/`. The results are
one table, `guide_gene_annotation.tsv`, and the
`guide_gene_annotation_columns.yaml` that
describes it, generated from the same definitions. See
[`references/outputs.md`](references/outputs.md) for the columns in full.

Seven things trip people up:

- **One row per site per gene.** A site in two overlapping genes gets two
  rows, so the row count is not the site count. Only the sites with a PAM are
  annotated, so a multi-mapper appears once per cuttable site.
- **`location`, `feature` and `promoter` are three axes, not one.**
  `location` is `genic` or `intergenic`, `feature` is `exon` or `intron`
  inside the reference transcript, and `promoter` is independent of both: a
  cut 300 bp into exon 1 is genic, exonic and in a promoter at once.
- **The table keys back to the alignment table** on `guide_id` plus the four
  site columns. The alias, accession and spacer are not repeated. A guide
  with two sites needs all five columns to join unambiguously.
- **`intergenic` means "not in a protein-coding gene"**, not "in no gene at
  all", whenever `--gene-types` is left at its default. A protospacer inside
  a lncRNA is `intergenic` with `other_gene_types` naming the type and
  `other_genes` its id, so the overlap is visible without deciding where the
  guide sits.
- **Distances are measured from the protospacer midpoint**, not from
  `target_site`. A dCas9 occupies the whole protospacer rather than acting at
  one base, so its centre is the honest reference; `--distance-from
  target_site` restores the nuclease convention. `cds_position` always uses
  `target_site`, because it asks where along the coding sequence Cas9 acts.
- **`distance_to_gene_start` tells a canonical promoter from an alternative
  one.** Both are flagged `promoter` and both name their transcripts, so the
  flag alone cannot separate them; the distance from the gene's annotated 5'
  end can.
- **A gene has as many promoters as transcript starts.** `promoter_enst`,
  `promoter_ensg` and `promoter_gene` name every transcript whose window the
  protospacer overlaps, semicolon-separated. `nearest_tss_*` stays singular,
  because the nearest is one.

## The reference transcript

A gene has many transcripts, so one is chosen for `enst`, `exon_number` and
the CDS position: **MANE Select, else Ensembl canonical, else the longest
CDS**. `feature` is called against it, and `exonic_any_transcript` says
whether any other transcript has an exon there, so a row called `intron` that
is someone else's exon is still visible.

## What is not here

**Off-target annotation.** Every site in the table is an exact match the
alignment found. Where else a guide might cut is a different question.

**Non-coding genes, enhancers and blacklists.** The store holds what
`--gene-types` kept.

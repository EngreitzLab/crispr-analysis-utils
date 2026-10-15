# Gene annotation

`cau guide-annotation` places every aligned guide in a gene. Alignment answers
where in the genome a guide is; this answers where in a gene: exon, intron or
intergenic, whether it falls in a promoter, which gene and transcript, how far
it is from the nearest protein-coding TSS, and a few properties of the
protospacer itself.

It reads the table [`cau guide-alignment`](guide-alignment.md) wrote, and a
store built once from GENCODE.

In Claude Code the `cau-guide-annotation` skill teaches Claude to run it and
read the output ([Agents and skills](agents-and-skills.md)).

## Quick start

Build the store once:

```bash
cau guide-annotation build-db gencode.v43.annotation.gtf.gz \
    --out gencode.v43.protein_coding.parquet
```

Then annotate:

```bash
cau guide-annotation run \
    --alignment results/alignments/guide_genome_annotation.tsv \
    --store gencode.v43.protein_coding.parquet \
    --promoters gencode.v43.protein_coding.TSS1bp.bed6 \
    --outdir results/
```

Nothing is written outside `--outdir`, and nothing beside your inputs.

## The store

Reading a full GENCODE takes about a minute, mostly because every `tag` on
every row has to be kept: a GENCODE transcript carries several, and a reader
that keeps only the last finds 4,290 of the 19,077 MANE Select transcripts.
The store is that work done once. It holds the `gene`, `transcript`, `exon`,
`CDS` and `UTR` rows of the gene types `--gene-types` kept, is about 13 MB for
protein-coding v43, and loads in about a second.

The GTF must be a full annotation. A genes-only GTF has no `exon` rows, so it
cannot say exon, intron or transcript, and `build-db` refuses it rather than
call everything intergenic.

## Promoters

`--promoters` takes a BED of transcription start sites, one per gene, which is
usually the file the rest of the lab already uses. The seventh column, if
present, is read as the gene id.

Without it, each gene's reference transcript start is used. That is a
reasonable default and not the same thing: a TSS file built from a different
release, or from a different transcript per gene, will disagree by a few bases
or a few kilobases.

`--promoter-window` (default 1000) is how far either side of a TSS counts as
its promoter.

## The reference transcript

A gene has many transcripts, so one is chosen for `enst`, `exon_number` and
the position along the CDS: **MANE Select, else Ensembl canonical, else the
longest CDS**, and the order is deterministic down to the transcript id, so a
rebuild chooses the same one.

`location` is called against that transcript. `exonic_any_transcript` records
whether any transcript of the gene has an exon there, so a row called `intron`
that is another transcript's exon is still visible.

## Location and promoter are separate columns

A cut 300 bp into exon 1 is both exonic and inside a promoter. A single
column would have to pick one and hide the other, so `location` is
`exon`/`intron`/`intergenic` and `promoter` is its own TRUE/FALSE.

With `--gene-types` at its default, `intergenic` means "not in a
protein-coding gene", not "in no gene at all".

## What is measured from where

Overlaps are decided on the **protospacer**, so a guide reaching into an exon
counts as exonic. Distances and the CDS position are measured from the **cut
site**. Every interval test ignores strand: a guide on one strand targets a
gene on either.

## Outputs

--8<-- "src/crispr_analysis_utils/skills/cau-guide-annotation/references/outputs.md:tables"

## Cost

A 43,000-site library against protein-coding v43 annotates in about a second,
once the store is built.

# MuData builder

`cau mudata-builder` assembles a sequenced channel's two count matrices into
one object. The four steps before it describe a library on paper;
[IGVF export](igvf-export.md) ends at `guide_metadata.tsv`. This is where the
design meets the data.

In Claude Code the `cau-mudata-builder` skill teaches Claude to run it and
read the output ([Agents and skills](agents-and-skills.md)).

## Quick start

```bash
cau mudata-builder run \
    --rna counts_unfiltered/adata.h5ad \
    --guide counts_unfiltered_modified/adata.h5ad \
    --guide-metadata results/igvf/guide_metadata.tsv \
    --store gencode.v43.annotation.parquet \
    --measurement-set IGVFDS8842PXCC \
    --gtf-accession IGVFFI9561BASO \
    --outdir results/
```

Nothing is written outside `--outdir`.

## What it assembles

| Where | What |
| --- | --- |
| `mod['gene'].X` | the gene counts, untouched, with the kb layers beside them |
| `mod['gene'].var` | `symbol`, `gene_type`, `gene_chr`, `gene_start`, `gene_end` |
| `mod['guide'].X` | the guide counts, untouched |
| `mod['guide'].var` | every column of `guide_metadata.tsv` |
| `obs` | `barcode`, `measurement_set`, and any `--covariates` given |
| `obsm['gene']`, `obsm['guide']` | which modality holds each cell |
| `uns['reference']` | the GTF and genome: name, md5 and accession |

The layout follows the IGVF CRISPR pipeline's own object, so the two can be
compared field by field.

## It assembles, and nothing else

No per-cell or per-feature statistics, no QC, no filtering, no doublet calls,
no guide assignment, no pairs to test. Every one of those is a judgement, and
a step that makes judgements while assembling hides them. They belong to the
steps after this one, where they can be seen and argued with.

## Cells are the union

A barcode seen in one modality and not the other is kept, and `obsm` records
which modality holds it. On one real channel that is 1,022,982 barcodes
against 207,274 seen by both: an intersection would discard four fifths of
the object before anyone had looked at it. Filtering comes later, once
somebody has decided what to filter on.

## Every gene type is annotated

The counts come from a transcriptome built on the whole GTF, so two thirds of
the features are not protein-coding:

```text
protein_coding           20,066
lncRNA                   18,874
processed_pseudogene     10,150
unprocessed_pseudogene    2,608
misc_RNA                  2,221
snRNA                     1,910
...                      62,757 in total
```

Annotating only the coding genes would leave 42,691 features nameless, so the
store is read with every gene type it holds, not only the ones guides are
annotated against.

## Two traps

**PAR genes appear twice.** 47 features carry a `_PAR_Y` suffix and have a
chrX twin with identical sequence, so kallisto cannot tell their reads apart
and counts them as two features. Both rows are carried through untouched.
Anything downstream that sums or de-duplicates them is making a choice; this
step does not.

**The ids carry GENCODE versions.** kb indexes genes as `ENSG00000290825.1`
while the store drops the version, and the two are joined on the store's own
rule: strip a trailing `.N`, leave `_PAR_Y` alone. A rule that split on the
first dot instead would hand the chrY copy the chrX gene's coordinates, and
nothing would complain.

## One GENCODE, checked

Genes and guides must be annotated against the same release, or the two
halves of the object describe different genomes. The run reads the record
beside `guide_metadata.tsv` and refuses a store other than the one that wrote
it. Where no record is found it says so rather than assuming.

`uns['reference']` then carries the GTF and the genome by **name, checksum
and accession**, so the object answers "which GENCODE, which genome" on its
own. A path is a claim; a checksum is evidence. Accessions come from
`--gtf-accession` and `--genome-accession` -- the CRISPR pipeline's inputs
carry the IGVF portal ids -- and fall back to the file's name, so the field
is never empty and never a guess. Only the name is stored, never the path:
the object is shared, and a path carries the filesystem it was built on.

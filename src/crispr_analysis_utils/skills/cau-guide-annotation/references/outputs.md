# The outputs of `cau guide-annotation`

<!-- --8<-- [start:tables] -->
`--outdir` gets three folders and nothing loose:

```text
annotations/     guide_gene_annotation.tsv
                 guide_gene_annotation_columns.yaml
run_metadata/    guide_annotation_<UTC timestamp>.json
logs/            guide_annotation.log
```

`guide_gene_annotation_columns.yaml` is written from the same definitions
the table is, so it is
always an accurate data dictionary of the folder it sits in.

## `guide_gene_annotation.tsv`

**One row per mapped site per overlapping gene.** A site in two genes gets
two rows; a site in none gets one row, `intergenic`, with the gene columns
empty. Only sites with a PAM are annotated.

| Column | Holds |
| --- | --- |
| `guide_id` | The aligned sequence. With the four site columns, the key back to `guide_genome_annotation.tsv`. |
| `chr`, `start`, `end`, `strand` | The protospacer, 0-based and half-open. |
| `ensg`, `gene_symbol` | The gene this row is about, without a version. |
| `enst` | Its reference transcript, when the protospacer hits one of that transcript's exons. |
| `location` | `genic` or `intergenic`, from whether a gene of an annotated type overlaps. |
| `feature` | `exon` or `intron` inside the reference transcript; empty when intergenic. |
| `exonic_any_transcript` | TRUE when any transcript of the gene has an exon here. |
| `region` | `CDS` or `UTR` for an exonic row. |
| `exon_number` | Which exon of the reference transcript, from its 5' end. |
| `cds_position`, `cds_fraction` | How far along the CDS the cut falls, in bases and as a fraction. |
| `distance_to_gene_start` | Signed bases from the anchor to this gene's annotated 5' end, in the gene's orientation. Near zero is the canonical promoter. |
| `sense` | Whether the protospacer reads on the gene's own strand. |
| `other_gene_types`, `other_genes` | Semicolon-separated types and ids of overlapping genes of other types, such as lncRNA. Empty when none. |
| `promoter` | TRUE within the promoter window of any protein-coding TSS. |
| `promoter_enst`, `promoter_ensg`, `promoter_gene` | Semicolon-separated transcripts whose promoter window the protospacer overlaps, and their genes. Empty when the TSS file carries no transcript id. |
| `distance_to_tss` | Signed bases from the anchor to the nearest TSS. The anchor is the protospacer midpoint unless `--distance-from target_site`. |
| `nearest_tss_gene`, `nearest_tss_ensg` | Whose TSS that is. |
| `gc_percent`, `max_polyt_run` | The protospacer's GC, and its longest T run. |

### Joining it back

The alias, accession and spacer are not repeated; they are in the alignment
table. Join on all five key columns, because a guide with two sites has two
rows in each table:

```python
annotation.merge(alignment, on=["guide_id", "chr", "start", "end", "strand"])
```

### What is measured from where

Overlaps are decided on the **protospacer**, so a guide reaching into an exon
counts as exonic. Distances and the CDS position are measured from the **cut
site**, which is where the break falls. Every interval test ignores strand: a
guide on one strand targets a gene on either.

## `run_metadata/guide_annotation_<timestamp>.json`

The parameters, the input checksums, the store's own provenance, the counts
and the timings. Nothing per-guide lives here.
<!-- --8<-- [end:tables] -->

## Reading a result

- **A `no_pam` or unmapped guide is simply absent.** Nothing can be said
  about where a guide falls if it does not map, so those rows are dropped
  rather than carried with empty columns.
- **`intergenic` is the right answer for an enhancer-targeting library.** A
  tiling screen puts most of its guides between genes on purpose.
- **`max_polyt_run` of 4 or more is a sequence warning, not an annotation.**
  A run of four T terminates Pol III, so the guide expresses poorly from a U6
  promoter however good its target is.

# The outputs of `cau guide-annotation`

<!-- --8<-- [start:tables] -->
`--outdir` gets three folders and nothing loose:

```text
annotations/     guides_gene_annotation.tsv  columns.yaml
run_metadata/    guide_annotation_<UTC timestamp>.json
logs/            guide_annotation.log
```

`columns.yaml` is written from the same definitions the table is, so it is
always an accurate data dictionary of the folder it sits in.

## `guides_gene_annotation.tsv`

**One row per mapped site per overlapping gene.** A site in two genes gets
two rows; a site in none gets one row, `intergenic`, with the gene columns
empty. Only sites with a PAM are annotated.

| Column | Holds |
| --- | --- |
| `guide_id` | The aligned sequence. With the four site columns, the key back to `guides_alignment_metadata.tsv`. |
| `chr`, `start`, `end`, `strand` | The protospacer, 0-based and half-open. |
| `ensg`, `gene_symbol` | The gene this row is about, without a version. |
| `enst` | Its reference transcript, when the protospacer hits one of that transcript's exons. |
| `location` | `exon`, `intron` or `intergenic`. |
| `exonic_any_transcript` | TRUE when any transcript of the gene has an exon here. |
| `region` | `CDS` or `UTR` for an exonic row. |
| `exon_number` | Which exon of the reference transcript, from its 5' end. |
| `cds_position`, `cds_fraction` | How far along the CDS the cut falls, in bases and as a fraction. |
| `sense` | Whether the protospacer reads on the gene's own strand. |
| `promoter` | TRUE within the promoter window of any protein-coding TSS. |
| `distance_to_tss` | Signed bases from the cut site to the nearest TSS. |
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

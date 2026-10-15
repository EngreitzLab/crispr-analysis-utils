# The output of `cau igvf-export`

<!-- --8<-- [start:tables] -->
`--outdir` gets three folders and nothing loose:

```text
igvf/            guide_metadata.tsv  guide_metadata_columns.yaml
run_metadata/    igvf_export_<UTC timestamp>.json
logs/            igvf_export.log
```

## `guide_metadata.tsv`

**One row per guide**, which the format requires: a `guide_id` must name one
guide and one spacer.

| Column | Holds |
| --- | --- |
| `guide_id` | The guide's name, from the input table's own name for it. |
| `spacer` | The targeting sequence, 5' to 3', in DNA bases. |
| `oligo` | The sequence ordered, without the PAM. |
| `n_oligos` | Entries sharing this spacer; above one the counts are a sum. |
| `targeting` | TRUE only when the guide targets a relevant locus. |
| `type` | `targeting`, `safe-targeting`, `non-targeting` or `positive control`. |
| `guide_chr`, `guide_start`, `guide_end` | The protospacer, 0-based and half-open, **without** the PAM. |
| `strand` | The strand targeted. |
| `pam` | The PAM as the genome holds it, 5' to 3' on that strand. |
| `genomic_element` | `promoter`, `enhancer`, `distal element`, `exon` or `intron`. |
| `intended_target_name` | Names the target; what it holds depends on the element. |
| `intended_target_chr/start/end` | The element's coordinates. |
| `putative_target_genes` | Every associated gene, comma-separated. |
| `reporter`, `imperfect` | Unused here; the format allows them. |

### What `intended_target_name` holds

| `genomic_element` | holds |
| --- | --- |
| `promoter` | the regulated gene's ENSEMBL id |
| `exon` | the ENSEMBL exon id |
| `intron` | the intron's coordinates |
| `enhancer`, `distal element` | the element's coordinates |

### A guide in two elements

Elements overlap, so a guide can sit in several, and the format still allows
one row. The element sharing the most bases with the protospacer wins, ties
broken by element id so a rerun gives the same answer. Nothing is lost:
`putative_target_genes` carries every gene of every element the guide hit.

## `run_metadata/igvf_export_<timestamp>.json`

The parameters, the input checksums, the store's provenance and the counts:
guides by `type`, elements by `genomic_element`, and how many guides are
missing a PAM, which should be zero for anything that maps.
<!-- --8<-- [end:tables] -->

## Before submitting

Three checks that catch a file which validates but is wrong:

- **`len(spacer) == guide_end - guide_start`.** A pipeline that strips a 5' G
  before aligning leaves a 19-nt spacer against 20-nt coordinates, and the
  two no longer describe the same guide.
- **`pam` is not `NGG`.** The pattern is not a PAM. Every row should hold the
  bases the genome has there.
- **`genomic_element` and `type` agree.** A `positive control` whose element
  is an `exon` contradicts itself.

# The outputs of `cau guide-alignment run`

<!-- --8<-- [start:tables] -->
`--outdir` gets three folders and nothing loose, so a run's output never sits
beside its inputs:

```text
alignments/      guides_alignment_metadata.tsv  guides_mapping.bed
                 guides_alignment_metadata_columns.yaml  qc.tsv
run_metadata/    guide_alignment_<UTC timestamp>.json
logs/            guide_alignment.log
```

`guides_alignment_metadata_columns.yaml` is written from the same
definitions the tables are, so it is
always an accurate data dictionary of the folder it sits in. Read it first.

## `guides_alignment_metadata.tsv`

**One row per library entry per site.** A guide matching three places gets
three rows; one matching nowhere gets a single row whose site columns are
empty. The guide's own columns repeat on each of its rows.

| Column | Holds |
| --- | --- |
| `guide_id` | The sequence that was aligned, ending in the PAM pattern. Identifies the guide. |
| `alias` | The name the input table gives it. |
| `spacer` | The protospacer: `guide_id` without any base 5' of it or the PAM. |
| `site_class` | `unique`, `multi` or `no_site`: how many places this guide can cut. |
| `n_matches` | How many places the protospacer occurs in the genome. |
| `n_cut_sites` | How many of those have a PAM. |
| `orientation` | `forward`, `reverse_complement`, `no_pam`, `no_match` or `not_checked`. |
| `orientation_evidence` | For `reverse_complement` only: `distinct` or `ambiguous`. |
| `chr`, `start`, `end`, `strand` | This row's match. 0-based, half-open, the protospacer alone. |
| `cut_site` | The boundary the cut falls on: between `cut_site - 1` and `cut_site`. |
| `pam_matches` | `TRUE` when the genomic PAM fits the pattern. |
| `mapped_sequence` | What the genome holds here, to compare against `guide_id`. |

### The two sequence columns

They line up base for base, so one diff answers three questions:

```text
guide_id         GAAAAGCCAACATGAATGCAGNGG     what was looked for
mapped_sequence   AAAAGCCAACATGAATGCAGTGG     what the genome holds
                 ^                     ^^^
                 the 5' G did not pair        the real PAM
```

`mapped_sequence` keeps the 5' base only where the genome has that same base,
so its length tells you whether the G is genomic or was added for the U6
promoter.

### `n_matches` is not `site_class`

A guide can occur many times and still be cuttable once. This is the most
common surprise:

```text
n_matches 89   n_cut_sites 1   site_class unique
```

The protospacer occurs 89 times; one of those has a PAM. Cas9 cuts there and
nowhere else. Every one of the 89 is a row in the table, 88 of them with
`pam_matches` `FALSE`.

## `guides_mapping.bed`

BED6, no header, one row per **cuttable** site, so a guide that cuts twice
gets two rows. `start`/`end` are the protospacer, without the PAM. The name
field is the `spacer`, which is the same on every row of one guide, so rows
group by name. Matches with no PAM are not here; they are in the TSV.

## `qc.tsv`

Three columns, `section`, `name`, `value`. The `class` and `orientation`
sections each sum to the number of entries. `scan` holds what the pass saw:
`n_matches`, `n_off_contig_end`, `n_no_pam`, `n_sites`.

## `run_metadata/guide_alignment_<timestamp>.json`

One file per invocation, named for the step and when it ran: the parameters,
the input checksums, the version, the counts and the timings. Nothing
per-guide lives here.
<!-- --8<-- [end:tables] -->

## Reading a result

- **A guide that matches nowhere is usually a fact about the library, not a
  failure.** Non-targeting controls are designed to be absent from the
  genome, and `no_match` is the correct answer for them.
- **`no_pam` means the protospacer is really there but has no PAM beside it.**
  Assemblies differ: a guide designed on hg19 can lose its PAM in hg38.
- **A `reverse_complement` verdict with `ambiguous` evidence is not proof.**
  For a guide with no base 5' of its spacer, a flipped entry and a chance
  `CCN` upstream look identical. Check the locus before acting.

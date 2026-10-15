---
name: cau-guide-alignment
description: >-
  Find where a CRISPR guide library maps in a genome with crispr-analysis-utils
  (the `cau guide-alignment` command, or
  `crispr_analysis_utils.guide_alignment.run`). Use when a user wants each
  guide's genomic site: which guides are unique, multi-mapping or absent;
  their protospacer coordinates and cut sites, for sceptre, element annotation
  or IGVF-style guide metadata; whether a guide table is reverse-complemented;
  or whether non-targeting controls really match nothing. Also for reading
  guide_genome_annotation.tsv, guide_mapping.bed, the columns file or
  qc.tsv,
  and for users of the earlier GEM-based workflow.
---

# Guide alignment (`cau guide-alignment`)

A guide is designed, not sequenced, so its site is an **exact** match of the
protospacer. The tool scans the reference once and answers three questions per
guide: where does the protospacer match, is there a PAM 3' of it, and what
sits 5' of it. There is no aligner and no index.

```bash
cau guide-alignment run \
    --guides guides.tsv --reference hg38.fa --outdir results/
```

A 44,000-guide library against hg38 takes about half a minute.

## Ask before running

The sequences must be the guide as ordered, 5' to 3'. Most tables are, but
check rather than assume:

- **What is in the sequence column?** The protospacer alone, or with a 5' G
  added for the U6 promoter? Pass `--spacer-length` for the protospacer
  length; anything 5' of it is treated as an added base whose mismatch is
  tolerated. If the table holds bare protospacers and the construct adds a G,
  pass `--add-leading-g`.
- **Is the PAM in the column?** It should not be. The tool appends `--pam`
  (default `NGG`) itself. One spacer in sixteen ends in `GG`, so a table that
  silently includes the PAM will look fine and be wrong by three bases.
- **Which contigs?** Without `--contigs` the whole FASTA is read, alt and
  random scaffolds included. For hg38 give the primary contigs.
- **Is the reference uncompressed?** A `.gz` FASTA is refused.

Never repair a sequence yourself. Whether a G was added, and which bases are
the PAM, is the library designer's answer, not something to infer.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--pam` | `NGG` | IUPAC, 3' of the protospacer |
| `--spacer-length` | `20` | The protospacer is the 3'-most this many bases |
| `--add-leading-g` | off | Put a G in front of every sequence |
| `--orientation-check` | `error` | `error`, `warn` or `off` |
| `--fix-orientation` | off | Align the reverse complement of flagged entries |
| `--contigs` | all | Comma-separated contig names |
| `--id-col`, `--spacer-col` | `0`, `1` | Position or name |
| `--sep`, `--header` | tab, off | `tab`, `comma`, or one character |

## Reading the result

`--outdir` gets `alignments/`, `run_metadata/` and `logs/`. Start with
`alignments/guide_genome_annotation_columns.yaml`: it is generated from
the same definitions as the
tables, so it always describes the folder it sits in. See
[`references/outputs.md`](references/outputs.md) for the columns in full.

Three things trip people up:

- **`n_matches` is not `site_class`.** A protospacer can occur 89 times and
  still be `unique`, because only one of those has a PAM. `n_cut_sites` is
  the one that matters for cutting; `n_matches` is the one that matters for
  specificity.
- **`no_match` is the right answer for a non-targeting control.** They are
  designed to be absent from the genome. A library's controls all coming back
  `no_match` is the check passing, not failing.
- **`guide_id` and `mapped_sequence` line up base for base.** Diff them to see
  whether the 5' base paired and what the real PAM is.

## Orientation

Both readings of every guide go into the same scan, so an entry is
`reverse_complement` only when its reverse complement has a PAM-bearing site
and the sequence as given does not. A `CCN` 5' of a match is **not** evidence:
it falls there by chance one time in sixteen.

A table reversed throughout stops the run; a handful of entries warns, because
one guide can look reversed when its PAM simply differs between assemblies.
`orientation_evidence` says whether the test was `distinct` or `ambiguous` —
ambiguous means the guide has nothing 5' of its spacer, so a flip and a chance
`CCN` are indistinguishable. Check the locus before believing it.

## What is not here

**Off-target search.** Matching is exact, so there are no mismatch counts. A
guide's off-targets are a different question from where it comes from.

## Coming from the GEM version

There is no index to build and no `--max-mismatches`, `--alt-pams` or
`--threads`; the files are renamed. See
[`references/migration.md`](references/migration.md).

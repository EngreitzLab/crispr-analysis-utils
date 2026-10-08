---
name: cau-accession
description: >-
  Give SpCas9 guides accession ids encoded from their own sequence, with
  crispr-analysis-utils (`import crispr_analysis_utils as cau`). Use when a
  user wants a lab-independent id for each guide (GU followed by base 32
  characters), wants guides keyed by sequence rather than by lab names, asks
  why one guide has different ids in two tables, wants the guide behind an
  accession, or hits a guide_accession error (no NGG PAM, no leading G, a base
  other than A, C, G or T, an empty or missing sequence).
---

# Guide accessions (`cau.accession`)

An accession encodes the guide itself, so the same guide gets the same
accession anywhere and two guides never share one.

```python
import crispr_analysis_utils as cau

cau.accession.guide_accession("GAAAAGCCAACATGAATGCAGNGG")  # 'GU602A170WJ'
cau.accession.guide_sequence("GU602A170WJ")  # 'GAAAAGCCAACATGAATGCAGNGG'
```

## Ask before accessioning anything

The guide must be written `G...NGG`: a 5' G, the protospacer, then the PAM
spelled `NGG`. Most tables hold neither end. **Ask; never repair a sequence.**

- **No `NGG` at the end?** Ask which bases are the PAM, then append `NGG`.
  Write `NGG` itself, never the genomic bases, so one guide keeps one
  accession whatever the genome holds at the N. Never read a sequence's own
  last three bases as a PAM: one spacer in sixteen ends in `GG`.
- **No G at the start?** Ask whether the guide carries its own G or is cloned
  with one added, then write the whole ordered sequence. A guide with an added
  G and the same guide without it are two RNAs, and get two accessions.
- **SpCas9 only.** Say so rather than pass another nuclease's PAM off as
  `NGG`.

## A library

```python
import pandas as pd

import crispr_analysis_utils as cau

guides = pd.read_csv(
    "guides.tsv", sep="\t", header=None, names=["name", "spacer"], dtype=str
)
# Only once the user has confirmed `spacer` is G + protospacer, with no PAM:
guides["accession"] = cau.accession.guide_accessions(guides["spacer"] + "NGG")
```

`dtype=str` makes an empty cell NaN, which raises an error naming its position
rather than quietly getting an accession.

## The encoding

Drop the PAM; pack the rest at two bits a base (A, C, G, T = 0, 1, 2, 3)
behind a leading 1 bit; write that in Crockford's base 32
(`0123456789ABCDEFGHJKMNPQRSTVWXYZ`, reading I and L as 1 and O as 0) after
`GU`. That is 9 characters for 20 to 22 bases before the PAM, 8 for 19.

Accessions are therefore not opaque: guides sharing their first bases share
their first characters, and anyone can decode one.

## Errors

| Error | Message | Cause |
| --- | --- | --- |
| `TypeError` | `A guide sequence must be a string, not float (nan).` | an empty cell or a number |
| `ValueError` | `A guide sequence must not be empty.` | an empty or blank string |
| `ValueError` | `... does not end with the PAM 'NGG'.` | no PAM, a genomic PAM, or another PAM |
| `ValueError` | `... does not start with G.` | the 5' G is missing |
| `ValueError` | `... holds 'U': a guide is written with A, C, G and T, ...` | RNA letters, IUPAC codes, gaps or inner spaces |
| `ValueError` | `... holds an N before the PAM ...` | an N in the protospacer |

`guide_accessions` prefixes each with `Sequence at position <i>:`, the 0-based
position in the input, not its pandas index label.

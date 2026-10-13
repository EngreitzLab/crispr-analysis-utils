---
name: cau-accession
description: >-
  Give CRISPR guides accession ids encoded from their protospacer, with
  crispr-analysis-utils (`import crispr_analysis_utils as cau`). Use when a
  user wants a lab-independent id for each guide (GU followed by base 32
  characters), wants to know whether two guides from different libraries hit
  the same target, asks why one guide has different ids in two tables, wants
  the protospacer behind an accession, or hits a guide_accession error (a base
  other than A, C, G or T, an empty or missing sequence).
---

# Guide accessions (`cau.accession`)

An accession encodes the **protospacer**, so the same target gets the same id
anywhere, whatever a lab adds at the 5' end and whatever it calls the guide.

```python
import crispr_analysis_utils as cau

cau.accession.guide_accession("AAAAGCCAACATGAATGCAG")  # 'GU102A170WJ'
cau.accession.spacer_of("GU102A170WJ")  # 'AAAAGCCAACATGAATGCAG'
```

Two different protospacers can never share an accession: the packing is
reversible, so `spacer_of` gives the sequence back.

## The protospacer, not the ordered sequence

This is the whole point, and it is what to explain if a user is surprised:

| | Answers |
| --- | --- |
| `accession` | which **target** — the same wherever that protospacer is used |
| `guide_alignment`'s `guide_id` | which **oligo** — the exact sequence ordered |

Two entries differing only 5' of the protospacer are two oligos and one
target, so they share an accession and differ in `guide_id`. In a real
library that happens whenever sub-libraries used different promoter
conventions: one ordered `GCGGG…` and another `GGCGGG…` for the same site.

## Ask before accessioning anything

The input must be the protospacer alone — no PAM, and nothing added 5' of it.
**Nothing in a sequence says where the protospacer ends**, so a sequence with
the PAM still attached is accepted and silently gets a different accession.
One spacer in sixteen ends in `GG`, so you cannot tell by looking.

- **Does the column include the PAM?** Ask. If it does, strip it before
  accessioning, and say how many bases you removed.
- **Does it carry a 5' G added for the U6 promoter?** Ask. If the library was
  designed at 20 nt and the sequences are 21, the first base is that G and is
  not part of the protospacer.
- If the user has already run `cau guide-alignment`, the `spacer` column of
  `guides_alignment_metadata.tsv` is the protospacer, and the `accession`
  column is already there.

Never guess which bases to strip.

## A whole library

```python
import pandas as pd

import crispr_analysis_utils as cau

guides = pd.read_csv(
    "guides.tsv", sep="\t", header=None, names=["name", "spacer"], dtype=str
)
guides["accession"] = cau.accession.guide_accessions(guides["spacer"])
```

`dtype=str` makes an empty cell NaN, which raises an error naming its position
rather than quietly getting an accession.

## How it is built

Pack the protospacer two bits per base (A, C, G, T = 0, 1, 2, 3) behind a
leading 1 bit, then write that number in Crockford's base 32
(`0123456789ABCDEFGHJKMNPQRSTVWXYZ`, reading I and L as 1 and O as 0) after
`GU`. A 20-nt protospacer gives 9 characters after the prefix, a 19-nt one 8.

Accessions are therefore not opaque: protospacers sharing their first bases
share their first characters, and anyone can decode one. That is the price of
never colliding.

## Errors

| Error | Message | Cause |
| --- | --- | --- |
| `TypeError` | `A protospacer must be a string, not float (nan).` | an empty cell or a number |
| `ValueError` | `A protospacer must not be empty.` | an empty or blank string |
| `ValueError` | `... holds 'U': a protospacer is written with A, C, G and T.` | RNA letters, N, IUPAC codes, gaps or inner spaces |

`guide_accessions` prefixes each with `Spacer at position <i>:`, the 0-based
position in the input, not its pandas index label.

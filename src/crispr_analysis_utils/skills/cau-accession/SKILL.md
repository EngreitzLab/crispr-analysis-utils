---
name: cau-accession
description: >-
  Give CRISPR guides stable accession ids computed from their sequence, with
  crispr-analysis-utils (`import crispr_analysis_utils as cau`). Use when a
  user wants a lab-independent id for each guide (GU followed by 12
  hexadecimal digits), wants guides keyed by sequence rather than by lab
  names, asks why one guide has different ids in two tables, or hits a
  guide_accession error (a base other than A, C, G or T, an empty or missing
  sequence, a collision). Also for code that used a create_accession_id.py
  script, whose ids this module does not reproduce.
---

# Guide accessions (`cau.accession`)

`cau.accession.guide_accession(sequence)` returns `GU` followed by 12
uppercase hexadecimal digits, computed from the sequence alone: the same
sequence gets the same accession on any machine, in any library, whatever the
guide is called. `cau.accession.guide_accessions(sequences)` does a whole
library and refuses collisions.

## Before running anything

- **Ask which sequence the table holds.** The accession is computed from the
  sequence that is aligned: 5' to 3', with the 5' G when the guide carries
  one, and without the PAM. A G added at cloning but missing from the table
  changes the accession, as does a PAM left in the table. Ask whether the
  cloning adds a G the table leaves out, and prepend it only with the user's
  agreement. Never strip a G from the table's sequence.
- **Don't reverse-complement.** Guides are written 5' to 3'; a sequence and
  its reverse complement are two guides with two accessions.
- **Do every library at once.** `guide_accessions` checks for collisions only
  among the sequences it is given. To keep a new library's accessions distinct
  from those already in use, pass both sets of sequences in one call.

## Usage

```python
import pandas as pd

import crispr_analysis_utils as cau

guides = pd.read_csv(
    "guides.tsv", sep="\t", header=None, names=["name", "sequence"], dtype=str
)
guides["accession"] = cau.accession.guide_accessions(guides["sequence"])

cau.accession.guide_accession("GAAAAGCCAACATGAATGCAG")  # 'GU5C664EE08629'
```

Read the table with `dtype=str`. An empty cell then becomes NaN, which raises
a `TypeError` naming its position rather than getting an accession.

## What it does

- Strips surrounding whitespace and uppercases. Anything else but A, C, G and
  T is an error: U, N, IUPAC codes, gaps, inner spaces.
- Takes the UUID version 5 of the result in `cau.accession.GUIDE_NAMESPACE`
  (`ee28d900-f51b-5e2b-9ca7-4c5808da45cc`) and keeps its first 12 hexadecimal
  digits, in uppercase, after `GU`. Outside Python:
  `"GU" + uuid5(namespace, sequence).hex[:12].upper()`.
- Repeats of one sequence share an accession. Two library entries with one
  sequence are one guide under two names.

A collision is rare: about 1 in 290,000 for 43,778 different guides, 1 in 560
for a million. If `guide_accessions` reports one, tell the user. Never edit an
accession by hand.

## Errors

| Error | Message | Cause |
| --- | --- | --- |
| `TypeError` | `A guide sequence must be a string, not float (nan).` | an empty cell or a number in the column |
| `ValueError` | `A guide sequence must not be empty.` | an empty or blank string |
| `ValueError` | `Guide sequence 'ACGU' holds 'U': only A, C, G and T are allowed.` | RNA letters, N, IUPAC codes, gaps or inner spaces |
| `ValueError` | `Two different sequences have the accession GU…: … and ….` | a collision |

`guide_accessions` prefixes the first three with `Sequence at position <i>:`,
the 0-based position in the input, not its pandas index label.

## Coming from create_accession_id.py

Lab repositories carry copies of a `create_accession_id.py` script. The copies
use different namespaces, id lengths and prefixes (`GU…` and `HGRMGU…`), and
take the id from a shifting slice of the UUID. This module does not reproduce
their ids: recompute every accession from the sequence with it.

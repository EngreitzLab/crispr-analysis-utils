---
name: cau-normalization
description: >-
  Normalize CRISPR screen count matrices with crispr-analysis-utils
  (`import crispr_analysis_utils as cau`). Use when a user wants counts per
  million (CPM) for guide- or gene-level counts held in a pandas DataFrame or a
  NumPy array, asks how to library-size normalize samples or guides, or hits a
  counts_per_million error (a zero-total margin, NaN or non-numeric counts, an
  invalid pseudocount). Also for code ported from the package's former R
  function counts_per_million(counts, pseudocount, margin).
---

# Normalization (`cau.counts_per_million`)

`counts_per_million` scales a count matrix so that every column (or every row)
sums to one million. It is the package's only normalization so far.

## Before running anything

- **Confirm the orientation.** The function cannot tell guides from samples.
  `axis=0`, the default, normalizes each *column*: the usual layout, with
  guides or genes in rows and samples in columns. `axis=1` normalizes each
  *row*. Ask when the layout is not obvious from the labels.
- **Ask whether a pseudocount is wanted.** The default is `0`. A pseudocount is
  added to every entry *before* scaling, so it changes every value; use one
  only when the user plans a log transform, and say which value you chose.
- **Look for empty samples.** A column (or row) summing to zero raises an
  error. Offer to drop it; with a pseudocount above zero no margin can be zero.

## Usage

```python
import pandas as pd

import crispr_analysis_utils as cau

counts = pd.DataFrame(
    {"sample_a": [100, 300], "sample_b": [50, 50]},
    index=["guide_1", "guide_2"],
)
cpm = cau.counts_per_million(counts)  # each column sums to 1e6
cpm_rows = cau.counts_per_million(counts, axis=1)  # each row sums to 1e6
cpm_pseudo = cau.counts_per_million(counts, pseudocount=0.5)
```

A DataFrame comes back as a float DataFrame with the same index and columns;
any other input (a NumPy array, nested lists) comes back as a NumPy array.

## Errors

Every failure is a `ValueError`:

| Message | Cause |
| --- | --- |
| `axis must be 0 for columns or 1 for rows.` | `axis` is neither 0 nor 1 |
| `pseudocount must be a single non-negative number.` | negative, NaN, infinite, boolean, text or a list |
| `counts must be numeric.` | a boolean or text column |
| `counts must not contain NaN.` | missing values: fill or drop them, with the user's agreement |
| `Cannot normalize a margin with total count zero.` | an all-zero column (`axis=0`) or row (`axis=1`) |

## Coming from the R package

The R function `counts_per_million(counts, pseudocount = 0, margin = 2)` was
removed when the package became Python-only. R's `margin = 2` (columns) is
`axis=0`, and `margin = 1` (rows) is `axis=1`. R returned a matrix; Python
returns a DataFrame for DataFrame input. Python also rejects an infinite
pseudocount, which R accepted and turned into an all-NaN result.

"""Normalization helpers for CRISPR count matrices."""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd


def counts_per_million(
    counts: pd.DataFrame | np.ndarray,
    pseudocount: float = 0.0,
    axis: Literal[0, 1] = 0,
) -> pd.DataFrame | np.ndarray:
    """Convert raw counts to counts per million.

    `axis=0` normalizes each column independently, which is the usual layout
    for CRISPR screens with guides/features in rows and samples in columns.
    `axis=1` normalizes each row independently.

    Parameters
    ----------
    counts
        Numeric count matrix as a pandas DataFrame or NumPy array.
    pseudocount
        Non-negative value added to every entry before normalization.
    axis
        Dimension to sum over before scaling. Use `0` for columns and `1` for
        rows.

    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        Counts scaled so each selected margin sums to one million.

    Examples
    --------
    >>> import pandas as pd
    >>> counts = pd.DataFrame({"sample_a": [100, 300], "sample_b": [50, 50]})
    >>> counts_per_million(counts)
       sample_a  sample_b
    0  250000.0  500000.0
    1  750000.0  500000.0
    """
    if axis not in (0, 1):
        raise ValueError("axis must be 0 for columns or 1 for rows.")
    if pseudocount < 0:
        raise ValueError("pseudocount must be non-negative.")

    if isinstance(counts, pd.DataFrame):
        values = counts.astype(float) + pseudocount
        totals = values.sum(axis=axis)
        if (totals == 0).any():
            raise ValueError("Cannot normalize a margin with total count zero.")
        return values.divide(totals, axis=1 if axis == 0 else 0) * 1_000_000

    values = np.asarray(counts, dtype=float) + pseudocount
    totals = values.sum(axis=axis, keepdims=True)
    if np.any(totals == 0):
        raise ValueError("Cannot normalize a margin with total count zero.")
    return values / totals * 1_000_000

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

    Raises
    ------
    ValueError
        If `axis` is not 0 or 1, `pseudocount` is not a single finite
        non-negative number, `counts` is not numeric (booleans and strings are
        rejected) or contains NaN, or a margin sums to zero.

    Notes
    -----
    This replaces the package's former R function
    `counts_per_million(counts, pseudocount, margin)`. R's `margin = 2`
    (columns) is `axis=0` here, and `margin = 1` (rows) is `axis=1`.

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
    if not _is_valid_pseudocount(pseudocount):
        raise ValueError("pseudocount must be a single non-negative number.")

    if isinstance(counts, pd.DataFrame):
        if not all(_is_numeric(dtype) for dtype in counts.dtypes):
            raise ValueError("counts must be numeric.")
        values = counts.astype(float) + pseudocount
        if values.isna().to_numpy().any():
            raise ValueError("counts must not contain NaN.")
        totals = values.sum(axis=axis)
        if (totals == 0).any():
            raise ValueError("Cannot normalize a margin with total count zero.")
        return values.divide(totals, axis=1 if axis == 0 else 0) * 1_000_000

    array = np.asarray(counts)
    if not _is_numeric(array.dtype):
        raise ValueError("counts must be numeric.")
    values = array.astype(float) + pseudocount
    if np.isnan(values).any():
        raise ValueError("counts must not contain NaN.")
    totals = values.sum(axis=axis, keepdims=True)
    if np.any(totals == 0):
        raise ValueError("Cannot normalize a margin with total count zero.")
    return values / totals * 1_000_000


def _is_valid_pseudocount(value: object) -> bool:
    """Whether `value` is a single finite non-negative real number."""
    if isinstance(value, (bool, np.bool_)):
        return False
    if not isinstance(value, (int, float, np.integer, np.floating)):
        return False
    return bool(np.isfinite(value)) and value >= 0


def _is_numeric(dtype: object) -> bool:
    """Whether `dtype` holds real numbers; booleans and complex do not count."""
    return (
        pd.api.types.is_numeric_dtype(dtype)
        and not pd.api.types.is_bool_dtype(dtype)
        and not pd.api.types.is_complex_dtype(dtype)
    )

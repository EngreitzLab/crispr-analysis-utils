import numpy as np
import pandas as pd
import pytest

from crispr_analysis_utils import counts_per_million


def test_counts_per_million_normalizes_columns():
    counts = pd.DataFrame(
        {"sample_a": [100, 300], "sample_b": [50, 50]},
        index=["guide_1", "guide_2"],
    )

    result = counts_per_million(counts)

    assert result.loc["guide_1", "sample_a"] == pytest.approx(250_000)
    assert result.loc["guide_2", "sample_a"] == pytest.approx(750_000)
    assert result["sample_a"].sum() == pytest.approx(1_000_000)
    assert result["sample_b"].sum() == pytest.approx(1_000_000)


def test_counts_per_million_normalizes_rows():
    counts = pd.DataFrame({"sample_a": [100, 300], "sample_b": [300, 100]})

    result = counts_per_million(counts, axis=1)

    assert result.iloc[0].tolist() == pytest.approx([250_000, 750_000])
    assert result.sum(axis=1).tolist() == pytest.approx([1_000_000, 1_000_000])


def test_counts_per_million_accepts_ndarray():
    counts = np.array([[100, 50], [300, 50]])

    result = counts_per_million(counts)

    assert isinstance(result, np.ndarray)
    assert result[:, 0].tolist() == pytest.approx([250_000, 750_000])
    assert result.sum(axis=0).tolist() == pytest.approx([1_000_000, 1_000_000])


def test_counts_per_million_adds_pseudocount_before_scaling():
    counts = pd.DataFrame({"sample_a": [0, 2]})

    result = counts_per_million(counts, pseudocount=1)

    assert result["sample_a"].tolist() == pytest.approx([250_000, 750_000])


def test_counts_per_million_rejects_zero_totals():
    counts = pd.DataFrame({"sample_a": [0, 0]})

    with pytest.raises(ValueError, match="total count zero"):
        counts_per_million(counts)


@pytest.mark.parametrize("axis", [-1, 2, "columns"])
def test_counts_per_million_rejects_bad_axis(axis):
    with pytest.raises(ValueError, match="axis must be 0"):
        counts_per_million(pd.DataFrame({"a": [1, 2]}), axis=axis)


@pytest.mark.parametrize(
    "pseudocount", [-1, float("nan"), float("inf"), [1, 2], True, "1", None]
)
def test_counts_per_million_rejects_invalid_pseudocount(pseudocount):
    with pytest.raises(ValueError, match="pseudocount must be a single"):
        counts_per_million(pd.DataFrame({"a": [1, 2]}), pseudocount=pseudocount)


@pytest.mark.parametrize(
    "counts",
    [
        pd.DataFrame({"a": [1.0, np.nan]}),
        np.array([[1.0, np.nan], [2.0, 3.0]]),
    ],
)
def test_counts_per_million_rejects_nan_counts(counts):
    with pytest.raises(ValueError, match="must not contain NaN"):
        counts_per_million(counts)


@pytest.mark.parametrize(
    "counts",
    [
        pd.DataFrame({"a": [True, False]}),
        pd.DataFrame({"a": [1, 2], "b": ["3", "4"]}),
        np.array([[True, False], [False, True]]),
        np.array([["1", "2"], ["3", "4"]]),
    ],
)
def test_counts_per_million_rejects_non_numeric_counts(counts):
    with pytest.raises(ValueError, match="counts must be numeric"):
        counts_per_million(counts)

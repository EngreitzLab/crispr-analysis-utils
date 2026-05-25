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


def test_counts_per_million_rejects_zero_totals():
    counts = pd.DataFrame({"sample_a": [0, 0]})

    with pytest.raises(ValueError, match="total count zero"):
        counts_per_million(counts)

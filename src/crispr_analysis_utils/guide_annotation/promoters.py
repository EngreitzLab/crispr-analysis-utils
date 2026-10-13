"""Transcription start sites, from a BED file or from the store itself."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

_BED6 = ("Chromosome", "Start", "End", "gene_symbol", "score", "Strand")
"""The six BED fields. A seventh column, if present, is taken as the gene id."""


class PromoterError(RuntimeError):
    """The promoter file is unusable, or shares no gene with the store."""


def read_promoters(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read a TSS BED file.

    Parameters
    ----------
    path
        A BED of transcription start sites, one per gene. Six columns at
        least; a seventh is read as the gene id, which is what lets a TSS be
        tied back to the store. A first line beginning with ``#`` is a header
        and is skipped.

    Returns
    -------
    pandas.DataFrame
        ``Chromosome``, ``Start``, ``End``, ``Strand``, ``gene_symbol`` and
        ``ensg``, with ``ensg`` empty when the file has no seventh column.

    Raises
    ------
    PromoterError
        The file is missing, empty, or has fewer than six columns.
    """
    path = Path(path)
    if not path.is_file():
        raise PromoterError(f"The promoter file {path} does not exist.")

    with open(path, encoding="utf-8") as handle:
        first = handle.readline()
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        skiprows=1 if first.startswith("#") else 0,
        dtype={0: str},
    )
    if frame.shape[1] < len(_BED6):
        raise PromoterError(
            f"{path} has {frame.shape[1]} columns; a TSS BED needs at least"
            f" {len(_BED6)}: {', '.join(_BED6)}."
        )
    if frame.empty:
        raise PromoterError(f"{path} holds no transcription start site.")

    names = list(_BED6) + ["ensg"] + [f"extra_{i}" for i in range(frame.shape[1] - 7)]
    frame.columns = names[: frame.shape[1]]
    if "ensg" not in frame:
        frame["ensg"] = ""
    return frame.assign(ensg=_unversioned(frame["ensg"]))[
        ["Chromosome", "Start", "End", "Strand", "gene_symbol", "ensg"]
    ]


def promoters_from_store(store, *, reference_only: bool = True) -> pd.DataFrame:
    """Derive one TSS per gene from the store, for when no BED is given.

    Parameters
    ----------
    store
        A loaded `Store`.
    reference_only
        Take the start of each gene's reference transcript. False takes the
        start of the gene itself.

    Returns
    -------
    pandas.DataFrame
        The same columns `read_promoters` returns, each TSS one base wide.
    """
    rows = store.rows("transcript", reference_only=reference_only)
    if rows.empty:
        rows = store.rows("gene")
    forward = rows["Strand"] == "+"
    start = rows["Start"].where(forward, rows["End"] - 1)
    return pd.DataFrame(
        {
            "Chromosome": rows["Chromosome"].to_numpy(),
            "Start": start.to_numpy(),
            "End": start.to_numpy() + 1,
            "Strand": rows["Strand"].to_numpy(),
            "gene_symbol": rows["gene_symbol"].to_numpy(),
            "ensg": rows["ensg"].to_numpy(),
        }
    )


def _unversioned(ids: pd.Series) -> pd.Series:
    """ENSG00000186092.7 -> ENSG00000186092."""
    return ids.astype(str).str.replace(r"\.\d+$", "", regex=True)

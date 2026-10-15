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
        A BED of transcription start sites, one row per site: a gene with
        several promoters needs one row each. Six columns at least; a seventh
        is read as the gene id, which is what lets a TSS be tied back to the
        store, and an eighth as the transcript id, which is what fills
        ``promoter_enst``. A first line beginning with ``#`` is a header and
        is skipped.

    Returns
    -------
    pandas.DataFrame
        ``Chromosome``, ``Start``, ``End``, ``Strand``, ``gene_symbol``,
        ``ensg`` and ``enst``, the last two empty when the file has no
        seventh or eighth column.

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
    try:
        frame = pd.read_csv(
            path,
            sep="\t",
            header=None,
            skiprows=1 if first.startswith("#") else 0,
            dtype={0: str},
        )
    except pd.errors.EmptyDataError as error:
        raise PromoterError(f"{path} holds no transcription start site.") from error
    if frame.shape[1] < len(_BED6):
        raise PromoterError(
            f"{path} has {frame.shape[1]} columns; a TSS BED needs at least"
            f" {len(_BED6)}: {', '.join(_BED6)}."
        )
    if frame.empty:
        raise PromoterError(f"{path} holds no transcription start site.")

    names = (
        list(_BED6)
        + ["ensg", "enst"]
        + [f"extra_{i}" for i in range(frame.shape[1] - 8)]
    )
    frame.columns = names[: frame.shape[1]]
    for column in ("ensg", "enst"):
        if column not in frame:
            frame[column] = ""
    return frame.assign(
        ensg=_unversioned(frame["ensg"]), enst=_unversioned(frame["enst"])
    )[["Chromosome", "Start", "End", "Strand", "gene_symbol", "ensg", "enst"]]


def promoters_from_store(store, *, reference_only: bool = False) -> pd.DataFrame:
    """Derive the TSS of every transcript from the store.

    A gene has as many promoters as it has distinct transcript starts, and a
    guide at an alternative promoter is invisible to a one-TSS-per-gene file,
    so every transcript contributes a site by default.

    Parameters
    ----------
    store
        A loaded `Store`.
    reference_only
        Take only each gene's reference transcript, one TSS per gene. This is
        the older behaviour and loses alternative promoters.

    Returns
    -------
    pandas.DataFrame
        The same columns `read_promoters` returns, each TSS one base wide,
        deduplicated so two transcripts sharing a start give one row per
        transcript but not two identical ones.
    """
    rows = store.rows("transcript", reference_only=reference_only)
    if rows.empty:
        rows = store.rows("gene")
    forward = rows["Strand"] == "+"
    start = rows["Start"].where(forward, rows["End"] - 1)
    frame = pd.DataFrame(
        {
            "Chromosome": rows["Chromosome"].to_numpy(),
            "Start": start.to_numpy(),
            "End": start.to_numpy() + 1,
            "Strand": rows["Strand"].to_numpy(),
            "gene_symbol": rows["gene_symbol"].to_numpy(),
            "ensg": rows["ensg"].to_numpy(),
            "enst": rows["enst"].to_numpy() if "enst" in rows else "",
        }
    )
    return frame.drop_duplicates().reset_index(drop=True)


def _unversioned(ids: pd.Series) -> pd.Series:
    """ENSG00000186092.7 -> ENSG00000186092."""
    return ids.astype(str).str.replace(r"\.\d+$", "", regex=True)

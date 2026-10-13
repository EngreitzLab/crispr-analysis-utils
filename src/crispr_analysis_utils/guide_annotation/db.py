"""The annotation store: a GENCODE GTF reduced to what annotating needs.

Reading a full GENCODE takes most of a minute, so a run reads this instead.
The store is parquet, which `pyranges1` loads in about a second.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .._version import __version__

FEATURES = ("gene", "transcript", "exon", "CDS", "UTR")
"""The feature rows the store keeps. Everything else is dropped."""

GENE_TYPES = ("protein_coding",)
"""The gene types kept by default."""

COLUMNS = (
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Feature",
    "ensg",
    "gene_symbol",
    "enst",
    "exon_number",
    "is_reference",
)
"""The store's columns. The first four are what pyranges needs."""

_METADATA_KEY = b"crispr_analysis_utils"


class StoreError(RuntimeError):
    """The store is missing, unreadable, or holds no usable annotation."""


@dataclass
class Store:
    """A loaded annotation store.

    Attributes
    ----------
    features
        Every kept row: genes, transcripts, exons, CDS and UTR.
    provenance
        What the store was built from, as `build_store` recorded it.
    """

    features: pd.DataFrame
    provenance: dict

    def rows(self, feature: str, *, reference_only: bool = False) -> pd.DataFrame:
        """The rows of one feature, optionally only the reference transcript's."""
        rows = self.features[self.features["Feature"] == feature]
        if reference_only:
            rows = rows[rows["is_reference"]]
        return rows

    @property
    def contigs(self) -> set[str]:
        """Every contig the store holds."""
        return set(self.features["Chromosome"].unique())


def build_store(
    gtf: str | os.PathLike[str],
    out: str | os.PathLike[str],
    *,
    gene_types: tuple[str, ...] = GENE_TYPES,
) -> Store:
    """Reduce a GENCODE GTF to an annotation store.

    Parameters
    ----------
    gtf
        A GENCODE GTF, gzipped or not. It must hold `transcript` and `exon`
        rows: a genes-only GTF cannot say exon, intron or transcript.
    out
        The parquet file to write.
    gene_types
        The gene types to keep.

    Returns
    -------
    Store
        The store just written.

    Raises
    ------
    StoreError
        The GTF holds no gene of those types, or no exons.
    """
    import pyranges1 as pr

    gtf = Path(gtf)
    if not gtf.is_file():
        raise FileNotFoundError(f"The GTF {gtf} does not exist.")

    # duplicate_attr is not optional: a GENCODE row carries several `tag`
    # keys, and without it the reader keeps one and drops most of MANE_Select.
    frame = pr.read_gtf(gtf, duplicate_attr=True)
    frame = frame.loc[
        frame["gene_type"].isin(gene_types) & frame["Feature"].isin(FEATURES)
    ]
    if frame.empty:
        raise StoreError(
            f"{gtf} holds no {' or '.join(gene_types)} gene. Is it the right GTF?"
        )
    if not (frame["Feature"] == "exon").any():
        raise StoreError(
            f"{gtf} holds no exon rows, so it cannot say exon, intron or"
            " transcript. A genes-only GTF is not enough; use the full"
            " GENCODE annotation."
        )

    tags = frame["tag"].fillna("")
    frame = frame.assign(
        ensg=_unversioned(frame["gene_id"]),
        enst=_unversioned(frame["transcript_id"]),
        gene_symbol=frame["gene_name"],
        mane_select=tags.str.contains("MANE_Select"),
        ensembl_canonical=tags.str.contains("Ensembl_canonical"),
    )

    reference = _reference_transcripts(frame)
    frame = frame.assign(
        is_reference=frame["ensg"].map(reference).eq(frame["enst"])
        & frame["enst"].notna()
    )

    provenance = {
        "gtf": str(gtf),
        "gene_types": list(gene_types),
        "cau_version": __version__,
        "counts": {
            feature: int((frame["Feature"] == feature).sum()) for feature in FEATURES
        },
        "genes": int(frame["ensg"][frame["Feature"] == "gene"].nunique()),
        "reference_transcripts": int(len(reference)),
        "mane_select": int(
            frame["mane_select"][frame["Feature"] == "transcript"].sum()
        ),
    }

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    # A plain DataFrame, so a built store and a loaded one are the same type.
    table = pd.DataFrame(frame[list(COLUMNS)]).reset_index(drop=True)
    _write_parquet(table, out, provenance)
    return Store(features=table, provenance=provenance)


def load_store(path: str | os.PathLike[str]) -> Store:
    """Load a store written by `build_store`.

    Parameters
    ----------
    path
        The parquet file.

    Returns
    -------
    Store
        Its rows and its provenance.

    Raises
    ------
    StoreError
        The file is missing, or was not written by `build_store`.
    """
    import pyarrow.parquet as pq

    path = Path(path)
    if not path.is_file():
        raise StoreError(
            f"The annotation store {path} does not exist. Build one with"
            " `cau guide-annotation build-db`."
        )
    table = pq.read_table(path)
    missing = set(COLUMNS) - set(table.column_names)
    if missing:
        raise StoreError(
            f"{path} is missing {', '.join(sorted(missing))}, so it was not"
            " written by this tool. Rebuild it."
        )
    raw = (table.schema.metadata or {}).get(_METADATA_KEY)
    return Store(
        features=table.to_pandas(),
        provenance=json.loads(raw) if raw else {},
    )


def _unversioned(ids: pd.Series) -> pd.Series:
    """ENSG00000186092.7 -> ENSG00000186092, leaving anything else alone."""
    return ids.str.replace(r"\.\d+$", "", regex=True)


def _reference_transcripts(frame: pd.DataFrame) -> pd.Series:
    """Each gene's reference transcript: MANE Select, else canonical, else longest CDS.

    One transcript has to be chosen for exon_number and position along the
    CDS to mean anything. The ordering is deterministic down to the id, so a
    rebuild picks the same one.
    """
    coding = frame[frame["Feature"] == "CDS"]
    cds_length = (coding["End"] - coding["Start"]).groupby(coding["enst"]).sum()

    transcripts = frame[frame["Feature"] == "transcript"].copy()
    transcripts["cds_length"] = transcripts["enst"].map(cds_length).fillna(0)
    transcripts = transcripts.sort_values(
        ["ensg", "mane_select", "ensembl_canonical", "cds_length", "enst"],
        ascending=[True, False, False, False, True],
    )
    return transcripts.groupby("ensg")["enst"].first()


def _write_parquet(table: pd.DataFrame, path: Path, provenance: dict) -> None:
    """Write the store, carrying its provenance in the parquet metadata."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    arrow = pa.Table.from_pandas(table, preserve_index=False)
    metadata = dict(arrow.schema.metadata or {})
    metadata[_METADATA_KEY] = json.dumps(provenance).encode("utf-8")
    pq.write_table(arrow.replace_schema_metadata(metadata), path, compression="zstd")

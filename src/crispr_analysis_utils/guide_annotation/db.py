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
"""The gene types a guide is annotated against.

Every call of genic, exon, intron and promoter follows these. Genes of other
types are kept as one row each, so a protospacer overlapping one can be
reported without that overlap deciding where the guide sits.
"""

COLUMNS = (
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Feature",
    "ensg",
    "gene_symbol",
    "enst",
    "exon_id",
    "exon_number",
    "is_reference",
    "gene_type",
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

    @property
    def gene_types(self) -> tuple[str, ...]:
        """The gene types this store's annotation is called against."""
        return tuple(self.provenance.get("gene_types") or GENE_TYPES)

    def rows(
        self, feature: str, *, reference_only: bool = False, primary: bool = True
    ) -> pd.DataFrame:
        """The rows of one feature.

        Parameters
        ----------
        feature
            One of `FEATURES`.
        reference_only
            Keep only the reference transcript's rows.
        primary
            Keep only genes of `gene_types`. On by default, so a store that
            also carries other gene types cannot widen what counts as genic.
        """
        rows = self.features[self.features["Feature"] == feature]
        if primary:
            rows = rows[rows["gene_type"].isin(self.gene_types)]
        if reference_only:
            rows = rows[rows["is_reference"]]
        return rows

    def other_genes(self) -> pd.DataFrame:
        """Gene rows whose type is not one the annotation is called against."""
        genes = self.features[self.features["Feature"] == "gene"]
        return genes[~genes["gene_type"].isin(self.gene_types)]

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
        The gene types a guide is annotated against: genic, exon, intron and
        promoter all follow these, and every feature row is kept for them.
        Genes of every other type are kept as one row each, so an overlap can
        be reported without it deciding where the guide sits.

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
    # Genes of the annotated types keep every feature row; genes of any other
    # type keep one row each, enough to report an overlap without that overlap
    # deciding where a guide sits.
    primary = frame["gene_type"].isin(gene_types)
    frame = frame.loc[
        (primary & frame["Feature"].isin(FEATURES))
        | (~primary & (frame["Feature"] == "gene"))
    ]
    kept_primary = frame["gene_type"].isin(gene_types)
    if not kept_primary.any():
        raise StoreError(
            f"{gtf} holds no {' or '.join(gene_types)} gene. Is it the right GTF?"
        )
    if not ((frame["Feature"] == "exon") & kept_primary).any():
        raise StoreError(
            f"{gtf} holds no exon rows, so it cannot say exon, intron or"
            " transcript. A genes-only GTF is not enough; use the full"
            " GENCODE annotation."
        )

    tags = frame["tag"].fillna("")
    frame = frame.assign(
        ensg=_unversioned(frame["gene_id"]),
        enst=_unversioned(frame["transcript_id"]),
        exon_id=_optional(frame, "exon_id"),
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
        "other_gene_types": sorted(
            frame["gene_type"][~frame["gene_type"].isin(gene_types)].unique().tolist()
        ),
        "cau_version": __version__,
        "counts": {
            feature: int((frame["Feature"] == feature).sum()) for feature in FEATURES
        },
        "genes": int(
            frame["ensg"][
                (frame["Feature"] == "gene") & frame["gene_type"].isin(gene_types)
            ].nunique()
        ),
        "other_genes": int(
            frame["ensg"][
                (frame["Feature"] == "gene") & ~frame["gene_type"].isin(gene_types)
            ].nunique()
        ),
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
        # gene_type arrived with the other-gene-type reporting, so a store
        # built before it is simply old rather than foreign.
        rebuild = (
            "Rebuild it with `cau guide-annotation build-db`."
            if missing == {"gene_type"}
            else "Rebuild it."
        )
        raise StoreError(
            f"{path} is missing {', '.join(sorted(missing))}, so it was not"
            f" written by this version of the tool. {rebuild}"
        )
    raw = (table.schema.metadata or {}).get(_METADATA_KEY)
    return Store(
        features=table.to_pandas(),
        provenance=json.loads(raw) if raw else {},
    )


def _optional(frame: pd.DataFrame, name: str) -> pd.Series:
    """An id column the GTF may not carry, without its version.

    GENCODE gives exons an `exon_id`; a GTF from elsewhere need not.
    """
    # A declared string dtype, because an all-missing object column loses its
    # type through parquet and comes back as floats.
    if name not in frame:
        return pd.Series(
            pd.array([None] * len(frame), dtype="string"), index=frame.index
        )
    return _unversioned(frame[name]).astype("string")


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

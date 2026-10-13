"""Place each aligned site in a gene: exon, intron, promoter or intergenic.

Overlaps are decided on the protospacer, so a guide reaching into an exon
counts as exonic. Distances and the position along a CDS are measured from
the cut site, which is where the break falls.

Every interval operation ignores strand: a guide on one strand targets a gene
on either.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .db import Store
from .promoters import PromoterError

EXON = "exon"
INTRON = "intron"
INTERGENIC = "intergenic"
LOCATIONS = (EXON, INTRON, INTERGENIC)
"""Where a site falls, in the order qc.tsv counts them."""

CDS = "CDS"
UTR = "UTR"

SENSE = "sense"
ANTISENSE = "antisense"

SITE_COLUMNS = ("chr", "start", "end", "strand", "cut_site")
"""What `annotate_sites` needs of every site."""

_IGNORE = "ignore"
"""Strand behaviour for every join: a guide targets a gene on either strand."""


class AnnotationError(RuntimeError):
    """The sites and the annotation cannot be compared at all."""


def annotate_sites(
    sites: pd.DataFrame,
    store: Store,
    promoters: pd.DataFrame,
    *,
    promoter_window: int = 1000,
) -> pd.DataFrame:
    """Annotate every site against the genes of the store.

    Parameters
    ----------
    sites
        One row per site, with `SITE_COLUMNS` at least. Any other column is
        carried through to the output untouched.
    store
        The annotation store, as `load_store` returns it.
    promoters
        Transcription start sites, as `read_promoters` returns them.
    promoter_window
        How far either side of a TSS counts as its promoter.

    Returns
    -------
    pandas.DataFrame
        One row per site per overlapping gene; a site in no gene keeps one
        row, `intergenic`, with the gene columns empty. The input's own
        columns come first.

    Raises
    ------
    AnnotationError
        No site contig appears in the store, so nothing could ever overlap.
    """
    import pyranges1 as pr

    missing = set(SITE_COLUMNS) - set(sites.columns)
    if missing:
        raise AnnotationError(f"The sites are missing {', '.join(sorted(missing))}.")
    if promoter_window < 0:
        raise ValueError(
            f"promoter_window must not be negative, got {promoter_window}."
        )

    sites = sites.reset_index(drop=True)
    _check_contigs(sites, store)

    spacers = pr.PyRanges(
        {
            "Chromosome": sites["chr"].to_numpy(),
            "Start": sites["start"].to_numpy(),
            "End": sites["end"].to_numpy(),
            "Strand": sites["strand"].to_numpy(),
            "_row": np.arange(len(sites)),
        }
    )
    cuts = pr.PyRanges(
        {
            "Chromosome": sites["chr"].to_numpy(),
            "Start": sites["cut_site"].to_numpy(),
            "End": sites["cut_site"].to_numpy() + 1,
            "Strand": sites["strand"].to_numpy(),
            "_row": np.arange(len(sites)),
        }
    )

    pairs = _gene_pairs(spacers, store, pr)
    pairs = pairs.merge(
        _exon_calls(spacers, store, pr), on=["_row", "ensg"], how="left"
    )
    pairs = pairs.merge(_cds_position(cuts, store, pr), on=["_row", "ensg"], how="left")

    pairs["location"] = np.where(
        pairs["ensg"].isna(),
        INTERGENIC,
        np.where(pairs["enst"].notna(), EXON, INTRON),
    )
    pairs["region"] = np.where(
        pairs["location"] == EXON,
        np.where(pairs["cds_position"].notna(), CDS, UTR),
        "",
    )
    pairs["exonic_any_transcript"] = pairs["exonic_any_transcript"].fillna(False)
    pairs["sense"] = np.where(
        pairs["gene_strand"].isna(),
        "",
        np.where(pairs["strand_site"] == pairs["gene_strand"], SENSE, ANTISENSE),
    )

    tss = _tss_columns(cuts, spacers, promoters, promoter_window, pr)
    frame = sites.join(pairs.set_index("_row").drop(columns=["strand_site"]))
    frame = frame.merge(tss, left_index=True, right_index=True, how="left")
    frame["gc_percent"] = _gc_percent(frame)
    frame["max_polyt_run"] = _max_polyt(frame)
    return frame.reset_index(drop=True)


def _check_contigs(sites: pd.DataFrame, store: Store) -> None:
    """Fail when no site contig is in the store, rather than call everything intergenic."""
    shared = set(sites["chr"].unique()) & store.contigs
    if not shared:
        raise AnnotationError(
            "No contig of the sites appears in the annotation store. The "
            f"sites use {_sample(sites['chr'])} and the store uses "
            f"{_sample(pd.Series(sorted(store.contigs)))}. Do they name "
            "contigs the same way?"
        )


def _sample(values: pd.Series, n: int = 3) -> str:
    """A few distinct values, for an error message."""
    names = list(dict.fromkeys(values.astype(str)))[:n]
    return ", ".join(names) + ("..." if len(names) == n else "")


def _gene_pairs(spacers, store: Store, pr) -> pd.DataFrame:
    """One row per site per overlapping gene, keeping sites that overlap none."""
    genes = pr.PyRanges(
        store.rows("gene")[
            ["Chromosome", "Start", "End", "Strand", "ensg", "gene_symbol"]
        ].reset_index(drop=True)
    )
    hits = spacers.join_overlaps(
        genes, join_type="left", strand_behavior=_IGNORE, suffix="_gene"
    )
    frame = pd.DataFrame(hits)
    return pd.DataFrame(
        {
            "_row": frame["_row"],
            "ensg": _joined(frame, "ensg", "_gene"),
            "gene_symbol": _joined(frame, "gene_symbol", "_gene"),
            "gene_strand": _joined(frame, "Strand", "_gene"),
            "strand_site": frame["Strand"],
        }
    )


def _exon_calls(spacers, store: Store, pr) -> pd.DataFrame:
    """Per site and gene: the reference exon hit, and whether any transcript has one."""
    reference = store.rows("exon", reference_only=True)
    columns = ["Chromosome", "Start", "End", "Strand", "ensg", "enst", "exon_number"]
    hits = spacers.join_overlaps(
        pr.PyRanges(reference[columns].reset_index(drop=True)),
        join_type="inner",
        strand_behavior=_IGNORE,
        suffix="_exon",
    )
    best = (
        pd.DataFrame(hits)[["_row", "ensg", "enst", "exon_number"]]
        .sort_values(["_row", "ensg", "exon_number"])
        .drop_duplicates(["_row", "ensg"])
    )

    every = store.rows("exon")
    any_hits = spacers.join_overlaps(
        pr.PyRanges(
            every[["Chromosome", "Start", "End", "Strand", "ensg"]].reset_index(
                drop=True
            )
        ),
        join_type="inner",
        strand_behavior=_IGNORE,
        suffix="_exon",
    )
    anywhere = (
        pd.DataFrame(any_hits)[["_row", "ensg"]]
        .drop_duplicates()
        .assign(exonic_any_transcript=True)
    )
    return best.merge(anywhere, on=["_row", "ensg"], how="outer")


def _cds_position(cuts, store: Store, pr) -> pd.DataFrame:
    """How far along the reference CDS the cut falls, and the CDS length."""
    blocks = _cds_blocks(store)
    if blocks.empty:
        return pd.DataFrame(columns=["_row", "ensg", "cds_position", "cds_fraction"])
    columns = ["Chromosome", "Start", "End", "Strand", "ensg", "_offset", "_total"]
    hits = pd.DataFrame(
        cuts.join_overlaps(
            pr.PyRanges(blocks[columns].reset_index(drop=True)),
            join_type="inner",
            strand_behavior=_IGNORE,
            suffix="_cds",
        )
    )
    if hits.empty:
        return pd.DataFrame(columns=["_row", "ensg", "cds_position", "cds_fraction"])

    forward = hits["Strand_cds"] == "+"
    within = np.where(
        forward,
        hits["Start"] - hits["Start_cds"],
        hits["End_cds"] - 1 - hits["Start"],
    )
    hits["cds_position"] = hits["_offset"] + within
    hits["cds_fraction"] = (hits["cds_position"] / hits["_total"]).round(4)
    return hits[["_row", "ensg", "cds_position", "cds_fraction"]].drop_duplicates(
        ["_row", "ensg"]
    )


def _cds_blocks(store: Store) -> pd.DataFrame:
    """The reference transcripts' CDS blocks, each with its offset along the CDS.

    Blocks are ordered 5' to 3', so the offset is how much coding sequence
    comes before the block in the transcript, not on the chromosome.
    """
    blocks = store.rows("CDS", reference_only=True).copy()
    if blocks.empty:
        return blocks
    blocks["_length"] = blocks["End"] - blocks["Start"]
    blocks["_key"] = np.where(
        blocks["Strand"] == "+", blocks["Start"], -blocks["Start"]
    )
    blocks = blocks.sort_values(["enst", "_key"])
    grouped = blocks.groupby("enst", sort=False)["_length"]
    blocks["_offset"] = grouped.cumsum() - blocks["_length"]
    blocks["_total"] = grouped.transform("sum")
    return blocks


def _tss_columns(
    cuts, spacers, promoters: pd.DataFrame, window: int, pr
) -> pd.DataFrame:
    """Promoter flag, and the signed distance to the nearest protein-coding TSS."""
    if promoters.empty:
        raise PromoterError("The promoter file holds no transcription start site.")

    tss = pr.PyRanges(
        promoters.rename(columns={"gene_symbol": "tss_gene", "ensg": "tss_ensg"})[
            ["Chromosome", "Start", "End", "Strand", "tss_gene", "tss_ensg"]
        ].reset_index(drop=True)
    )
    windows = pr.PyRanges(
        {
            "Chromosome": promoters["Chromosome"].to_numpy(),
            "Start": np.maximum(promoters["Start"].to_numpy() - window, 0),
            "End": promoters["End"].to_numpy() + window,
            "Strand": promoters["Strand"].to_numpy(),
        }
    )
    inside = pd.DataFrame(
        spacers.join_overlaps(
            windows, join_type="inner", strand_behavior=_IGNORE, suffix="_promoter"
        )
    )
    promoter = pd.Series(False, index=range(len(cuts)), name="promoter")
    if not inside.empty:
        promoter.iloc[np.unique(inside["_row"].to_numpy())] = True

    near = pd.DataFrame(
        cuts.nearest_ranges(
            tss,
            strand_behavior=_IGNORE,
            suffix="_tss",
            dist_col="_distance",
            ties="first",
        )
    )
    near = near.drop_duplicates("_row").set_index("_row")
    rows = range(len(promoter))
    # Signed in the gene's orientation: negative upstream of the TSS.
    toward = near["Start"] - _joined(near, "Start", "_tss")
    signed = np.where(_joined(near, "Strand", "_tss") == "+", toward, -toward)
    return pd.DataFrame(
        {
            "promoter": promoter.to_numpy(),
            "distance_to_tss": pd.Series(signed, index=near.index).reindex(rows),
            "nearest_tss_gene": _joined(near, "tss_gene", "_tss").reindex(rows),
            "nearest_tss_ensg": _joined(near, "tss_ensg", "_tss").reindex(rows),
        }
    )


def _joined(frame: pd.DataFrame, name: str, suffix: str) -> pd.Series:
    """A column of the joined-in table, suffixed or not.

    `join_overlaps` suffixes only the names that clash; `nearest_ranges`
    suffixes all of them.
    """
    if name + suffix in frame.columns:
        return frame[name + suffix]
    return frame[name]


def _gc_percent(frame: pd.DataFrame) -> pd.Series:
    """The spacer's GC percentage, to one decimal."""
    if "spacer" not in frame:
        return pd.Series(np.nan, index=frame.index)
    spacer = frame["spacer"].fillna("").astype(str).str.upper()
    length = spacer.str.len()
    gc = spacer.str.count("[GC]")
    return (100 * gc / length.where(length > 0)).round(1)


def _max_polyt(frame: pd.DataFrame) -> pd.Series:
    """The longest run of T in the spacer; 4 or more terminates Pol III."""
    if "spacer" not in frame:
        return pd.Series(pd.NA, index=frame.index, dtype="Int64")
    spacer = frame["spacer"].fillna("").astype(str).str.upper()
    runs = spacer.str.findall("T+").map(
        lambda hits: max((len(h) for h in hits), default=0)
    )
    return runs.astype("Int64")

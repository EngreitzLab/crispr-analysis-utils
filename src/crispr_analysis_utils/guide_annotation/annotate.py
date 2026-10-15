"""Place each aligned site in a gene: exon, intron, promoter or intergenic.

Overlaps are decided on the protospacer, so a guide reaching into an exon
counts as exonic. Distances and the position along a CDS are measured from
the target site, which is where the break falls.

Every interval operation ignores strand: a guide on one strand targets a gene
on either.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .db import Store
from .promoters import PromoterError

GENIC = "genic"
INTERGENIC = "intergenic"
EXON = "exon"
INTRON = "intron"
DELIMITER = ";"
"""Separates the ids in a list-valued column.

A semicolon rather than a comma, because a comma already separates ids inside
the IGVF format's own ``putative_target_genes``.
"""
LOCATIONS = (GENIC, INTERGENIC)
FEATURES = (EXON, INTRON)
"""Where a site falls, in the order qc.tsv counts them."""

CDS = "CDS"
UTR = "UTR"

SENSE = "sense"
ANTISENSE = "antisense"

SITE_COLUMNS = ("chr", "start", "end", "strand", "target_site")
"""What `annotate_sites` needs of every site."""

MIDPOINT = "midpoint"
TARGET_SITE = "target_site"
DISTANCE_FROM = (MIDPOINT, TARGET_SITE)
"""Where a distance is measured from.

The midpoint of the protospacer by default. A dCas9 occupies the whole
protospacer rather than acting at one base, so its centre is the more honest
reference for how far a guide sits from a TSS. `target_site` measures from
the boundary 3 bases 5' of the PAM instead, which is where a nuclease cuts.
"""

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
    distance_from: str = MIDPOINT,
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
    if distance_from not in DISTANCE_FROM:
        raise ValueError(
            f"distance_from must be one of {', '.join(DISTANCE_FROM)}, got"
            f" {distance_from!r}."
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
    # The CDS position asks where along the coding sequence Cas9 acts, so it
    # stays on the target site whatever distances are measured from.
    targets = pr.PyRanges(
        {
            "Chromosome": sites["chr"].to_numpy(),
            "Start": sites["target_site"].to_numpy(),
            "End": sites["target_site"].to_numpy() + 1,
            "Strand": sites["strand"].to_numpy(),
            "_row": np.arange(len(sites)),
        }
    )

    # Where distances are measured from. The midpoint rounds down, so a
    # protospacer of even length takes the base left of centre.
    anchor_at = (
        (sites["start"].to_numpy() + sites["end"].to_numpy()) // 2
        if distance_from == MIDPOINT
        else sites["target_site"].to_numpy()
    )
    anchors = pr.PyRanges(
        {
            "Chromosome": sites["chr"].to_numpy(),
            "Start": anchor_at,
            "End": anchor_at + 1,
            "Strand": sites["strand"].to_numpy(),
            "_row": np.arange(len(sites)),
        }
    )

    pairs = _gene_pairs(spacers, store, pr)
    pairs = pairs.merge(
        _exon_calls(spacers, store, pr), on=["_row", "ensg"], how="left"
    )
    pairs = pairs.merge(
        _cds_position(targets, store, pr), on=["_row", "ensg"], how="left"
    )

    # Two axes, not one. Whether a gene is there at all, and where inside it:
    # a guide can be genic and in a promoter and in an exon at once, which a
    # single column could not say.
    pairs["location"] = np.where(pairs["ensg"].isna(), INTERGENIC, GENIC)
    pairs["feature"] = np.where(
        pairs["ensg"].isna(),
        "",
        np.where(pairs["enst"].notna(), EXON, INTRON),
    )
    pairs["region"] = np.where(
        pairs["feature"] == EXON,
        np.where(pairs["cds_position"].notna(), CDS, UTR),
        "",
    )
    pairs["exonic_any_transcript"] = pairs["exonic_any_transcript"].fillna(False)
    # How far into the gene the target site falls, from the gene's own 5'
    # end. Near zero is the annotated start; a large one is deeper in, which
    # is what distinguishes a canonical promoter from an alternative one.
    five_prime = np.where(
        pairs["gene_strand"] == "+", pairs["gene_start"], pairs["gene_end"] - 1
    )
    toward = anchor_at[pairs["_row"].to_numpy()] - five_prime
    pairs["distance_to_gene_start"] = np.where(
        pairs["gene_strand"].isna(),
        np.nan,
        np.where(pairs["gene_strand"] == "+", toward, -toward),
    )

    pairs["sense"] = np.where(
        pairs["gene_strand"].isna(),
        "",
        np.where(pairs["strand_site"] == pairs["gene_strand"], SENSE, ANTISENSE),
    )

    tss = _tss_columns(anchors, spacers, promoters, promoter_window, pr)
    others = _other_gene_overlaps(spacers, store, pr, len(sites))
    frame = sites.join(
        pairs.set_index("_row").drop(columns=["strand_site", "gene_start", "gene_end"])
    )
    frame = frame.merge(tss, left_index=True, right_index=True, how="left")
    frame = frame.merge(others, left_index=True, right_index=True, how="left")
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
            "gene_start": _joined(frame, "Start", "_gene"),
            "gene_end": _joined(frame, "End", "_gene"),
            "strand_site": frame["Strand"],
        }
    )


def _other_gene_overlaps(spacers, store: Store, pr, n_rows: int) -> pd.DataFrame:
    """Per site: the genes of other types the protospacer overlaps.

    Site-level, not per gene: an overlap with a lncRNA says the guide sits in
    ambiguous territory whichever protein-coding gene the row is about.
    """
    empty = pd.DataFrame(
        {
            "other_gene_types": pd.Series([""] * n_rows, dtype="string"),
            "other_genes": pd.Series([""] * n_rows, dtype="string"),
        },
        index=range(n_rows),
    )
    others = store.other_genes()
    if others.empty:
        return empty

    hits = spacers.join_overlaps(
        pr.PyRanges(
            others[
                ["Chromosome", "Start", "End", "Strand", "ensg", "gene_type"]
            ].reset_index(drop=True)
        ),
        join_type="inner",
        strand_behavior=_IGNORE,
        suffix="_other",
    )
    frame = pd.DataFrame(hits)
    if frame.empty:
        return empty
    frame = pd.DataFrame(
        {
            "_row": frame["_row"],
            "ensg": _joined(frame, "ensg", "_other"),
            "gene_type": _joined(frame, "gene_type", "_other"),
        }
    )
    grouped = frame.groupby("_row").agg(
        other_gene_types=("gene_type", _joined_uniques),
        other_genes=("ensg", _joined_uniques),
    )
    rows = range(n_rows)
    return pd.DataFrame(
        {
            name: grouped[name].reindex(rows).fillna("").astype("string")
            for name in ("other_gene_types", "other_genes")
        },
        index=rows,
    )


def _joined_uniques(values: pd.Series) -> str:
    """The distinct values, sorted, as one delimited string."""
    return DELIMITER.join(sorted({str(v) for v in values if pd.notna(v) and str(v)}))


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
    # exon_number comes off the GTF as a string, where "10" sorts before "2".
    best = pd.DataFrame(hits)[["_row", "ensg", "enst", "exon_number"]].assign(
        _order=lambda frame: pd.to_numeric(frame["exon_number"], errors="coerce")
    )
    best = (
        best.sort_values(["_row", "ensg", "_order"])
        .drop_duplicates(["_row", "ensg"])
        .drop(columns="_order")
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


def _cds_position(targets, store: Store, pr) -> pd.DataFrame:
    """How far along the reference CDS the target site falls, and its length."""
    blocks = _cds_blocks(store)
    if blocks.empty:
        return pd.DataFrame(columns=["_row", "ensg", "cds_position", "cds_fraction"])
    columns = ["Chromosome", "Start", "End", "Strand", "ensg", "_offset", "_total"]
    hits = pd.DataFrame(
        targets.join_overlaps(
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
    anchors, spacers, promoters: pd.DataFrame, window: int, pr
) -> pd.DataFrame:
    """Promoter flag, and the signed distance to the nearest protein-coding TSS."""
    if promoters.empty:
        raise PromoterError("The promoter file holds no transcription start site.")

    tss = pr.PyRanges(
        promoters.rename(columns={"gene_symbol": "tss_gene", "ensg": "tss_ensg"})[
            ["Chromosome", "Start", "End", "Strand", "tss_gene", "tss_ensg"]
        ].reset_index(drop=True)
    )
    # The windows carry their own ids, so an overlap says which transcripts'
    # promoters a guide sits in, not merely that it sits in one.
    windows = pr.PyRanges(
        {
            "Chromosome": promoters["Chromosome"].to_numpy(),
            "Start": np.maximum(promoters["Start"].to_numpy() - window, 0),
            "End": promoters["End"].to_numpy() + window,
            "Strand": promoters["Strand"].to_numpy(),
            "w_gene": promoters["gene_symbol"].to_numpy(),
            "w_ensg": promoters["ensg"].to_numpy(),
            "w_enst": (
                promoters["enst"].to_numpy()
                if "enst" in promoters
                else np.full(len(promoters), "")
            ),
        }
    )
    inside = pd.DataFrame(
        spacers.join_overlaps(
            windows, join_type="inner", strand_behavior=_IGNORE, suffix="_promoter"
        )
    )
    rows = range(len(anchors))
    promoter = pd.Series(False, index=rows, name="promoter")
    hits = {
        name: pd.Series([""] * len(anchors), index=rows, dtype="string")
        for name in ("promoter_enst", "promoter_ensg", "promoter_gene")
    }
    if not inside.empty:
        promoter.iloc[np.unique(inside["_row"].to_numpy())] = True
        grouped = inside.groupby("_row").agg(
            promoter_enst=("w_enst", _joined_uniques),
            promoter_ensg=("w_ensg", _joined_uniques),
            promoter_gene=("w_gene", _joined_uniques),
        )
        hits = {
            name: grouped[name].reindex(rows).fillna("").astype("string")
            for name in hits
        }

    near = pd.DataFrame(
        anchors.nearest_ranges(
            tss,
            strand_behavior=_IGNORE,
            suffix="_tss",
            dist_col="_distance",
            ties="first",
        )
    )
    near = near.drop_duplicates("_row").set_index("_row")
    # Signed in the gene's orientation: negative upstream of the TSS.
    toward = near["Start"] - _joined(near, "Start", "_tss")
    signed = np.where(_joined(near, "Strand", "_tss") == "+", toward, -toward)
    return pd.DataFrame(
        {
            "promoter": promoter.to_numpy(),
            **hits,
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

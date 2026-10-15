"""The annotation table, and the data dictionary that describes it."""

from __future__ import annotations

import os
from collections import Counter
from pathlib import Path

import pandas as pd

from .annotate import (
    ANTISENSE,
    CDS,
    FEATURES,
    INTERGENIC,
    LOCATIONS,
    SENSE,
    UTR,
)

GROUPS = (
    ("which site this row is about, and the key back to the alignment table", 5),
    ("the gene it falls in; empty when it falls in none", 12),
    ("the genes of other types it overlaps", 2),
    ("the promoters it sits in, and the nearest one", 7),
    ("the protospacer itself", 2),
)
"""The column groups, in order, and how many columns each has."""

COLUMNS = {
    # which site, and the key back to the alignment table
    "guide_id": (
        "The sequence that was aligned, ending in the PAM pattern. Join on it"
        " and the four site columns for the alias, accession and spacer, which"
        " are the same on every row of a guide and so are not repeated here."
    ),
    "chr": "Contig of this site.",
    "start": "First base of the protospacer, 0-based.",
    "end": "One past its last base, so the protospacer is [start, end).",
    "strand": (
        "Strand the protospacer reads 5' to 3' on. With the three columns"
        " above it says which of a multi-mapping guide's sites this row is."
    ),
    # the gene
    "ensg": "Ensembl gene id, without its version. Empty when intergenic.",
    "gene_symbol": "The gene's name.",
    "enst": (
        "The gene's reference transcript, whose exons and CDS this row is"
        " called against: MANE Select, else Ensembl canonical, else the"
        " longest CDS. Empty unless the protospacer overlaps one of its exons."
    ),
    "location": (
        "genic when the protospacer overlaps a gene of an annotated type,"
        " intergenic when it overlaps none. Independent of feature and"
        " promoter: a guide can be all three at once."
    ),
    "feature": (
        "exon or intron, from whether the protospacer overlaps an exon of the"
        " reference transcript. Empty on an intergenic row."
    ),
    "exonic_any_transcript": (
        "TRUE when any transcript of this gene has an exon here, so a row"
        " called intron is another transcript's exon."
    ),
    "region": "CDS or UTR for an exonic row, which is what decides a frameshift.",
    "exon_number": "Which exon of the reference transcript, counting from its 5' end.",
    "cds_position": (
        "How many coding bases precede the target site in the reference"
        " transcript,"
        " so a small number is an early truncation."
    ),
    "cds_fraction": "cds_position as a fraction of the whole CDS, 0 to 1.",
    "distance_to_gene_start": (
        "Bases from the target site to this gene's annotated 5' end, signed in"
        " the gene's orientation: negative upstream of it, positive inside."
        " Near zero is the canonical promoter; a large positive distance in a"
        " promoter window is an alternative one."
    ),
    "sense": (
        "sense when the protospacer reads on the gene's own strand,"
        " antisense otherwise."
    ),
    # the genes of other types
    "other_gene_types": (
        "Semicolon-separated types of the genes the protospacer overlaps that"
        " are not of an annotated type, such as lncRNA. Empty when none. A"
        " guide here sits where a read cannot be assigned unambiguously."
    ),
    "other_genes": "Semicolon-separated Ensembl ids of those genes.",
    # the promoters
    "promoter": (
        "TRUE when the protospacer falls within the promoter window of any"
        " protein-coding TSS, whichever gene it belongs to."
    ),
    "distance_to_tss": (
        "Bases from the target site (end - 3 on +, start + 3 on -) to the"
        " nearest protein-coding TSS, signed in that gene's orientation:"
        " negative upstream, positive downstream."
    ),
    "nearest_tss_gene": "The gene whose TSS that is; not always this row's gene.",
    "nearest_tss_ensg": "That gene's Ensembl id.",
    "promoter_enst": (
        "Semicolon-separated transcripts whose promoter window the"
        " protospacer overlaps. Empty when the TSS file carries no transcript"
        " id, which a BED of six or seven columns does not."
    ),
    "promoter_ensg": "Semicolon-separated gene ids of those transcripts.",
    "promoter_gene": "Semicolon-separated gene symbols of those transcripts.",
    # the protospacer
    "gc_percent": "The protospacer's GC content, as a percentage.",
    "max_polyt_run": (
        "Longest run of T in the protospacer. Four or more terminates Pol III,"
        " so the guide expresses poorly from a U6 promoter."
    ),
}
"""Every column of guide_gene_annotation.tsv, and what it holds."""

TABLE = "guide_gene_annotation.tsv"
"""The one results table."""


def annotation_columns() -> list[str]:
    """The columns of the annotation table, in order."""
    return list(COLUMNS)


def write_annotation_tsv(path: str | os.PathLike[str], annotated: pd.DataFrame) -> None:
    """Write the annotation table.

    Parameters
    ----------
    path
        The file to write.
    annotated
        What `annotate_sites` returned. Columns it does not carry are written
        empty, so the table always has the same shape.
    """
    frame = annotated.reindex(columns=annotation_columns())
    for column in ("promoter", "exonic_any_transcript"):
        frame[column] = frame[column].map({True: "TRUE", False: "FALSE"})
    for column in (
        "exon_number",
        "cds_position",
        "distance_to_tss",
        "distance_to_gene_start",
        "max_polyt_run",
    ):
        frame[column] = _whole_numbers(frame[column])
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def write_columns_yaml(path: str | os.PathLike[str]) -> None:
    """Write guide_gene_annotation_columns.yaml: what every column holds.

    Parameters
    ----------
    path
        The file to write.

    Notes
    -----
    Generated from `COLUMNS`, so it cannot drift from the table itself.
    """
    lines = ["# What each column of this folder's output holds.", "", f"{TABLE}:"]
    columns = iter(COLUMNS.items())
    for heading, count in GROUPS:
        lines.append(f"  # {heading}")
        for _ in range(count):
            column, description = next(columns)
            lines.append(f'  {column}: "{_one_line(description)}"')
    lines.append("")
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def counts(annotated: pd.DataFrame) -> dict:
    """Count the rows by location, region, promoter and sense.

    Parameters
    ----------
    annotated
        What `annotate_sites` returned.

    Returns
    -------
    dict
        One section per grouping, plus the totals. `location` sums to the
        number of rows, which is more than the number of sites when a site
        falls in two genes.
    """
    location = Counter(annotated["location"])
    region = Counter(annotated.loc[annotated["region"] != "", "region"])
    sense = Counter(annotated.loc[annotated["sense"] != "", "sense"])
    return {
        "rows": len(annotated),
        "sites": int(
            annotated[["chr", "start", "end", "strand"]].drop_duplicates().shape[0]
        ),
        "guides": int(annotated["guide_id"].nunique())
        if "guide_id" in annotated
        else 0,
        "genes": int(annotated["ensg"].nunique()),
        "location": {name: int(location.get(name, 0)) for name in LOCATIONS},
        "feature": {
            name: int(Counter(annotated["feature"]).get(name, 0)) for name in FEATURES
        },
        "region": {name: int(region.get(name, 0)) for name in (CDS, UTR)},
        "promoter": int(annotated["promoter"].sum()),
        "sense": {name: int(sense.get(name, 0)) for name in (SENSE, ANTISENSE)},
        "intergenic": int((annotated["location"] == INTERGENIC).sum()),
    }


def _whole_numbers(values: pd.Series) -> pd.Series:
    """Write a column of whole numbers without a trailing .0."""
    return values.astype("Float64").astype("Int64")


def _one_line(text: str) -> str:
    """One line, with any quote escaped, so the YAML stays valid."""
    return " ".join(text.split()).replace('"', '\\"')

"""The annotation table, and the data dictionary that describes it."""

from __future__ import annotations

import os
from collections import Counter
from pathlib import Path

import pandas as pd

from .annotate import ANTISENSE, CDS, INTERGENIC, LOCATIONS, SENSE, UTR

GROUPS = (
    ("the key back to guides_alignment_metadata.tsv, which holds the rest", 5),
    ("the gene it falls in; empty when it falls in none", 10),
    ("the nearest protein-coding promoter", 4),
    ("the protospacer itself", 2),
)
"""The column groups, in order, and how many columns each has."""

COLUMNS = {
    # the key back to the alignment table
    "guide_id": (
        "The sequence that was aligned, ending in the PAM pattern. Join on it"
        " and the four site columns to pick up the alias, accession, spacer"
        " and the rest; a guide with two sites needs all five to be unambiguous."
    ),
    "chr": "Contig of this site.",
    "start": "First base of the protospacer, 0-based.",
    "end": "One past its last base, so the protospacer is [start, end).",
    "strand": "Strand the protospacer reads 5' to 3' on.",
    # the gene
    "ensg": "Ensembl gene id, without its version. Empty when intergenic.",
    "gene_symbol": "The gene's name.",
    "enst": (
        "The gene's reference transcript, whose exons and CDS this row is"
        " called against: MANE Select, else Ensembl canonical, else the"
        " longest CDS. Empty unless the protospacer overlaps one of its exons."
    ),
    "location": (
        "exon, intron or intergenic, from whether the protospacer overlaps an"
        " exon of the reference transcript. Independent of promoter."
    ),
    "exonic_any_transcript": (
        "TRUE when any transcript of this gene has an exon here, so a row"
        " called intron is another transcript's exon."
    ),
    "region": "CDS or UTR for an exonic row, which is what decides a frameshift.",
    "exon_number": "Which exon of the reference transcript, counting from its 5' end.",
    "cds_position": (
        "How many coding bases precede the cut in the reference transcript,"
        " so a small number is an early truncation."
    ),
    "cds_fraction": "cds_position as a fraction of the whole CDS, 0 to 1.",
    "sense": (
        "sense when the protospacer reads on the gene's own strand,"
        " antisense otherwise."
    ),
    # the nearest promoter
    "promoter": (
        "TRUE when the protospacer falls within the promoter window of any"
        " protein-coding TSS, whichever gene it belongs to."
    ),
    "distance_to_tss": (
        "Bases from the cut site (end - 3 on +, start + 3 on -) to the"
        " nearest protein-coding TSS, signed in that gene's orientation:"
        " negative upstream, positive downstream."
    ),
    "nearest_tss_gene": "The gene whose TSS that is; not always this row's gene.",
    "nearest_tss_ensg": "That gene's Ensembl id.",
    # the protospacer
    "gc_percent": "The protospacer's GC content, as a percentage.",
    "max_polyt_run": (
        "Longest run of T in the protospacer. Four or more terminates Pol III,"
        " so the guide expresses poorly from a U6 promoter."
    ),
}
"""Every column of guides_gene_annotation.tsv, and what it holds."""

TABLE = "guides_gene_annotation.tsv"
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
    for column in ("exon_number", "cds_position", "distance_to_tss", "max_polyt_run"):
        frame[column] = _whole_numbers(frame[column])
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def write_columns_yaml(path: str | os.PathLike[str]) -> None:
    """Write columns.yaml: what every column of the table holds.

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

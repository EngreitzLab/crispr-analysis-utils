"""The elements a screen targets, and the aligned sites tested against them.

Every interval here is 0-based half-open, the way `element_mapping` reads and
writes one, so a test asserts the constant it reads.

The elements:

==== ====== ========= ========= ============================================
id   contig span      category  note
==== ====== ========= ========= ============================================
E1   chr1   1000-1100 enhancer  overlaps E2 in 1050-1100
E2   chr1   1050-1200 promoter
E3   chr2   2000-2100 enhancer  alone on its contig, so its edges are
                                unambiguous
==== ====== ========= ========= ============================================
"""

import pandas as pd

from crispr_analysis_utils.element_mapping import SITE_COLUMNS

LENGTH = 20
"""A protospacer's length, so a site is named by its start alone."""

ELEMENTS = (
    ("chr1", 1000, 1100, "E1", "enhancer"),
    ("chr1", 1050, 1200, "E2", "promoter"),
    ("chr2", 2000, 2100, "E3", "enhancer"),
)
"""The screen's elements, as the fields of a BED5."""


def bed_lines(rows=ELEMENTS, columns=5):
    """`rows` as BED lines, keeping the first `columns` fields of each."""
    return ["\t".join(str(field) for field in row[:columns]) for row in rows]


def write_bed(path, rows=ELEMENTS, *, columns=5, header=None):
    """Write the elements as a BED and return the path."""
    lines = bed_lines(rows, columns)
    if header is not None:
        lines.insert(0, header)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def site(guide_id, chrom, start, strand="+", *, length=LENGTH, pam_matches=True):
    """One row of guides_alignment_metadata.tsv, as `guide_alignment` writes it."""
    return {
        "guide_id": guide_id,
        "chr": chrom,
        "start": start,
        "end": start + length,
        "strand": strand,
        "pam_matches": pam_matches,
    }


def unmapped(guide_id, *, pam_matches=True):
    """A library entry that matched nowhere: a row with no coordinates."""
    return {
        "guide_id": guide_id,
        "chr": "",
        "start": pd.NA,
        "end": pd.NA,
        "strand": "",
        "pam_matches": pam_matches,
    }


def sites_frame(*rows):
    """The rows as a guides_alignment_metadata table."""
    return pd.DataFrame(list(rows))


def write_alignment(path, rows, *, as_text=True):
    """The same table on disk, with the flag spelled as a written one has it."""
    frame = sites_frame(*rows)
    if as_text:
        frame["pam_matches"] = frame["pam_matches"].map({True: "TRUE", False: "FALSE"})
    frame.to_csv(path, sep="\t", index=False)
    return path


def annotation_frame(*pairs):
    """guides_gene_annotation.tsv: the site columns and a distance for each."""
    frame = sites_frame(*(row for row, _ in pairs))[list(SITE_COLUMNS)].copy()
    frame["distance_to_tss"] = [distance for _, distance in pairs]
    return frame

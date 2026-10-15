"""The elements an IGVF export names, and the library aimed at them.

Every interval here is 0-based half-open, as is every coordinate
`igvf_export` reads or writes. The elements sit among the genes of
`guide_annotation_fixtures`, so what the genome says about each one is
worked out by hand from the intervals there.

==== ====== =========== ============== ========================================
id   contig span        the design     the genome, at WINDOW
==== ====== =========== ============== ========================================
P1   chr1   950-1050    promoter       promoter: it holds AAA's TSS at 1000
A2   chr1   1520-1540   enhancer       exon: A2's exon, in the reference intron
I1   chr1   1700-1800   enhancer       intron: inside AAA, in no exon
X1   chr1   2850-2950   promoter       exon: 1549 bases from the nearest TSS
D1   chr1   20000-20100 enhancer       distal element
Z1   chr1   30000-30100 promoter       distal element: nowhere near a TSS
V1   chr2   2950-3050   e2g_validated  enhancer, though it holds CCC's TSS
O2   chr2   20000-20200 enhancer       distal element; overlaps O1 in 20100-20200
O1   chr2   20100-20400 enhancer       distal element
==== ====== =========== ============== ========================================

A guide is named by what its element is, so a failure names the rule it
broke.
"""

import pandas as pd

from crispr_analysis_utils.element_mapping import map_to_elements

WINDOW = 100
"""A narrow promoter window: every exon of the fixture genes lies within the
2000 of `PROMOTER_WINDOW`, so at the default none of them could be an exon."""

LENGTH = 20
"""A protospacer's length, so a site is named by its start alone."""

ELEMENTS = (
    ("chr1", 950, 1050, "P1", "promoter"),
    ("chr1", 1520, 1540, "A2", "enhancer"),
    ("chr1", 1700, 1800, "I1", "enhancer"),
    ("chr1", 2850, 2950, "X1", "promoter"),
    ("chr1", 20000, 20100, "D1", "enhancer"),
    ("chr1", 30000, 30100, "Z1", "promoter"),
    ("chr2", 2950, 3050, "V1", "e2g_validated"),
    ("chr2", 20000, 20200, "O2", "enhancer"),
    ("chr2", 20100, 20400, "O1", "enhancer"),
)
"""The screen's elements, as the fields of a BED5."""

GUIDES = (
    # alias, contig, start, strand, the PAM the genome holds there
    ("promoter", "chr1", 980, "+", "TGG"),
    ("other_exon", "chr1", 1525, "+", "AGG"),
    ("intron", "chr1", 1710, "+", "AGG"),
    ("exon", "chr1", 2860, "+", "CGG"),
    ("distal", "chr1", 20010, "+", "TGG"),
    ("mislabel", "chr1", 30010, "-", "TGG"),
    ("validated", "chr2", 2960, "-", "GGG"),
    # Inside O2 and O1 both: 20 bases of O2 against 10 of O1.
    ("overlap", "chr2", 20090, "+", "TGG"),
    # Inside both again, 20 bases of each, so only the id can decide.
    ("tie", "chr2", 20100, "+", "TGG"),
    # Maps, and in no element.
    ("safe", "chr1", 50000, "+", "TGG"),
    # In the library and nowhere in the genome, so it has no site at all.
    ("non_targeting", "", None, "", ""),
)
"""Every library entry, in the order the table reports them."""

ALIASES = tuple(row[0] for row in GUIDES)
"""Each entry's name, which is the `guide_id` of the exported table."""

GENE_PAIRS = (
    # A bidirectional promoter: both genes were measured, not inferred.
    ("P1", "ENSG00000000001"),
    ("P1", "ENSG00000000002"),
    ("O2", "ENSG00000000003"),
    ("O2", "ENSG00000000005"),
    ("V1", "ENSG00000000003"),
)
"""Measured element-to-gene pairs, as a ``--genes`` table holds them."""


def spacer(alias):
    """A distinct 20-mer for one entry, so no two share a `guide_id`."""
    index = ALIASES.index(alias)
    tail = "".join("ACGT"[index // 4**place % 4] for place in reversed(range(3)))
    return "ACGTACGTACGTACGTA" + tail


def sequence_id(alias):
    """An entry's `guide_id`: what was aligned, ending in the PAM pattern."""
    return spacer(alias) + "NGG"


def alignment_frame(guides=GUIDES):
    """guide_genome_annotation.tsv, one row per entry.

    Every guide here matches at most once, so the per-site columns do not
    repeat.
    """
    rows = []
    for alias, chrom, start, strand, pam in guides:
        rows.append(
            {
                "guide_id": sequence_id(alias),
                "alias": alias,
                "spacer": spacer(alias),
                "chr": chrom,
                "start": pd.NA if start is None else start,
                "end": pd.NA if start is None else start + LENGTH,
                "strand": strand,
                "pam_matches": bool(pam),
                "mapped_sequence": spacer(alias) + pam if pam else "",
            }
        )
    return pd.DataFrame(rows)


def sites_frame(guides=GUIDES):
    """The entries that matched somewhere, as `element_mapping` takes them."""
    frame = alignment_frame(guides)
    return frame[frame["chr"] != ""].reset_index(drop=True)


def mapping_frame(elements, guides=GUIDES):
    """guide_element_mapping.tsv, made by the step that writes it."""
    return map_to_elements(
        sites_frame(guides),
        elements,
        guides=pd.Series([sequence_id(alias) for alias, *_ in guides]),
    )


def genes_map(pairs=GENE_PAIRS):
    """The measured pairs as `classify_elements` takes them."""
    measured = {}
    for element_id, ensg in pairs:
        measured.setdefault(element_id, []).append(ensg)
    return measured


def write_alignment(path, guides=GUIDES):
    """Write the alignment table and return the path."""
    alignment_frame(guides).to_csv(path, sep="\t", index=False)
    return path


def write_genes(path, pairs=GENE_PAIRS, columns=("element_id", "ensembl_id")):
    """Write the measured pairs as a ``--genes`` table and return the path."""
    pd.DataFrame(list(pairs), columns=list(columns)).to_csv(path, sep="\t", index=False)
    return path

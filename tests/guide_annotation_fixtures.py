"""A hand-built GENCODE-shaped GTF, and the sites annotated against it.

Every interval here is 0-based half-open, the way pyranges reads it back, so
a test asserts the constant it reads. `write_gtf` is the one place that turns
a start into the 1-based inclusive coordinate a GTF holds.

The genes:

===== ====== ====== =========== ==========================================
gene  contig strand span        reference transcript
===== ====== ====== =========== ==========================================
AAA   chr1   +      1000-3000   A1, by MANE_Select
BBB   chr1   -      2500-4500   B1, by Ensembl_canonical over B2's longer CDS
CCC   chr2   -      1000-3000   C2, by the longest CDS
EEE   chr2   +      5000-6000   E1, by the lowest id of two equal CDS
DDD   chr1   +      7000-8000   lncRNA, so no store keeps it by default
===== ====== ====== =========== ==========================================

AAA and BBB overlap in 2500-3000, so a site there falls in two genes.
"""

import pandas as pd

FIELD = "\t"

SPACER = "ACGTACGTACGTACGTACGT"
"""A 20-mer: 50.0 percent GC, longest run of T is 1."""

GENES = (
    {
        "ensg": "ENSG00000000001",
        "symbol": "AAA",
        "gene_type": "protein_coding",
        "chrom": "chr1",
        "strand": "+",
        "span": (1000, 3000),
        "transcripts": (
            {
                # MANE_Select sits between two other tags: a reader that keeps
                # one tag per row loses it.
                "enst": "ENST00000000011",
                "span": (1000, 3000),
                "tags": ("basic", "MANE_Select", "CCDS"),
                "exons": ((1000, 1200), (2000, 2200), (2800, 3000)),
                "cds": ((1100, 1200), (2000, 2200), (2800, 2900)),
                "utr": ((1000, 1100), (2900, 3000)),
            },
            {
                # Canonical, and its second exon falls in A1's intron.
                "enst": "ENST00000000012",
                "span": (1000, 3000),
                "tags": ("basic", "Ensembl_canonical"),
                "exons": ((1000, 1200), (1500, 1600), (2800, 3000)),
                "cds": ((1100, 1200), (1500, 1600)),
                "utr": ((1000, 1100), (2800, 3000)),
            },
        ),
    },
    {
        "ensg": "ENSG00000000002",
        "symbol": "BBB",
        "gene_type": "protein_coding",
        "chrom": "chr1",
        "strand": "-",
        "span": (2500, 4500),
        "transcripts": (
            {
                "enst": "ENST00000000021",
                "span": (2500, 4500),
                "tags": ("basic", "Ensembl_canonical"),
                "exons": ((2500, 2700), (3500, 3600), (4300, 4500)),
                "cds": ((2600, 2650), (3500, 3600)),
                "utr": ((2500, 2600), (2650, 2700), (4300, 4500)),
            },
            {
                # 400 coding bases against B1's 150, and still not the
                # reference: canonical outranks length.
                "enst": "ENST00000000022",
                "span": (2500, 4500),
                "tags": ("basic",),
                "exons": ((2500, 2700), (3000, 3400), (4300, 4500)),
                "cds": ((3000, 3400),),
                "utr": ((2500, 2700), (4300, 4500)),
            },
        ),
    },
    {
        "ensg": "ENSG00000000003",
        "symbol": "CCC",
        "gene_type": "protein_coding",
        "chrom": "chr2",
        "strand": "-",
        "span": (1000, 3000),
        "transcripts": (
            {
                "enst": "ENST00000000031",
                "span": (1000, 3000),
                "tags": ("basic",),
                "exons": ((1000, 1500), (2500, 3000)),
                "cds": ((2600, 2700),),
                "utr": ((1000, 1500), (2500, 2600), (2700, 3000)),
            },
            {
                # 300 coding bases, the longest of the three, and no tag.
                "enst": "ENST00000000032",
                "span": (1000, 3000),
                "tags": ("basic",),
                "exons": ((1000, 1500), (2000, 2200), (2500, 3000)),
                "cds": ((2600, 2900),),
                "utr": ((1000, 1500), (2000, 2200), (2500, 2600), (2900, 3000)),
            },
            {
                "enst": "ENST00000000033",
                "span": (1000, 3000),
                "tags": ("basic",),
                "exons": ((1000, 1500), (2500, 3000)),
                "cds": ((1200, 1300),),
                "utr": ((1000, 1200), (1300, 1500), (2500, 3000)),
            },
        ),
    },
    {
        "ensg": "ENSG00000000005",
        "symbol": "EEE",
        "gene_type": "protein_coding",
        "chrom": "chr2",
        "strand": "+",
        "span": (5000, 6000),
        "transcripts": (
            {
                # Equal CDS lengths, so the lower id breaks the tie.
                "enst": "ENST00000000051",
                "span": (5000, 5500),
                "tags": ("basic",),
                "exons": ((5000, 5500),),
                "cds": ((5100, 5200),),
                "utr": ((5000, 5100), (5200, 5500)),
            },
            {
                "enst": "ENST00000000052",
                "span": (5500, 6000),
                "tags": ("basic",),
                "exons": ((5500, 6000),),
                "cds": ((5600, 5700),),
                "utr": ((5500, 5600), (5700, 6000)),
            },
        ),
    },
    {
        "ensg": "ENSG00000000004",
        "symbol": "DDD",
        "gene_type": "lncRNA",
        "chrom": "chr1",
        "strand": "+",
        "span": (7000, 8000),
        "transcripts": (
            {
                "enst": "ENST00000000041",
                "span": (7000, 8000),
                "tags": ("basic",),
                "exons": ((7000, 8000),),
                "cds": (),
                "utr": ((7000, 8000),),
            },
        ),
    },
)

REFERENCE = {
    "ENSG00000000001": "ENST00000000011",
    "ENSG00000000002": "ENST00000000021",
    "ENSG00000000003": "ENST00000000032",
    "ENSG00000000005": "ENST00000000051",
}
"""Each gene's reference transcript, worked out by hand from `GENES`."""

TSS = {
    "ENSG00000000001": 1000,
    "ENSG00000000002": 4499,
    "ENSG00000000003": 2999,
    "ENSG00000000005": 5000,
}
"""The reference transcript's start, which is what `promoters_from_store` takes."""


def _attributes(pairs):
    """GTF attributes, as `key "value";` in the order given."""
    return " ".join(f'{key} "{value}";' for key, value in pairs)


def _exon_numbers(gene, transcript):
    """Each exon's number, counting from the transcript's own 5' end."""
    exons = list(transcript["exons"])
    order = exons if gene["strand"] == "+" else exons[::-1]
    return {block: number for number, block in enumerate(order, start=1)}


def _containing(block, numbers):
    """The number of the exon a CDS or UTR block sits in."""
    for exon, number in numbers.items():
        if exon[0] <= block[0] and block[1] <= exon[1]:
            return number
    raise AssertionError(f"{block} is in no exon of {sorted(numbers)}")


def _row(gene, start, end, feature, attributes, frame="."):
    """One GTF line, turning the 0-based half-open interval into 1-based inclusive."""
    return FIELD.join(
        [
            gene["chrom"],
            "TEST",
            feature,
            str(start + 1),
            str(end),
            ".",
            gene["strand"],
            frame,
            attributes,
        ]
    )


def gtf_lines(genes=GENES):
    """Every GTF line of `genes`, in GENCODE's gene-then-transcript order."""
    lines = []
    for gene in genes:
        shared = [
            ("gene_id", f"{gene['ensg']}.5"),
            ("gene_type", gene["gene_type"]),
            ("gene_name", gene["symbol"]),
        ]
        lines.append(_row(gene, *gene["span"], "gene", _attributes(shared)))
        for transcript in gene["transcripts"]:
            numbers = _exon_numbers(gene, transcript)
            named = [*shared, ("transcript_id", f"{transcript['enst']}.3")]
            tags = [("tag", tag) for tag in transcript["tags"]]
            lines.append(
                _row(
                    gene,
                    *transcript["span"],
                    "transcript",
                    _attributes(named + tags),
                )
            )
            for feature, blocks in (
                ("exon", transcript["exons"]),
                ("CDS", transcript["cds"]),
                ("UTR", transcript["utr"]),
            ):
                for block in blocks:
                    number = (
                        numbers[block]
                        if feature == "exon"
                        else _containing(block, numbers)
                    )
                    lines.append(
                        _row(
                            gene,
                            *block,
                            feature,
                            _attributes([*named, ("exon_number", number)]),
                            frame="0" if feature == "CDS" else ".",
                        )
                    )
    return lines


def write_gtf(path, genes=GENES):
    """Write `genes` as a GTF and return the path."""
    text = "\n".join(["#!genome-build TEST", *gtf_lines(genes)]) + "\n"
    path.write_text(text, encoding="utf-8")
    return path


def write_genes_only_gtf(path):
    """A GTF with gene rows and nothing else, which cannot say exon or intron."""
    lines = [line for line in gtf_lines() if FIELD.join(["TEST", "gene"]) in line]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def site(guide_id, chrom, start, strand, *, spacer=SPACER, length=20, pam_matches=True):
    """One aligned site, with the cut where `guide_alignment` puts it."""
    end = start + length
    return {
        "guide_id": guide_id,
        "chr": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "cut_site": end - 3 if strand == "+" else start + 3,
        "spacer": spacer,
        "pam_matches": pam_matches,
    }


def sites_frame(*rows):
    """The sites as `annotate_sites` takes them."""
    return pd.DataFrame(list(rows))

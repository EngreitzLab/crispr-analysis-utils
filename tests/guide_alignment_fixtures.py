"""Hand-built reference and SAM files for the guide_alignment tests."""

import pysam

from crispr_analysis_utils.guide_alignment.iupac import reverse_complement

SPACER_A = "ACGTTGCAAGCTTCGATCGA"
SPACER_C = "CTTGACCGATAGCATGCAAT"
SPACER_G = "GATCCGTAGCTAGGCTTACG"


class Contig:
    """A reference contig built piece by piece; `add` returns where a piece starts."""

    def __init__(self):
        self.sequence = ""

    def add(self, piece):
        start = len(self.sequence)
        self.sequence += piece
        return start

    def add_minus(self, piece):
        """Add a piece given in guide orientation for a minus-strand site."""
        return self.add(reverse_complement(piece))


def write_fasta(path, contigs, width=60):
    """Write {name: sequence} as a FASTA with `width`-base lines."""
    with open(path, "w", encoding="ascii") as handle:
        for name, sequence in contigs.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")
    return path


def hit(name, chrom, pos, strand, length, *, nm=None, secondary=False, cigar=None):
    """A mapped record: 0-based `pos`, flag 0 or 16, plus 256 if secondary."""
    return {
        "name": name,
        "chrom": chrom,
        "pos": pos,
        "strand": strand,
        "length": length,
        "nm": nm,
        "secondary": secondary,
        "cigar": cigar,
    }


def unmapped(name):
    """An unmapped record (flag 4)."""
    return {"name": name, "chrom": None}


def write_sam(path, references, records):
    """Write records against {name: length} @SQ lines; BAM if `path` ends in .bam."""
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": name, "LN": length} for name, length in references.items()],
    }
    mode = "wb" if str(path).endswith(".bam") else "w"
    with pysam.AlignmentFile(str(path), mode, header=header) as out:
        for record in records:
            segment = pysam.AlignedSegment(out.header)
            segment.query_name = record["name"]
            if record["chrom"] is None:
                segment.flag = 4
                segment.reference_id = -1
                segment.reference_start = -1
            else:
                segment.flag = (16 if record["strand"] == "-" else 0) | (
                    256 if record["secondary"] else 0
                )
                segment.reference_id = out.get_tid(record["chrom"])
                segment.reference_start = record["pos"]
                segment.mapping_quality = 0 if record["secondary"] else 60
                segment.cigarstring = record["cigar"] or f"{record['length']}M"
                if record["nm"] is not None:
                    segment.set_tag("NM", record["nm"])
            out.write(segment)
    return path

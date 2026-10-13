"""Hand-built reference pieces for the guide_alignment tests."""

from crispr_analysis_utils.guide_alignment.iupac import reverse_complement

SPACER_A = "ACGTTGCAAGCTTCGATCGA"
SPACER_C = "CTTGACCGATAGCATGCAAT"
SPACER_G = "GATCCGTAGCTAGGCTTACG"
SPACER_T = "TTGCCAGTACGATCAGGTCA"


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

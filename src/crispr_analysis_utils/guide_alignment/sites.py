"""One genomic site of a guide, as the scan reports it."""

from __future__ import annotations

from dataclasses import dataclass

NO_PAM_CLASS = "none"
"""The `Site.pam_class` of a site whose genomic PAM matches no pattern."""


@dataclass(frozen=True)
class Site:
    """One genomic site of a guide, as verified against the reference.

    Attributes
    ----------
    guide_id
        The guide's id.
    chr
        The contig.
    start
        The first base of the protospacer, 0-based.
    end
        One past the last base of the protospacer: the protospacer is
        ``[start, end)``.
    strand
        ``+`` or ``-``: the strand the protospacer reads 5' to 3' on.
    pam
        The genomic PAM, uppercased, in guide orientation.
    pam_class
        The first pattern the genomic PAM matches, trying the PAM and then any
        alternatives in their order, or ``none``. A genomic ``N`` matches only
        an ``N`` of the pattern.
    n_mismatches
        The mismatches between the spacer and the genomic protospacer. A
        genomic ``N`` counts as one; the added G never counts.
    mismatch_positions
        The 1-based spacer positions of those mismatches, counted from the
        guide's 5' end.
    mapped_sequence
        What was aligned here, 5' to 3': the base 5' of the protospacer where
        the genome has that same base, the protospacer, then the genomic PAM.
    genomic_protospacer
        The reference bases of the protospacer, uppercased, in guide
        orientation.
    genomic_5p_base
        The genomic base under the sequence's base 5' of the protospacer, in
        guide orientation; None when it has none.
    """

    guide_id: str
    chr: str
    start: int
    end: int
    strand: str
    pam: str
    pam_class: str
    n_mismatches: int
    mismatch_positions: tuple[int, ...]
    mapped_sequence: str
    genomic_protospacer: str
    genomic_5p_base: str | None

    @property
    def paired_5p(self) -> bool | None:
        """Whether the guide's base 5' of the spacer pairs with the genome.

        None when there is nothing to ask: a guide with no 5' base, or a
        match so close to a contig end that the base under it was never
        read. A 5' base that pairs is recognised like any other, which is
        why the BED's interval covers it and the declared protospacer does
        not move.
        """
        if self.genomic_5p_base is None:
            return None
        return len(self.mapped_sequence) > len(self.genomic_protospacer) + len(self.pam)

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
        The guide's own bases that pair with the genome here, 5' to 3': the
        spacer, keeping an added 5' G only when the genomic base under it is
        also G. The PAM is genomic, so it is not part of it.
    genomic_protospacer
        The reference bases of the protospacer, uppercased, in guide
        orientation.
    added_g_base
        For a guide with bases 5' of its spacer, the genomic bases under
        them, in guide orientation; None otherwise.
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
    added_g_base: str | None

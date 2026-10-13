"""The genome scan: every exact site of every guide, read from the reference.

A guide is designed, not sequenced, so its own site is an exact match of the
spacer. One pass of the compiled core finds every such match on both strands
and reads the bases around it, which answers all three questions at once:
where the spacer matches, what sits 3' of it (the PAM) and what sits 5' of it
(the base an added G lands on).

Orientation falls out of the same pass, because both directions go into the
same table: each guide contributes its spacer and the spacer of its reverse
complement. An entry is `reverse_complement` only when that second sequence
matches somewhere with a PAM and the first does not. Reading `CCN` 5' of a
match is not the test -- it happens by chance one time in sixteen -- and for a
guide carrying an added G the reversed spacer is not even the reverse
complement of the spacer, since the trim moves to the other end.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Collection, Sequence
from dataclasses import dataclass, field

from .iupac import pattern_matches, reverse_complement, validate_pams
from .library import Guide
from .orientation import (
    AMBIGUOUS,
    DISTINCT,
    FORWARD,
    NO_MATCH,
    NO_PAM,
    REVERSE_COMPLEMENT,
)
from .sites import NO_PAM_CLASS, Site

logger = logging.getLogger(__name__)


@dataclass
class ScanResult:
    """Every guide's exact sites, its orientation, and what was set aside.

    Attributes
    ----------
    sites
        The sites of each guide row, in reference order.
    verdicts
        Each guide row's orientation verdict.
    evidence
        For a ``reverse_complement`` row, whether the reversed reading was a
        distinct sequence or only the spacer's reverse complement; None for
        every other row.
    n_matches
        Exact spacer matches found, before any were set aside.
    n_outside_contigs
        Matches dropped because `contigs` did not list their contig. Always 0:
        the scan is told which contigs to read.
    n_off_contig_end
        Matches whose PAM would run off the end of their contig.
    n_no_pam
        Matches kept with a genomic PAM that matches no pattern.
    n_sites
        Matches kept.
    """

    sites: dict[int, list[Site]] = field(default_factory=dict)
    verdicts: dict[int, str] = field(default_factory=dict)
    evidence: dict[int, str | None] = field(default_factory=dict)
    n_matches: int = 0
    n_outside_contigs: int = 0
    n_off_contig_end: int = 0
    n_no_pam: int = 0
    n_sites: int = 0

    def counts(self) -> dict[str, int]:
        """The counts, for qc.tsv and the step's metadata."""
        return {
            "n_matches": self.n_matches,
            "n_outside_contigs": self.n_outside_contigs,
            "n_off_contig_end": self.n_off_contig_end,
            "n_no_pam": self.n_no_pam,
            "n_sites": self.n_sites,
        }


def scan_guides(
    guides: Sequence[Guide],
    reference: str | os.PathLike[str],
    *,
    pam: str = "NGG",
    contigs: Collection[str] | None = None,
) -> ScanResult:
    """Find every exact site of every guide in `reference`.

    Parameters
    ----------
    guides
        The library, as returned by `read_guides`.
    reference
        The reference FASTA. It is read as text, so it must not be compressed.
    pam
        The PAM pattern, in IUPAC codes, 3' of the protospacer.
    contigs
        If given, read only these contigs.

    Returns
    -------
    ScanResult
        The sites of each guide row and its orientation verdict.

    Raises
    ------
    ValueError
        If `pam` is not IUPAC codes, or the guides' spacers differ in length.
    OSError
        If `reference` cannot be read.

    Notes
    -----
    Each distinct spacer is scanned once and its matches are given to every
    guide that uses it, so two entries differing only in an added G cost one
    search between them.
    """
    from .. import _core

    patterns = validate_pams([pam])
    primary = patterns[0]
    if not guides:
        return ScanResult()

    rows_of: dict[str, list[int]] = {}
    flipped_rows: dict[str, list[int]] = {}
    for row, guide in enumerate(guides):
        rows_of.setdefault(guide.spacer, []).append(row)
        # Not reverse_complement(spacer): the trim moves to the other end, so
        # a guide with an added G reverses to a different 20-mer entirely.
        flipped = reverse_complement(guide.aligned_sequence)[-len(guide.spacer) :]
        flipped_rows.setdefault(flipped, []).append(row)
    spacers = list(rows_of)
    reversed_spacers = [s for s in flipped_rows if s not in rows_of]
    lengths = {len(spacer) for spacer in spacers}
    if len(lengths) > 1:
        raise ValueError(
            "Every guide must have the same spacer length; found "
            f"{', '.join(str(length) for length in sorted(lengths))}."
        )

    # Both directions in one pass: the second block of spacers is only there
    # to answer the orientation question, and never yields a site.
    queries = spacers + reversed_spacers
    upstream_width = max((guide.n_extra for guide in guides), default=0)
    listed = None if contigs is None else list(contigs)
    matches = _core.scan_fasta(
        os.fspath(reference), queries, primary, upstream_width, listed
    )

    result = ScanResult()
    reversed_has_pam: set[int] = set()
    for index, chrom, start, strand, genomic_pam, upstream, matched in matches:
        if index >= len(spacers):
            if matched:
                reversed_has_pam.update(flipped_rows[queries[index]])
            continue
        result.n_matches += 1
        if len(genomic_pam) != len(primary):
            result.n_off_contig_end += 1
            continue
        if not matched:
            result.n_no_pam += 1
        result.n_sites += 1
        for row in rows_of[spacers[index]]:
            guide = guides[row]
            result.sites.setdefault(row, []).append(
                _site(guide, chrom, start, strand, genomic_pam, upstream, primary)
            )

    for row, guide in enumerate(guides):
        verdict = _verdict(result.sites.get(row, ()), row in reversed_has_pam)
        result.verdicts[row] = verdict
        result.evidence[row] = (
            None
            if verdict != REVERSE_COMPLEMENT
            else (DISTINCT if guide.n_extra else AMBIGUOUS)
        )
    logger.info("Genome scan: %s", result.counts())
    return result


def _site(
    guide: Guide,
    chrom: str,
    start: int,
    strand: str,
    genomic_pam: str,
    upstream: str,
    pam: str,
) -> Site:
    """One match, as a `Site`. An exact match has no mismatches by definition."""
    extra = guide.aligned_sequence[: guide.n_extra]
    genomic_extra = upstream[-guide.n_extra :] if guide.n_extra else None
    return Site(
        guide_id=guide.guide_id,
        chr=chrom,
        start=start,
        end=start + len(guide.spacer),
        strand=strand,
        pam=genomic_pam,
        pam_class=pam if pattern_matches(genomic_pam, pam) else NO_PAM_CLASS,
        n_mismatches=0,
        mismatch_positions=(),
        mapped_sequence=(_kept(extra, genomic_extra) + guide.spacer + genomic_pam),
        genomic_protospacer=guide.spacer,
        genomic_5p_base=genomic_extra,
    )


def _kept(extra: str, genomic: str | None) -> str:
    """The guide's bases 5' of the spacer that pair with the genome."""
    if not extra or genomic is None:
        return ""
    return "".join(
        base for base, found in zip(extra, genomic, strict=True) if base == found
    )


def _verdict(sites: Sequence[Site], reversed_has_pam: bool) -> str:
    """One guide's orientation, from both directions' matches."""
    if any(site.pam_class != NO_PAM_CLASS for site in sites):
        return FORWARD
    if reversed_has_pam:
        return REVERSE_COMPLEMENT
    return NO_PAM if sites else NO_MATCH

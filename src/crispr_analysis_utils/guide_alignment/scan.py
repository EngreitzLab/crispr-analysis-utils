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

from .iupac import reverse_complement, validate_pams
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
    matches
        How many places each guide row's protospacer occurs, counted whether
        or not the site was kept. Under `pam_only` it is the only record of
        the sites that were not.
    n_matches
        Exact spacer matches found, before any were set aside.
    n_outside_contigs
        Matches dropped because `contigs` did not list their contig. Always 0:
        the scan is told which contigs to read.
    n_off_contig_end
        Matches whose PAM would run off the end of their contig.
    n_no_pam
        Matches kept with a genomic PAM that matches no pattern. Under
        `pam_only` these are counted and not written, so the table holds
        ``n_sites - n_no_pam`` rows.
    n_sites
        Matches kept.
    """

    sites: dict[int, list[Site]] = field(default_factory=dict)
    verdicts: dict[int, str] = field(default_factory=dict)
    evidence: dict[int, str | None] = field(default_factory=dict)
    matches: dict[int, int] = field(default_factory=dict)
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
    pam_only: bool = False,
    tolerate_5p: int = 0,
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
    pam_only
        Keep only the matches whose genomic PAM matches `pam`. The core
        drops the rest as it finds them, so they never cross into Python
        and nothing holds them. Every match is still counted, in
        `ScanResult.matches` and `n_matches`: a short spacer matches
        everywhere, and the places it cannot be cut at are the bulk of
        them, but how many there were is worth knowing.
    tolerate_5p
        How many bases at the protospacer's PAM-distal end need not pair with
        the genome. The exact search is that much shorter, while the
        protospacer and its coordinates keep their declared length: 1 finds
        the site of a guide whose designer wrote a G at position 1 where the
        genome holds another base. What the genome holds there is reported in
        `Site.genomic_spacer_5p` and counted nowhere, as the added G's base is.

    Returns
    -------
    ScanResult
        The sites of each guide row and its orientation verdict.

    Raises
    ------
    ValueError
        If `pam` is not IUPAC codes, the guides' spacers differ in length, or
        `tolerate_5p` is negative or leaves no bases to search.
    OSError
        If `reference` cannot be read.

    Notes
    -----
    Each distinct search core is scanned once and its matches are given to
    every guide that uses it, so two entries differing only in an added G cost
    one search between them.

    `pam_only` shrinks what is held and what is written, not what is
    searched: every match is found and counted either way. The counters
    come from the scan rather than from the matches returned, which is what
    lets them report what they never handed back.
    """
    from .. import _core

    patterns = validate_pams([pam])
    primary = patterns[0]
    if not guides:
        return ScanResult()
    if (
        isinstance(tolerate_5p, bool)
        or not isinstance(tolerate_5p, int)
        or tolerate_5p < 0
    ):
        raise ValueError("tolerate_5p must be a non-negative integer.")

    # The core is what is searched: the spacer, less the PAM-distal bases the
    # run tolerates. The spacer itself still decides the coordinates.
    rows_of: dict[str, list[int]] = {}
    flipped_rows: dict[str, list[int]] = {}
    for row, guide in enumerate(guides):
        rows_of.setdefault(guide.spacer[tolerate_5p:], []).append(row)
        # Not reverse_complement(spacer): the trim moves to the other end, so
        # a guide with an added G reverses to a different 20-mer entirely.
        flipped = reverse_complement(guide.aligned_sequence)[-len(guide.spacer) :]
        flipped_rows.setdefault(flipped[tolerate_5p:], []).append(row)
    cores = list(rows_of)
    reversed_cores = [s for s in flipped_rows if s not in rows_of]
    lengths = {len(core) for core in cores}
    if len(lengths) > 1:
        raise ValueError(
            "Every guide must have the same spacer length; found "
            f"{', '.join(str(length + tolerate_5p) for length in sorted(lengths))}."
        )
    if min(lengths) < 1:
        raise ValueError(
            f"tolerate_5p={tolerate_5p} leaves no bases of a "
            f"{tolerate_5p}-base protospacer to search."
        )

    # Both directions in one pass: the second block of cores is only there
    # to answer the orientation question, and never yields a site.
    queries = cores + reversed_cores
    upstream_width = tolerate_5p + max((guide.n_extra for guide in guides), default=0)
    listed = None if contigs is None else list(contigs)
    matches, counts = _core.scan_fasta(
        os.fspath(reference), queries, primary, upstream_width, listed, pam_only
    )

    result = ScanResult()
    result.matches = dict.fromkeys(range(len(guides)), 0)
    # The counts come from the scan itself rather than from the matches, so
    # they say the same thing whether or not `pam_only` returned them all.
    # The reversed readings are not sites and are not counted here.
    for index, (found, off_end, no_pam) in enumerate(counts[: len(cores)]):
        result.n_matches += found
        result.n_off_contig_end += off_end
        result.n_no_pam += no_pam
        result.n_sites += found - off_end
        for row in rows_of[cores[index]]:
            result.matches[row] += found - off_end

    reversed_has_pam: set[int] = set()
    for index, chrom, start, strand, genomic_pam, upstream, matched in matches:
        # One guide's reversed reading can be another's own core, which is
        # then searched once, in the first block: gather it from either.
        if matched:
            reversed_has_pam.update(flipped_rows.get(queries[index], ()))
        if index >= len(cores):
            continue
        if len(genomic_pam) != len(primary):
            continue
        for row in rows_of[cores[index]]:
            guide = guides[row]
            site = _site(
                guide,
                chrom,
                start,
                strand,
                genomic_pam,
                upstream,
                primary if matched else NO_PAM_CLASS,
                tolerate_5p,
            )
            if site is not None:
                result.sites.setdefault(row, []).append(site)

    for row, guide in enumerate(guides):
        verdict = _verdict(
            result.sites.get(row, ()), row in reversed_has_pam, result.matches[row]
        )
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
    pam_class: str,
    tolerate_5p: int = 0,
) -> Site | None:
    """One match, as a `Site`, or None when no whole protospacer fits here.

    `pam_class` is the core's verdict on `genomic_pam`, which is not decided
    again here.

    The core matched exactly, so the genome can differ only at the
    `tolerate_5p` bases the search left out, which it supplies here. `start`
    is the core's left end, and the tolerated bases lie 5' of it: on ``+``
    that is the lower coordinate, on ``-`` the higher one.
    """
    needed = tolerate_5p + guide.n_extra
    flank = upstream[len(upstream) - needed :] if needed else ""
    if len(flank) < tolerate_5p:
        # The protospacer would run off the contig, so there is no site here.
        return None
    genomic_5p = flank[len(flank) - tolerate_5p :] if tolerate_5p else ""
    before = flank[: len(flank) - tolerate_5p]
    extra = guide.aligned_sequence[: guide.n_extra]
    # The core returns an empty flank when it would run off the contig, so a
    # guide with a 5' base matching at base 0 has nothing under that base.
    genomic_extra = (
        before[-guide.n_extra :]
        if guide.n_extra and len(before) >= guide.n_extra
        else None
    )
    genomic_protospacer = genomic_5p + guide.spacer[tolerate_5p:]
    # The core matched exactly, so only the tolerated positions can differ.
    mismatches = tuple(
        position + 1
        for position, (want, found) in enumerate(
            zip(guide.spacer[:tolerate_5p], genomic_5p, strict=True)
        )
        if want != found
    )
    first = start - tolerate_5p if strand == "+" else start
    return Site(
        guide_id=guide.guide_id,
        chr=chrom,
        start=first,
        end=first + len(guide.spacer),
        strand=strand,
        pam=genomic_pam,
        pam_class=pam_class,
        n_mismatches=len(mismatches),
        mismatch_positions=mismatches,
        mapped_sequence=(
            _kept(extra, genomic_extra) + genomic_protospacer + genomic_pam
        ),
        genomic_protospacer=genomic_protospacer,
        genomic_5p_base=genomic_extra,
        genomic_spacer_5p=genomic_5p,
    )


def _kept(extra: str, genomic: str | None) -> str:
    """The guide's bases 5' of the spacer that pair with the genome."""
    if not extra or genomic is None:
        return ""
    return "".join(
        base for base, found in zip(extra, genomic, strict=True) if base == found
    )


def _verdict(sites: Sequence[Site], reversed_has_pam: bool, n_matches: int) -> str:
    """One guide's orientation, from both directions' matches.

    The last branch counts matches rather than kept sites: under `pam_only`
    a guide that matched everywhere without a PAM keeps none of them, and
    reading the empty list would call it unmatched.
    """
    if any(site.pam_class != NO_PAM_CLASS for site in sites):
        return FORWARD
    if reversed_has_pam:
        return REVERSE_COMPLEMENT
    return NO_PAM if n_matches else NO_MATCH

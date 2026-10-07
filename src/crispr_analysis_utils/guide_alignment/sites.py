"""Pass 2 of the guide alignment: the genomic sites of each guide.

Every alignment of a ``[G] + spacer + PAM`` read (see `write_site_reads`) is a
candidate site. Its protospacer, PAM and leading-G base are read from the
reference in guide orientation, and its mismatches recomputed against the
spacer; the aligner's own ``NM`` is only compared. Coordinates are 0-based
and half-open on the reference's forward strand, and the protospacer excludes
the PAM and the G.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from ._alignment import (
    as_paths,
    check_cigar,
    check_contigs,
    check_header,
    check_max_mismatches,
    fetch_oriented,
    guide_row,
    is_unplaced,
    open_alignments,
    open_reference,
    strand_of,
)
from .iupac import pattern_matches, validate_pams
from .library import Guide

if TYPE_CHECKING:
    import pysam

logger = logging.getLogger(__name__)

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
        The first pattern the genomic PAM matches, trying the primary PAM and
        then the alternatives in their order, or ``none``. A genomic ``N``
        matches only an ``N`` of the pattern.
    n_mismatches
        The mismatches between the spacer and the genomic protospacer. A
        genomic ``N`` counts as one; the added G never counts.
    mismatch_positions
        The 1-based spacer positions of those mismatches, counted from the
        guide's 5' end.
    genomic_protospacer
        The reference bases of the protospacer, uppercased, in guide
        orientation.
    leading_g_base
        For a guide with an added G, the genomic base under that G, in guide
        orientation; None otherwise.
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
    genomic_protospacer: str
    leading_g_base: str | None


@dataclass
class SiteResult:
    """The kept sites of each guide, and counts of what was set aside.

    Attributes
    ----------
    sites
        The kept sites, by guide row (the guide's position in the library).
        Only rows with at least one kept site appear, and only guides without
        `Guide.duplicate_of` have rows. Each list is sorted by contig in
        reference order, then by start, then by strand.
    n_records
        SAM records read.
    n_unmapped
        Records of unmapped reads, skipped.
    n_outside_contigs
        Alignments on a contig that the `contigs` filter leaves out, skipped.
    n_off_contig_end
        Alignments whose G, protospacer or PAM runs past the end of the
        contig, skipped.
    n_over_max_mismatches
        Alignments with more spacer mismatches than `max_mismatches`, dropped.
    n_repeated
        Alignments at a site already kept for the same guide, merged into it.
    n_nm_disagreements
        Alignments whose ``NM`` tag differs from the recomputed total: the
        spacer mismatches, plus each PAM position whose code is not A, C, G
        or T or differs from the genomic base, plus a mismatched added G.
        Alignments without an ``NM`` tag are not compared.
    """

    sites: dict[int, list[Site]] = field(default_factory=dict)
    n_records: int = 0
    n_unmapped: int = 0
    n_outside_contigs: int = 0
    n_off_contig_end: int = 0
    n_over_max_mismatches: int = 0
    n_repeated: int = 0
    n_nm_disagreements: int = 0

    def counts(self) -> dict[str, int]:
        """The counters by attribute name, plus ``n_sites``, the kept sites."""
        return {
            "n_records": self.n_records,
            "n_unmapped": self.n_unmapped,
            "n_outside_contigs": self.n_outside_contigs,
            "n_off_contig_end": self.n_off_contig_end,
            "n_over_max_mismatches": self.n_over_max_mismatches,
            "n_repeated": self.n_repeated,
            "n_nm_disagreements": self.n_nm_disagreements,
            "n_sites": sum(len(sites) for sites in self.sites.values()),
        }


def read_sites(
    sam: str | os.PathLike[str] | Iterable[str | os.PathLike[str]],
    guides: Sequence[Guide],
    reference: str | os.PathLike[str],
    *,
    pam: str = "NGG",
    alt_pams: Iterable[str] = ("NAG", "NGA"),
    max_mismatches: int = 3,
    contigs: Collection[str] | None = None,
    fasta_index: str | os.PathLike[str] | None = None,
) -> SiteResult:
    """Read the site-search alignments and verify each site against the reference.

    The records are streamed; only kept sites are held in memory.

    Parameters
    ----------
    sam
        One SAM or BAM file, or several (such as the alignments of the reads
        with and without an added G). A SAM may be gzip-compressed. Read names
        must be ``<row>:<pam>`` as `write_site_reads` writes them.
    guides
        The library the reads were written from, as returned by `read_guides`.
    reference
        The reference FASTA the aligner's index was built from.
    pam
        The primary PAM pattern, in IUPAC codes.
    alt_pams
        The alternative PAM patterns, in IUPAC codes, tried in order after
        `pam` to classify a genomic PAM.
    max_mismatches
        Keep a site only if its spacer has at most this many mismatches.
    contigs
        If given, keep only alignments on these contigs.
    fasta_index
        The FASTA index to use, built there if it does not exist. By default
        the index is ``<reference>.fai``, which pysam builds next to the FASTA
        if it does not exist.

    Returns
    -------
    SiteResult
        The kept sites by guide row, with the counts of skipped, dropped and
        merged alignments. A site found by the reads of several PAMs, or
        reported twice, is kept once: sites are merged on contig, start and
        strand.

    Raises
    ------
    ValueError
        If the PAM patterns fail `validate_pams`, `max_mismatches` is negative,
        `contigs` names a contig the reference lacks, a file's @SQ names and
        lengths differ from the reference's, a read name is not
        ``<row>:<pam>`` for a non-duplicate row and one of the patterns, or a
        mapped record is not a gapless, unclipped alignment of the whole read.
    """
    patterns = validate_pams([pam, *alt_pams])
    check_max_mismatches(max_mismatches)
    result = SiteResult()
    found: dict[int, dict[tuple[str, int, str], Site]] = {}
    with open_reference(reference, fasta_index) as fasta:
        allowed = check_contigs(contigs, fasta)
        lengths = dict(zip(fasta.references, fasta.lengths, strict=True))
        for path in as_paths(sam):
            with open_alignments(path) as alignments:
                check_header(alignments, fasta, path, reference)
                for record in alignments.fetch(until_eof=True):
                    _read_record(
                        record,
                        guides=guides,
                        patterns=patterns,
                        max_mismatches=max_mismatches,
                        allowed=allowed,
                        lengths=lengths,
                        fasta=fasta,
                        result=result,
                        found=found,
                    )
        order = {name: index for index, name in enumerate(fasta.references)}
    for row in sorted(found):
        result.sites[row] = sorted(
            found[row].values(),
            key=lambda site: (order[site.chr], site.start, site.strand),
        )
    logger.info("Site search: %s", result.counts())
    if result.n_nm_disagreements:
        logger.warning(
            "%d alignments have an NM tag that differs from the recomputed mismatches.",
            result.n_nm_disagreements,
        )
    return result


def _read_record(
    record: pysam.AlignedSegment,
    *,
    guides: Sequence[Guide],
    patterns: tuple[str, ...],
    max_mismatches: int,
    allowed: frozenset[str] | None,
    lengths: dict[str, int],
    fasta: pysam.FastaFile,
    result: SiteResult,
    found: dict[int, dict[tuple[str, int, str], Site]],
) -> None:
    """Verify one record, then keep its site or count why it was set aside."""
    result.n_records += 1
    if is_unplaced(record):
        result.n_unmapped += 1
        return
    row, read_pam = _parse_read_name(record.query_name, guides, patterns)
    chrom = record.reference_name
    if allowed is not None and chrom not in allowed:
        result.n_outside_contigs += 1
        return

    guide = guides[row]
    n_g = int(guide.leading_g)
    n_spacer = len(guide.spacer)
    n_pam = len(read_pam)
    check_cigar(record, n_g + n_spacer + n_pam)
    position = record.reference_start
    stop = position + n_g + n_spacer + n_pam
    if position < 0 or stop > lengths[chrom]:
        result.n_off_contig_end += 1
        return

    strand = strand_of(record)
    bases = fetch_oriented(fasta, chrom, position, stop, strand, lengths[chrom])
    leading_g_base = bases[:n_g] or None
    protospacer = bases[n_g : n_g + n_spacer]
    genomic_pam = bases[n_g + n_spacer :]
    positions = tuple(
        index
        for index, (genomic, base) in enumerate(
            zip(protospacer, guide.spacer, strict=True), start=1
        )
        if genomic != base
    )

    if record.has_tag("NM"):
        recomputed = (
            len(positions)
            + _pam_errors(genomic_pam, read_pam)
            + int(leading_g_base not in (None, "G"))
        )
        if record.get_tag("NM") != recomputed:
            result.n_nm_disagreements += 1

    if len(positions) > max_mismatches:
        result.n_over_max_mismatches += 1
        return

    start = position + n_g if strand == "+" else position + n_pam
    key = (chrom, start, strand)
    guide_sites = found.setdefault(row, {})
    if key in guide_sites:
        result.n_repeated += 1
        return
    guide_sites[key] = Site(
        guide_id=guide.guide_id,
        chr=chrom,
        start=start,
        end=start + n_spacer,
        strand=strand,
        pam=genomic_pam,
        pam_class=next(
            (pattern for pattern in patterns if pattern_matches(genomic_pam, pattern)),
            NO_PAM_CLASS,
        ),
        n_mismatches=len(positions),
        mismatch_positions=positions,
        genomic_protospacer=protospacer,
        leading_g_base=leading_g_base,
    )


def _parse_read_name(
    name: str, guides: Sequence[Guide], patterns: tuple[str, ...]
) -> tuple[int, str]:
    """The guide row and the PAM pattern of a ``<row>:<pam>`` read name."""
    row_text, colon, read_pam = name.partition(":")
    if not colon or read_pam not in patterns:
        raise ValueError(
            f"Read name {name!r} is not <row>:<pam> with one of the PAM patterns "
            f"{', '.join(patterns)}."
        )
    return guide_row(row_text, guides, name), read_pam


def _pam_errors(genomic_pam: str, pattern: str) -> int:
    """PAM mismatches as the aligner counts them: every non-ACGT code is one."""
    return sum(
        code not in "ACGT" or base != code
        for base, code in zip(genomic_pam, pattern, strict=True)
    )

"""Pass 1 of the guide alignment: is the library given 5' to 3'?

The bare spacers (see `write_spacer_reads`) are aligned with no mismatch
allowed, once as given and once reverse-complemented. For each perfect hit the
flanks are read from the reference in read orientation, and the verdict says
where a PAM was found: ``pam_3p`` is what a correctly oriented guide gives.

`resolve_orientation` then gives each guide its verdict from both directions:
an entry is ``reverse_complement`` only when its reverse complement has a
perfect hit with a PAM and the sequence as given does not. A table can mix
orientations, so this is per guide; `reversed_rows` names the ones to flip.
Coordinates are 0-based and half-open on the reference's forward strand.
"""

from __future__ import annotations

import logging
import os
from collections import Counter
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from ._alignment import (
    as_paths,
    check_cigar,
    check_contigs,
    check_header,
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

PAM_3P = "pam_3p"
"""One direction's result: a perfect hit has a PAM on its 3' side."""
NO_PAM = "no_pam"
"""The guide has perfect hits, and none has a PAM 3' of it."""
NO_PERFECT_HIT = "no_perfect_hit"
"""The guide has no perfect hit."""
NOT_CHECKED = "not_checked"
"""The orientation check did not run."""

FORWARD = "forward"
"""The sequence as given has a perfect hit with a PAM 3' of it."""
REVERSE_COMPLEMENT = "reverse_complement"
"""Only the reverse complement does: this entry is written 3' to 5'."""
VERDICTS = (FORWARD, REVERSE_COMPLEMENT, NO_PAM, NO_PERFECT_HIT, NOT_CHECKED)
"""Every orientation verdict, in report order."""


@dataclass(frozen=True)
class PerfectHit:
    """A genomic copy of a spacer with no mismatch, and its flanks.

    Attributes
    ----------
    chr
        The contig.
    start
        The first base of the protospacer, 0-based.
    end
        One past the last base of the protospacer.
    strand
        ``+`` or ``-``: the strand the protospacer reads 5' to 3' on.
    flank_5p
        The PAM-length stretch of reference 5' of the protospacer, uppercased,
        in guide orientation. Shorter when the contig ends first.
    flank_3p
        The same on the protospacer's 3' side.
    """

    chr: str
    start: int
    end: int
    strand: str
    flank_5p: str
    flank_3p: str


@dataclass
class OrientationResult:
    """The perfect hits and the orientation verdict of each guide.

    Attributes
    ----------
    hits
        The perfect hits, by guide row (the guide's position in the library).
        Only rows with a perfect hit appear. Each list is sorted by contig in
        reference order, then by start, then by strand.
    verdicts
        The verdict of every guide without `Guide.duplicate_of`, by row: one
        of ``ok``, ``reversed``, ``no_pam`` and ``no_perfect_hit``.
    n_records
        SAM records read.
    n_unmapped
        Records of unmapped reads, skipped.
    n_imperfect
        Alignments whose reference bases differ from the spacer, set aside.
    n_outside_contigs
        Alignments on a contig that the `contigs` filter leaves out, skipped.
    n_off_contig_end
        Alignments that run past the end of the contig, skipped.
    n_repeated
        Alignments at a hit already recorded for the same guide.
    """

    hits: dict[int, list[PerfectHit]] = field(default_factory=dict)
    verdicts: dict[int, str] = field(default_factory=dict)
    n_records: int = 0
    n_unmapped: int = 0
    n_imperfect: int = 0
    n_outside_contigs: int = 0
    n_off_contig_end: int = 0
    n_repeated: int = 0

    def counts(self) -> dict[str, int]:
        """The counters by attribute name, plus ``n_hits``, the perfect hits."""
        return {
            "n_records": self.n_records,
            "n_unmapped": self.n_unmapped,
            "n_imperfect": self.n_imperfect,
            "n_outside_contigs": self.n_outside_contigs,
            "n_off_contig_end": self.n_off_contig_end,
            "n_repeated": self.n_repeated,
            "n_hits": sum(len(hits) for hits in self.hits.values()),
        }


def check_orientation(
    sam: str | os.PathLike[str] | Iterable[str | os.PathLike[str]],
    guides: Sequence[Guide],
    reference: str | os.PathLike[str],
    *,
    pam: str = "NGG",
    alt_pams: Iterable[str] = ("NAG", "NGA"),
    contigs: Collection[str] | None = None,
    fasta_index: str | os.PathLike[str] | None = None,
) -> OrientationResult:
    """Read the orientation-check alignments and give each guide a verdict.

    Parameters
    ----------
    sam
        One SAM or BAM file, or several; a SAM may be gzip-compressed. Read
        names must be guide rows, as `write_spacer_reads` writes them.
    guides
        The library the reads were written from, as returned by `read_guides`.
    reference
        The reference FASTA the aligner's index was built from.
    pam
        The primary PAM pattern, in IUPAC codes.
    alt_pams
        The alternative PAM patterns, in IUPAC codes.
    contigs
        If given, use only alignments on these contigs.
    fasta_index
        The FASTA index to use, built there if it does not exist. By default
        the index is ``<reference>.fai``, which pysam builds next to the FASTA
        if it does not exist.

    Returns
    -------
    OrientationResult
        The perfect hits and the verdict of each guide. A guide is ``ok`` if
        any perfect hit has a PAM (`pam` or one of `alt_pams`) on its 3' side;
        else ``reversed`` if any has the reverse complement of a PAM on its 5'
        side; else ``no_pam``. A guide without a perfect hit is
        ``no_perfect_hit``. A flank cut short by a contig end matches no PAM.

    Raises
    ------
    ValueError
        If the PAM patterns fail `validate_pams`, `contigs` names a contig the
        reference lacks, a file's @SQ names and lengths differ from the
        reference's, a read name is not the row of a non-duplicate guide, or
        a mapped record is not a gapless, unclipped alignment of the whole
        spacer.
    """
    patterns = validate_pams([pam, *alt_pams])
    flank_length = len(patterns[0])
    result = OrientationResult()
    found: dict[int, dict[tuple[str, int, str], PerfectHit]] = {}
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
                        flank_length=flank_length,
                        allowed=allowed,
                        lengths=lengths,
                        fasta=fasta,
                        result=result,
                        found=found,
                    )
        order = {name: index for index, name in enumerate(fasta.references)}
    for row in sorted(found):
        result.hits[row] = sorted(
            found[row].values(),
            key=lambda hit: (order[hit.chr], hit.start, hit.strand),
        )

    for row, guide in enumerate(guides):
        if guide.duplicate_of is None:
            result.verdicts[row] = _verdict(result.hits.get(row, []), patterns)
    logger.info(
        "Orientation check: %s; verdicts %s",
        result.counts(),
        dict(Counter(result.verdicts.values())),
    )
    return result


def resolve_orientation(
    forward: Mapping[int, str], reverse: Mapping[int, str]
) -> dict[int, str]:
    """Each guide's orientation, from its result in both directions.

    A guide is only called ``reverse_complement`` when its reverse complement
    has a perfect hit with a PAM and the sequence as given does not. Nothing
    is inferred from a flank: the other direction was aligned, so there is no
    need to guess at it.

    Parameters
    ----------
    forward
        Each guide row's result for the sequences as given, as in
        `OrientationResult.verdicts`.
    reverse
        The same for the reverse-complemented sequences.

    Returns
    -------
    dict
        The verdict of every row in either mapping: ``forward``,
        ``reverse_complement``, ``no_pam`` or ``no_perfect_hit``.
    """
    verdicts = {}
    for row in sorted(set(forward) | set(reverse)):
        as_given = forward.get(row, NO_PERFECT_HIT)
        flipped = reverse.get(row, NO_PERFECT_HIT)
        if as_given == PAM_3P:
            verdicts[row] = FORWARD
        elif flipped == PAM_3P:
            verdicts[row] = REVERSE_COMPLEMENT
        elif NO_PAM in (as_given, flipped):
            verdicts[row] = NO_PAM
        else:
            verdicts[row] = NO_PERFECT_HIT
    return verdicts


def reversed_rows(verdicts: Mapping[int, str]) -> list[int]:
    """The rows whose entry is written 3' to 5', in order.

    Parameters
    ----------
    verdicts
        The verdicts of `resolve_orientation`.

    Returns
    -------
    list of int
        Every row whose verdict is ``reverse_complement``.
    """
    return [
        row
        for row, verdict in sorted(verdicts.items())
        if verdict == REVERSE_COMPLEMENT
    ]


def _read_record(
    record: pysam.AlignedSegment,
    *,
    guides: Sequence[Guide],
    flank_length: int,
    allowed: frozenset[str] | None,
    lengths: dict[str, int],
    fasta: pysam.FastaFile,
    result: OrientationResult,
    found: dict[int, dict[tuple[str, int, str], PerfectHit]],
) -> None:
    """Verify one record, then record its hit or count why it was set aside."""
    result.n_records += 1
    if is_unplaced(record):
        result.n_unmapped += 1
        return
    row = guide_row(record.query_name, guides, record.query_name)
    chrom = record.reference_name
    if allowed is not None and chrom not in allowed:
        result.n_outside_contigs += 1
        return

    spacer = guides[row].spacer
    check_cigar(record, len(spacer))
    start = record.reference_start
    end = start + len(spacer)
    length = lengths[chrom]
    if start < 0 or end > length:
        result.n_off_contig_end += 1
        return
    strand = strand_of(record)
    if fetch_oriented(fasta, chrom, start, end, strand, length) != spacer:
        result.n_imperfect += 1
        return

    left = fetch_oriented(fasta, chrom, start - flank_length, start, strand, length)
    right = fetch_oriented(fasta, chrom, end, end + flank_length, strand, length)
    key = (chrom, start, strand)
    guide_hits = found.setdefault(row, {})
    if key in guide_hits:
        result.n_repeated += 1
        return
    guide_hits[key] = PerfectHit(
        chr=chrom,
        start=start,
        end=end,
        strand=strand,
        flank_5p=left if strand == "+" else right,
        flank_3p=right if strand == "+" else left,
    )


def _verdict(hits: Sequence[PerfectHit], patterns: Sequence[str]) -> str:
    """One direction's result for a guide, from its perfect hits."""
    if not hits:
        return NO_PERFECT_HIT
    if any(pattern_matches(hit.flank_3p, p) for hit in hits for p in patterns):
        return PAM_3P
    return NO_PAM

"""Shared helpers for reading the aligner's SAM output against the reference."""

from __future__ import annotations

import os
from collections.abc import Collection, Iterable, Iterator, Sequence
from contextlib import contextmanager
from typing import TYPE_CHECKING

from .iupac import reverse_complement
from .library import Guide, _listing

if TYPE_CHECKING:
    import pysam

_ALIGNED_OPS = frozenset({0, 7, 8})  # CIGAR M, = and X


def as_paths(
    paths: str | os.PathLike[str] | Iterable[str | os.PathLike[str]],
) -> list[str]:
    """One path, or several, as a list of strings."""
    if isinstance(paths, (str, os.PathLike)):
        return [os.fspath(paths)]
    return [os.fspath(path) for path in paths]


@contextmanager
def open_reference(
    reference: str | os.PathLike[str],
    fasta_index: str | os.PathLike[str] | None = None,
) -> Iterator[pysam.FastaFile]:
    """Open the reference FASTA, building `fasta_index` first if it is missing."""
    import pysam

    reference = os.fspath(reference)
    if fasta_index is None:
        fasta = pysam.FastaFile(reference)
    else:
        fasta_index = os.fspath(fasta_index)
        if not os.path.exists(fasta_index):
            pysam.faidx(reference, "--fai-idx", fasta_index)
        fasta = pysam.FastaFile(reference, filepath_index=fasta_index)
    try:
        yield fasta
    finally:
        fasta.close()


@contextmanager
def open_alignments(path: str) -> Iterator[pysam.AlignmentFile]:
    """Open a SAM (plain or gzip-compressed) or BAM file for streaming."""
    import pysam

    with pysam.AlignmentFile(path, "r", check_sq=False) as alignments:
        yield alignments


def check_header(
    alignments: pysam.AlignmentFile,
    fasta: pysam.FastaFile,
    path: str,
    reference: str | os.PathLike[str],
) -> None:
    """Raise ValueError unless the @SQ names and lengths equal the reference's."""
    sam = dict(zip(alignments.references, alignments.lengths, strict=True))
    if not sam:
        raise ValueError(
            f"{path} has no @SQ header lines to check against the reference."
        )
    genome = dict(zip(fasta.references, fasta.lengths, strict=True))
    if sam == genome:
        return
    problems = []
    only_sam = [name for name in sam if name not in genome]
    only_genome = [name for name in genome if name not in sam]
    differ = [
        f"{name} ({sam[name]} vs {genome[name]})"
        for name in sam
        if name in genome and sam[name] != genome[name]
    ]
    if only_sam:
        problems.append(f"only in the SAM: {_listing(only_sam)}")
    if only_genome:
        problems.append(f"only in the reference: {_listing(only_genome)}")
    if differ:
        problems.append(f"lengths differ: {_listing(differ)}")
    raise ValueError(
        f"The @SQ lines of {path} do not match the reference "
        f"{os.fspath(reference)}; {'; '.join(problems)}."
    )


def check_contigs(
    contigs: Collection[str] | None, fasta: pysam.FastaFile
) -> frozenset[str] | None:
    """The contig filter as a set, after checking every name is in the reference."""
    if contigs is None:
        return None
    if isinstance(contigs, str):
        raise TypeError("contigs must be a collection of contig names, not a string.")
    wanted = frozenset(contigs)
    if not wanted:
        raise ValueError("contigs must name at least one contig.")
    unknown = sorted(wanted - set(fasta.references))
    if unknown:
        raise ValueError(
            f"contigs names contigs that are not in the reference: {_listing(unknown)}."
        )
    return wanted


def check_max_mismatches(max_mismatches: int) -> None:
    """Raise ValueError unless `max_mismatches` is a non-negative integer."""
    if (
        isinstance(max_mismatches, bool)
        or not isinstance(max_mismatches, int)
        or max_mismatches < 0
    ):
        raise ValueError("max_mismatches must be a non-negative integer.")


def guide_row(text: str, guides: Sequence[Guide], name: str) -> int:
    """The guide row a read name carries, checked against the library."""
    if not (text.isascii() and text.isdigit()):
        raise ValueError(f"Read name {name!r} does not carry a guide row.")
    row = int(text)
    if row >= len(guides):
        raise ValueError(
            f"Read {name!r} names guide row {row}, but the library has "
            f"{len(guides)} guides."
        )
    if guides[row].duplicate_of is not None:
        raise ValueError(
            f"Read {name!r} names guide row {row}, a duplicate of "
            f"{guides[row].duplicate_of!r}: reads are written only for the first "
            "guide with each spacer."
        )
    return row


def is_unplaced(record: pysam.AlignedSegment) -> bool:
    """Whether a record is unmapped (flag 4) or names no reference."""
    return record.is_unmapped or record.reference_id < 0


def check_cigar(record: pysam.AlignedSegment, expected: int) -> None:
    """Raise ValueError unless the record aligns all `expected` bases, gapless."""
    cigar = record.cigartuples or []
    if (
        not cigar
        or any(op not in _ALIGNED_OPS for op, _ in cigar)
        or sum(length for _, length in cigar) != expected
    ):
        raise ValueError(
            f"Read {record.query_name!r} at {record.reference_name}:"
            f"{record.reference_start + 1} has CIGAR {record.cigarstring or '*'}; "
            f"expected one gapless, unclipped alignment of all {expected} bases "
            f"(such as {expected}M)."
        )


def strand_of(record: pysam.AlignedSegment) -> str:
    """``-`` for a reverse-strand record, ``+`` otherwise."""
    return "-" if record.is_reverse else "+"


def fetch_oriented(
    fasta: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    contig_length: int,
) -> str:
    """The reference bases of [start, end), uppercased, in guide orientation.

    The interval is first clipped to the contig, so a flank that runs past
    either end comes back shorter. Minus-strand bases are reverse-complemented.
    """
    start = max(start, 0)
    end = min(end, contig_length)
    if start >= end:
        return ""
    sequence = fasta.fetch(chrom, start, end).upper()
    return reverse_complement(sequence) if strand == "-" else sequence

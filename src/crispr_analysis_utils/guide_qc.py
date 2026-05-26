"""Guide-level QC helpers."""

from __future__ import annotations

from pathlib import Path
import gzip
import re
from collections import defaultdict

import pandas as pd

CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
MD_PATTERN = re.compile(r"(\d+|\^[A-Za-z]+|[A-Za-z])")
_DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

_CIGAR_MATCH = {0, 7, 8}  # M, =, X
_CIGAR_INS = 1            # I
_CIGAR_DEL = 2            # D
_CIGAR_SOFT = 4           # S

DEFAULT_PRIMARY_CONTIGS = None


def guides_to_fastq(
    guides_input: str | Path | pd.DataFrame,
    output_path: str | Path,
    *,
    quality_char: str = "I",
    pam: str = "",
    add_leading_g: bool = False,
    id_col: int | str = 0,
    sequence_col: int | str = 1,
) -> Path:
    """Write guide sequences into a synthetic FASTQ file.

    Parameters
    ----------
    guides_input
        Either:
        - Path to a tab-separated file
        - pandas DataFrame
    output_path
        FASTQ output path. Use `.gz` suffix for gzip-compressed output.
    quality_char
        Single character repeated to sequence length for quality scores.
        Defaults to `I` (Q = 40 in Phred+33 encoding).
    pam
        Optional PAM sequence appended to each guide sequence.
    add_leading_g
        If `True`, prepend `G` only when the sequence does not already start
        with `G`.
    id_col
        Guide id column in `guides_input` (integer position or column name).
    sequence_col
        Guide sequence column in `guides_input` (integer position or column
        name).

    Returns
    -------
    pathlib.Path
        Path to the written FASTQ file.
    """
    if isinstance(guides_input, pd.DataFrame):
        guides = guides_input
    else:
        guides = pd.read_csv(guides_input, sep="\t", header=None)

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    open_fn = gzip.open if path.suffix == ".gz" else open
    with open_fn(path, "wt", encoding="utf-8") as handle:
        for _, row in guides.iterrows():
            if isinstance(id_col, int):
                record_id = str(row.iloc[id_col])
            else:
                record_id = str(row[id_col])

            if isinstance(sequence_col, int):
                sequence = str(row.iloc[sequence_col]).strip().upper()
            else:
                sequence = str(row[sequence_col]).strip().upper()
            if add_leading_g and not sequence.startswith("G"):
                sequence = f"G{sequence}"
            sequence = f"{sequence}{pam}"
            quality = quality_char * len(sequence)
            handle.write(f"@{record_id}\n{sequence}\n+\n{quality}\n")

    return path


def _mismatch_positions_from_md(md_tag: str, cigartuples: list[tuple[int, int]]) -> set[int]:
    """Infer mismatch query positions from MD + CIGAR."""
    mismatches: set[int] = set()
    tokens = MD_PATTERN.findall(md_tag)
    token_i = 0
    qpos = 0

    for op, length in cigartuples:
        if op in (_CIGAR_MATCH - {8}):  # M and = (X handled directly from CIGAR)
            remaining = length
            while remaining > 0 and token_i < len(tokens):
                tok = tokens[token_i]
                if tok.isdigit():
                    n = int(tok)
                    step = min(n, remaining)
                    qpos += step
                    remaining -= step
                    if n > step:
                        tokens[token_i] = str(n - step)
                    else:
                        token_i += 1
                elif tok.startswith("^"):
                    token_i += 1
                else:
                    mismatches.add(qpos)
                    qpos += 1
                    remaining -= 1
                    token_i += 1
            qpos += remaining
        elif op == 8:  # X
            for i in range(length):
                mismatches.add(qpos + i)
            qpos += length
        elif op in (1, 4):  # I, S
            qpos += length
    return mismatches


def _evaluate_alignment_layout(
    cigartuples: list[tuple[int, int]],
    query_length: int,
    mismatch_positions: set[int],
    *,
    pam: str = "NGG",
    is_reverse_strand: bool = False,
    allow_leading_g_softclip: bool = True,
) -> tuple[bool, str]:
    """Evaluate whether alignment satisfies guide QC layout rules."""
    if query_length <= len(pam):
        return False, "query_too_short_for_pam"

    query_ops: list[int | None] = [None] * query_length
    qpos = 0
    deletion_query_positions: set[int] = set()

    for op, length in cigartuples:
        if op in _CIGAR_MATCH | {_CIGAR_INS, _CIGAR_SOFT}:
            for _ in range(length):
                if qpos < query_length:
                    query_ops[qpos] = op
                qpos += 1
        elif op == _CIGAR_DEL:
            deletion_query_positions.add(qpos)

    if is_reverse_strand:
        pam_start = 0
        pam_indices = range(0, len(pam))
    else:
        pam_start = query_length - len(pam)
        pam_indices = range(pam_start, query_length)

    if not all(query_ops[i] in _CIGAR_MATCH for i in pam_indices):
        return False, "pam_not_fully_aligned"

    mismatch_check_start = 1 if pam and pam[0].upper() == "N" else 0
    for i in range(mismatch_check_start, len(pam)):
        if (pam_start + i) in mismatch_positions:
            return False, "pam_gg_mismatch"

    if is_reverse_strand:
        spacer_start = len(pam)
        spacer_end = query_length - (1 if allow_leading_g_softclip else 0)
    else:
        spacer_start = 1 if allow_leading_g_softclip else 0
        spacer_end = pam_start

    for i in range(spacer_start, spacer_end):
        if query_ops[i] == _CIGAR_INS:
            return False, "spacer_insertion"

    for q in deletion_query_positions:
        if spacer_start <= q < spacer_end:
            return False, "spacer_deletion"

    for i, op in enumerate(query_ops):
        if op == _CIGAR_SOFT:
            allow_idx = 0 if not is_reverse_strand else (query_length - 1)
            if not (allow_leading_g_softclip and i == allow_idx):
                return False, "softclip_not_allowed"

    return True, "valid"


def _query_length_from_cigartuples(cigartuples: list[tuple[int, int]]) -> int:
    """Infer query length from CIGAR tuples.

    Counts query-consuming operations: M, I, S, =, X.
    """
    query_consuming = {0, 1, 4, 7, 8}
    return sum(length for op, length in cigartuples if op in query_consuming)


def _revcomp(sequence: str) -> str:
    """Reverse-complement DNA sequence."""
    return sequence.translate(_DNA_COMPLEMENT)[::-1]


def filter_guide_alignments(
    input_sam: str | Path,
    output_unique_sam: str | Path | None = None,
    output_multi_sam: str | Path | None = None,
    *,
    output_invalid_tsv: str | Path = "auto",
    output_valid_bed: str | Path = "auto",
    output_discarded_tsv: str | Path = "auto",
    output_unmapped_tsv: str | Path = "auto",
    output_guide_log_tsv: str | Path = "auto",
    output_summary_tsv: str | Path = "auto",
    alias_by_guide_id: dict[str, str] | None = None,
    pam: str = "NGG",
    allow_leading_g_softclip: bool = True,
    primary_contigs: set[str] | list[str] | tuple[str, ...] | None = None,
    chromsizes: str | Path | pd.DataFrame | None = None,
) -> dict[str, int]:
    """Filter SAM/BAM guide alignments and split valid unique vs multi-mapping.

    Rules:
    - must map to a primary contig
    - PAM must be aligned; mismatches/indels/soft clips not allowed in PAM
      except PAM's first base (`N`)
    - spacer cannot contain insertions or deletions
    - leading G soft-clipping is allowed by default
    """
    try:
        import pysam
    except ImportError as exc:
        raise ImportError("pysam is required for SAM/BAM filtering.") from exc

    primary = _resolve_allowed_contigs(primary_contigs=primary_contigs, chromsizes=chromsizes)
    input_path = Path(input_sam)
    mode = "rb" if input_path.suffix == ".bam" else "r"

    valid_by_guide: dict[str, list] = {}
    invalid_rows: list[tuple[str, str, str]] = []
    discarded_rows: list[tuple[str, str, int, int, int, str, str, int, int, str, str]] = []
    unmapped_rows: list[tuple[str, int]] = []
    valid_bed_rows: list[tuple[str, int, int, str, int, str, int, int, str]] = []
    guide_stats = defaultdict(
        lambda: {"n_aligned": 0, "n_valid": 0, "n_discarded": 0, "n_not_mapped": 0}
    )
    qname_to_sequence: dict[str, str] = {}

    with pysam.AlignmentFile(str(input_path), mode) as in_sam:
        header = in_sam.header
        for aln in in_sam.fetch(until_eof=True):
            qname = aln.query_name
            seq = aln.query_sequence
            if seq is not None and seq != "*":
                seq_u = seq.upper()
                qname_to_sequence[qname] = _revcomp(seq_u) if aln.is_reverse else seq_u
            recovered_seq = qname_to_sequence.get(qname)
            guide_id = recovered_seq if recovered_seq is not None else qname

            if aln.is_unmapped:
                invalid_rows.append((guide_id, ".", "unmapped"))
                unmapped_rows.append((guide_id, aln.flag))
                guide_stats[guide_id]["n_not_mapped"] = 1
                continue

            contig = aln.reference_name
            guide_stats[guide_id]["n_aligned"] += 1
            if primary is not None and contig not in primary:
                invalid_rows.append((guide_id, contig, "non_primary_contig"))
                guide_stats[guide_id]["n_discarded"] += 1
                discarded_rows.append(
                    (
                        guide_id,
                        contig,
                        int(aln.reference_start) + 1,
                        int(aln.flag),
                        int(aln.mapping_quality),
                        "-" if aln.is_reverse else "+",
                        aln.cigarstring or "*",
                        int(aln.get_tag("NM")) if aln.has_tag("NM") else -1,
                        int(aln.get_tag("AS")) if aln.has_tag("AS") else -1,
                        aln.get_tag("MD") if aln.has_tag("MD") else "",
                        "non_primary_contig",
                    )
                )
                continue

            md_tag = aln.get_tag("MD") if aln.has_tag("MD") else ""
            cigartuples = aln.cigartuples or []
            query_len = aln.query_length or _query_length_from_cigartuples(cigartuples)
            mismatch_positions = _mismatch_positions_from_md(md_tag, cigartuples)
            is_valid, reason = _evaluate_alignment_layout(
                cigartuples,
                query_len,
                mismatch_positions,
                pam=pam,
                is_reverse_strand=aln.is_reverse,
                allow_leading_g_softclip=allow_leading_g_softclip,
            )
            if is_valid:
                valid_by_guide.setdefault(guide_id, []).append(aln)
                guide_stats[guide_id]["n_valid"] += 1
                nm_tag = int(aln.get_tag("NM")) if aln.has_tag("NM") else -1
                as_tag = int(aln.get_tag("AS")) if aln.has_tag("AS") else -1
                bed_span = _protospacer_bed_span(
                    aln,
                    query_len=query_len,
                    pam_len=len(pam),
                    allow_leading_g_softclip=allow_leading_g_softclip,
                )
                if bed_span is not None:
                    valid_bed_rows.append(
                        (
                            contig,
                            bed_span[0],
                            bed_span[1],
                            guide_id,
                            int(aln.mapping_quality),
                            "-" if aln.is_reverse else "+",
                            nm_tag,
                            as_tag,
                            alias_by_guide_id.get(guide_id, guide_id)
                            if alias_by_guide_id is not None
                            else guide_id,
                        )
                    )
            else:
                if reason in {"pam_not_fully_aligned", "query_too_short_for_pam", "softclip_not_allowed"}:
                    normalized_reason = "discarded_tail_unaligned"
                elif reason == "pam_gg_mismatch":
                    normalized_reason = "discarded_tail_mismatch"
                elif reason in {"spacer_insertion", "spacer_deletion"}:
                    normalized_reason = "discarded_protospacer_indel"
                else:
                    normalized_reason = reason
                invalid_rows.append((guide_id, contig, reason))
                guide_stats[guide_id]["n_discarded"] += 1
                discarded_rows.append(
                    (
                        guide_id,
                        contig,
                        int(aln.reference_start) + 1,
                        int(aln.flag),
                        int(aln.mapping_quality),
                        "-" if aln.is_reverse else "+",
                        aln.cigarstring or "*",
                        int(aln.get_tag("NM")) if aln.has_tag("NM") else -1,
                        int(aln.get_tag("AS")) if aln.has_tag("AS") else -1,
                        aln.get_tag("MD") if aln.has_tag("MD") else "",
                        normalized_reason,
                    )
                )

    base_output_dir = Path(".")
    if output_unique_sam is not None:
        unique_path = Path(output_unique_sam)
        base_output_dir = unique_path.parent
    elif output_multi_sam is not None:
        multi_path = Path(output_multi_sam)
        base_output_dir = multi_path.parent

    if output_unique_sam is not None and output_multi_sam is not None:
        unique_path = Path(output_unique_sam)
        multi_path = Path(output_multi_sam)
        unique_path.parent.mkdir(parents=True, exist_ok=True)
        multi_path.parent.mkdir(parents=True, exist_ok=True)

        with (
            pysam.AlignmentFile(str(unique_path), "w", header=header) as unique_sam,
            pysam.AlignmentFile(str(multi_path), "w", header=header) as multi_sam,
        ):
            for alignments in valid_by_guide.values():
                target = unique_sam if len(alignments) == 1 else multi_sam
                for aln in alignments:
                    target.write(aln)

    if output_invalid_tsv == "auto":
        invalid_path = base_output_dir / "invalid_alignments.tsv"
    else:
        invalid_path = Path(output_invalid_tsv)
    invalid_path.parent.mkdir(parents=True, exist_ok=True)
    with invalid_path.open("w", encoding="utf-8") as handle:
        handle.write("guide_id\tcontig\treason\n")
        for row in invalid_rows:
            handle.write("\t".join(row) + "\n")

    if output_valid_bed == "auto":
        valid_bed_path = base_output_dir / "valid_alignments.bed"
    else:
        valid_bed_path = Path(output_valid_bed)
    valid_bed_path.parent.mkdir(parents=True, exist_ok=True)
    with valid_bed_path.open("w", encoding="utf-8") as handle:
        for row in valid_bed_rows:
            handle.write("\t".join(map(str, row)) + "\n")

    if output_discarded_tsv == "auto":
        discarded_path = base_output_dir / "discarded_alignments.tsv"
    else:
        discarded_path = Path(output_discarded_tsv)
    discarded_path.parent.mkdir(parents=True, exist_ok=True)
    with discarded_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "read_name\tchromosome\tpos1\tflag\tmapq\tstrand\tcigar\tNM\tAS\tMD\treason\n"
        )
        for row in discarded_rows:
            handle.write("\t".join(map(str, row)) + "\n")

    if output_unmapped_tsv == "auto":
        unmapped_path = base_output_dir / "unmapped.tsv"
    else:
        unmapped_path = Path(output_unmapped_tsv)
    unmapped_path.parent.mkdir(parents=True, exist_ok=True)
    with unmapped_path.open("w", encoding="utf-8") as handle:
        handle.write("read_name\tflag\n")
        for row in unmapped_rows:
            handle.write("\t".join(map(str, row)) + "\n")

    if output_guide_log_tsv == "auto":
        guide_log_path = base_output_dir / "guide_alignment_log.tsv"
    else:
        guide_log_path = Path(output_guide_log_tsv)
    guide_log_path.parent.mkdir(parents=True, exist_ok=True)
    with guide_log_path.open("w", encoding="utf-8") as handle:
        handle.write("guide_id\tn_aligned\tn_valid\tn_discarded\tn_not_mapped\n")
        for guide_id in sorted(guide_stats):
            s = guide_stats[guide_id]
            handle.write(
                f"{guide_id}\t{s['n_aligned']}\t{s['n_valid']}\t{s['n_discarded']}\t{s['n_not_mapped']}\n"
            )

    n_guides_unique_valid = 0
    n_guides_multi_valid = 0
    n_guides_aligned_none_valid = 0
    n_guides_unmapped = 0
    n_guides_one_valid_plus_invalid = 0

    for s in guide_stats.values():
        if s["n_not_mapped"] == 1 and s["n_aligned"] == 0:
            n_guides_unmapped += 1
            continue
        if s["n_valid"] == 0 and s["n_aligned"] > 0:
            n_guides_aligned_none_valid += 1
            continue
        if s["n_valid"] == 1:
            n_guides_unique_valid += 1
            if s["n_discarded"] > 0:
                n_guides_one_valid_plus_invalid += 1
            continue
        if s["n_valid"] > 1:
            n_guides_multi_valid += 1

    if output_summary_tsv == "auto":
        summary_path = base_output_dir / "alignment_summary.tsv"
    else:
        summary_path = Path(output_summary_tsv)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("metric\tcount\n")
        handle.write(f"guides_unique_valid\t{n_guides_unique_valid}\n")
        handle.write(f"guides_multi_valid\t{n_guides_multi_valid}\n")
        handle.write(f"guides_aligned_none_valid\t{n_guides_aligned_none_valid}\n")
        handle.write(f"guides_unmapped\t{n_guides_unmapped}\n")
        handle.write(f"guides_one_valid_plus_invalid\t{n_guides_one_valid_plus_invalid}\n")

    n_valid_guides = len(valid_by_guide)
    n_unique_guides = sum(1 for v in valid_by_guide.values() if len(v) == 1)
    n_multi_guides = sum(1 for v in valid_by_guide.values() if len(v) > 1)
    n_valid_alignments = sum(len(v) for v in valid_by_guide.values())

    return {
        "valid_guides": n_valid_guides,
        "unique_guides": n_unique_guides,
        "multi_guides": n_multi_guides,
        "valid_alignments": n_valid_alignments,
        "invalid_alignments": len(invalid_rows),
        "guides_unique_valid": n_guides_unique_valid,
        "guides_multi_valid": n_guides_multi_valid,
        "guides_aligned_none_valid": n_guides_aligned_none_valid,
        "guides_unmapped": n_guides_unmapped,
        "guides_one_valid_plus_invalid": n_guides_one_valid_plus_invalid,
    }


def _protospacer_bed_span(
    aln,
    *,
    query_len: int,
    pam_len: int,
    allow_leading_g_softclip: bool,
) -> tuple[int, int] | None:
    """Map protospacer query region (no PAM) to genomic BED span."""
    if query_len <= pam_len:
        return None

    query_ops = [None] * query_len
    qpos = 0
    for op, length in (aln.cigartuples or []):
        if op in _CIGAR_MATCH | {_CIGAR_INS, _CIGAR_SOFT}:
            for _ in range(length):
                if qpos < query_len:
                    query_ops[qpos] = op
                qpos += 1

    if aln.is_reverse:
        protospacer_start = pam_len
        protospacer_end = query_len
        if allow_leading_g_softclip and query_ops and query_ops[-1] == _CIGAR_SOFT:
            protospacer_end -= 1
    else:
        protospacer_start = 0
        if allow_leading_g_softclip and query_ops and query_ops[0] == _CIGAR_SOFT:
            protospacer_start = 1
        protospacer_end = query_len - pam_len
    if protospacer_start >= protospacer_end:
        return None

    q2r = {q: r for q, r in aln.get_aligned_pairs(matches_only=True)}
    ref_positions = [q2r[i] for i in range(protospacer_start, protospacer_end) if i in q2r]
    if not ref_positions:
        return None
    return min(ref_positions), max(ref_positions) + 1


def _resolve_allowed_contigs(
    *,
    primary_contigs: set[str] | list[str] | tuple[str, ...] | None,
    chromsizes: str | Path | pd.DataFrame | None,
) -> set[str] | None:
    """Resolve allowed contigs from direct list/set or chromsizes input.

    If both are None, returns None (no contig filtering).
    """
    if primary_contigs is not None:
        return {str(x) for x in primary_contigs}

    if chromsizes is None:
        return None

    if isinstance(chromsizes, pd.DataFrame):
        if chromsizes.shape[1] == 0:
            return set()
        return set(chromsizes.iloc[:, 0].astype(str).tolist())

    path = Path(chromsizes)
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] == 0:
        return set()
    return set(df.iloc[:, 0].astype(str).tolist())

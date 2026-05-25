#!/usr/bin/env python3
"""End-to-end guide alignment and QC pipeline.

Steps:
1. Build guide FASTQ from TSV/DataFrame-like input
2. Build GEM index
3. Map guides with GEM
4. Filter alignments into QC outputs
"""

from __future__ import annotations

import argparse
from pathlib import Path

import crispr_analysis_utils as cau


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run guide alignment QC pipeline.")

    p.add_argument("--guides-tsv", required=True, help="Input guide TSV (id, sequence).")
    p.add_argument("--reference-fasta", required=True, help="Reference FASTA path.")
    p.add_argument("--chromsizes", default=None, help="Chromsizes path for allowed contigs.")
    p.add_argument("--outdir", required=True, help="Output directory.")

    p.add_argument("--threads", type=int, default=8, help="Threads for GEM steps.")
    p.add_argument("--pam", default="NGG", help="PAM string (default: NGG).")
    p.add_argument("--add-leading-g", action="store_true", help="Prepend G if missing.")
    p.add_argument("--mapping-mode", default="sensitive", help="GEM mapping mode.")
    p.add_argument(
        "--allow-leading-g-softclip",
        action="store_true",
        default=True,
        help="Allow leading G soft-clipping in filtering (default: enabled).",
    )

    return p


def main() -> None:
    args = build_parser().parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    guides_fastq = outdir / "guides_input.fastq"
    gem_index_prefix = outdir / "genome_index"
    gem_index = outdir / "genome_index.gem"
    mapped_sam = outdir / "guides_mapped.sam"
    valid_unique_sam = outdir / "guides_valid_unique.sam"
    valid_multi_sam = outdir / "guides_valid_multi.sam"

    cau.guide_qc.guides_to_fastq(
        args.guides_tsv,
        guides_fastq,
        pam=args.pam,
        add_leading_g=args.add_leading_g,
    )

    cau.gem_mapper.build_gem_index(
        args.reference_fasta,
        gem_index_prefix,
        threads=args.threads,
    )

    cau.gem_mapper.map_guides_with_gem(
        gem_index,
        guides_fastq,
        mapped_sam,
        threads=args.threads,
        mapping_mode=args.mapping_mode,
        sam_compact=False,
    )

    summary = cau.guide_qc.filter_guide_alignments(
        mapped_sam,
        valid_unique_sam,
        valid_multi_sam,
        pam=args.pam,
        chromsizes=args.chromsizes,
        allow_leading_g_softclip=args.allow_leading_g_softclip,
    )
    print(summary)


if __name__ == "__main__":
    main()

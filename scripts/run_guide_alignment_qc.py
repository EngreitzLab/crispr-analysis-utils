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
from datetime import datetime
import logging
from pathlib import Path
import shlex
import time

import crispr_analysis_utils as cau
import pandas as pd


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

    gem_index_dir = outdir / "gem_index"
    alignment_dir = outdir / "alignment_outputs"
    filtering_dir = outdir / "guide_alignments_outputs"
    for p in (gem_index_dir, alignment_dir, filtering_dir):
        p.mkdir(parents=True, exist_ok=True)

    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "pipeline.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="w"),
            logging.StreamHandler(),
        ],
    )

    started = time.time()
    start_iso = datetime.now().isoformat(timespec="seconds")
    logging.info("Pipeline start: %s", start_iso)
    logging.info("Output directory: %s", outdir)
    logging.info("Log file: %s", log_path)

    guides_fastq = alignment_dir / "guides_input.fastq"
    gem_index_prefix = gem_index_dir / "genome_index"
    gem_index = gem_index_dir / "genome_index.gem"
    mapped_sam = alignment_dir / "guides_mapped.sam"

    logging.info("Step 1/4: Build guide FASTQ")
    guides_df = pd.read_csv(args.guides_tsv, sep="\t", header=None)
    alias_by_sequence = {
        str(row.iloc[1]).strip().upper(): str(row.iloc[0])
        for _, row in guides_df.iterrows()
    }
    if guides_fastq.exists():
        logging.info("FASTQ exists, skipping: %s", guides_fastq)
    else:
    cau.guide_qc.guides_to_fastq(
        guides_df,
        guides_fastq,
        pam=args.pam,
        add_leading_g=args.add_leading_g,
        id_col=0,
        sequence_col=1,
    )

    logging.info("Step 2/4: Build GEM index")
    index_cmd_parts = [
        "gem-indexer",
        "-i",
        str(args.reference_fasta),
        "-o",
        str(gem_index_prefix),
        "-t",
        str(args.threads),
    ]
    index_cmd = " ".join(shlex.quote(x) for x in index_cmd_parts)
    (gem_index_dir / "index_command.sh").write_text(index_cmd + "\n", encoding="utf-8")
    (gem_index_dir / "index_inputs.txt").write_text(
        f"reference_fasta={Path(args.reference_fasta).resolve()}\n",
        encoding="utf-8",
    )
    if gem_index.exists():
        logging.info("GEM index exists, skipping re-index: %s", gem_index)
    else:
        cau.gem_mapper.build_gem_index(
            args.reference_fasta,
            gem_index_prefix,
            threads=args.threads,
            log_path=gem_index_dir / "gem_index.log",
        )

    logging.info("Step 3/4: Map guides with GEM")
    if mapped_sam.exists():
        logging.info("Mapped SAM exists, skipping alignment: %s", mapped_sam)
    else:
        cau.gem_mapper.map_guides_with_gem(
            gem_index,
            guides_fastq,
            mapped_sam,
            threads=args.threads,
            mapping_mode=args.mapping_mode,
            sam_compact=False,
            log_path=alignment_dir / "guides_mapped.log",
        )

    logging.info("Step 4/4: Filter alignments and compute QC outputs")
    summary_tsv = filtering_dir / "alignment_summary.tsv"
    if summary_tsv.exists():
        logging.info("Filtering summary exists, skipping filtering: %s", summary_tsv)
        summary = {"status": "skipped_filtering", "alignment_summary_tsv": str(summary_tsv)}
    else:
        summary = cau.guide_qc.filter_guide_alignments(
            mapped_sam,
            None,
            None,
            pam=args.pam,
            chromsizes=args.chromsizes,
            allow_leading_g_softclip=args.allow_leading_g_softclip,
            alias_by_guide_id=alias_by_sequence,
            output_valid_bed=filtering_dir / "valid_alignments.bed",
            output_discarded_tsv=filtering_dir / "discarded_alignments.tsv",
            output_unmapped_tsv=filtering_dir / "unmapped.tsv",
            output_guide_log_tsv=filtering_dir / "guide_alignment_log.tsv",
            output_invalid_tsv=filtering_dir / "invalid_alignments.tsv",
            output_summary_tsv=summary_tsv,
        )
    finished = time.time()
    finish_iso = datetime.now().isoformat(timespec="seconds")
    logging.info("Pipeline finish: %s", finish_iso)
    logging.info("Elapsed seconds: %.2f", finished - started)
    logging.info("Summary: %s", summary)
    print(summary)


if __name__ == "__main__":
    main()

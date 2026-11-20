"""``cau pairs``: the discovery and positive-control pairs a screen tests."""

from __future__ import annotations

import argparse

from . import CommandError


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Declare the ``run`` action and its arguments."""
    actions = parser.add_subparsers(dest="action", metavar="ACTION", required=True)

    run_parser = actions.add_parser(
        "run",
        help="Build the guide-gene pairs from the IGVF per-guide table.",
        description=(
            "Pair every positive-control guide with the genes the design names"
            " for it, and every targeting guide with each readout gene whose"
            " TSS lies near its element. --outdir gets pairs/ with"
            " pairs_to_test.csv, its columns YAML, qc.tsv and"
            " panel_genes_without_pairs.tsv, plus run_metadata/ and logs/."
        ),
    )
    run_parser.add_argument(
        "--guide-metadata",
        required=True,
        help="guide_metadata.tsv from `cau igvf-export run`.",
    )
    run_parser.add_argument(
        "--genes",
        required=True,
        help=(
            "The readout panel: a .csv, .tsv or .txt (optionally .gz) with a"
            " column of ENSEMBL gene ids."
        ),
    )
    run_parser.add_argument(
        "--gene-column",
        default="gene_id",
        help="The panel's column of ENSEMBL gene ids (default: gene_id).",
    )
    run_parser.add_argument(
        "--store",
        required=True,
        help="The annotation store the guide table was built with.",
    )
    run_parser.add_argument(
        "--outdir",
        required=True,
        help="Where the table goes. Nothing is written outside it.",
    )
    run_parser.add_argument(
        "--max-distance",
        type=int,
        default=1_000_000,
        help=(
            "Largest distance, in bases, from an element's midpoint to a"
            " readout gene's TSS for a discovery pair (default: 1000000)."
        ),
    )


def run(args: argparse.Namespace) -> int:
    """Run the ``run`` action."""
    from ..guide_annotation.db import StoreError
    from ..pairs import DISCOVERY, POSITIVE_CONTROL, PairsError
    from ..pairs import run as build

    try:
        result = build(
            args.guide_metadata,
            args.genes,
            args.store,
            args.outdir,
            max_distance=args.max_distance,
            gene_column=args.gene_column,
        )
    except (FileNotFoundError, PairsError, StoreError, ValueError, OSError) as error:
        raise CommandError(str(error)) from error

    found = result.metadata["counts"]
    for pair_type in (POSITIVE_CONTROL, DISCOVERY):
        part = found[pair_type]
        print(
            f"{pair_type}: {part['pairs']} pairs, {part['guides']} guides,"
            f" {part['elements']} elements, {part['genes']} genes"
        )
    print(
        f"{found['panel_genes_without_pairs']} of {found['panel_genes']} panel genes have no pair"
    )
    print(f"Outputs in {result.outdir}")
    return 0

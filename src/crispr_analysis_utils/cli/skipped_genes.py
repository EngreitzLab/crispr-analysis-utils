"""``cau skipped-genes``: how many genes a pair's link passes over."""

from __future__ import annotations

import argparse

from . import CommandError


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Declare the ``run`` action and its arguments."""
    actions = parser.add_subparsers(dest="action", metavar="ACTION", required=True)

    run_parser = actions.add_parser(
        "run",
        help="Classify every pair by the genes its link passes over.",
        description=(
            "For each pair, rank every annotated gene on the element's"
            " chromosome by distance to the element's midpoint, and classify"
            " it by how many nearer genes on the same side are skipped, and"
            " whether those are in the readout panel. --outdir gets"
            " skipped_genes/ with skipped_genes.csv, its columns YAML and"
            " qc.tsv, plus run_metadata/ and logs/."
        ),
    )
    run_parser.add_argument(
        "--pairs",
        required=True,
        help="pairs_to_test.csv from `cau pairs run`.",
    )
    run_parser.add_argument(
        "--store",
        required=True,
        help="The annotation store the pairs were built with.",
    )
    run_parser.add_argument(
        "--genes",
        required=True,
        help="The readout panel, the same file given to `cau pairs run --genes`.",
    )
    run_parser.add_argument(
        "--gene-column",
        default="gene_id",
        help="The panel's column of ENSEMBL gene ids (default: gene_id).",
    )
    run_parser.add_argument(
        "--outdir",
        required=True,
        help="Where the table goes. Nothing is written outside it.",
    )


def run(args: argparse.Namespace) -> int:
    """Run the ``run`` action."""
    from ..guide_annotation.db import StoreError
    from ..pairs import PairsError
    from ..skipped_genes import CATEGORIES, SkippedGenesError
    from ..skipped_genes import run as classify

    try:
        result = classify(
            args.pairs,
            args.store,
            args.genes,
            args.outdir,
            gene_column=args.gene_column,
        )
    except (
        FileNotFoundError,
        PairsError,
        SkippedGenesError,
        StoreError,
        ValueError,
        OSError,
    ) as error:
        raise CommandError(str(error)) from error

    found = result.metadata["counts"]
    for category in CATEGORIES:
        print(f"{category}: {found[category]} pairs")
    print(f"Outputs in {result.outdir}")
    return 0

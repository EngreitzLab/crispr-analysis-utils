"""``cau element-mapping``: map a guide library to the elements a screen targets."""

from __future__ import annotations

import argparse

from . import CommandError


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Declare the ``run`` action and its arguments."""
    actions = parser.add_subparsers(dest="action", metavar="ACTION", required=True)

    run_parser = actions.add_parser(
        "run",
        help="Map every guide to the elements it targets.",
        description=(
            "Classify every guide against the elements a screen targets. "
            "--outdir gets elements/ with guide_element_mapping.tsv and "
            "columns.yaml, run_metadata/ and logs/."
        ),
    )
    run_parser.add_argument(
        "--alignment",
        required=True,
        help="guides_alignment_metadata.tsv from `cau guide-alignment run`.",
    )
    run_parser.add_argument(
        "--elements",
        required=True,
        help="A BED of the regions the screen targets.",
    )
    run_parser.add_argument(
        "--outdir",
        required=True,
        help="Where the table goes. Nothing is written outside it.",
    )
    run_parser.add_argument(
        "--annotation",
        default=None,
        help=(
            "guides_gene_annotation.tsv from `cau guide-annotation run`. "
            "Without it no site can be called safe_targeting."
        ),
    )
    run_parser.add_argument(
        "--safe-distance",
        type=int,
        default=1_000_000,
        help=(
            "How far from a protein-coding TSS a site must be to count as "
            "safe_targeting (default: 1000000)."
        ),
    )


def run(args: argparse.Namespace) -> int:
    """Run the ``run`` action."""
    from ..element_mapping import CLASSES, ElementError
    from ..element_mapping import run as map_elements

    try:
        result = map_elements(
            args.alignment,
            args.elements,
            args.outdir,
            annotation=args.annotation,
            safe_distance=args.safe_distance,
        )
    except (FileNotFoundError, ElementError, ValueError, OSError, KeyError) as error:
        raise CommandError(str(error)) from error

    found = result.metadata["counts"]
    print(
        f"{found['guides']} guides, {found['sites']} sites: "
        + ", ".join(f"{found['classes'][name]} {name}" for name in CLASSES)
    )
    print(f"{found['elements_with_a_guide']} elements have at least one guide")
    print(f"Outputs in {result.outdir}")
    return 0

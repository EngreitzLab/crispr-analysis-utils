"""``cau igvf-export``: write the guide design in the IGVF per-guide format."""

from __future__ import annotations

import argparse

from . import CommandError


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Declare the ``run`` action and its arguments."""
    actions = parser.add_subparsers(dest="action", metavar="ACTION", required=True)

    run_parser = actions.add_parser(
        "run",
        help="Write the IGVF per-guide metadata table.",
        description=(
            "Assemble the IGVF per-guide metadata table from the alignment, "
            "the element mapping and the annotation store. --outdir gets "
            "igvf/ with guide_metadata.tsv and columns.yaml, run_metadata/ "
            "and logs/."
        ),
    )
    run_parser.add_argument(
        "--alignment",
        required=True,
        help="guides_alignment_metadata.tsv from `cau guide-alignment run`.",
    )
    run_parser.add_argument(
        "--mapping",
        required=True,
        help="guide_element_mapping.tsv from `cau element-mapping run`.",
    )
    run_parser.add_argument(
        "--elements", required=True, help="The elements BED the mapping used."
    )
    run_parser.add_argument(
        "--store", required=True, help="The annotation store, from `build-db`."
    )
    run_parser.add_argument(
        "--outdir",
        required=True,
        help="Where the table goes. Nothing is written outside it.",
    )
    run_parser.add_argument(
        "--promoters",
        default=None,
        help="A TSS BED. Without one, the store's reference transcript starts.",
    )
    run_parser.add_argument(
        "--genes",
        default=None,
        help=(
            "A table of measured element-to-gene pairs: an element id column "
            "and a column of ENSEMBL gene ids."
        ),
    )
    run_parser.add_argument(
        "--promoter-window",
        type=int,
        default=2000,
        help=(
            "Within this of a TSS an element is a promoter (default: 2000, "
            "the distance the format uses to define a distal element)."
        ),
    )


def run(args: argparse.Namespace) -> int:
    """Run the ``run`` action."""
    from ..element_mapping import ElementError
    from ..guide_annotation import PromoterError, StoreError
    from ..igvf_export import TYPES, ExportError
    from ..igvf_export import run as export

    try:
        result = export(
            args.alignment,
            args.mapping,
            args.elements,
            args.store,
            args.outdir,
            promoters=args.promoters,
            genes=args.genes,
            promoter_window=args.promoter_window,
        )
    except (
        FileNotFoundError,
        ExportError,
        ElementError,
        PromoterError,
        StoreError,
        ValueError,
        OSError,
        KeyError,
    ) as error:
        raise CommandError(str(error)) from error

    found = result.metadata["counts"]
    print(
        f"{found['guides']} guides: "
        + ", ".join(f"{found['type'][name]} {name}" for name in TYPES)
    )
    print(
        "elements: "
        + ", ".join(
            f"{count} {name}"
            for name, count in found["genomic_element"].items()
            if count
        )
    )
    print(f"Outputs in {result.outdir}")
    return 0

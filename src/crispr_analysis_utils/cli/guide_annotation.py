"""``cau guide-annotation``: place every aligned site in a gene.

``build-db`` reduces a GENCODE GTF to an annotation store, once. ``run`` reads
that store and the alignment table and writes one annotated table.
"""

from __future__ import annotations

import argparse

from . import CommandError


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Declare the ``build-db`` and ``run`` actions and their arguments."""
    actions = parser.add_subparsers(dest="action", metavar="ACTION", required=True)

    build = actions.add_parser(
        "build-db",
        help="Reduce a GENCODE GTF to an annotation store.",
        description=(
            "Read a GENCODE GTF and write the store a run reads. Reading the "
            "GTF takes about a minute; the store is read in about a second. "
            "The GTF must hold transcript and exon rows: a genes-only GTF "
            "cannot say exon, intron or transcript."
        ),
    )
    build.add_argument("gtf", help="The GENCODE GTF, gzipped or not.")
    build.add_argument(
        "--out",
        required=True,
        help="The store to write, a parquet file.",
    )
    build.add_argument(
        "--gene-types",
        default="protein_coding",
        help="The gene types to keep, comma-separated (default: protein_coding).",
    )

    annotate = actions.add_parser(
        "run",
        help="Annotate every aligned site against the genes of a store.",
        description=(
            "Place every mapped site in a gene. --outdir gets annotations/ "
            "with guides_gene_annotation.tsv and guides_gene_annotation_columns.yaml, run_metadata/ "
            "with this invocation's record, and logs/."
        ),
    )
    annotate.add_argument(
        "--alignment",
        required=True,
        help="guides_alignment_metadata.tsv from `cau guide-alignment run`.",
    )
    annotate.add_argument(
        "--store",
        required=True,
        help="The annotation store, from `build-db`.",
    )
    annotate.add_argument(
        "--outdir",
        required=True,
        help="Where the table goes. Nothing is written outside it.",
    )
    annotate.add_argument(
        "--promoters",
        default=None,
        help=(
            "A TSS BED, one site per gene. Without one, each gene's reference "
            "transcript start is used."
        ),
    )
    annotate.add_argument(
        "--promoter-window",
        type=int,
        default=1000,
        help="How far either side of a TSS counts as its promoter (default: 1000).",
    )


def run(args: argparse.Namespace) -> int:
    """Run the chosen action."""
    if args.action == "build-db":
        return _build_db(args)
    return _annotate(args)


def _build_db(args: argparse.Namespace) -> int:
    """Build the annotation store."""
    from ..guide_annotation import StoreError, build_store

    gene_types = tuple(
        name.strip() for name in args.gene_types.split(",") if name.strip()
    )
    if not gene_types:
        raise CommandError("--gene-types was given no gene type.")
    try:
        store = build_store(args.gtf, args.out, gene_types=gene_types)
    except (FileNotFoundError, StoreError, OSError, ValueError) as error:
        raise CommandError(str(error)) from error

    counts = store.provenance["counts"]
    print(
        f"{store.provenance['genes']} genes, "
        + ", ".join(f"{counts[name]} {name}" for name in ("transcript", "exon", "CDS"))
    )
    print(f"Store in {args.out}")
    return 0


def _annotate(args: argparse.Namespace) -> int:
    """Annotate the aligned sites."""
    from ..guide_annotation import AnnotationError, PromoterError, StoreError
    from ..guide_annotation import run as annotate
    from ..guide_annotation.annotate import LOCATIONS

    try:
        result = annotate(
            args.alignment,
            args.store,
            args.outdir,
            promoters=args.promoters,
            promoter_window=args.promoter_window,
        )
    except (
        FileNotFoundError,
        AnnotationError,
        PromoterError,
        StoreError,
        ValueError,
        OSError,
    ) as error:
        raise CommandError(str(error)) from error

    found = result.metadata["counts"]
    print(
        f"{found['rows']} rows for {found['sites']} sites in {found['genes']} genes: "
        + ", ".join(f"{found['location'][name]} {name}" for name in LOCATIONS)
        + f", {found['promoter']} in a promoter"
    )
    print(f"Outputs in {result.outdir}")
    return 0

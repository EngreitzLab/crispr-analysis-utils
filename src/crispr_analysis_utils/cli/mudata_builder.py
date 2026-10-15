"""``cau mudata-builder``: assemble one channel's matrices into a MuData."""

from __future__ import annotations

import argparse
from pathlib import Path

from . import CommandError


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Declare the ``run`` action and its arguments."""
    actions = parser.add_subparsers(dest="action", metavar="ACTION", required=True)
    run_parser = actions.add_parser(
        "run",
        help="Assemble a gene and a guide count matrix into one MuData.",
        description=(
            "Assembly only: nothing is computed from the counts, no cell is "
            "dropped and no guide is called. Writes mudata/mudata.h5mu, "
            "run_metadata/ and logs/."
        ),
    )
    run_parser.add_argument("--rna", required=True, help="The gene counts, as h5ad.")
    run_parser.add_argument("--guide", required=True, help="The guide counts, as h5ad.")
    run_parser.add_argument(
        "--guide-metadata",
        required=True,
        help="guide_metadata.tsv from `cau igvf-export run`.",
    )
    run_parser.add_argument(
        "--store",
        required=True,
        help=(
            "The annotation store, which must be the one the guide metadata "
            "was built against. Every gene type in it is used, so a lncRNA "
            "feature is named like any other."
        ),
    )
    run_parser.add_argument(
        "--measurement-set",
        required=True,
        help=(
            "The channel these matrices are. The cell index becomes "
            "{barcode}_{measurement_set}, so channels can be concatenated."
        ),
    )
    run_parser.add_argument("--outdir", required=True, help="Where the object goes.")
    run_parser.add_argument(
        "--genome", help="The reference FASTA, for its checksum. Nothing is read."
    )
    run_parser.add_argument(
        "--gtf-accession",
        help=(
            "The GTF's IGVF portal id. Without it the object records the "
            "file's name instead."
        ),
    )
    run_parser.add_argument("--genome-accession", help="The genome's IGVF portal id.")
    run_parser.add_argument(
        "--guide-annotation",
        help=(
            "guide_gene_annotation.tsv. Carried whole in guide.uns: it is "
            "one row per guide per site per gene, so it cannot go in var."
        ),
    )
    run_parser.add_argument(
        "--element-mapping",
        help="guide_element_mapping.tsv. Carried whole in guide.uns.",
    )
    run_parser.add_argument(
        "--channel",
        action="append",
        default=[],
        metavar="NAME=VALUE",
        help=(
            "What this channel is, repeatable: --channel replicate=1 "
            "--channel day=4 --channel chip=A. One value for the whole "
            "channel, joined to every cell's obs and kept in uns."
        ),
    )
    run_parser.add_argument(
        "--readout-genes",
        help=(
            "The panel a targeted assay reads out, as a table with a column "
            "of ENSEMBL ids. Adds a boolean readout column to gene.var; "
            "without it the column is absent."
        ),
    )
    run_parser.add_argument(
        "--covariates",
        help=(
            "Per-cell covariates, barcode in the first column. Carried into "
            "the root obs unchanged."
        ),
    )


def _channel(pairs: list[str]) -> dict[str, str]:
    """NAME=VALUE arguments as a mapping, refusing one without the equals."""
    out = {}
    for pair in pairs:
        name, sep, value = pair.partition("=")
        if not sep or not name:
            raise CommandError(f"--channel takes NAME=VALUE, got {pair!r}.")
        out[name] = value
    return out


def run(args: argparse.Namespace) -> int:
    """Run the action."""
    from ..mudata_builder import BuildError, run

    try:
        result = run(
            args.rna,
            args.guide,
            args.guide_metadata,
            args.store,
            args.outdir,
            measurement_set=args.measurement_set,
            genome=args.genome,
            gtf_accession=args.gtf_accession,
            genome_accession=args.genome_accession,
            covariates=args.covariates,
            metadata=_channel(args.channel),
            readout_genes=args.readout_genes,
            guide_annotation=args.guide_annotation,
            element_mapping=args.element_mapping,
        )
    except (BuildError, FileNotFoundError, OSError) as error:
        raise CommandError(str(error)) from error

    found = result.metadata["counts"]
    print(
        f"{found['cells']} cells ({found['cells_both']} in both, "
        f"{found['cells_gene_only']} gene only, "
        f"{found['cells_guide_only']} guide only): "
        f"{found['genes']} genes, {found['guides']} guides"
    )
    if found["genes_unannotated"] or found["guides_unannotated"]:
        print(
            f"Unannotated: {found['genes_unannotated']} genes, "
            f"{found['guides_unannotated']} guides"
        )
    for name, rows in found.get("uns_tables", {}).items():
        print(f"guide.uns[{name!r}]: {rows} rows")
    print(f"Outputs in {Path(result.outdir)}")
    return 0

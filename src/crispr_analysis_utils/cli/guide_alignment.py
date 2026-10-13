"""``cau guide-alignment``: find every guide's genomic site.

``cau guide-alignment run --guides GUIDES --reference FASTA --outdir DIR``
scans the reference once and writes the tables (see
``crispr_analysis_utils.guide_alignment.run``). There is no index to build:
a guide is designed, not sequenced, so its site is an exact match.
"""

from __future__ import annotations

import argparse
from collections import Counter

from . import CommandError


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Declare the ``run`` action and its arguments."""
    actions = parser.add_subparsers(dest="action", metavar="ACTION", required=True)

    search = actions.add_parser(
        "run",
        help="Find every guide's genomic site.",
        description=(
            "Scan a reference for every guide's exact site. --outdir gets "
            "alignments/ with the results (guides_metadata.tsv, sites.tsv, "
            "orientation.tsv, guides.bed, cut_sites.bed, qc.tsv), "
            "run_metadata/ with this invocation's record, and logs/."
        ),
    )
    search.add_argument(
        "--guides",
        required=True,
        help="The guide table: a name column and a sequence column.",
    )
    search.add_argument(
        "--reference",
        required=True,
        help="The reference FASTA, uncompressed.",
    )
    search.add_argument(
        "--outdir",
        required=True,
        help="Where the tables go. Nothing is written outside it.",
    )
    search.add_argument(
        "--pam",
        default="NGG",
        help="The PAM, IUPAC codes, 3' of the protospacer (default: NGG).",
    )
    search.add_argument(
        "--spacer-length",
        type=int,
        default=20,
        help="The protospacer length, 3'-most in each sequence (default: 20).",
    )
    search.add_argument(
        "--add-leading-g",
        action="store_true",
        help=(
            "Put a 5' G in front of every sequence, for a table that gives "
            "bare protospacers while the construct adds one."
        ),
    )
    search.add_argument(
        "--orientation-check",
        choices=("warn", "error", "off"),
        default="error",
        help=(
            "Check that the table reads 5' to 3', by holding both readings of "
            "every guide: stop when the reversed ones outnumber the rest, only "
            "warn, or skip the check (default: error)."
        ),
    )
    search.add_argument(
        "--fix-orientation",
        action="store_true",
        help=(
            "Align the reverse complement of the entries that read 3' to 5', "
            "rather than leaving them as given."
        ),
    )
    search.add_argument(
        "--contigs",
        default=None,
        help="Read only these contigs, comma-separated (default: all of them).",
    )
    search.add_argument(
        "--id-col",
        default="0",
        help="The name column: a 0-based position, or a name (default: 0).",
    )
    search.add_argument(
        "--spacer-col",
        default="1",
        help="The sequence column: a 0-based position, or a name (default: 1).",
    )
    search.add_argument(
        "--sep",
        type=_separator,
        default="\t",
        help=(
            "The guide table's field separator: tab, comma, or one character "
            "(default: tab)."
        ),
    )
    search.add_argument(
        "--header",
        action="store_true",
        help="The guide table's first line holds column names.",
    )


def run(args: argparse.Namespace) -> int:
    """Run the ``run`` action."""
    from ..guide_alignment import OrientationError
    from ..guide_alignment import run as align
    from ..guide_alignment.summary import CLASSES

    try:
        result = align(
            args.guides,
            args.reference,
            args.outdir,
            pam=args.pam,
            spacer_length=args.spacer_length,
            add_leading_g=args.add_leading_g,
            orientation_check=args.orientation_check,
            fix_orientation=args.fix_orientation,
            contigs=_contigs(args.contigs),
            id_col=_column(args.id_col),
            spacer_col=_column(args.spacer_col),
            sep=args.sep,
            header=args.header,
        )
    except (FileNotFoundError, ValueError, OrientationError, OSError) as error:
        raise CommandError(str(error)) from error

    classes = Counter(summary.site_class for summary in result.summaries)
    print(
        f"{len(result.summaries)} guides: "
        + ", ".join(f"{classes[name]} {name}" for name in CLASSES)
    )
    print(f"Outputs in {result.outdir}")
    return 0


def _contigs(value: str | None) -> list[str] | None:
    """A comma-separated contig list, or None for every contig."""
    if value is None:
        return None
    names = [name.strip() for name in value.split(",") if name.strip()]
    if not names:
        raise CommandError("--contigs was given no contig names.")
    return names


def _separator(value: str) -> str:
    """The field separator: a name, an escape, or a single character."""
    names = {"tab": "\t", "comma": ",", "semicolon": ";", "space": " ", "\\t": "\t"}
    if value.lower() in names:
        return names[value.lower()]
    if len(value) != 1:
        raise argparse.ArgumentTypeError(
            f"--sep takes tab, comma, or one character, not {value!r}."
        )
    return value


def _column(value: str) -> int | str:
    """A column given as a 0-based position, or as a name."""
    try:
        return int(value)
    except ValueError:
        return value

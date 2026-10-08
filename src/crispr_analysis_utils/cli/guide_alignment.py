"""``cau guide-alignment``: index a reference, or align a guide library to it.

``cau guide-alignment index REFERENCE PREFIX`` builds the GEM index.
``cau guide-alignment run`` maps a guide library with it and writes the
output tables (see ``crispr_analysis_utils.guide_alignment.run``).
"""

from __future__ import annotations

import argparse
import logging
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager

from . import CommandError

_SEPARATORS = {"tab": "\t", "\\t": "\t", "comma": ","}


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Declare the ``index`` and ``run`` actions and their arguments."""
    actions = parser.add_subparsers(dest="action", metavar="ACTION", required=True)

    index = actions.add_parser(
        "index",
        help="Build a GEM index from a reference FASTA.",
        description=(
            "Build a GEM index from an uncompressed reference FASTA. Writes "
            "PREFIX.gem and PREFIX.info, gem-indexer's output in PREFIX.log, "
            "and PREFIX.gem.json, which records the reference's MD5 checksum "
            "so that `run` can check it is given the same reference."
        ),
    )
    index.add_argument("reference", help="The reference FASTA, uncompressed.")
    index.add_argument(
        "prefix", help="The index path without .gem. Missing folders are created."
    )
    index.add_argument(
        "--threads",
        type=_positive_int,
        default=8,
        help="gem-indexer threads (default: 8).",
    )

    run = actions.add_parser(
        "run",
        help="Align a guide library to its genomic sites and classify each guide.",
        description=(
            "Align a guide library with GEM3, check every hit against the "
            "reference, and write guides.tsv, sites.tsv, orientation.tsv, "
            "guides.bed, cut_sites.bed, summary.tsv and run.json into --outdir."
        ),
    )
    inputs = run.add_argument_group("inputs and outputs")
    inputs.add_argument(
        "--guides",
        required=True,
        help="The guide table: a guide id and its sequence, 5' to 3', per line.",
    )
    inputs.add_argument(
        "--reference",
        required=True,
        help="The reference FASTA the index was built from, uncompressed.",
    )
    inputs.add_argument(
        "--index", required=True, help="The GEM index: PREFIX.gem or PREFIX."
    )
    inputs.add_argument(
        "--outdir",
        required=True,
        help="The output folder, created if missing. Nothing is written elsewhere.",
    )

    search = run.add_argument_group("search")
    search.add_argument(
        "--max-mismatches",
        type=_non_negative_int,
        default=3,
        help="Spacer mismatches a site may have (default: 3).",
    )
    search.add_argument(
        "--pam",
        default="NGG",
        help="The primary PAM, IUPAC codes, 3' of the protospacer (default: NGG).",
    )
    search.add_argument(
        "--alt-pams",
        type=_pam_list,
        default=["NAG", "NGA"],
        help=(
            'Alternative PAMs, comma-separated; "" for none. They are counted '
            "but never decide a guide's class (default: NAG,NGA)."
        ),
    )
    search.add_argument(
        "--spacer-length",
        type=_positive_int,
        default=20,
        help="The spacer is the 3'-most bases of each sequence (default: 20).",
    )
    search.add_argument(
        "--add-leading-g",
        action="store_true",
        help=(
            "Align a spacer that does not start with G with an extra 5' G, "
            "whose genomic base is reported but never counted as a mismatch."
        ),
    )
    search.add_argument(
        "--orientation-check",
        choices=("warn", "error", "off"),
        default="warn",
        help=(
            "Check that the spacers are given 5' to 3', from their exact hits: "
            "warn, stop with an error when reversed guides outnumber ok ones, "
            "or skip the check (default: warn)."
        ),
    )
    search.add_argument(
        "--contigs",
        type=_name_list,
        default=None,
        help="Keep only hits on these contigs, comma-separated (default: all).",
    )
    search.add_argument(
        "--threads",
        type=_positive_int,
        default=8,
        help="gem-mapper threads (default: 8).",
    )

    table = run.add_argument_group("guide table")
    table.add_argument(
        "--id-col",
        type=_column,
        default=0,
        help="The guide id column: a 0-based position, or a name (default: 0).",
    )
    table.add_argument(
        "--spacer-col",
        type=_column,
        default=1,
        help="The sequence column: a 0-based position, or a name (default: 1).",
    )
    table.add_argument(
        "--sep",
        type=_separator,
        default="\t",
        help="The field separator: tab, comma, or one character (default: tab).",
    )
    table.add_argument(
        "--header",
        action="store_true",
        help="The first line holds column names, as in a Cell Ranger CSV.",
    )


def run(args: argparse.Namespace) -> int:
    """Run the ``index`` or ``run`` action."""
    from ..guide_alignment import GemError, OrientationError, build_index
    from ..guide_alignment import run as align
    from ..guide_alignment.summary import CLASSES

    expected = (FileNotFoundError, ValueError, GemError, OrientationError)
    with _console_logging():
        try:
            if args.action == "index":
                index = build_index(args.reference, args.prefix, threads=args.threads)
                print(index)
                return 0
            result = align(
                args.guides,
                args.reference,
                args.index,
                args.outdir,
                max_mismatches=args.max_mismatches,
                pam=args.pam,
                alt_pams=args.alt_pams,
                spacer_length=args.spacer_length,
                add_leading_g=args.add_leading_g,
                orientation_check=args.orientation_check,
                contigs=args.contigs,
                threads=args.threads,
                id_col=args.id_col,
                spacer_col=args.spacer_col,
                sep=args.sep,
                header=args.header,
            )
        except expected as error:
            raise CommandError(str(error)) from error
    classes = Counter(summary.alignment_class for summary in result.summaries)
    print(
        f"{len(result.summaries)} guides: "
        + ", ".join(f"{classes[name]} {name}" for name in CLASSES)
    )
    print(f"Outputs in {result.outdir}")
    return 0


@contextmanager
def _console_logging() -> Iterator[None]:
    """Show the package's log records, from INFO up, on standard error."""
    package = logging.getLogger("crispr_analysis_utils")
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S")
    )
    level = package.level
    package.setLevel(logging.INFO)
    package.addHandler(handler)
    try:
        yield
    finally:
        package.removeHandler(handler)
        package.setLevel(level)


def _positive_int(text: str) -> int:
    value = _integer(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"expected a positive integer, got {text!r}")
    return value


def _non_negative_int(text: str) -> int:
    value = _integer(text)
    if value < 0:
        raise argparse.ArgumentTypeError(
            f"expected a non-negative integer, got {text!r}"
        )
    return value


def _integer(text: str) -> int:
    try:
        return int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected an integer, got {text!r}") from None


def _pam_list(text: str) -> list[str]:
    """Comma-separated PAMs; an empty string for none."""
    if not text.strip():
        return []
    return _items(text, "PAM")


def _name_list(text: str) -> list[str]:
    """Comma-separated contig names, at least one."""
    return _items(text, "contig name")


def _items(text: str, what: str) -> list[str]:
    items = [item.strip() for item in text.split(",")]
    if not all(items):
        raise argparse.ArgumentTypeError(f"empty {what} in {text!r}")
    return items


def _column(text: str) -> int | str:
    """A 0-based column position when all digits, else a column name."""
    return int(text) if text.isascii() and text.isdigit() else text


def _separator(text: str) -> str:
    separator = _SEPARATORS.get(text.lower(), text)
    if len(separator) != 1:
        raise argparse.ArgumentTypeError(
            f"expected tab, comma or one character, got {text!r}"
        )
    return separator

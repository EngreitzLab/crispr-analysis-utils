"""``cau inference-input``: one MuData an inference reads, in the pipeline's layout."""

from __future__ import annotations

import argparse

from . import CommandError


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Declare the ``run`` action and its arguments."""
    actions = parser.add_subparsers(dest="action", metavar="ACTION", required=True)

    run_parser = actions.add_parser(
        "run",
        help="Join QC'd channels with guide calls into one inference object.",
        description=(
            "Concatenate per-channel MuData files that carry guide calls in"
            " mod['guide'].layers['guide_assignment'], key the guides and"
            " bucket the controls as the IGVF CRISPR pipeline does, and attach"
            " the pairs to test. --outdir gets inference_input/<name>.h5mu,"
            " plus run_metadata/ and logs/."
        ),
    )
    run_parser.add_argument(
        "--inputs",
        required=True,
        nargs="+",
        help="QC'd per-channel .h5mu files with guide calls, one object's worth.",
    )
    run_parser.add_argument(
        "--pairs",
        required=True,
        help="pairs_to_test.csv from `cau pairs run`.",
    )
    run_parser.add_argument(
        "--outdir",
        required=True,
        help="Where the object goes. Nothing is written outside it.",
    )
    run_parser.add_argument(
        "--name",
        required=True,
        help="The object's file stem, such as day4.",
    )
    run_parser.add_argument(
        "--capture-method",
        required=True,
        help="How guides were captured, as the pipeline records it (CROP-seq, direct capture).",
    )
    run_parser.add_argument(
        "--moi",
        choices=("high", "low"),
        default="high",
        help="Multiplicity of infection (default: high).",
    )
    run_parser.add_argument(
        "--min-guides",
        type=int,
        default=None,
        help="Keep cells with at least this many assigned guides (default: keep all).",
    )
    run_parser.add_argument(
        "--max-guides",
        type=int,
        default=None,
        help="Keep cells with at most this many assigned guides (default: keep all).",
    )


def run(args: argparse.Namespace) -> int:
    """Run the ``run`` action."""
    from ..inference_input import InputError, _report_lines
    from ..inference_input import run as build

    try:
        result = build(
            args.inputs,
            args.pairs,
            args.outdir,
            name=args.name,
            moi=args.moi,
            capture_method=args.capture_method,
            min_guides=args.min_guides,
            max_guides=args.max_guides,
        )
    except (FileNotFoundError, InputError, ValueError, OSError) as error:
        raise CommandError(str(error)) from error

    for line in _report_lines(result.report):
        print(line)
    print(f"Wrote {result.path}")
    return 0

"""``cau guide-alignment``: argument parsing, and the ``run`` action end to end."""

import json

import pytest
from guide_alignment_fixtures import SPACER_A, SPACER_C, Contig, write_fasta

from crispr_analysis_utils import cli
from crispr_analysis_utils.guide_alignment.iupac import reverse_complement

REQUIRED = [
    "--guides",
    "guides.tsv",
    "--reference",
    "ref.fa",
    "--outdir",
    "out",
]


def parse(*argv):
    return cli.build_parser().parse_args(["guide-alignment", *argv])


def plant(tmp_path, *, reversed_guides=False):
    """A one-contig reference with SPACER_A in it, and a guide table beside it."""
    contig = Contig()
    contig.add("ACGTTTGACCATTTACGAGT")
    start = contig.add(SPACER_A + "TGG")
    contig.add("TTGCAATCCGGATTACATGA")
    (tmp_path / "ref").mkdir(exist_ok=True)
    # A second, empty contig so a test can check that a comma-separated list
    # with spaces around it is split properly.
    fasta = write_fasta(
        tmp_path / "ref" / "genome.fa",
        {"chr1": contig.sequence, "chr2": "ACGT" * 10},
    )
    sequence = reverse_complement(SPACER_A) if reversed_guides else SPACER_A
    guides = tmp_path / "guides.csv"
    guides.write_text(f"id,sequence\ng1,{sequence}\n")
    return fasta, guides, start


def test_guide_alignment_is_registered():
    assert cli.COMMANDS["guide-alignment"][0] == "guide_alignment"


def test_parse_run_defaults():
    args = parse("run", *REQUIRED)
    assert (args.command, args.action) == ("guide-alignment", "run")
    assert (args.guides, args.reference, args.outdir) == (
        "guides.tsv",
        "ref.fa",
        "out",
    )
    assert args.pam == "NGG"
    assert args.spacer_length == 20
    assert args.add_leading_g is False
    assert args.orientation_check == "error"
    assert args.fix_orientation is False
    assert args.contigs is None
    # The columns stay strings here: `run` turns a digit into a position.
    assert (args.id_col, args.spacer_col, args.sep, args.header) == (
        "0",
        "1",
        "\t",
        False,
    )


def test_parse_run_every_flag():
    args = parse(
        "run",
        *REQUIRED,
        "--pam",
        "NRG",
        "--spacer-length",
        "19",
        "--add-leading-g",
        "--orientation-check",
        "warn",
        "--fix-orientation",
        "--contigs",
        "chr1, chr2",
        "--id-col",
        "id",
        "--spacer-col",
        "2",
        "--sep",
        ",",
        "--header",
    )
    assert args.pam == "NRG"
    assert args.spacer_length == 19
    assert args.add_leading_g is True
    assert args.orientation_check == "warn"
    assert args.fix_orientation is True
    assert args.contigs == "chr1, chr2"
    assert (args.id_col, args.spacer_col, args.sep, args.header) == (
        "id",
        "2",
        ",",
        True,
    )


@pytest.mark.parametrize(
    "argv",
    [
        ["guide-alignment"],
        ["guide-alignment", "align"],
        # The GEM index is gone, and so is its action.
        ["guide-alignment", "index", "ref.fa", "idx/genome"],
        ["guide-alignment", "run", *REQUIRED[:-2]],
        ["guide-alignment", "run", *REQUIRED, "--spacer-length", "twenty"],
        ["guide-alignment", "run", *REQUIRED, "--orientation-check", "maybe"],
        # Options the Rust scan does away with.
        ["guide-alignment", "run", *REQUIRED, "--index", "idx/genome.gem"],
        ["guide-alignment", "run", *REQUIRED, "--threads", "4"],
        ["guide-alignment", "run", *REQUIRED, "--max-mismatches", "2"],
        ["guide-alignment", "run", *REQUIRED, "--alt-pams", "NAG"],
    ],
)
def test_bad_arguments_are_usage_errors(argv):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(argv)
    assert exit_info.value.code == 2


def test_run_end_to_end(tmp_path, capsys):
    fasta, guides, start = plant(tmp_path)

    status = cli.main(
        [
            "guide-alignment",
            "run",
            "--guides",
            str(guides),
            "--reference",
            str(fasta),
            "--outdir",
            "out",
            "--sep",
            ",",
            "--header",
            "--id-col",
            "id",
            "--spacer-col",
            "sequence",
            "--contigs",
            "chr1, chr2",
        ]
    )

    assert status == 0
    captured = capsys.readouterr()
    assert captured.out.splitlines() == [
        "1 guides: 1 unique, 0 multi, 0 no_site",
        "Outputs in out",
    ]
    # Nothing is logged to the console: the run's log is a file in outdir.
    assert captured.err == ""
    assert "Done in" in (tmp_path / "out" / "logs" / "guide_alignment.log").read_text()
    rows = (
        (tmp_path / "out" / "alignments" / "guide_mapping.bed").read_text().splitlines()
    )
    # The BED's name field is the spacer, not the alias the table gives it.
    assert rows == [f"chr1\t{start}\t{start + 20}\t{SPACER_A}\t0\t+"]
    records = (tmp_path / "out" / "run_metadata").glob("guide_alignment_*.json")
    state = json.loads(sorted(records)[-1].read_text())
    # --contigs is split and stripped, and a digit column stays a name here.
    assert state["parameters"]["contigs"] == ["chr1", "chr2"]
    assert state["parameters"]["id_col"] == "id"
    assert state["parameters"]["sep"] == ","


def test_column_positions_are_read_as_numbers(tmp_path, capsys):
    fasta, _, start = plant(tmp_path)
    guides = tmp_path / "guides.tsv"
    guides.write_text(f"note\t{SPACER_A}\tg1\n")

    status = cli.main(
        [
            "guide-alignment",
            "run",
            "--guides",
            str(guides),
            "--reference",
            str(fasta),
            "--outdir",
            "out",
            "--id-col",
            "2",
            "--spacer-col",
            "1",
        ]
    )

    assert status == 0
    assert capsys.readouterr().out.startswith("1 guides: 1 unique")
    rows = (
        (tmp_path / "out" / "alignments" / "guide_mapping.bed").read_text().splitlines()
    )
    # The BED's name field is the spacer, not the alias the table gives it.
    assert rows == [f"chr1\t{start}\t{start + 20}\t{SPACER_A}\t0\t+"]


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["--contigs", ""], "--contigs was given no contig names"),
        (["--pam", "NGZ"], "not IUPAC nucleotide codes"),
        (["--spacer-length", "0"], "positive integer"),
    ],
)
def test_run_errors_are_reported_without_a_traceback(tmp_path, capsys, argv, message):
    fasta, guides, _ = plant(tmp_path)

    with pytest.raises(SystemExit) as exit_info:
        cli.main(
            [
                "guide-alignment",
                "run",
                "--guides",
                str(guides),
                "--reference",
                str(fasta),
                "--outdir",
                "out",
                "--sep",
                ",",
                "--header",
                "--id-col",
                "id",
                "--spacer-col",
                "sequence",
                *argv,
            ]
        )

    assert exit_info.value.code == 1
    error = capsys.readouterr().err
    assert error.startswith("cau: error: ")
    assert message in error
    assert "Traceback" not in error


def test_a_missing_reference_is_reported_without_a_traceback(tmp_path, capsys):
    (tmp_path / "guides.tsv").write_text(f"g1\t{SPACER_A}\n")

    with pytest.raises(SystemExit) as exit_info:
        cli.main(
            [
                "guide-alignment",
                "run",
                "--guides",
                "guides.tsv",
                "--reference",
                "missing.fa",
                "--outdir",
                "out",
            ]
        )

    assert exit_info.value.code == 1
    error = capsys.readouterr().err
    assert "cau: error: The reference" in error
    assert "Traceback" not in error


def test_a_reversed_library_is_reported_without_a_traceback(tmp_path, capsys):
    fasta, guides, _ = plant(tmp_path, reversed_guides=True)

    with pytest.raises(SystemExit) as exit_info:
        cli.main(
            [
                "guide-alignment",
                "run",
                "--guides",
                str(guides),
                "--reference",
                str(fasta),
                "--outdir",
                "out",
                "--sep",
                ",",
                "--header",
                "--id-col",
                "id",
                "--spacer-col",
                "sequence",
            ]
        )

    assert exit_info.value.code == 1
    error = capsys.readouterr().err
    assert "cau: error: 1 of 1 entries read 3' to 5'" in error
    assert "--fix-orientation" in error
    assert "Traceback" not in error


def test_fix_orientation_turns_a_reversed_entry_round(tmp_path, capsys):
    fasta, guides, start = plant(tmp_path, reversed_guides=True)

    status = cli.main(
        [
            "guide-alignment",
            "run",
            "--guides",
            str(guides),
            "--reference",
            str(fasta),
            "--outdir",
            "out",
            "--sep",
            ",",
            "--header",
            "--id-col",
            "id",
            "--spacer-col",
            "sequence",
            "--orientation-check",
            "warn",
            "--fix-orientation",
        ]
    )

    assert status == 0
    assert capsys.readouterr().out.startswith("1 guides: 1 unique")
    rows = (
        (tmp_path / "out" / "alignments" / "guide_mapping.bed").read_text().splitlines()
    )
    # The BED's name field is the spacer, not the alias the table gives it.
    assert rows == [f"chr1\t{start}\t{start + 20}\t{SPACER_A}\t0\t+"]


def test_the_sep_is_taken_literally(tmp_path, capsys):
    fasta, _, _ = plant(tmp_path)
    guides = tmp_path / "guides.psv"
    guides.write_text(f"id|sequence\ng1|{SPACER_A}\ng2|{SPACER_C}\n")

    status = cli.main(
        [
            "guide-alignment",
            "run",
            "--guides",
            str(guides),
            "--reference",
            str(fasta),
            "--outdir",
            "out",
            "--sep",
            "|",
            "--header",
            "--id-col",
            "id",
            "--spacer-col",
            "sequence",
        ]
    )

    assert status == 0
    assert capsys.readouterr().out.startswith("2 guides: 1 unique")

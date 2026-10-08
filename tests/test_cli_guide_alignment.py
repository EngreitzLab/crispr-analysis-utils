"""``cau guide-alignment``: argument parsing, and both actions on the fake GEM."""

import pytest
from guide_alignment_fixtures import SPACER_A, Contig, write_fasta

from crispr_analysis_utils import cli

REQUIRED = [
    "--guides",
    "guides.tsv",
    "--reference",
    "ref.fa",
    "--index",
    "idx/genome.gem",
    "--outdir",
    "out",
]


def parse(*argv):
    return cli.build_parser().parse_args(["guide-alignment", *argv])


def test_guide_alignment_is_registered():
    assert cli.COMMANDS["guide-alignment"][0] == "guide_alignment"


def test_parse_index():
    args = parse("index", "ref.fa", "idx/genome")
    assert (args.command, args.action) == ("guide-alignment", "index")
    assert (args.reference, args.prefix, args.threads) == ("ref.fa", "idx/genome", 8)
    assert parse("index", "ref.fa", "idx/genome", "--threads", "2").threads == 2


def test_parse_run_defaults():
    args = parse("run", *REQUIRED)
    assert (args.command, args.action) == ("guide-alignment", "run")
    assert (args.guides, args.reference, args.index, args.outdir) == (
        "guides.tsv",
        "ref.fa",
        "idx/genome.gem",
        "out",
    )
    assert args.max_mismatches == 3
    assert args.pam == "NGG"
    assert args.alt_pams == ["NAG", "NGA"]
    assert args.spacer_length == 20
    assert args.add_leading_g is False
    assert args.orientation_check == "warn"
    assert args.contigs is None
    assert args.threads == 8
    assert (args.id_col, args.spacer_col, args.sep, args.header) == (0, 1, "\t", False)


def test_parse_run_every_flag():
    args = parse(
        "run",
        *REQUIRED,
        "--max-mismatches",
        "2",
        "--pam",
        "NRG",
        "--alt-pams",
        "",
        "--spacer-length",
        "19",
        "--add-leading-g",
        "--orientation-check",
        "error",
        "--contigs",
        "chr1, chr2",
        "--threads",
        "4",
        "--id-col",
        "id",
        "--spacer-col",
        "2",
        "--sep",
        ",",
        "--header",
    )
    assert args.max_mismatches == 2
    assert args.pam == "NRG"
    assert args.alt_pams == []
    assert args.spacer_length == 19
    assert args.add_leading_g is True
    assert args.orientation_check == "error"
    assert args.contigs == ["chr1", "chr2"]
    assert args.threads == 4
    assert (args.id_col, args.spacer_col, args.sep, args.header) == ("id", 2, ",", True)


@pytest.mark.parametrize("given", ["tab", "TAB", "\\t", "\t"])
def test_parse_run_tab_separator(given):
    assert parse("run", *REQUIRED, "--sep", given).sep == "\t"


@pytest.mark.parametrize(
    "argv",
    [
        ["guide-alignment"],
        ["guide-alignment", "align"],
        ["guide-alignment", "index", "ref.fa"],
        ["guide-alignment", "run", *REQUIRED[:-2]],
        ["guide-alignment", "run", *REQUIRED, "--alt-pams", "NAG,,NGA"],
        ["guide-alignment", "run", *REQUIRED, "--contigs", ""],
        ["guide-alignment", "run", *REQUIRED, "--threads", "0"],
        ["guide-alignment", "run", *REQUIRED, "--max-mismatches", "-1"],
        ["guide-alignment", "run", *REQUIRED, "--spacer-length", "twenty"],
        ["guide-alignment", "run", *REQUIRED, "--orientation-check", "maybe"],
        ["guide-alignment", "run", *REQUIRED, "--sep", ";;"],
    ],
)
def test_bad_arguments_are_usage_errors(argv):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(argv)
    assert exit_info.value.code == 2


def test_index_then_run(tmp_path, fake_gem, capsys):
    contig = Contig()
    contig.add("ACGTTTGACCATTTACGAGT")
    start = contig.add(SPACER_A + "TGG")
    contig.add("TTGCAATCCGGATTACATGA")
    (tmp_path / "ref").mkdir()
    fasta = write_fasta(tmp_path / "ref" / "genome.fa", {"chr1": contig.sequence})
    guides = tmp_path / "guides.csv"
    guides.write_text(f"id,sequence\ng1,{SPACER_A}\n")

    assert cli.main(["guide-alignment", "index", str(fasta), "idx/genome"]) == 0
    index = tmp_path / "idx" / "genome.gem"
    assert capsys.readouterr().out.strip() == str(index.resolve())
    assert index.with_name("genome.gem.json").is_file()

    status = cli.main(
        [
            "guide-alignment",
            "run",
            "--guides",
            str(guides),
            "--reference",
            str(fasta),
            "--index",
            "idx/genome",
            "--outdir",
            "out",
            "--sep",
            "comma",
            "--header",
            "--id-col",
            "id",
            "--spacer-col",
            "sequence",
            "--threads",
            "1",
        ]
    )
    assert status == 0
    captured = capsys.readouterr()
    assert captured.out.splitlines() == [
        "1 guides: 1 unique, 0 multi, 0 off_target_only, 0 no_site",
        "Outputs in out",
    ]
    assert "Done in" in captured.err
    rows = (tmp_path / "out" / "guides.bed").read_text().splitlines()
    assert rows == [f"chr1\t{start}\t{start + 20}\tg1\t0\t+"]


def test_run_errors_are_reported_without_a_traceback(tmp_path, fake_gem, capsys):
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
                "--index",
                "missing",
                "--outdir",
                "out",
            ]
        )
    assert exit_info.value.code == 1
    error = capsys.readouterr().err
    assert "cau: error: The reference" in error
    assert "Traceback" not in error

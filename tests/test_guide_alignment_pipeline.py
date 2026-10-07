"""The guide alignment pipeline, run on a small reference with the fake GEM."""

import csv
import json
import logging
import random

import fake_gem as fake_gem_tools
import pandas as pd
import pytest
from guide_alignment_fixtures import SPACER_A, SPACER_C, SPACER_G, Contig, write_fasta

from crispr_analysis_utils.guide_alignment import (
    GemError,
    OrientationError,
    build_index,
    run,
)
from crispr_analysis_utils.guide_alignment.iupac import reverse_complement
from crispr_analysis_utils.guide_alignment.pipeline import (
    PASS_ORIENTATION,
    PASS_SITES,
)

SPACER_NT = "TTAGCGCATAGCGTACGTAA"
GUIDES = [
    ("g_a", SPACER_A),
    ("g_c", SPACER_C),
    ("g_g", SPACER_G),
    ("g_nt", SPACER_NT),
    ("g_dup", SPACER_A),
]


def mutate(spacer, position):
    """Change the base at this 1-based position."""
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    return spacer[: position - 1] + swap[spacer[position - 1]] + spacer[position:]


def build_reference(path):
    """Two contigs: g_a once (+), g_c once (-) plus a 1-mismatch NAG site,
    g_g once on each contig. Every site sits next to TTT (or AAA), so no
    flank of a reverse-complemented spacer reads as a PAM.
    """
    rng = random.Random(7)

    def filler(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    sites = {}
    chr1 = Contig()
    chr1.add(filler(40) + "TTT")
    sites["g_a"] = ("chr1", chr1.add(SPACER_A + "TGG"), "+")
    chr1.add(filler(40))
    sites["g_c"] = ("chr1", chr1.add_minus(SPACER_C + "AGG") + 3, "-")
    chr1.add("TTT" + filler(40) + "TTT")
    sites["g_c_nag"] = ("chr1", chr1.add(mutate(SPACER_C, 10) + "TAG"), "+")
    chr1.add(filler(40) + "TTT")
    sites["g_g_1"] = ("chr1", chr1.add(SPACER_G + "CGG"), "+")
    chr1.add(filler(40))
    chr2 = Contig()
    chr2.add(filler(30) + "TTT")
    sites["g_g_2"] = ("chr2", chr2.add(SPACER_G + "AGG"), "+")
    chr2.add(filler(30))
    path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(path, {"chr1": chr1.sequence, "chr2": chr2.sequence})
    return sites


def write_guides(path, guides=GUIDES):
    path.write_text("".join(f"{gid}\t{spacer}\n" for gid, spacer in guides))
    return path


@pytest.fixture
def setup(tmp_path, fake_gem):
    reference = tmp_path / "ref" / "genome.fa"
    sites = build_reference(reference)
    index = build_index(reference, tmp_path / "index" / "genome")
    guides = write_guides(tmp_path / "guides.tsv")
    return {
        "reference": reference,
        "index": index,
        "guides": guides,
        "sites": sites,
        "calls": fake_gem,
        "outdir": tmp_path / "out",
    }


def go(setup, **kwargs):
    arguments = {"threads": 1, **kwargs}
    guides = arguments.pop("guides", setup["guides"])
    return run(guides, setup["reference"], setup["index"], setup["outdir"], **arguments)


def table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def by_id(rows):
    return {row["guide_id"]: row for row in rows}


def mapper_calls(setup):
    return [
        call for call in fake_gem_tools.calls(setup["calls"]) if call[0] == "gem-mapper"
    ]


def run_json(setup):
    return json.loads((setup["outdir"] / "run.json").read_text())


def test_run_writes_every_output(setup):
    result = go(setup)
    out = setup["outdir"]
    for name in (
        "guides.tsv",
        "sites.tsv",
        "orientation.tsv",
        "guides.bed",
        "cut_sites.bed",
        "summary.tsv",
        "run.json",
        "logs/guide_alignment.log",
        "logs/gem-mapper.orientation.log",
        "logs/gem-mapper.sites.log",
        "alignments/orientation.fastq",
        "alignments/orientation.sam.gz",
        "alignments/sites.fastq",
        "alignments/sites.sam.gz",
    ):
        assert (out / name).is_file(), name
    assert not (out / "alignments" / "sites_leading_g.fastq").exists()

    guides = by_id(table(out / "guides.tsv"))
    chrom, a_start, _ = setup["sites"]["g_a"]
    assert guides["g_a"]["alignment_class"] == "unique"
    assert (guides["g_a"]["guide_chr"], guides["g_a"]["guide_start"]) == (
        chrom,
        str(a_start),
    )
    assert guides["g_a"]["strand"] == "+"
    assert guides["g_a"]["pam"] == "TGG"
    assert guides["g_a"]["cut_site"] == str(a_start + 17)
    _, c_start, _ = setup["sites"]["g_c"]
    assert guides["g_c"]["alignment_class"] == "unique"
    assert (guides["g_c"]["guide_start"], guides["g_c"]["strand"]) == (
        str(c_start),
        "-",
    )
    assert guides["g_c"]["cut_site"] == str(c_start + 3)
    assert guides["g_c"]["n_alt_pam_1mm"] == "1"
    assert guides["g_g"]["alignment_class"] == "multi"
    assert guides["g_g"]["n_ngg_0mm"] == "2"
    assert guides["g_nt"]["alignment_class"] == "no_site"
    assert guides["g_dup"]["duplicate_of"] == "g_a"
    assert guides["g_dup"]["guide_start"] == str(a_start)
    assert {g: row["orientation"] for g, row in guides.items()} == {
        "g_a": "ok",
        "g_c": "ok",
        "g_g": "ok",
        "g_nt": "no_perfect_hit",
        "g_dup": "ok",
    }
    assert [s.alignment_class for s in result.summaries] == [
        "unique",
        "unique",
        "multi",
        "no_site",
        "unique",
    ]

    bed = (out / "cut_sites.bed").read_text().splitlines()
    assert f"chr1\t{a_start + 17}\t{a_start + 18}\tg_a\t0\t+" in bed
    assert len(bed) == 3

    summary = {
        (r["section"], r["name"]): r["value"] for r in table(out / "summary.tsv")
    }
    assert summary[("class", "unique")] == "3"
    assert summary[("site_search", "n_out_of_range")] == "0"
    assert summary[("parameter", "sep")] == "\\t"
    assert summary[("parameter", "gem_version")] == fake_gem_tools.VERSION

    record = run_json(setup)
    assert record["complete"] is True and record["error"] is None
    assert record["gem_version"] == fake_gem_tools.VERSION
    assert set(record["passes"]) == {PASS_ORIENTATION, PASS_SITES}
    assert record["passes"][PASS_SITES]["reused"] is False
    assert record["passes"][PASS_SITES]["n_reads"] == 4 * 3
    assert (
        record["inputs"]["reference"]["md5"]
        == record["inputs"]["index"]["reference_md5"]
    )
    assert record["classes"] == {
        "unique": 3,
        "multi": 1,
        "off_target_only": 0,
        "no_site": 1,
    }
    assert record["alignments_per_guide"]["g_g"]["orientation"] == 2
    assert "g_dup" not in record["alignments_per_guide"]
    assert result.passes.keys() == record["passes"].keys()


def test_run_uses_the_tuned_budget_for_each_pass(setup):
    go(setup, add_leading_g=True)
    budgets = {
        call[call.index("-i") + 1].rsplit("/", 1)[-1]: [
            arg
            for arg in call
            if arg.startswith(("--alignment-max-error", "--complete-"))
        ]
        for call in mapper_calls(setup)
    }
    assert budgets == {
        "orientation.fastq": [
            "--alignment-max-error=0",
            "--complete-search-error=0",
            "--complete-strata-after-best=0",
        ],
        "sites.fastq": [
            "--alignment-max-error=3",
            "--complete-search-error=4",
            "--complete-strata-after-best=4",
        ],
        "sites_leading_g.fastq": [
            "--alignment-max-error=4",
            "--complete-search-error=5",
            "--complete-strata-after-best=5",
        ],
    }
    sites = table(setup["outdir"] / "sites.tsv")
    g_bases = {
        (row["guide_id"], row["start"]): row["leading_g_base"]
        for row in sites
        if row["pam_class"] == "NGG"
    }
    _, a_start, _ = setup["sites"]["g_a"]
    _, c_start, _ = setup["sites"]["g_c"]
    # TTT precedes g_a on +; TTT follows g_c's site on +, so A in guide orientation.
    assert g_bases[("g_a", str(a_start))] == "T"
    assert g_bases[("g_c", str(c_start))] == "A"
    assert all(base == "" for (guide, _), base in g_bases.items() if guide == "g_g")


def test_rerun_reuses_every_pass(setup):
    first = go(setup)
    n_calls = len(mapper_calls(setup))
    second = go(setup)
    assert len(mapper_calls(setup)) == n_calls
    record = run_json(setup)
    assert record["complete"] is True
    assert all(entry["reused"] for entry in record["passes"].values())
    assert [s.alignment_class for s in second.summaries] == [
        s.alignment_class for s in first.summaries
    ]
    assert record["alignments_per_guide"] == {
        gid: entry for gid, entry in record["alignments_per_guide"].items()
    }
    log = (setup["outdir"] / "logs" / "guide_alignment.log").read_text()
    assert log.count("Guide alignment into") == 2
    assert "reusing" in log


def test_rerun_maps_again_only_the_pass_whose_options_changed(setup):
    go(setup)
    go(setup, max_mismatches=2)
    passes = run_json(setup)["passes"]
    assert passes[PASS_ORIENTATION]["reused"] is True
    assert passes[PASS_SITES]["reused"] is False
    assert "--alignment-max-error=2" in passes[PASS_SITES]["argv"]


def test_rerun_maps_again_when_the_alignments_changed(setup):
    go(setup)
    sam = setup["outdir"] / "alignments" / "sites.sam.gz"
    sam.write_bytes(sam.read_bytes() + b"\n")
    go(setup)
    passes = run_json(setup)["passes"]
    assert passes[PASS_ORIENTATION]["reused"] is True
    assert passes[PASS_SITES]["reused"] is False


def test_failed_run_keeps_the_earlier_passes_reusable(setup):
    go(setup)
    n_calls = len(mapper_calls(setup))
    with pytest.raises(ValueError, match="chrUn"):
        go(setup, contigs=["chrUn"])
    record = run_json(setup)
    assert record["complete"] is False
    assert record["error"].startswith("ValueError")
    assert all(entry["carried_over"] for entry in record["passes"].values())
    go(setup)
    assert len(mapper_calls(setup)) == n_calls
    record = run_json(setup)
    assert all(entry["reused"] for entry in record["passes"].values())
    assert not any("carried_over" in entry for entry in record["passes"].values())


def test_contigs_keep_only_their_hits(setup):
    go(setup, contigs=["chr1"])
    guides = by_id(table(setup["outdir"] / "guides.tsv"))
    assert guides["g_g"]["alignment_class"] == "unique"
    assert guides["g_g"]["guide_chr"] == "chr1"


def test_reversed_library_stops_the_run_when_the_check_is_error(setup, tmp_path):
    reversed_guides = write_guides(
        tmp_path / "reversed.tsv",
        [(gid, reverse_complement(spacer)) for gid, spacer in GUIDES[:3]],
    )
    with pytest.raises(OrientationError, match="3 against 0"):
        go(setup, guides=reversed_guides, orientation_check="error")
    out = setup["outdir"]
    verdicts = {
        row["guide_id"]: row["verdict"] for row in table(out / "orientation.tsv")
    }
    assert verdicts == {"g_a": "reversed", "g_c": "reversed", "g_g": "reversed"}
    assert not (out / "alignments" / "sites.fastq").exists()
    assert not (out / "guides.tsv").exists()
    record = run_json(setup)
    assert record["complete"] is False
    assert record["error"].startswith("OrientationError")


def test_reversed_library_is_only_logged_when_the_check_is_warn(
    setup, tmp_path, caplog
):
    reversed_guides = write_guides(
        tmp_path / "reversed.tsv",
        [(gid, reverse_complement(spacer)) for gid, spacer in GUIDES[:3]],
    )
    with caplog.at_level(logging.WARNING, logger="crispr_analysis_utils"):
        result = go(setup, guides=reversed_guides)
    assert "Reversed guides outnumber ok ones" in caplog.text
    assert {s.orientation for s in result.summaries} == {"reversed"}
    assert (setup["outdir"] / "guides.tsv").is_file()


def test_orientation_check_off_skips_pass_1(setup):
    go(setup)
    result = go(setup, orientation_check="off")
    out = setup["outdir"]
    assert not (out / "orientation.tsv").exists()
    assert not (out / "alignments" / "orientation.sam.gz").exists()
    assert result.orientation is None
    assert {s.orientation for s in result.summaries} == {"not_checked"}
    record = run_json(setup)
    assert set(record["passes"]) == {PASS_SITES}
    # One record per PAM read, all at g_a's one site.
    assert record["alignments_per_guide"]["g_a"] == {"sites": 3}


def test_reference_must_match_the_index_sidecar(setup, tmp_path):
    other = tmp_path / "other" / "genome.fa"
    other.parent.mkdir()
    other.write_text(setup["reference"].read_text().replace("TTT", "TTA", 1))
    with pytest.raises(ValueError, match="was built from a reference with MD5"):
        run(setup["guides"], other, setup["index"], setup["outdir"], threads=1)
    assert mapper_calls(setup) == []


def test_index_without_a_sidecar_is_used_with_a_warning(setup, caplog):
    setup["index"].with_name("genome.gem.json").unlink()
    with caplog.at_level(logging.WARNING, logger="crispr_analysis_utils"):
        go(setup)
    assert "has no sidecar" in caplog.text
    assert run_json(setup)["inputs"]["index"]["reference_md5"] is None


def test_index_can_be_given_by_its_prefix(setup):
    result = run(
        setup["guides"],
        setup["reference"],
        setup["index"].with_suffix(""),
        setup["outdir"],
        threads=1,
    )
    assert result.summaries[0].alignment_class == "unique"


def test_nothing_is_written_outside_outdir(setup, tmp_path):
    def snapshot():
        return {
            path.relative_to(tmp_path)
            for path in tmp_path.rglob("*")
            if not path.is_relative_to(setup["outdir"])
            and path.name != "fake-gem-calls.jsonl"
        }

    before = snapshot()
    go(setup)
    assert snapshot() == before
    # No .fai next to the reference: the pipeline built its own in outdir.
    assert (setup["outdir"] / "genome.fa.fai").is_file()
    assert run_json(setup)["inputs"]["fasta_index"]["built_in_outdir"] is True


def test_the_reference_index_is_used_when_it_exists(setup):
    import pysam

    pysam.faidx(str(setup["reference"]))
    go(setup)
    assert not (setup["outdir"] / "genome.fa.fai").exists()
    assert run_json(setup)["inputs"]["fasta_index"]["built_in_outdir"] is False


def test_records_outside_their_contig_are_dropped_and_counted(
    setup, monkeypatch, tmp_path
):
    clean = run(
        setup["guides"],
        setup["reference"],
        setup["index"],
        tmp_path / "clean",
        threads=1,
    )
    monkeypatch.setenv("FAKE_GEM_OUT_OF_RANGE", "1")
    go(setup)
    assert table(setup["outdir"] / "guides.tsv") == table(clean.outdir / "guides.tsv")
    assert table(setup["outdir"] / "sites.tsv") == table(clean.outdir / "sites.tsv")
    summary = {
        (r["section"], r["name"]): r["value"]
        for r in table(setup["outdir"] / "summary.tsv")
    }
    assert int(summary[("site_search", "n_out_of_range")]) > 0
    assert int(summary[("orientation_check", "n_out_of_range")]) > 0


def test_a_gem_crash_names_the_fastq_and_the_guide(setup, monkeypatch):
    monkeypatch.setenv("FAKE_GEM_CRASH", "1")
    with pytest.raises(GemError) as error:
        go(setup, add_leading_g=True, orientation_check="off")
    message = str(error.value)
    # The fake crashes on the first read of the first FASTQ mapped: g_g's.
    assert "sites.fastq" in message
    assert "read '2:NGG' (guide 'g_g')" in message
    assert "Signal raised" in message
    assert not (setup["outdir"] / "alignments" / "sites.sam.gz").exists()
    record = run_json(setup)
    assert record["complete"] is False
    assert record["error"].startswith("GemError")


def test_a_gem_crash_on_a_leading_g_read_says_so(setup, monkeypatch, tmp_path):
    guides = write_guides(tmp_path / "a_only.tsv", [("g_a", SPACER_A)])
    monkeypatch.setenv("FAKE_GEM_CRASH", "1")
    with pytest.raises(
        GemError, match=r"read '0:NGG' \(guide 'g_a', with the added G\)"
    ) as error:
        go(setup, guides=guides, add_leading_g=True, orientation_check="off")
    assert "sites_leading_g.fastq" in str(error.value)


def test_a_dataframe_library_is_accepted(setup):
    frame = pd.DataFrame(GUIDES, columns=["id", "spacer"])
    result = go(setup, guides=frame, id_col="id", spacer_col="spacer")
    assert [s.guide_id for s in result.summaries] == [gid for gid, _ in GUIDES]
    assert run_json(setup)["inputs"]["guides"] == {"path": None, "md5": None}


@pytest.mark.parametrize(
    ("kwargs", "error"),
    [
        ({"orientation_check": "maybe"}, ValueError),
        ({"alt_pams": "NAG"}, TypeError),
        ({"alt_pams": ["NAGG"]}, ValueError),
        ({"threads": 0}, ValueError),
        ({"contigs": "chr1"}, TypeError),
        ({"max_mismatches": -1}, ValueError),
    ],
)
def test_run_rejects_bad_parameters(setup, kwargs, error):
    with pytest.raises(error):
        go(setup, **kwargs)
    assert mapper_calls(setup) == []
    assert not setup["outdir"].exists()


def test_missing_inputs_are_reported(setup, tmp_path):
    with pytest.raises(FileNotFoundError, match="GEM index"):
        run(setup["guides"], setup["reference"], tmp_path / "nope.gem", setup["outdir"])
    with pytest.raises(FileNotFoundError, match="reference"):
        run(setup["guides"], tmp_path / "nope.fa", setup["index"], setup["outdir"])

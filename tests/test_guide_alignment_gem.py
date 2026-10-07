import gzip
import json
import os

import fake_gem as fake_gem_tools
import pysam
import pytest
from guide_alignment_fixtures import SPACER_A, SPACER_C, Contig, write_fasta

from crispr_analysis_utils.guide_alignment import gem
from crispr_analysis_utils.guide_alignment._files import file_md5

BASE = [
    "--mapping-mode=customed",
    "--alignment-model=hamming",
    "--clipping=none",
    "--alignment-local=never",
    "--sam-compact=false",
    "--max-reported-matches=all",
    "--alignment-max-bandwidth=0",
]


def budget(errors, search):
    return [
        f"--alignment-max-error={errors}",
        f"--complete-search-error={search}",
        f"--complete-strata-after-best={search}",
    ]


@pytest.mark.parametrize(
    ("kwargs", "expected"),
    [
        # Pass 1: exact matches of the bare spacer.
        ({"max_mismatches": 0}, budget(0, 0)),
        # Pass 2, NGG/NAG/NGA (one N), with and without the added G.
        ({"max_mismatches": 3, "n_ambiguous": 1}, budget(3, 4)),
        ({"max_mismatches": 3, "n_ambiguous": 1, "leading_g": True}, budget(4, 5)),
        # A PAM such as NRG holds two non-ACGT codes.
        ({"max_mismatches": 3, "n_ambiguous": 2}, budget(3, 5)),
        ({"max_mismatches": 0, "n_ambiguous": 1, "leading_g": True}, budget(1, 2)),
    ],
)
def test_mapper_argv_is_pinned(kwargs, expected):
    argv = gem.mapper_argv("idx/genome.gem", "reads.fastq", threads=8, **kwargs)
    assert argv == [
        "gem-mapper",
        "-I",
        "idx/genome.gem",
        "-i",
        "reads.fastq",
        "-t",
        "8",
        *BASE,
        *expected,
    ]


def test_mapper_options_are_the_argv_without_inputs():
    options = gem.mapper_options(max_mismatches=2, n_ambiguous=1, leading_g=True)
    assert options == [*BASE, *budget(3, 4)]
    assert (
        gem.mapper_argv(
            "i.gem", "r.fq", max_mismatches=2, n_ambiguous=1, leading_g=True
        )[7:]
        == options
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_mismatches": -1},
        {"max_mismatches": 1.5},
        {"max_mismatches": True},
        {"max_mismatches": 3, "n_ambiguous": -1},
        {"max_mismatches": 3, "threads": 0},
    ],
)
def test_mapper_argv_rejects_bad_budgets(kwargs):
    with pytest.raises(ValueError):
        gem.mapper_argv("i.gem", "r.fq", **kwargs)


@pytest.mark.parametrize("given", ["idx/genome", "idx/genome.gem"])
def test_index_file_accepts_the_prefix_or_the_index(given):
    assert gem.index_file(given) == gem.index_file("idx/genome.gem")
    assert str(gem.index_file(given)).endswith("idx/genome.gem")
    assert str(gem.sidecar_path(given)).endswith("idx/genome.gem.json")


def test_index_file_keeps_dots_in_the_prefix():
    assert gem.index_file("hg38.p14").name == "hg38.p14.gem"


def reference(tmp_path):
    contig = Contig()
    contig.add("ACGTACGTAC")
    plus = contig.add(SPACER_A + "TGG")
    contig.add("CCCCCCCCCC")
    minus = contig.add_minus(SPACER_C + "AGG")
    contig.add("TTTTTTTTTT")
    (tmp_path / "ref").mkdir(exist_ok=True)
    path = write_fasta(tmp_path / "ref" / "genome.fa", {"chr1": contig.sequence})
    return path, plus, minus


def write_reads(path, reads):
    with open(path, "w", encoding="ascii") as handle:
        for name, sequence in reads:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    return path


@pytest.fixture
def indexed(tmp_path, fake_gem):
    fasta, plus, minus = reference(tmp_path)
    index = gem.build_index(fasta, tmp_path / "index" / "genome", threads=2)
    reads = write_reads(
        tmp_path / "reads.fastq",
        [("0:NGG", SPACER_A + "NGG"), ("1:NGG", SPACER_C + "NGG"), ("2:NGG", "A" * 23)],
    )
    return {
        "fasta": fasta,
        "index": index,
        "reads": reads,
        "plus": plus,
        "minus": minus,
    }


def test_build_index_writes_the_index_and_its_sidecar(tmp_path, fake_gem, indexed):
    index = indexed["index"]
    assert index == (tmp_path / "index" / "genome.gem").resolve()
    assert index.is_file()
    sidecar = json.loads(gem.sidecar_path(index).read_text())
    assert sidecar == gem.read_sidecar(tmp_path / "index" / "genome")
    assert sidecar["reference"] == str(indexed["fasta"].resolve())
    assert sidecar["reference_md5"] == file_md5(indexed["fasta"])
    assert sidecar["reference_size"] == indexed["fasta"].stat().st_size
    assert sidecar["gem_version"] == fake_gem_tools.VERSION
    ((tool, *args),) = fake_gem_tools.calls(fake_gem)
    assert tool == "gem-indexer"
    assert args[:4] == [
        "-i",
        str(indexed["fasta"].resolve()),
        "-o",
        str(index.with_suffix("")),
    ]
    assert args[-2:] == ["-t", "2"]
    scratch = args[args.index("--tmp-folder") + 1]
    # gem-indexer joins the folder and a file name with no separator.
    assert scratch.endswith(os.sep)
    assert os.path.dirname(scratch.rstrip(os.sep)) == str(index.parent)
    assert not os.path.exists(scratch)
    assert (index.parent / "genome.log").is_file()


def test_build_index_drops_a_trailing_gem_from_the_prefix(tmp_path, fake_gem):
    fasta, _, _ = reference(tmp_path)
    index = gem.build_index(fasta, tmp_path / "idx" / "genome.gem")
    assert index.name == "genome.gem"
    ((_, *args),) = fake_gem_tools.calls(fake_gem)
    assert "-t" not in args


def test_build_index_failure_removes_the_partial_index(tmp_path, fake_gem, monkeypatch):
    fasta, _, _ = reference(tmp_path)
    monkeypatch.setenv("FAKE_GEM_INDEX_FAIL", "1")
    with pytest.raises(gem.GemError, match="gem-indexer failed") as error:
        gem.build_index(fasta, tmp_path / "idx" / "genome")
    assert "GEM::FatalError" in str(error.value)
    assert not (tmp_path / "idx" / "genome.gem").exists()
    assert not (tmp_path / "idx" / "genome.gem.json").exists()


def test_build_index_rejects_a_compressed_reference(tmp_path, fake_gem):
    fasta, _, _ = reference(tmp_path)
    packed = tmp_path / "genome.fa.gz"
    packed.write_bytes(gzip.compress(fasta.read_bytes()))
    with pytest.raises(ValueError, match="gzip-compressed"):
        gem.build_index(packed, tmp_path / "idx" / "genome")
    with pytest.raises(FileNotFoundError):
        gem.build_index(tmp_path / "missing.fa", tmp_path / "idx" / "genome")
    assert fake_gem_tools.calls(fake_gem) == []


def test_read_sidecar_is_none_without_one(tmp_path):
    assert gem.read_sidecar(tmp_path / "genome") is None
    (tmp_path / "genome.gem.json").write_text("[1, 2]")
    with pytest.raises(ValueError, match="not a JSON object"):
        gem.read_sidecar(tmp_path / "genome")


def test_gem_version(fake_gem):
    assert gem.gem_version() == fake_gem_tools.VERSION
    with pytest.raises(gem.GemError, match="was not found"):
        gem.gem_version("no-such-gem-mapper")


def test_map_reads_writes_the_alignments(tmp_path, fake_gem, indexed):
    sam = tmp_path / "out" / "reads.sam.gz"
    log = tmp_path / "logs" / "map.log"
    result = gem.map_reads(
        indexed["index"],
        indexed["reads"],
        sam,
        log=log,
        max_mismatches=0,
        n_ambiguous=1,
        threads=3,
    )
    assert result.sam == sam and sam.is_file()
    assert not sam.with_name(sam.name + ".partial").exists()
    assert result.counts() == {"n_records": 3, "n_unmapped": 1, "n_out_of_range": 0}
    assert result.mapped_per_read == {"0:NGG": 1, "1:NGG": 1}
    expected = gem.mapper_argv(
        indexed["index"],
        indexed["reads"].resolve(),
        max_mismatches=0,
        n_ambiguous=1,
        threads=3,
    )
    assert list(result.argv) == expected
    assert fake_gem_tools.calls(fake_gem)[-1] == expected
    assert log.read_text().startswith("$ gem-mapper -I ")
    with pysam.AlignmentFile(str(sam), "r", check_sq=False) as alignments:
        records = list(alignments.fetch(until_eof=True))
    placed = {
        (r.query_name, r.reference_start, r.is_reverse)
        for r in records
        if not r.is_unmapped
    }
    assert placed == {
        ("0:NGG", indexed["plus"], False),
        ("1:NGG", indexed["minus"], True),
    }


def test_map_reads_writes_plain_sam_too(tmp_path, fake_gem, indexed):
    sam = tmp_path / "reads.sam"
    gem.map_reads(
        indexed["index"],
        indexed["reads"],
        sam,
        log=tmp_path / "map.log",
        max_mismatches=0,
        n_ambiguous=1,
    )
    assert sam.read_text().startswith("@HD\t")


def test_map_reads_drops_records_outside_their_contig(
    tmp_path, fake_gem, indexed, monkeypatch
):
    monkeypatch.setenv("FAKE_GEM_OUT_OF_RANGE", "1")
    sam = tmp_path / "reads.sam.gz"
    result = gem.map_reads(
        indexed["index"],
        indexed["reads"],
        sam,
        log=tmp_path / "map.log",
        max_mismatches=0,
        n_ambiguous=1,
    )
    # Two mapped reads, each with a record past its contig and one underflowed.
    assert result.counts() == {"n_records": 7, "n_unmapped": 1, "n_out_of_range": 4}
    assert result.mapped_per_read == {"0:NGG": 1, "1:NGG": 1}
    with pysam.AlignmentFile(str(sam), "r", check_sq=False) as alignments:
        assert len(list(alignments.fetch(until_eof=True))) == 3


def test_map_reads_reports_a_crash_with_the_fastq_and_the_read(
    tmp_path, fake_gem, indexed, monkeypatch
):
    monkeypatch.setenv("FAKE_GEM_CRASH", "1")
    sam = tmp_path / "reads.sam.gz"
    log = tmp_path / "map.log"
    with pytest.raises(gem.GemError) as error:
        gem.map_reads(
            indexed["index"], indexed["reads"], sam, log=log, max_mismatches=0
        )
    message = str(error.value)
    assert str(indexed["reads"].resolve()) in message
    assert "'0:NGG'" in message
    assert "Signal raised" in message
    assert str(log) in message
    assert "first contig" in message
    assert error.value.read_name == "0:NGG"
    assert error.value.log == log
    assert not sam.exists()
    assert not sam.with_name(sam.name + ".partial").exists()


def test_map_reads_treats_a_signal_as_a_crash_even_with_exit_status_0(
    tmp_path, fake_gem, indexed, monkeypatch
):
    monkeypatch.setenv("FAKE_GEM_CRASH", "1")
    monkeypatch.setenv("FAKE_GEM_SIGNAL_EXIT", "0")
    with pytest.raises(gem.GemError, match="Signal raised"):
        gem.map_reads(
            indexed["index"],
            indexed["reads"],
            tmp_path / "r.sam",
            log=tmp_path / "m.log",
            max_mismatches=0,
        )
    assert not (tmp_path / "r.sam").exists()


def test_map_reads_reports_a_failure_with_the_end_of_the_log(
    tmp_path, fake_gem, indexed
):
    broken = tmp_path / "broken.gem"
    broken.write_text(str(tmp_path / "missing.fa"))
    with pytest.raises(gem.GemError, match="exit status 1") as error:
        gem.map_reads(
            broken,
            indexed["reads"],
            tmp_path / "r.sam",
            log=tmp_path / "m.log",
            max_mismatches=0,
        )
    assert "FileNotFoundError" in str(error.value)
    assert "first contig" not in str(error.value)
    assert not (tmp_path / "r.sam").exists()


def test_map_reads_checks_its_inputs(tmp_path, fake_gem, indexed):
    with pytest.raises(FileNotFoundError, match="index"):
        gem.map_reads(
            tmp_path / "nope",
            indexed["reads"],
            tmp_path / "r.sam",
            log=tmp_path / "m.log",
            max_mismatches=0,
        )
    with pytest.raises(FileNotFoundError, match="FASTQ"):
        gem.map_reads(
            indexed["index"],
            tmp_path / "nope.fq",
            tmp_path / "r.sam",
            log=tmp_path / "m.log",
            max_mismatches=0,
        )
    with pytest.raises(gem.GemError, match="was not found"):
        gem.map_reads(
            indexed["index"],
            indexed["reads"],
            tmp_path / "r.sam",
            log=tmp_path / "m.log",
            max_mismatches=0,
            gem_mapper="no-such-gem-mapper",
        )
    assert not (tmp_path / "r.sam").exists()


def test_filter_sam_rejects_a_record_on_an_unlisted_contig():
    lines = [
        b"@SQ\tSN:chr1\tLN:100\n",
        b"r\t0\tchr9\t1\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n",
    ]
    with open(os.devnull, "wb") as sink, pytest.raises(ValueError, match="chr9"):
        gem._filter_sam(lines, sink)


@pytest.mark.parametrize(
    ("position", "cigar", "kept"),
    [
        (b"1", b"5M", True),
        (b"96", b"5M", True),
        (b"97", b"5M", False),
        (b"0", b"5M", False),
        (b"-3", b"5M", False),
        (b"18446744073709551594", b"5M", False),
        (b"1e5", b"5M", False),
        (b"95", b"2M1D3M", True),
        (b"96", b"2M1D3M", False),
        (b"96", b"2M1I2M", True),
    ],
)
def test_filter_sam_keeps_records_whose_span_fits(tmp_path, position, cigar, kept):
    out = tmp_path / "out.sam"
    record = b"r\t0\tchr1\t" + position + b"\t60\t" + cigar + b"\t*\t0\t0\t*\t*\n"
    with open(out, "wb") as handle:
        counts = gem._filter_sam([b"@SQ\tSN:chr1\tLN:100\n", record], handle)
    assert counts.n_out_of_range == (0 if kept else 1)
    assert (record in out.read_bytes()) == kept

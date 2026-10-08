import gzip
import logging

import pandas as pd
import pytest

from crispr_analysis_utils.guide_alignment import (
    Guide,
    read_guides,
    write_site_reads,
    write_spacer_reads,
)

SPACER_A = "ACGTTGCAAGCTTCGATCGA"
SPACER_C = "CTTGACCGATAGCATGCAAT"
SPACER_G = "GATCCGTAGCTAGGCTTACG"
PAMS = ["NGG", "NAG", "NGA"]


def write_table(path, rows, *, header=None, sep="\t"):
    lines = [sep.join(header)] if header else []
    lines += [sep.join(row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def read_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as handle:
        lines = handle.read().splitlines()
    return [tuple(lines[i : i + 4]) for i in range(0, len(lines), 4)]


def test_read_guides_keeps_ids_as_strings(tmp_path):
    table = write_table(
        tmp_path / "guides.tsv",
        [("0001", SPACER_A), ("NA", SPACER_C), ("1e3", SPACER_G)],
    )

    guides = read_guides(table)

    assert [guide.guide_id for guide in guides] == ["0001", "NA", "1e3"]


def test_read_guides_uppercases_and_keeps_the_3prime_spacer(tmp_path):
    table = write_table(
        tmp_path / "guides.tsv",
        [("g1", "tt" + SPACER_A.lower()), ("g2", f" {SPACER_C} ")],
    )

    guides = read_guides(table)

    assert guides == [
        Guide("g1", "TT" + SPACER_A, SPACER_A, False, None),
        Guide("g2", SPACER_C, SPACER_C, False, None),
    ]


def test_read_guides_spacer_length(tmp_path):
    table = write_table(tmp_path / "guides.tsv", [("g1", SPACER_A)])

    (guide,) = read_guides(table, spacer_length=19)

    assert guide.spacer == SPACER_A[1:]
    assert guide.input_sequence == SPACER_A


def test_read_guides_with_a_header_by_name_or_position(tmp_path):
    table = write_table(
        tmp_path / "guides.csv",
        [("note", SPACER_A, "g1")],
        header=("comment", "protospacer", "name"),
        sep=",",
    )

    by_name = read_guides(
        table, id_col="name", spacer_col="protospacer", sep=",", header=True
    )
    by_position = read_guides(table, id_col=2, spacer_col=1, sep=",", header=True)

    assert by_name == by_position == [Guide("g1", SPACER_A, SPACER_A, False, None)]


def test_read_guides_from_a_dataframe():
    frame = pd.DataFrame({"id": [7, 8], "seq": [SPACER_A, SPACER_C]})

    guides = read_guides(frame, id_col="id", spacer_col="seq")

    assert [guide.guide_id for guide in guides] == ["7", "8"]
    assert [guide.spacer for guide in guides] == [SPACER_A, SPACER_C]


def test_add_leading_g_marks_spacers_that_do_not_start_with_g():
    frame = pd.DataFrame(
        [("a", SPACER_A), ("g", SPACER_G), ("extra_g", "G" + SPACER_C)]
    )

    with_g = read_guides(frame, add_leading_g=True)
    without_g = read_guides(frame)

    # The 21-base input starts with G, but its 20-base spacer does not.
    assert [guide.leading_g for guide in with_g] == [True, False, True]
    assert not any(guide.leading_g for guide in without_g)


def test_duplicate_spacers_point_at_the_first_guide(caplog):
    frame = pd.DataFrame(
        [("a", SPACER_A), ("b", SPACER_C), ("c", "TT" + SPACER_A), ("d", SPACER_A)]
    )

    with caplog.at_level(logging.WARNING):
        guides = read_guides(frame)

    assert [guide.duplicate_of for guide in guides] == [None, None, "a", "a"]
    assert [guide.input_sequence for guide in guides][2] == "TT" + SPACER_A
    assert "2 guides" in caplog.text
    assert "duplicate_of" in caplog.text


def test_read_guides_logs_nothing_without_duplicates(caplog):
    with caplog.at_level(logging.WARNING):
        read_guides(pd.DataFrame([("a", SPACER_A), ("b", SPACER_C)]))

    assert not caplog.records


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([("g1", SPACER_A), ("g1", SPACER_C)], r"unique.*'g1' \(2 times\)"),
        ([("g1", SPACER_A[:-1] + "N")], r"only A, C, G and T.*'g1' \('N'\)"),
        ([("g1", SPACER_A[:-1] + "U")], r"only A, C, G and T.*'g1' \('U'\)"),
        ([("g1", SPACER_A[:19])], r"at least 20 bases.*'g1' \(19 bases\)"),
        ([("g1", SPACER_A), ("", SPACER_C)], r"must not be empty.*data rows 2"),
    ],
)
def test_read_guides_rejects(tmp_path, rows, message):
    table = write_table(tmp_path / "guides.tsv", rows)

    with pytest.raises(ValueError, match=message):
        read_guides(table)


def test_read_guides_lists_only_the_first_offenders():
    frame = pd.DataFrame([(f"g{i}", "ACGT") for i in range(8)])

    with pytest.raises(ValueError, match="and 3 more"):
        read_guides(frame)


def test_read_guides_rejects_missing_values_in_a_dataframe():
    frame = pd.DataFrame({"id": ["a", None], "seq": [SPACER_A, SPACER_C]})

    with pytest.raises(ValueError, match="Missing guide id.*data rows 2"):
        read_guides(frame)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"id_col": 5}, "out of range"),
        ({"spacer_col": "sequence"}, "not a column"),
        ({"spacer_length": 0}, "positive integer"),
    ],
)
def test_read_guides_rejects_bad_arguments(tmp_path, kwargs, message):
    table = write_table(tmp_path / "guides.tsv", [("g1", SPACER_A)])

    with pytest.raises(ValueError, match=message):
        read_guides(table, **kwargs)


def test_read_guides_rejects_an_empty_table(tmp_path):
    empty = tmp_path / "empty.tsv"
    empty.write_text("", encoding="utf-8")
    header_only = write_table(tmp_path / "header.tsv", [], header=("id", "seq"))

    with pytest.raises(ValueError, match="No guides"):
        read_guides(empty)
    with pytest.raises(ValueError, match="No guides"):
        read_guides(header_only, header=True)


def test_write_spacer_reads_writes_each_spacer_once(tmp_path):
    guides = read_guides(
        pd.DataFrame(
            [("a", SPACER_A), ("b", SPACER_C), ("c", SPACER_A), ("d", SPACER_G)]
        ),
        add_leading_g=True,
    )

    n_reads = write_spacer_reads(guides, tmp_path / "spacers.fastq")

    assert n_reads == 3
    assert read_fastq(tmp_path / "spacers.fastq") == [
        ("@0", SPACER_A, "+", "I" * 20),
        ("@1", SPACER_C, "+", "I" * 20),
        ("@3", SPACER_G, "+", "I" * 20),
    ]


def test_write_site_reads_splits_guides_by_leading_g(tmp_path):
    guides = read_guides(
        pd.DataFrame(
            [("a", SPACER_A), ("g", SPACER_G), ("dup", SPACER_A), ("c", SPACER_C)]
        ),
        add_leading_g=True,
    )

    n_with = write_site_reads(
        guides, tmp_path / "with_g.fastq", pams=PAMS, leading_g=True
    )
    n_without = write_site_reads(
        guides, tmp_path / "without_g.fastq", pams=PAMS, leading_g=False
    )

    assert (n_with, n_without) == (6, 3)
    with_g = read_fastq(tmp_path / "with_g.fastq")
    assert [read[0] for read in with_g] == [
        "@0:NGG",
        "@0:NAG",
        "@0:NGA",
        "@3:NGG",
        "@3:NAG",
        "@3:NGA",
    ]
    assert with_g[0] == ("@0:NGG", "G" + SPACER_A + "NGG", "+", "I" * 24)
    assert with_g[5][1] == "G" + SPACER_C + "NGA"
    assert read_fastq(tmp_path / "without_g.fastq") == [
        ("@1:NGG", SPACER_G + "NGG", "+", "I" * 23),
        ("@1:NAG", SPACER_G + "NAG", "+", "I" * 23),
        ("@1:NGA", SPACER_G + "NGA", "+", "I" * 23),
    ]


def test_write_site_reads_uppercases_the_pams(tmp_path):
    guides = read_guides(pd.DataFrame([("a", SPACER_A)]))

    write_site_reads(guides, tmp_path / "reads.fastq", pams=["nrg"], leading_g=False)

    assert read_fastq(tmp_path / "reads.fastq") == [
        ("@0:NRG", SPACER_A + "NRG", "+", "I" * 23)
    ]


def test_write_site_reads_rejects_bad_pams(tmp_path):
    guides = read_guides(pd.DataFrame([("a", SPACER_A)]))

    with pytest.raises(ValueError, match="same length"):
        write_site_reads(
            guides, tmp_path / "reads.fastq", pams=["NGG", "NG"], leading_g=False
        )


def test_reads_are_gzipped_reproducibly(tmp_path):
    guides = read_guides(pd.DataFrame([("a", SPACER_A), ("b", SPACER_C)]))
    first = tmp_path / "first" / "reads.fastq.gz"
    second = tmp_path / "second" / "reads.fastq.gz"
    first.parent.mkdir()
    second.parent.mkdir()

    write_site_reads(guides, first, pams=PAMS, leading_g=False)
    write_site_reads(guides, second, pams=PAMS, leading_g=False)
    write_spacer_reads(guides, tmp_path / "spacers.fastq.gz")

    assert first.read_bytes()[:2] == b"\x1f\x8b"
    assert first.read_bytes() == second.read_bytes()
    assert len(read_fastq(first)) == 6
    assert read_fastq(tmp_path / "spacers.fastq.gz")[1] == (
        "@1",
        SPACER_C,
        "+",
        "I" * 20,
    )


@pytest.mark.parametrize("guide_id", ["a\tb", "a\nb", "a\rb"])
def test_read_guides_rejects_ids_that_would_break_the_tables(guide_id):
    frame = pd.DataFrame([(guide_id, SPACER_A)])

    with pytest.raises(ValueError, match="tabs or line breaks"):
        read_guides(frame)

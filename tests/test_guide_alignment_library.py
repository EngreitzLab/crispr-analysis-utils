"""`read_guides`: reading, validating and de-duplicating a guide table."""

import logging

import pandas as pd
import pytest
from guide_alignment_fixtures import SPACER_A, SPACER_C, SPACER_G

from crispr_analysis_utils.guide_alignment import Guide, read_guides


def write_table(path, rows, *, header=None, sep="\t"):
    lines = [sep.join(header)] if header else []
    lines += [sep.join(row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_read_guides_keeps_aliases_as_strings(tmp_path):
    table = write_table(
        tmp_path / "guides.tsv",
        [("0001", SPACER_A), ("NA", SPACER_C), ("1e3", SPACER_G)],
    )

    guides = read_guides(table)

    assert [guide.alias for guide in guides] == ["0001", "NA", "1e3"]


def test_read_guides_uppercases_and_keeps_the_3prime_spacer(tmp_path):
    table = write_table(
        tmp_path / "guides.tsv",
        [("g1", "tt" + SPACER_A.lower()), ("g2", f" {SPACER_C} ")],
    )

    guides = read_guides(table)

    assert guides == [
        Guide(
            guide_id="TT" + SPACER_A + "NGG",
            alias="g1",
            input_sequence="TT" + SPACER_A,
            aligned_sequence="TT" + SPACER_A,
            spacer=SPACER_A,
            n_extra=2,
            duplicate_of=None,
        ),
        Guide(
            guide_id=SPACER_C + "NGG",
            alias="g2",
            input_sequence=SPACER_C,
            aligned_sequence=SPACER_C,
            spacer=SPACER_C,
            n_extra=0,
            duplicate_of=None,
        ),
    ]


def test_the_pam_completes_the_guide_id():
    frame = pd.DataFrame([("a", SPACER_A)])

    (default,) = read_guides(frame)
    (other,) = read_guides(frame, pam="NRG")
    (with_g,) = read_guides(frame, add_leading_g=True, pam="TTTV")

    assert default.guide_id == SPACER_A + "NGG"
    assert other.guide_id == SPACER_A + "NRG"
    assert with_g.guide_id == "G" + SPACER_A + "TTTV"


def test_read_guides_spacer_length(tmp_path):
    table = write_table(tmp_path / "guides.tsv", [("g1", SPACER_A)])

    (guide,) = read_guides(table, spacer_length=19)

    assert guide.spacer == SPACER_A[1:]
    assert guide.input_sequence == SPACER_A
    assert guide.aligned_sequence == SPACER_A
    assert guide.n_extra == 1


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

    assert by_name == by_position
    assert by_name == [
        Guide(
            guide_id=SPACER_A + "NGG",
            alias="g1",
            input_sequence=SPACER_A,
            aligned_sequence=SPACER_A,
            spacer=SPACER_A,
            n_extra=0,
            duplicate_of=None,
        )
    ]


def test_read_guides_from_a_dataframe():
    frame = pd.DataFrame({"id": [7, 8], "seq": [SPACER_A, SPACER_C]})

    guides = read_guides(frame, id_col="id", spacer_col="seq")

    assert [guide.alias for guide in guides] == ["7", "8"]
    assert [guide.spacer for guide in guides] == [SPACER_A, SPACER_C]


def test_add_leading_g_prepends_a_g_to_every_sequence():
    frame = pd.DataFrame(
        [("a", SPACER_A), ("g", SPACER_G), ("extra_g", "G" + SPACER_C)]
    )

    with_g = read_guides(frame, add_leading_g=True)
    without_g = read_guides(frame)

    # The G goes in front whatever the first base already is: SPACER_G starts
    # with a G and still gets one.
    assert [guide.aligned_sequence for guide in with_g] == [
        "G" + SPACER_A,
        "G" + SPACER_G,
        "GG" + SPACER_C,
    ]
    assert [guide.input_sequence for guide in with_g] == [
        SPACER_A,
        SPACER_G,
        "G" + SPACER_C,
    ]
    assert [guide.aligned_sequence for guide in without_g] == [
        SPACER_A,
        SPACER_G,
        "G" + SPACER_C,
    ]


def test_g_added_comes_from_the_length_not_from_the_first_base():
    frame = pd.DataFrame(
        [("a", SPACER_A), ("g", SPACER_G), ("extra_g", "G" + SPACER_C)]
    )

    with_g = read_guides(frame, add_leading_g=True)
    without_g = read_guides(frame)

    # SPACER_G starts with a G but is only 20 bases long, so nothing sits 5'
    # of its spacer; the 21-base entry carries one without being told to.
    assert [guide.n_extra for guide in without_g] == [0, 0, 1]
    assert [guide.g_added for guide in without_g] == [False, False, True]
    assert [guide.n_extra for guide in with_g] == [1, 1, 2]
    assert [guide.g_added for guide in with_g] == [True, True, True]


def test_duplicates_are_keyed_on_the_whole_sequence(caplog):
    frame = pd.DataFrame(
        [("a", SPACER_A), ("b", SPACER_C), ("c", "TT" + SPACER_A), ("d", SPACER_A)]
    )

    with caplog.at_level(logging.WARNING):
        guides = read_guides(frame)

    # 'c' shares a spacer with 'a' but not a sequence, so it is its own guide.
    assert [guide.duplicate_of for guide in guides] == [None, None, None, "a"]
    assert guides[2].input_sequence == "TT" + SPACER_A
    assert guides[2].spacer == guides[0].spacer
    assert guides[2].guide_id != guides[0].guide_id
    assert "1 entries" in caplog.text
    assert "'d' (same sequence as 'a')" in caplog.text


def test_a_leading_g_makes_two_guides_not_a_duplicate(caplog):
    frame = pd.DataFrame([("bare", SPACER_A), ("with_g", "G" + SPACER_A)])

    with caplog.at_level(logging.WARNING):
        plain = read_guides(frame)
        prepended = read_guides(frame, add_leading_g=True)

    assert [guide.duplicate_of for guide in plain] == [None, None]
    assert [guide.guide_id for guide in plain] == [
        SPACER_A + "NGG",
        "G" + SPACER_A + "NGG",
    ]
    # Adding a G to both keeps them apart as well.
    assert [guide.duplicate_of for guide in prepended] == [None, None]
    assert [guide.guide_id for guide in prepended] == [
        "G" + SPACER_A + "NGG",
        "GG" + SPACER_A + "NGG",
    ]
    assert not caplog.records


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

    with pytest.raises(ValueError, match="Missing guide alias.*data rows 2"):
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


@pytest.mark.parametrize("alias", ["a\tb", "a\nb", "a\rb"])
def test_read_guides_rejects_aliases_that_would_break_the_tables(alias):
    frame = pd.DataFrame([(alias, SPACER_A)])

    with pytest.raises(ValueError, match="tabs or line breaks"):
        read_guides(frame)

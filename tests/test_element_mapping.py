"""Mapping a guide library to the elements a screen targets.

Every class asserted here is worked out by hand from the intervals in
`element_mapping_fixtures`, which are 0-based half-open, as is every
coordinate `map_to_elements` reads or writes.
"""

import pandas as pd
import pytest
import yaml
from element_mapping_fixtures import (
    LENGTH,
    annotation_frame,
    site,
    sites_frame,
    write_bed,
)
from pandas.testing import assert_frame_equal

from crispr_analysis_utils.element_mapping import (
    CLASSES,
    COLUMNS,
    NON_TARGETING,
    SAFE_TARGETING,
    SITE_COLUMNS,
    TABLE,
    TARGETING,
    ElementError,
    counts,
    map_to_elements,
    read_elements,
    write_columns_yaml,
    write_mapping_tsv,
)

LIBRARY = ("in_one", "in_two", "split", "safe", "never")
"""Every guide of the library, in the order the table reports them."""

SITES = (
    # Inside E1, which no other element overlaps there.
    site("in_one", "chr1", 1000),
    # Where E1 and E2 overlap.
    site("in_two", "chr1", 1060),
    # One guide, two sites: this one is inside E3.
    site("split", "chr2", 2000),
    # ... and this one is in no element at all.
    site("split", "chr1", 5000),
    # Maps, far from every element, and on the other strand.
    site("safe", "chr1", 6000, "-"),
    # "never" has no site, so it is in LIBRARY and not here.
)


@pytest.fixture
def elements(tmp_path):
    """The fixture elements, read back as `map_to_elements` takes them."""
    return read_elements(write_bed(tmp_path / "elements.bed"))


@pytest.fixture
def mapping(elements):
    return map_to_elements(sites_frame(*SITES), elements, guides=pd.Series(LIBRARY))


def rows_for(mapping, guide_id):
    """Every row of one guide, in the order `map_to_elements` returned them."""
    return mapping[mapping["guide_id"] == guide_id]


def one(mapping, guide_id):
    """The one row of a guide with a single site in at most one element."""
    rows = rows_for(mapping, guide_id)
    assert len(rows) == 1
    return rows.iloc[0]


# --- the elements BED -----------------------------------------------------


def test_the_frame_has_the_five_columns_map_to_elements_reads(elements):
    assert list(elements.columns) == [
        "Chromosome",
        "Start",
        "End",
        "element_id",
        "element_category",
    ]
    assert len(elements) == 3


def test_a_bed3_derives_an_id_whose_start_is_one_based(tmp_path):
    frame = read_elements(write_bed(tmp_path / "e.bed", columns=3))
    assert frame["element_id"].tolist() == [
        "chr1:1001-1100",
        "chr1:1051-1200",
        "chr2:2001-2100",
    ]
    assert frame["element_category"].tolist() == ["", "", ""]


def test_a_bed4_keeps_the_files_own_name(tmp_path):
    frame = read_elements(write_bed(tmp_path / "e.bed", columns=4))
    assert frame["element_id"].tolist() == ["E1", "E2", "E3"]
    assert frame["element_category"].tolist() == ["", "", ""]


def test_a_name_that_reads_like_a_locus_is_not_re_derived(tmp_path):
    # An id written elsewhere keeps its spelling, whatever it says.
    rows = [("chr1", 1000, 1100, "chr1:1-10")]
    frame = read_elements(write_bed(tmp_path / "e.bed", rows, columns=4))
    assert frame["element_id"].tolist() == ["chr1:1-10"]


def test_a_bed5_takes_the_fifth_column_as_the_category(elements):
    assert elements["element_category"].tolist() == [
        "enhancer",
        "promoter",
        "enhancer",
    ]


def test_a_sixth_column_is_dropped(tmp_path):
    rows = [("chr1", 1000, 1100, "E1", "enhancer", "+")]
    frame = read_elements(write_bed(tmp_path / "e.bed", rows, columns=6))
    assert list(frame.columns)[-1] == "element_category"
    assert frame["element_category"].tolist() == ["enhancer"]


@pytest.mark.parametrize(
    "header",
    ["# a comment", 'track name="elements"', "browser position chr1:1-2000"],
)
def test_a_comment_track_or_browser_line_is_skipped(tmp_path, header):
    frame = read_elements(write_bed(tmp_path / "e.bed", header=header))
    assert frame["element_id"].tolist() == ["E1", "E2", "E3"]


def test_a_missing_elements_file_is_refused(tmp_path):
    with pytest.raises(ElementError, match="nowhere.bed does not exist"):
        read_elements(tmp_path / "nowhere.bed")


@pytest.mark.parametrize("text", ["", "# only a header\n"])
def test_an_elements_file_with_no_element_is_refused(tmp_path, text):
    path = tmp_path / "e.bed"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ElementError, match="holds no element"):
        read_elements(path)


def test_fewer_than_three_columns_is_refused(tmp_path):
    path = tmp_path / "e.bed"
    path.write_text("chr1\t1000\n", encoding="utf-8")
    with pytest.raises(ElementError, match="has 2 columns; a BED needs at least three"):
        read_elements(path)


# --- the three classes ----------------------------------------------------


def test_a_site_inside_an_element_is_targeting(mapping):
    row = one(mapping, "in_one")
    assert row["target_class"] == TARGETING
    assert (row["element_id"], row["element_category"]) == ("E1", "enhancer")
    assert (row["chr"], row["start"], row["end"]) == ("chr1", 1000, 1020)


def test_a_site_in_no_element_is_safe_targeting(mapping):
    row = one(mapping, "safe")
    assert row["target_class"] == SAFE_TARGETING
    assert (row["element_id"], row["element_category"]) == ("", "")
    # The site itself is still reported, since it is what makes the guide safe.
    assert (row["chr"], row["start"], row["end"], row["strand"]) == (
        "chr1",
        6000,
        6020,
        "-",
    )


def test_a_guide_that_matches_nowhere_is_one_non_targeting_row(mapping):
    row = one(mapping, "never")
    assert row["target_class"] == NON_TARGETING
    assert (row["chr"], row["strand"]) == ("", "")
    assert (row["element_id"], row["element_category"]) == ("", "")
    assert pd.isna(row["start"]) and pd.isna(row["end"])


def test_every_row_carries_one_of_the_three_classes(mapping):
    assert set(mapping["target_class"]) == set(CLASSES)


def test_without_a_library_only_the_guides_with_sites_are_classified(elements):
    frame = map_to_elements(sites_frame(*SITES), elements)
    assert set(frame["guide_id"]) == {"in_one", "in_two", "split", "safe"}
    assert NON_TARGETING not in set(frame["target_class"])


def test_a_library_with_no_site_at_all_is_all_non_targeting(elements):
    frame = map_to_elements(
        sites_frame(*SITES).iloc[:0], elements, guides=pd.Series(LIBRARY)
    )
    assert frame["guide_id"].tolist() == list(LIBRARY)
    assert frame["target_class"].tolist() == [NON_TARGETING] * len(LIBRARY)


# --- per site, not per guide ----------------------------------------------


def test_a_guide_whose_sites_disagree_reports_each_of_them(mapping):
    # The behaviour the class exists for: one site of "split" is inside E3 and
    # the other is in no element, so the guide is neither wholly targeting nor
    # wholly safe. It gets one row of each.
    rows = rows_for(mapping, "split")
    assert len(rows) == 2
    assert dict(zip(rows["chr"], rows["target_class"], strict=True)) == {
        "chr2": TARGETING,
        "chr1": SAFE_TARGETING,
    }
    assert sorted(rows["element_id"]) == ["", "E3"]


def test_a_guide_whose_sites_disagree_counts_under_both_classes(mapping):
    found = counts(mapping)
    assert found["guides"] == len(LIBRARY)
    # Every guide is counted once per class it has a site in, and "split" has
    # two, so the per-class guide counts over-count the library by exactly one.
    assert sum(found["guides_by_class"].values()) == len(LIBRARY) + 1
    assert found["guides_by_class"] == {
        TARGETING: 3,
        SAFE_TARGETING: 2,
        NON_TARGETING: 1,
    }


# --- overlapping elements -------------------------------------------------


def test_a_site_inside_two_elements_gets_a_row_for_each(mapping):
    rows = rows_for(mapping, "in_two")
    assert rows["element_id"].tolist() == ["E1", "E2"]
    assert rows["element_category"].tolist() == ["enhancer", "promoter"]
    assert set(rows["target_class"]) == {TARGETING}
    # One site, reported twice: both rows hold the same coordinates.
    assert len(rows[["chr", "start", "end", "strand"]].drop_duplicates()) == 1


# --- strand ---------------------------------------------------------------


@pytest.mark.parametrize("strand", ["+", "-"])
def test_the_overlap_ignores_the_strand_the_guide_reads_on(elements, strand):
    # A strand-aware join would keep one of these two and drop the other.
    frame = map_to_elements(sites_frame(site("g", "chr1", 1000, strand)), elements)
    assert frame["target_class"].tolist() == [TARGETING]
    assert frame["element_id"].tolist() == ["E1"]


def test_an_element_with_a_strand_of_its_own_still_ignores_it(tmp_path):
    # A BED6 is read to five columns, so a '+' element and a '-' guide meet.
    rows = [("chr1", 1000, 1100, "E1", "enhancer", "+")]
    elements = read_elements(write_bed(tmp_path / "e.bed", rows, columns=6))
    frame = map_to_elements(sites_frame(site("g", "chr1", 1000, "-")), elements)
    assert frame["target_class"].tolist() == [TARGETING]
    assert frame["element_id"].tolist() == ["E1"]


# --- the edges of an element ----------------------------------------------


@pytest.mark.parametrize(
    ("start", "target_class"),
    [
        # E3 is chr2 2000-2100, and a protospacer is 20 bases.
        pytest.param(2000 - LENGTH, SAFE_TARGETING, id="ends_where_it_starts"),
        pytest.param(2000 - LENGTH + 1, TARGETING, id="one_base_past_the_start"),
        pytest.param(2100, SAFE_TARGETING, id="starts_where_it_ends"),
        pytest.param(2100 - 1, TARGETING, id="one_base_before_the_end"),
    ],
)
def test_the_overlap_is_half_open_at_both_edges(elements, start, target_class):
    frame = map_to_elements(sites_frame(site("g", "chr2", start)), elements)
    assert frame["target_class"].tolist() == [target_class]


# --- the distance to the nearest TSS --------------------------------------


DISTANCES = annotation_frame(
    # A site in two genes is two rows of the annotation, and the distance is
    # the site's, so it is the same on both.
    (site("in_one", "chr1", 1000), 12),
    (site("in_one", "chr1", 1000), 12),
    (site("in_two", "chr1", 1060), -300),
    (site("split", "chr2", 2000), 4000),
    (site("split", "chr1", 5000), -50000),
    # Two decoys, so the join is shown to be on all five site columns: the
    # "safe" site reads on -, and "decoy" is another guide on in_one's bases.
    (site("safe", "chr1", 6000, "+"), 99),
    (site("decoy", "chr1", 1000), 999),
)


@pytest.fixture
def with_distances(elements):
    return map_to_elements(
        sites_frame(*SITES),
        elements,
        guides=pd.Series(LIBRARY),
        distances=DISTANCES,
    )


def test_without_an_annotation_every_distance_is_empty(mapping):
    assert mapping["distance_to_tss"].isna().all()


def test_the_distance_is_carried_across_for_a_targeting_site_too(with_distances):
    # 12 is in_one's own, not the 999 of the decoy that shares its bases.
    assert one(with_distances, "in_one")["distance_to_tss"] == 12
    assert set(rows_for(with_distances, "in_two")["distance_to_tss"]) == {-300}


def test_each_site_of_a_guide_keeps_its_own_distance(with_distances):
    rows = rows_for(with_distances, "split")
    assert dict(zip(rows["chr"], rows["distance_to_tss"], strict=True)) == {
        "chr2": 4000,
        "chr1": -50000,
    }


@pytest.mark.parametrize("guide_id", ["safe", "never"])
def test_a_site_the_annotation_never_saw_has_no_distance(with_distances, guide_id):
    # "safe" reads on -, and the annotation holds those bases on +, so the
    # strand is part of the key and the 99 beside it is not this site's.
    assert pd.isna(one(with_distances, guide_id)["distance_to_tss"])


def test_a_repeated_annotation_row_does_not_repeat_the_site(with_distances, mapping):
    assert len(with_distances) == len(mapping)


def test_the_classes_do_not_depend_on_the_annotation(with_distances, mapping):
    assert_frame_equal(
        with_distances.drop(columns=["distance_to_tss"]),
        mapping.drop(columns=["distance_to_tss"]),
    )


def test_an_annotation_without_the_column_leaves_the_distance_empty(elements, mapping):
    frame = map_to_elements(
        sites_frame(*SITES),
        elements,
        guides=pd.Series(LIBRARY),
        distances=DISTANCES.drop(columns=["distance_to_tss"]),
    )
    assert frame["distance_to_tss"].isna().all()
    assert_frame_equal(frame, mapping)


# --- when the sites and the elements cannot be compared -------------------


def test_sites_on_contigs_no_element_names_are_refused(elements):
    sites = sites_frame(site("g", "1", 1000), site("h", "2", 1000))
    with pytest.raises(ElementError, match="No contig of the sites appears"):
        map_to_elements(sites, elements)


def test_the_contig_message_shows_both_namings(elements):
    with pytest.raises(
        ElementError, match="sites use 1 and the elements use chr1, chr2"
    ):
        map_to_elements(sites_frame(site("g", "1", 1000)), elements)


def test_one_shared_contig_is_enough(elements):
    sites = sites_frame(site("here", "chr1", 1000), site("away", "chrZ", 10))
    frame = map_to_elements(sites, elements)
    assert frame["target_class"].tolist() == [TARGETING, SAFE_TARGETING]


@pytest.mark.parametrize("missing", SITE_COLUMNS)
def test_a_missing_site_column_is_named(elements, missing):
    sites = sites_frame(site("g", "chr1", 1000)).drop(columns=[missing])
    with pytest.raises(ElementError, match=f"sites are missing .*{missing}"):
        map_to_elements(sites, elements)


# --- the counts -----------------------------------------------------------


def test_the_class_counts_are_over_sites_and_sum_to_them(mapping):
    found = counts(mapping)
    assert found["rows"] == len(mapping) == 7
    # Five mapped sites. The guide that matches nowhere has none, so it is
    # counted as a class but not as a site.
    assert found["sites"] == 5
    assert found["classes"] == {TARGETING: 3, SAFE_TARGETING: 2, NON_TARGETING: 1}
    assert (
        sum(found["classes"].values())
        == found["sites"] + found["classes"][NON_TARGETING]
    )


def test_elements_with_a_guide_counts_distinct_ids(mapping):
    # E1 twice, E2 and E3 once each; the empty id of a safe site counts nowhere.
    assert counts(mapping)["elements_with_a_guide"] == 3


def test_two_guides_in_one_element_count_it_once(elements):
    sites = sites_frame(site("a", "chr2", 2000), site("b", "chr2", 2050))
    found = counts(map_to_elements(sites, elements))
    assert (found["rows"], found["sites"]) == (2, 2)
    assert found["elements_with_a_guide"] == 1


# --- the table and its data dictionary ------------------------------------


def test_the_site_columns_lead_the_table():
    assert list(COLUMNS)[: len(SITE_COLUMNS)] == list(SITE_COLUMNS)


def test_the_column_order_is_exactly_the_data_dictionary(mapping, tmp_path):
    assert list(mapping.columns) == list(COLUMNS)
    write_mapping_tsv(tmp_path / "m.tsv", mapping)
    header = (tmp_path / "m.tsv").read_text().splitlines()[0]
    assert header.split("\t") == list(COLUMNS)


def test_whole_numbers_are_written_without_a_trailing_zero(tmp_path):
    # A merge against an absent row leaves float64 behind, and 1000.0 is not a
    # coordinate.
    frame = pd.DataFrame(
        {
            "guide_id": ["a", "b"],
            "chr": ["chr1", ""],
            "start": [1000.0, float("nan")],
            "end": [1020.0, float("nan")],
            "strand": ["+", ""],
            "element_id": ["E1", ""],
            "element_category": ["enhancer", ""],
            "target_class": [TARGETING, NON_TARGETING],
            "distance_to_tss": [47.0, float("nan")],
        }
    )
    write_mapping_tsv(tmp_path / "m.tsv", frame)
    lines = (tmp_path / "m.tsv").read_text().splitlines()
    assert lines[1] == "a\tchr1\t1000\t1020\t+\tE1\tenhancer\ttargeting\t47"
    assert lines[2] == "b\t\t\t\t\t\t\tnon_targeting\t"


def test_a_written_mapping_holds_no_decimal_point(with_distances, tmp_path):
    write_mapping_tsv(tmp_path / "m.tsv", with_distances)
    text = (tmp_path / "m.tsv").read_text()
    assert ".0" not in text
    cells = dict(zip(COLUMNS, text.splitlines()[1].split("\t"), strict=True))
    assert cells["guide_id"] == "in_one"
    assert (cells["start"], cells["end"]) == ("1000", "1020")
    assert cells["distance_to_tss"] == "12"


def test_a_non_targeting_row_is_written_with_empty_cells(mapping, tmp_path):
    write_mapping_tsv(tmp_path / "m.tsv", mapping)
    lines = (tmp_path / "m.tsv").read_text().splitlines()
    assert lines[-1] == "never\t\t\t\t\t\t\tnon_targeting\t"


def test_a_column_map_to_elements_never_made_is_written_empty(tmp_path):
    bare = pd.DataFrame({"guide_id": ["a"], "target_class": [NON_TARGETING]})
    write_mapping_tsv(tmp_path / "bare.tsv", bare)
    lines = (tmp_path / "bare.tsv").read_text().splitlines()
    assert lines[0].split("\t") == list(COLUMNS)
    assert lines[1] == "a\t\t\t\t\t\t\tnon_targeting\t"


def test_columns_yaml_parses_and_documents_the_table(mapping, tmp_path):
    write_mapping_tsv(tmp_path / "m.tsv", mapping)
    write_columns_yaml(tmp_path / "columns.yaml")
    parsed = yaml.safe_load((tmp_path / "columns.yaml").read_text())
    header = (tmp_path / "m.tsv").read_text().splitlines()[0].split("\t")
    assert list(parsed) == [TABLE]
    assert list(parsed[TABLE]) == list(COLUMNS) == header
    assert all(parsed[TABLE].values())


def test_columns_yaml_describes_no_column_the_table_lacks(tmp_path):
    write_columns_yaml(tmp_path / "columns.yaml")
    parsed = yaml.safe_load((tmp_path / "columns.yaml").read_text())
    assert set(parsed[TABLE]) == set(COLUMNS)


def test_every_description_is_folded_onto_one_line(tmp_path):
    write_columns_yaml(tmp_path / "columns.yaml")
    parsed = yaml.safe_load((tmp_path / "columns.yaml").read_text())
    assert all("\n" not in text for text in parsed[TABLE].values())

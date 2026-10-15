"""The guide design in the IGVF per-guide metadata format.

What a guide targets is decided from the genome, not from the design's own
label for the element, so every class asserted here is worked out by hand
from the intervals in `igvf_export_fixtures` and `guide_annotation_fixtures`.
Both are 0-based half-open, as is every coordinate the format holds.
"""

import pandas as pd
import pytest
import yaml
from element_mapping_fixtures import write_bed
from guide_annotation_fixtures import EXON_IDS, write_gtf
from igvf_export_fixtures import (
    ELEMENTS,
    GUIDES,
    LENGTH,
    WINDOW,
    alignment_frame,
    genes_map,
    mapping_frame,
    sequence_id,
    spacer,
)

from crispr_analysis_utils.element_mapping import (
    NON_TARGETING,
    SAFE_TARGETING,
    TARGETING,
    read_elements,
)
from crispr_analysis_utils.guide_annotation import build_store, promoters_from_store
from crispr_analysis_utils.igvf_export import (
    COLUMNS,
    NON,
    POSITIVE,
    POSITIVE_CATEGORIES,
    PROMOTER_WINDOW,
    SAFE,
    TABLE,
    TYPES,
    VALIDATED,
    ExportError,
    build_table,
    classify_elements,
    counts,
    write_columns_yaml,
    write_table,
)
from crispr_analysis_utils.igvf_export import (
    ELEMENTS as ELEMENT_KINDS,
)
from crispr_analysis_utils.igvf_export import (
    TARGETING as TARGETING_TYPE,
)

AAA = "ENSG00000000001"
BBB = "ENSG00000000002"
CCC = "ENSG00000000003"
EEE = "ENSG00000000005"
"""The fixture's genes, by the ids `guide_annotation_fixtures` gives them."""


@pytest.fixture
def store(tmp_path):
    return build_store(write_gtf(tmp_path / "test.gtf"), tmp_path / "store.parquet")


@pytest.fixture
def tss(store):
    return promoters_from_store(store)


@pytest.fixture
def elements(tmp_path):
    return read_elements(write_bed(tmp_path / "elements.bed", ELEMENTS))


@pytest.fixture
def info(elements, store, tss):
    """Every element, classified from the genome, with its measured genes."""
    return classify_elements(
        elements, store, tss, promoter_window=WINDOW, genes=genes_map()
    )


@pytest.fixture
def alignment():
    return alignment_frame()


@pytest.fixture
def mapping(elements):
    return mapping_frame(elements)


@pytest.fixture
def table(alignment, mapping, info):
    return build_table(alignment, mapping, info)


def element(info, element_id):
    """The one row of an element."""
    rows = info[info["element_id"] == element_id]
    assert len(rows) == 1
    return rows.iloc[0]


def guide(table, alias):
    """The one row of a library entry."""
    rows = table[table["guide_id"] == alias]
    assert len(rows) == 1
    return rows.iloc[0]


# --- what the genome says an element is -----------------------------------


def test_every_element_is_classified_once_and_keeps_its_id(info):
    assert info["element_id"].tolist() == [row[3] for row in ELEMENTS]
    assert not info["element_id"].duplicated().any()


def test_a_repeated_element_is_classified_once(elements, store, tss):
    doubled = pd.concat([elements, elements], ignore_index=True)
    info = classify_elements(doubled, store, tss, promoter_window=WINDOW)
    assert info["element_id"].tolist() == [row[3] for row in ELEMENTS]


def test_classifying_no_element_at_all_is_refused(elements, store, tss):
    with pytest.raises(ExportError, match="no elements to classify"):
        classify_elements(elements.iloc[:0], store, tss)


@pytest.mark.parametrize(
    ("element_id", "genomic_element"),
    [
        pytest.param("P1", "promoter", id="over_a_tss"),
        pytest.param("A2", "exon", id="an_exon_of_another_transcript"),
        pytest.param("I1", "intron", id="in_a_gene_and_in_no_exon"),
        pytest.param("X1", "exon", id="in_an_exon_far_from_a_tss"),
        pytest.param("D1", "distal element", id="in_no_gene"),
        pytest.param("Z1", "distal element", id="in_no_gene_either"),
        pytest.param("V1", "enhancer", id="validated"),
        pytest.param("O2", "distal element", id="overlapping_its_neighbour"),
        pytest.param("O1", "distal element", id="overlapping_its_neighbour_too"),
    ],
)
def test_each_element_is_what_the_genome_says_it_is(info, element_id, genomic_element):
    assert element(info, element_id)["genomic_element"] == genomic_element


def test_every_class_is_one_the_format_allows(info):
    assert set(info["genomic_element"]) <= set(ELEMENT_KINDS)


@pytest.mark.parametrize("element_id", ["X1", "Z1"])
def test_an_element_the_design_calls_a_promoter_need_not_be_one(info, element_id):
    # The rule exists for this case: the design named these promoters, and
    # the nearest TSS is 1549 bases from X1 and 10 kb from Z1. What they are
    # is read off the genome instead.
    assert element(info, element_id)["element_category"] == "promoter"
    assert element(info, element_id)["genomic_element"] != "promoter"


def test_a_validated_element_is_an_enhancer_even_over_a_tss(info):
    # V1 holds CCC's TSS, so the genome alone would call it a promoter. That
    # it regulates a gene is measured, and no annotation can say it.
    assert element(info, "V1")["element_category"] == VALIDATED
    assert element(info, "V1")["genomic_element"] == "enhancer"


def test_without_a_window_the_formats_two_kilobases_apply(elements, store, tss):
    # X1 is 1549 bases from BBB's TSS, so the default window reaches it.
    wide = classify_elements(elements, store, tss)
    assert element(wide, "X1")["genomic_element"] == "promoter"
    assert element(wide, "X1")["intended_target_name"] == BBB
    # Z1 is 10 kb from the nearest TSS, so the window does not reach it.
    assert element(wide, "Z1")["genomic_element"] == "distal element"


def test_the_window_defaults_to_the_distance_the_format_uses():
    assert PROMOTER_WINDOW == 2000


# --- what the element targets ---------------------------------------------


@pytest.mark.parametrize(
    ("element_id", "intended_target_name"),
    [
        pytest.param("P1", AAA, id="a_promoter_names_its_gene"),
        pytest.param("X1", EXON_IDS["chr1", (2800, 3000)], id="an_exon_names_itself"),
        pytest.param(
            "A2", EXON_IDS["chr1", (1500, 1600)], id="of_whichever_transcript"
        ),
        pytest.param("I1", "chr1:1200-2000", id="an_intron_names_the_gap"),
        pytest.param("V1", "chr2:2950-3050", id="an_enhancer_names_itself"),
        pytest.param("D1", "chr1:20000-20100", id="a_distal_element_too"),
    ],
)
def test_the_target_is_named_by_what_the_element_is(
    info, element_id, intended_target_name
):
    assert element(info, element_id)["intended_target_name"] == intended_target_name


@pytest.mark.parametrize("element_id", ["A2", "X1"])
def test_an_exon_is_named_by_an_ensembl_exon_id(info, element_id):
    name = element(info, element_id)["intended_target_name"]
    assert name.startswith("ENSE")
    # The store strips the version the GTF writes, so the id is stable.
    assert "." not in name


def test_an_intron_is_named_by_the_gap_between_the_reference_exons(info):
    # AAA's reference transcript has exons at 1000-1200 and 2000-2200, and I1
    # sits at 1700-1800, inside the gap.
    assert element(info, "I1")["intended_target_name"] == "chr1:1200-2000"


@pytest.mark.parametrize("element_id", [row[3] for row in ELEMENTS])
def test_the_targets_coordinates_are_the_elements_own(info, element_id):
    row = element(info, element_id)
    assert (row["intended_target_chr"], row["intended_target_start"]) == (
        row["Chromosome"],
        row["Start"],
    )
    assert row["intended_target_end"] == row["End"]


# --- the genes an element is associated with ------------------------------


def test_a_bidirectional_promoter_keeps_both_genes(info):
    # Both genes were measured. The genome offers one TSS per element, so
    # only a measured pair can name a second gene.
    row = element(info, "P1")
    assert row["putative_target_genes"] == f"{AAA},{BBB}"
    assert row["intended_target_name"] == AAA


def test_a_promoter_without_a_measured_gene_takes_the_nearest(elements, store, tss):
    info = classify_elements(elements, store, tss, promoter_window=WINDOW)
    assert element(info, "P1")["putative_target_genes"] == AAA


def test_an_enhancer_keeps_its_genes_and_is_still_named_by_its_coordinates(info):
    assert element(info, "V1")["putative_target_genes"] == CCC
    assert element(info, "V1")["intended_target_name"] == "chr2:2950-3050"


@pytest.mark.parametrize("element_id", ["A2", "I1", "X1", "D1", "Z1", "O1"])
def test_an_element_with_no_gene_has_an_empty_cell(info, element_id):
    assert element(info, element_id)["putative_target_genes"] == ""


# --- one row per library entry --------------------------------------------


def test_the_column_order_is_exactly_the_data_dictionary(table):
    assert list(table.columns) == list(COLUMNS)


def test_there_is_one_row_per_entry_in_library_order(table):
    assert table["guide_id"].tolist() == [alias for alias, *_ in GUIDES]


def test_the_guide_id_is_the_name_the_table_gave_the_entry(table, alignment):
    # guide_id in the alignment is the sequence that was aligned; the format
    # wants a unique name, which is the library's alias for the entry.
    assert set(table["guide_id"]) == set(alignment["alias"])
    assert not set(table["guide_id"]) & set(alignment["guide_id"])


def test_the_spacer_is_the_entrys_own(table):
    assert guide(table, "promoter")["spacer"] == spacer("promoter")
    assert not table["spacer"].duplicated().any()


def test_two_entries_sharing_an_alias_are_refused(alignment, mapping, info):
    doubled = alignment.copy()
    doubled.loc[1, "alias"] = doubled.loc[0, "alias"]
    with pytest.raises(ExportError, match="share an alias"):
        build_table(doubled, mapping, info)


# --- the type and the targeting flag --------------------------------------


@pytest.mark.parametrize(
    ("alias", "expected"),
    [
        pytest.param("promoter", POSITIVE, id="a_promoter_the_genome_agrees_with"),
        pytest.param("validated", POSITIVE, id="a_validated_enhancer"),
        pytest.param("exon", TARGETING_TYPE, id="a_design_promoter_that_is_an_exon"),
        pytest.param("mislabel", TARGETING_TYPE, id="a_design_promoter_far_from_a_tss"),
        pytest.param("intron", TARGETING_TYPE, id="an_enhancer_in_an_intron"),
        pytest.param("distal", TARGETING_TYPE, id="a_distal_element"),
        pytest.param("safe", SAFE, id="a_site_in_no_element"),
        pytest.param("non_targeting", NON, id="no_site_at_all"),
    ],
)
def test_each_entry_gets_the_type_its_element_earns(table, alias, expected):
    assert guide(table, alias)["type"] == expected


def test_a_design_promoter_the_genome_calls_an_exon_is_not_a_positive_control(table):
    # X1's category is one of POSITIVE_CATEGORIES, and that is not enough:
    # the genome has to agree the element is a promoter.
    assert "promoter" in POSITIVE_CATEGORIES
    assert guide(table, "exon")["genomic_element"] == "exon"
    assert guide(table, "exon")["type"] == TARGETING_TYPE


def test_without_a_positive_category_nothing_is_a_positive_control(
    alignment, mapping, info
):
    plain = build_table(alignment, mapping, info, positive_categories=())
    assert POSITIVE not in set(plain["type"])
    assert guide(plain, "promoter")["type"] == TARGETING_TYPE


def test_every_type_is_one_the_format_allows(table):
    assert set(table["type"]) <= set(TYPES)


@pytest.mark.parametrize(
    ("alias", "target_class", "targeting"),
    [
        pytest.param("promoter", TARGETING, "TRUE", id="a_positive_control_targets"),
        pytest.param("distal", TARGETING, "TRUE", id="a_targeting_guide"),
        pytest.param("safe", SAFE_TARGETING, "FALSE", id="a_safe_targeting_control"),
        pytest.param("non_targeting", NON_TARGETING, "FALSE", id="a_non_targeting_one"),
    ],
)
def test_only_a_targeting_guide_is_flagged_targeting(
    table, mapping, alias, target_class, targeting
):
    rows = mapping[mapping["guide_id"] == sequence_id(alias)]
    assert set(rows["target_class"]) == {target_class}
    assert guide(table, alias)["targeting"] == targeting


# --- the site, and the PAM the genome holds -------------------------------


def test_the_site_is_the_protospacer_without_the_pam(table):
    row = guide(table, "promoter")
    assert (row["guide_chr"], row["guide_start"], row["guide_end"]) == (
        "chr1",
        980,
        1000,
    )
    assert row["guide_end"] - row["guide_start"] == LENGTH == len(row["spacer"])
    assert row["strand"] == "+"


@pytest.mark.parametrize(
    ("alias", "pam"),
    [("promoter", "TGG"), ("exon", "CGG"), ("validated", "GGG"), ("safe", "TGG")],
)
def test_the_pam_is_the_bases_the_genome_holds(table, alias, pam):
    assert guide(table, alias)["pam"] == pam


def test_the_pam_is_never_the_pattern(table):
    # mapped_sequence ends in the real PAM, which is why it is read from
    # there rather than from guide_id.
    assert "NGG" not in set(table["pam"].dropna())


def test_a_non_targeting_entry_has_no_site_and_no_element(table):
    row = guide(table, "non_targeting")
    assert row["type"] == NON
    for column in ("guide_chr", "guide_start", "guide_end", "strand", "pam"):
        assert pd.isna(row[column])
    for column in ("genomic_element", "intended_target_name", "intended_target_chr"):
        assert pd.isna(row[column])


def test_a_safe_targeting_entry_keeps_its_site_and_has_no_element(table):
    row = guide(table, "safe")
    assert (row["guide_chr"], row["guide_start"], row["guide_end"]) == (
        "chr1",
        50000,
        50020,
    )
    assert row["pam"] == "TGG"
    assert pd.isna(row["genomic_element"]) and pd.isna(row["intended_target_name"])


# --- a guide inside two elements ------------------------------------------


def test_a_guide_in_two_elements_is_one_row_naming_the_longer_overlap(table, mapping):
    # The mapping reports the guide once per element; the format allows one
    # row, and O2 shares 20 bases of the protospacer against O1's 10.
    rows = mapping[mapping["guide_id"] == sequence_id("overlap")]
    assert sorted(rows["element_id"]) == ["O1", "O2"]
    assert len(table[table["guide_id"] == "overlap"]) == 1
    assert guide(table, "overlap")["intended_target_name"] == "chr2:20000-20200"


def test_an_equal_overlap_is_broken_by_the_element_id(table):
    # 20 bases of each, so only the id can decide, and O1 sorts before O2.
    assert guide(table, "tie")["intended_target_name"] == "chr2:20100-20400"


def test_the_chosen_element_does_not_depend_on_the_rows_order(alignment, mapping, info):
    reversed_rows = mapping.iloc[::-1].reset_index(drop=True)
    again = build_table(alignment, reversed_rows, info)
    assert guide(again, "overlap")["intended_target_name"] == "chr2:20000-20200"
    assert guide(again, "tie")["intended_target_name"] == "chr2:20100-20400"


def test_the_chosen_element_still_carries_every_gene(table):
    assert guide(table, "overlap")["putative_target_genes"] == f"{CCC},{EEE}"


# --- the table on disk ----------------------------------------------------


def test_the_written_header_is_the_data_dictionary(table, tmp_path):
    write_table(tmp_path / "igvf.tsv", table)
    header = (tmp_path / "igvf.tsv").read_text().splitlines()[0]
    assert header.split("\t") == list(COLUMNS)


def test_a_written_coordinate_has_no_decimal_point(table, tmp_path):
    write_table(tmp_path / "igvf.tsv", table)
    text = (tmp_path / "igvf.tsv").read_text()
    assert ".0" not in text
    cells = dict(zip(COLUMNS, text.splitlines()[1].split("\t"), strict=True))
    assert (cells["guide_start"], cells["guide_end"]) == ("980", "1000")
    assert (cells["intended_target_start"], cells["intended_target_end"]) == (
        "950",
        "1050",
    )


def test_a_field_that_does_not_apply_is_an_empty_cell(table, tmp_path):
    write_table(tmp_path / "igvf.tsv", table)
    lines = (tmp_path / "igvf.tsv").read_text().splitlines()
    # guide_id, spacer, oligo, n_oligos, targeting, type: the fixture's
    # alignment predates oligo and n_oligos, so both are empty.
    empty = "\t" * (len(COLUMNS) - 6)
    assert lines[-1] == (
        f"non_targeting\t{spacer('non_targeting')}\t\t\tFALSE\t{NON}{empty}"
    )
    # The reporter and the imperfect guide are unused, so every row ends empty.
    assert all(line.endswith("\t\t") for line in lines[1:])


def test_a_column_the_table_lacks_is_written_empty(tmp_path):
    bare = pd.DataFrame({"guide_id": ["a"], "type": [NON]})
    write_table(tmp_path / "bare.tsv", bare)
    lines = (tmp_path / "bare.tsv").read_text().splitlines()
    assert lines[0].split("\t") == list(COLUMNS)
    assert lines[1].split("\t") == ["a", "", "", "", "", NON] + [""] * (
        len(COLUMNS) - 6
    )


def test_columns_yaml_parses_and_documents_the_table(table, tmp_path):
    write_table(tmp_path / "igvf.tsv", table)
    write_columns_yaml(tmp_path / "columns.yaml")
    parsed = yaml.safe_load((tmp_path / "columns.yaml").read_text())
    header = (tmp_path / "igvf.tsv").read_text().splitlines()[0].split("\t")
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


# --- the counts -----------------------------------------------------------


def test_the_counts_are_over_the_table_and_the_elements(table, info):
    assert counts(table, info) == {
        "guides": len(GUIDES),
        "type": {TARGETING_TYPE: 7, SAFE: 1, NON: 1, POSITIVE: 2},
        "genomic_element": {
            "promoter": 1,
            "enhancer": 1,
            "distal element": 4,
            "exon": 2,
            "intron": 1,
        },
        "elements": len(ELEMENTS),
        # P1, O2 and V1 were measured; the rest have no gene.
        "with_a_gene": 3,
        "missing_pam": 0,
    }


def test_the_type_counts_sum_to_the_library(table, info):
    found = counts(table, info)
    assert sum(found["type"].values()) == found["guides"] == len(table)


def test_a_targeting_guide_without_a_pam_is_counted(alignment, mapping, info):
    # A guide that matches twice has no one mapped_sequence, so the column is
    # empty and the PAM cannot be read off it.
    without = alignment.copy()
    without.loc[without["alias"] == "distal", "mapped_sequence"] = pd.NA
    table = build_table(without, mapping, info)
    assert pd.isna(guide(table, "distal")["pam"])
    assert counts(table, info)["missing_pam"] == 1


def test_a_measured_pair_matches_across_a_one_base_start_difference():
    """chr:start-end is browser notation, so 1-based; a pairs table built
    from the BED's own columns is 0-based."""
    from crispr_analysis_utils.igvf_export import _measured_genes

    known = {"chr11:35618822-35619124": ["ENSG00000179431"]}
    # The element as a browser would name it, one base along.
    assert _measured_genes(known, "chr11:35618823-35619124") == ["ENSG00000179431"]
    # The reverse does not arise from this convention and is not tried, so a
    # neighbouring element one base below cannot be matched by accident.
    assert _measured_genes(known, "chr11:35618821-35619124") == []
    # An exact match still wins, and an unrelated element still finds nothing.
    assert _measured_genes(known, "chr11:35618822-35619124") == ["ENSG00000179431"]
    assert _measured_genes(known, "chr11:99999999-99999999") == []
    assert _measured_genes(known, "not_coordinates") == []

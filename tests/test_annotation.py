"""Placing each aligned site in a gene: exon, intron, promoter or intergenic.

Every number asserted here is worked out by hand from the intervals in
`guide_annotation_fixtures`, which are 0-based half-open, as is every
coordinate `annotate_sites` reads or writes.
"""

import pandas as pd
import pytest
from guide_annotation_fixtures import site, sites_frame, write_gtf
from pandas.testing import assert_frame_equal

from crispr_analysis_utils.guide_annotation import (
    AnnotationError,
    PromoterError,
    annotate_sites,
    build_store,
    promoters_from_store,
)

WINDOW = 100
"""A narrow promoter window, so the fixture's genes stay 1000 bases apart."""

SITES = (
    # In AAA's exon 2, which is wholly CDS.
    site("exon_plus", "chr1", 2030, "+"),
    # The same bases read the other way: the gene is on +, the guide on -.
    site("exon_minus", "chr1", 2030, "-"),
    # In A2's extra exon, which falls inside the reference transcript's intron.
    site("other_transcript_exon", "chr1", 1520, "+"),
    # An intron of every transcript of AAA.
    site("deep_intron", "chr1", 1700, "+"),
    # Where AAA and BBB overlap: AAA's intron and BBB's exon 3.
    site("two_genes", "chr1", 2600, "+"),
    # In BBB's 5'-most CDS block, which is the highest-coordinate one.
    site("minus_cds_first_block", "chr1", 3540, "-"),
    # AAA's exon 1, in its 5' UTR, and inside AAA's promoter window.
    site("utr_and_promoter", "chr1", 1030, "+"),
    # Far from every gene.
    site("intergenic", "chr1", 20000, "+"),
    # Downstream of BBB's TSS at 4499, which is on the minus strand.
    site("downstream_of_minus_tss", "chr1", 4390, "+"),
    # Upstream of AAA's TSS at 1000, and so in no gene but in a promoter.
    site("upstream_of_plus_tss", "chr1", 960, "+"),
    # Inside DDD, which is lncRNA: no store annotates against it, so the row
    # is intergenic and the overlap is reported on its own columns.
    site("in_a_lncrna", "chr1", 7500, "+"),
    # EEE's second transcript starts at 5500; its first starts at 5000, so
    # this sits in the alternative promoter and no other.
    site("alternative_promoter", "chr2", 5480, "+"),
)


@pytest.fixture
def store(tmp_path):
    return build_store(write_gtf(tmp_path / "test.gtf"), tmp_path / "store.parquet")


@pytest.fixture
def tss(store):
    return promoters_from_store(store)


@pytest.fixture
def annotated(store, tss):
    return annotate_sites(sites_frame(*SITES), store, tss, promoter_window=WINDOW)


def rows_for(annotated, guide_id):
    """Every row of one site, in the order `annotate_sites` returned them."""
    return annotated[annotated["guide_id"] == guide_id]


def one(annotated, guide_id):
    """The one row of a site that falls in at most one gene."""
    rows = rows_for(annotated, guide_id)
    assert len(rows) == 1
    return rows.iloc[0]


# --- what comes back ------------------------------------------------------


def test_the_input_columns_come_first_and_untouched(annotated):
    carried = list(SITES[0])
    assert list(annotated.columns)[: len(carried)] == carried
    # Every input row comes back as it went in, repeated where a site falls in
    # two genes, and in the order it was given.
    assert_frame_equal(
        annotated[carried].drop_duplicates(ignore_index=True), sites_frame(*SITES)
    )
    assert annotated["target_site"].tolist()[:2] == [2047, 2033]


def test_one_row_per_site_except_where_two_genes_overlap(annotated):
    assert len(annotated) == len(SITES) + 1
    assert len(rows_for(annotated, "two_genes")) == 2


# --- exon, intron, intergenic ---------------------------------------------


@pytest.mark.parametrize(
    ("guide_id", "location", "feature"),
    [
        ("exon_plus", "genic", "exon"),
        ("exon_minus", "genic", "exon"),
        ("other_transcript_exon", "genic", "intron"),
        ("deep_intron", "genic", "intron"),
        ("minus_cds_first_block", "genic", "exon"),
        ("utr_and_promoter", "genic", "exon"),
        ("intergenic", "intergenic", ""),
        ("upstream_of_plus_tss", "intergenic", ""),
    ],
)
def test_where_each_site_falls(annotated, guide_id, location, feature):
    row = one(annotated, guide_id)
    assert (row["location"], row["feature"]) == (location, feature)


def test_an_exonic_row_names_the_reference_transcript_and_its_exon(annotated):
    row = one(annotated, "exon_plus")
    assert (row["ensg"], row["gene_symbol"]) == ("ENSG00000000001", "AAA")
    assert row["enst"] == "ENST00000000011"
    assert row["exon_number"] == "2"


def test_an_intronic_row_names_the_gene_but_no_transcript(annotated):
    row = one(annotated, "deep_intron")
    assert row["gene_symbol"] == "AAA"
    assert pd.isna(row["enst"])
    assert pd.isna(row["exon_number"])


def test_a_site_in_no_gene_gets_one_row_with_the_gene_columns_empty(annotated):
    row = one(annotated, "intergenic")
    assert row["location"] == "intergenic"
    assert all(pd.isna(row[name]) for name in ("ensg", "gene_symbol", "enst"))
    assert (row["region"], row["sense"]) == ("", "")
    assert row["exonic_any_transcript"] is False


def test_a_site_in_two_genes_gets_a_row_for_each(annotated):
    rows = rows_for(annotated, "two_genes")
    assert rows["gene_symbol"].tolist() == ["AAA", "BBB"]
    # The same bases are AAA's intron and BBB's third exon.
    assert rows["feature"].tolist() == ["intron", "exon"]
    assert rows["location"].tolist() == ["genic", "genic"]
    assert rows["enst"].tolist()[1] == "ENST00000000021"
    assert rows["exon_number"].tolist()[1] == "3"


# --- another transcript's exon --------------------------------------------


def test_a_non_reference_exon_is_flagged_without_changing_the_feature(annotated):
    row = one(annotated, "other_transcript_exon")
    # A2 has an exon at 1500-1600; the reference transcript A1 has an intron.
    assert row["feature"] == "intron"
    assert row["exonic_any_transcript"] is True
    assert pd.isna(row["enst"])
    assert row["region"] == ""


def test_an_intron_of_every_transcript_is_not_flagged(annotated):
    assert one(annotated, "deep_intron")["exonic_any_transcript"] is False


def test_an_exon_of_the_reference_transcript_is_flagged_too(annotated):
    assert one(annotated, "exon_plus")["exonic_any_transcript"] is True


# --- strand ---------------------------------------------------------------


def test_a_guide_on_the_other_strand_is_still_in_the_gene(annotated):
    # AAA reads on +, this guide on -. Every join ignores strand, so the site
    # is annotated; a same-strand join would call it intergenic.
    row = one(annotated, "exon_minus")
    assert (row["gene_symbol"], row["feature"]) == ("AAA", "exon")
    assert row["enst"] == "ENST00000000011"


def test_a_minus_strand_gene_annotates_a_plus_strand_guide(annotated):
    # The other direction: BBB reads on -, this guide on +.
    row = rows_for(annotated, "two_genes").iloc[1]
    assert (row["gene_symbol"], row["feature"]) == ("BBB", "exon")


@pytest.mark.parametrize(
    ("guide_id", "sense"),
    [
        ("exon_plus", "sense"),
        ("exon_minus", "antisense"),
        ("minus_cds_first_block", "sense"),
        ("deep_intron", "sense"),
    ],
)
def test_sense_follows_the_site_strand_against_the_genes(annotated, guide_id, sense):
    assert one(annotated, guide_id)["sense"] == sense


def test_one_site_is_sense_to_one_gene_and_antisense_to_the_other(annotated):
    rows = rows_for(annotated, "two_genes")
    assert rows["sense"].tolist() == ["sense", "antisense"]


# --- along the CDS --------------------------------------------------------


def test_the_cut_is_counted_from_the_start_codon_on_a_plus_strand_gene(annotated):
    # AAA's CDS blocks are 1100-1200, 2000-2200 and 2800-2900: 400 bases. The
    # cut at 2047 is 47 into the second block, which 100 bases precede.
    row = one(annotated, "exon_plus")
    assert row["cds_position"] == 147
    assert row["cds_fraction"] == pytest.approx(0.3675)
    assert row["region"] == "CDS"


def test_the_cut_of_a_minus_strand_guide_is_counted_from_its_own_position(annotated):
    # The same protospacer, cut at 2033 instead: 33 into the second block.
    assert one(annotated, "exon_minus")["cds_position"] == 133


def test_on_a_minus_strand_gene_the_highest_block_comes_first(annotated):
    # BBB's CDS blocks are 2600-2650 and 3500-3600: 150 bases, and 3500-3600
    # is the 5' one. The cut at 3543 is 3600 - 1 - 3543 = 56 bases in.
    row = one(annotated, "minus_cds_first_block")
    assert row["cds_position"] == 56
    assert row["cds_fraction"] == pytest.approx(0.3733)


def test_the_lower_block_of_a_minus_strand_gene_carries_the_offset(annotated):
    # The cut at 2617 is 2650 - 1 - 2617 = 32 into the second block, which the
    # first block's 100 bases precede.
    row = rows_for(annotated, "two_genes").iloc[1]
    assert row["cds_position"] == 132
    assert row["cds_fraction"] == pytest.approx(0.88)


def test_an_exon_whose_cut_misses_the_cds_is_utr(annotated):
    # The cut at 1047 is in AAA's exon 1 but 5' of the start codon at 1100.
    row = one(annotated, "utr_and_promoter")
    assert (row["feature"], row["region"], row["exon_number"]) == ("exon", "UTR", "1")
    assert pd.isna(row["cds_position"])


@pytest.mark.parametrize(
    "guide_id", ["deep_intron", "other_transcript_exon", "intergenic"]
)
def test_region_is_empty_off_an_exon(annotated, guide_id):
    row = one(annotated, guide_id)
    assert row["region"] == ""
    assert pd.isna(row["cds_position"]) and pd.isna(row["cds_fraction"])


# --- the nearest promoter -------------------------------------------------


@pytest.mark.parametrize(
    ("guide_id", "distance", "gene"),
    [
        # AAA's TSS is at 1000 and reads on +, so a later cut is downstream.
        # Measured from the protospacer midpoint, which is the default.
        ("utr_and_promoter", 40, "AAA"),
        ("exon_plus", 1040, "AAA"),
        # Before the TSS of a + gene: upstream, so negative.
        ("upstream_of_plus_tss", -30, "AAA"),
        # BBB's TSS is at 4499 and reads on -, so an earlier site is downstream.
        ("downstream_of_minus_tss", 99, "BBB"),
        ("minus_cds_first_block", 949, "BBB"),
        # After the TSS of a - gene: upstream, so negative.
        ("intergenic", -15511, "BBB"),
    ],
)
def test_the_distance_to_the_tss_is_signed_in_the_genes_orientation(
    annotated, guide_id, distance, gene
):
    row = one(annotated, guide_id)
    assert row["distance_to_tss"] == distance
    assert row["nearest_tss_gene"] == gene


def test_the_nearest_tss_is_not_always_this_rows_gene(annotated):
    # The cut at 2617 is in BBB, but AAA's TSS at 1000 is the nearer one.
    row = rows_for(annotated, "two_genes").iloc[1]
    assert row["gene_symbol"] == "BBB"
    assert (row["nearest_tss_gene"], row["nearest_tss_ensg"]) == (
        "AAA",
        "ENSG00000000001",
    )


@pytest.mark.parametrize(
    ("guide_id", "promoter"),
    [
        ("utr_and_promoter", True),
        ("upstream_of_plus_tss", True),
        ("downstream_of_minus_tss", True),
        ("exon_plus", False),
        ("deep_intron", False),
        ("intergenic", False),
    ],
)
def test_the_promoter_flag_is_the_window_around_any_tss(annotated, guide_id, promoter):
    assert bool(one(annotated, guide_id)["promoter"]) is promoter


def test_a_site_can_be_an_exon_and_a_promoter_at_once(annotated):
    row = one(annotated, "utr_and_promoter")
    assert row["feature"] == "exon"
    assert bool(row["promoter"]) is True


def test_a_site_can_be_a_promoter_and_in_no_gene(annotated):
    row = one(annotated, "upstream_of_plus_tss")
    assert row["location"] == "intergenic"
    assert bool(row["promoter"]) is True


def test_a_wider_window_takes_in_more_sites(store, tss):
    narrow = annotate_sites(sites_frame(*SITES), store, tss, promoter_window=WINDOW)
    wide = annotate_sites(sites_frame(*SITES), store, tss, promoter_window=2000)
    assert wide["promoter"].sum() > narrow["promoter"].sum()


def test_a_zero_window_is_the_tss_base_itself(store, tss):
    sites = sites_frame(
        site("on_the_tss", "chr1", 995, "+"), site("past_it", "chr1", 1005, "+")
    )
    frame = annotate_sites(sites, store, tss, promoter_window=0)
    assert frame["promoter"].tolist() == [True, False]


# --- the protospacer itself -----------------------------------------------


@pytest.mark.parametrize(
    ("spacer", "gc_percent", "polyt"),
    [
        ("ACGTACGTACGTACGTACGT", 50.0, 1),
        ("TTTTACGCGCGCGCGCGCGC", 75.0, 4),
        ("AAAAAAAAAAAAAAAAAAAA", 0.0, 0),
        ("GCGCGCGCGCGCGCGCGCGC", 100.0, 0),
        ("ATTTTTTTTTTTTTTTTTTA", 0.0, 18),
        ("acgtacgtacgtacgtacgt", 50.0, 1),
        ("TTACGTTTTACGTTTACGTT", 30.0, 4),
    ],
)
def test_gc_and_the_longest_run_of_t(store, tss, spacer, gc_percent, polyt):
    frame = annotate_sites(
        sites_frame(site("g", "chr1", 2030, "+", spacer=spacer)),
        store,
        tss,
        promoter_window=WINDOW,
    )
    assert frame.loc[0, "gc_percent"] == pytest.approx(gc_percent)
    assert frame.loc[0, "max_polyt_run"] == polyt


def test_without_a_spacer_column_both_are_empty(store, tss):
    sites = sites_frame(site("g", "chr1", 2030, "+")).drop(columns=["spacer"])
    frame = annotate_sites(sites, store, tss, promoter_window=WINDOW)
    assert pd.isna(frame.loc[0, "gc_percent"])
    assert pd.isna(frame.loc[0, "max_polyt_run"])


# --- when the sites and the annotation cannot be compared -----------------


def test_sites_on_contigs_the_store_never_heard_of_are_refused(store, tss):
    sites = sites_frame(site("g", "1", 2030, "+"), site("h", "2", 2030, "+"))
    with pytest.raises(AnnotationError, match="No contig of the sites appears"):
        annotate_sites(sites, store, tss)


def test_the_contig_message_shows_both_namings(store, tss):
    sites = sites_frame(site("g", "1", 2030, "+"))
    with pytest.raises(AnnotationError, match="sites use 1 and the store uses chr1"):
        annotate_sites(sites, store, tss)


def test_one_shared_contig_is_enough(store, tss):
    sites = sites_frame(site("here", "chr1", 2030, "+"), site("away", "chrZ", 10, "+"))
    frame = annotate_sites(sites, store, tss, promoter_window=WINDOW)
    assert frame["location"].tolist() == ["genic", "intergenic"]


@pytest.mark.parametrize("missing", ["chr", "start", "end", "strand", "target_site"])
def test_a_missing_site_column_is_named(store, tss, missing):
    sites = sites_frame(site("g", "chr1", 2030, "+")).drop(columns=[missing])
    with pytest.raises(AnnotationError, match=f"sites are missing .*{missing}"):
        annotate_sites(sites, store, tss)


def test_a_negative_promoter_window_is_refused(store, tss):
    with pytest.raises(ValueError, match="must not be negative, got -1"):
        annotate_sites(sites_frame(*SITES), store, tss, promoter_window=-1)


def test_a_promoter_table_with_no_tss_is_refused(store, tss):
    with pytest.raises(PromoterError, match="holds no transcription start site"):
        annotate_sites(sites_frame(*SITES), store, tss.iloc[:0])


def test_a_lncrna_overlap_is_reported_without_making_the_row_genic(annotated):
    row = one(annotated, "in_a_lncrna")
    assert (row["location"], row["feature"]) == ("intergenic", "")
    assert pd.isna(row["ensg"])
    assert row["other_gene_types"] == "lncRNA"
    assert row["other_genes"] == "ENSG00000000004"


def test_a_row_over_no_other_gene_type_leaves_those_columns_empty(annotated):
    row = one(annotated, "exon_plus")
    assert (row["other_gene_types"], row["other_genes"]) == ("", "")


def test_an_alternative_promoter_is_found_and_names_its_transcript(annotated):
    row = one(annotated, "alternative_promoter")
    assert row["promoter"]
    # The one-TSS-per-gene file would have had only ENST00000000051 at 5000.
    assert row["promoter_enst"] == "ENST00000000052"
    assert row["promoter_ensg"] == "ENSG00000000005"
    assert row["promoter_gene"] == "EEE"


def test_a_site_in_no_promoter_leaves_the_promoter_ids_empty(annotated):
    row = one(annotated, "deep_intron")
    assert not row["promoter"]
    assert (row["promoter_enst"], row["promoter_ensg"], row["promoter_gene"]) == (
        "",
        "",
        "",
    )


def test_the_distance_to_the_gene_start_is_signed_in_the_genes_orientation(annotated):
    # AAA runs + from 1000; this cut is in its exon 2, well inside.
    plus = one(annotated, "exon_plus")
    assert plus["distance_to_gene_start"] > 0

    # upstream_of_plus_tss sits before AAA entirely, so it is in no gene and
    # the column is empty rather than negative.
    outside = one(annotated, "upstream_of_plus_tss")
    assert pd.isna(outside["distance_to_gene_start"])


def test_a_minus_strand_gene_measures_from_its_own_5_prime_end(annotated, store):
    row = one(annotated, "minus_cds_first_block")
    gene = store.rows("gene")
    bbb = gene[gene["gene_symbol"] == row["gene_symbol"]].iloc[0]
    assert bbb["Strand"] == "-"
    # Measured from End - 1, the gene's 5'-most base, counting inwards.
    midpoint = (row["start"] + row["end"]) // 2
    assert row["distance_to_gene_start"] == bbb["End"] - 1 - midpoint


def test_the_gene_start_distance_separates_canonical_from_alternative(annotated):
    # EEE starts at 5000; the alternative promoter sits 480 bases in, which is
    # what tells the two apart when both are flagged promoter.
    alt = one(annotated, "alternative_promoter")
    assert alt["promoter"]
    # EEE starts at 5000 and the protospacer midpoint is 5490.
    assert alt["distance_to_gene_start"] == 490


@pytest.mark.parametrize(
    ("distance_from", "expected"),
    [("midpoint", 40), ("target_site", 47)],
)
def test_distance_from_chooses_the_reference_point(store, tss, distance_from, expected):
    frame = annotate_sites(
        sites_frame(*SITES),
        store,
        tss,
        promoter_window=WINDOW,
        distance_from=distance_from,
    )
    assert one(frame, "utr_and_promoter")["distance_to_tss"] == expected


def test_an_unknown_distance_reference_is_refused(store, tss):
    with pytest.raises(ValueError, match="distance_from must be one of"):
        annotate_sites(sites_frame(*SITES), store, tss, distance_from="pam")

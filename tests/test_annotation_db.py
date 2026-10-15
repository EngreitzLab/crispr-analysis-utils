"""The annotation store, and the transcription start sites a run measures from.

The store is a GENCODE GTF reduced to the rows annotating needs. Everything
here is built from `guide_annotation_fixtures`, so no GENCODE is downloaded
and every coordinate is one a test can state.
"""

import pandas as pd
import pytest
from guide_annotation_fixtures import (
    GENES,
    REFERENCE,
    TSS,
    write_genes_only_gtf,
    write_gtf,
)
from pandas.testing import assert_frame_equal

from crispr_analysis_utils.guide_annotation import (
    PromoterError,
    Store,
    StoreError,
    build_store,
    load_store,
    promoters_from_store,
    read_promoters,
)
from crispr_analysis_utils.guide_annotation.db import COLUMNS, FEATURES


@pytest.fixture
def gtf(tmp_path):
    return write_gtf(tmp_path / "test.gtf")


@pytest.fixture
def store(gtf, tmp_path):
    return build_store(gtf, tmp_path / "store.parquet")


def references(store):
    """Each gene's reference transcript, as the store marks it."""
    rows = store.rows("transcript", reference_only=True)
    return dict(zip(rows["ensg"], rows["enst"], strict=True))


# --- what build_store keeps -----------------------------------------------


def test_the_store_holds_the_five_features_and_their_columns(store):
    assert list(store.features.columns) == list(COLUMNS)
    assert set(store.features["Feature"].unique()) == set(FEATURES)


def test_only_the_asked_for_gene_types_are_annotated_against(store):
    # DDD is lncRNA. It is kept, so an overlap with it can be reported, but
    # only as a gene row: nothing of it can decide where a guide sits.
    assert sorted(store.rows("gene")["gene_symbol"].unique()) == [
        "AAA",
        "BBB",
        "CCC",
        "EEE",
    ]
    assert sorted(store.other_genes()["gene_symbol"].unique()) == ["DDD"]
    assert store.provenance["genes"] == 4
    assert store.provenance["other_genes"] == 1
    assert store.contigs == {"chr1", "chr2"}


def test_a_gene_of_another_type_keeps_no_feature_but_its_own_row(store):
    ddd = store.features[store.features["gene_symbol"] == "DDD"]
    assert ddd["Feature"].tolist() == ["gene"]
    # So it can never reach the exon, intron or promoter calls.
    for feature in ("transcript", "exon", "CDS", "UTR"):
        assert "DDD" not in set(store.rows(feature, primary=False)["gene_symbol"])


def test_another_gene_type_swaps_which_genes_are_annotated_against(gtf, tmp_path):
    store = build_store(gtf, tmp_path / "lnc.parquet", gene_types=("lncRNA",))
    assert sorted(store.rows("gene")["gene_symbol"].unique()) == ["DDD"]
    assert sorted(store.other_genes()["gene_symbol"].unique()) == [
        "AAA",
        "BBB",
        "CCC",
        "EEE",
    ]
    assert store.provenance["gene_types"] == ["lncRNA"]


def test_the_gtfs_1_based_rows_come_back_0_based_half_open(store):
    gene = store.rows("gene")
    aaa = gene[gene["gene_symbol"] == "AAA"].iloc[0]
    # The GTF line reads 1001 3000; the fixture declares [1000, 3000).
    assert (aaa["Start"], aaa["End"]) == (1000, 3000)
    assert (aaa["Chromosome"], aaa["Strand"]) == ("chr1", "+")


def test_gene_and_transcript_ids_lose_their_version(store):
    assert not store.features["ensg"].str.contains(r"\.", regex=True).any()
    assert set(store.rows("transcript")["enst"]) == {
        transcript["enst"]
        for gene in GENES
        if gene["gene_type"] == "protein_coding"
        for transcript in gene["transcripts"]
    }


def test_the_counts_in_the_provenance_are_the_rows_kept(store):
    counts = store.provenance["counts"]
    for feature, count in counts.items():
        assert count == int((store.features["Feature"] == feature).sum())
    # Five genes, because the lncRNA keeps its one row; every other feature
    # count is of the annotated types alone.
    assert counts == {"gene": 5, "transcript": 9, "exon": 21, "CDS": 13, "UTR": 23}


def test_exons_are_numbered_from_the_transcripts_own_5_prime_end(store):
    exons = store.rows("exon", reference_only=True)
    b1 = exons[exons["enst"] == "ENST00000000021"]
    # BBB reads on the minus strand, so its highest-coordinate exon is exon 1.
    assert dict(zip(b1["Start"], b1["exon_number"], strict=True)) == {
        4300: "1",
        3500: "2",
        2500: "3",
    }


# --- choosing the reference transcript ------------------------------------


@pytest.mark.parametrize(("ensg", "enst"), sorted(REFERENCE.items()))
def test_each_gene_gets_the_reference_transcript_worked_out_by_hand(store, ensg, enst):
    assert references(store)[ensg] == enst


def test_mane_select_is_found_among_several_tags(gtf, store):
    # The row carries basic, MANE_Select and CCDS; a reader that keeps one tag
    # per row drops the middle one and the gene falls through to canonical.
    line = next(
        line
        for line in gtf.read_text().splitlines()
        if "ENST00000000011" in line and "\ttranscript\t" in line
    )
    assert line.count('tag "') == 3
    assert store.provenance["mane_select"] == 1
    assert references(store)["ENSG00000000001"] == "ENST00000000011"


def test_ensembl_canonical_beats_a_longer_cds(store):
    # B2 has 400 coding bases against B1's 150, and is still not the reference.
    cds = store.rows("CDS")
    lengths = (cds["End"] - cds["Start"]).groupby(cds["enst"]).sum()
    assert lengths["ENST00000000022"] > lengths["ENST00000000021"]
    assert references(store)["ENSG00000000002"] == "ENST00000000021"


def test_with_no_tag_at_all_the_longest_cds_wins(store):
    assert references(store)["ENSG00000000003"] == "ENST00000000032"


def test_equal_cds_lengths_are_broken_by_the_lower_id(store):
    assert references(store)["ENSG00000000005"] == "ENST00000000051"


def test_a_rebuild_from_a_reordered_gtf_picks_the_same_transcripts(tmp_path):
    reordered = [
        {**gene, "transcripts": tuple(reversed(gene["transcripts"]))}
        for gene in reversed(GENES)
    ]
    first = build_store(write_gtf(tmp_path / "a.gtf"), tmp_path / "a.parquet")
    second = build_store(
        write_gtf(tmp_path / "b.gtf", reordered), tmp_path / "b.parquet"
    )
    assert references(second) == references(first) == REFERENCE


def test_only_the_reference_transcripts_own_rows_are_reference(store):
    reference = references(store)
    rows = store.features[store.features["enst"].notna()]
    expected = rows["ensg"].map(reference).eq(rows["enst"])
    assert rows["is_reference"].tolist() == expected.tolist()
    assert not store.rows("gene")["is_reference"].any()


# --- when the GTF cannot answer the question ------------------------------


def test_a_genes_only_gtf_is_refused_for_having_no_exons(tmp_path):
    genes_only = write_genes_only_gtf(tmp_path / "genes.gtf")
    with pytest.raises(StoreError, match="holds no exon rows"):
        build_store(genes_only, tmp_path / "store.parquet")


def test_the_genes_only_message_says_what_to_use_instead(tmp_path):
    genes_only = write_genes_only_gtf(tmp_path / "genes.gtf")
    with pytest.raises(StoreError, match="use the full GENCODE annotation"):
        build_store(genes_only, tmp_path / "store.parquet")


def test_a_gene_type_no_gene_has_is_refused(gtf, tmp_path):
    with pytest.raises(StoreError, match="holds no miRNA gene"):
        build_store(gtf, tmp_path / "store.parquet", gene_types=("miRNA",))


def test_a_missing_gtf_is_a_file_not_found_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="does not exist"):
        build_store(tmp_path / "nowhere.gtf", tmp_path / "store.parquet")


def test_nothing_is_written_when_the_gtf_is_refused(gtf, tmp_path):
    with pytest.raises(StoreError):
        build_store(gtf, tmp_path / "store.parquet", gene_types=("miRNA",))
    assert not (tmp_path / "store.parquet").exists()


# --- the round trip -------------------------------------------------------


def test_build_then_load_gives_back_the_same_rows(store, tmp_path):
    loaded = load_store(tmp_path / "store.parquet")
    # build_store hands back pyranges' own DataFrame subclass and load_store a
    # plain one, so the comparison is of the rows, not of the class.
    assert_frame_equal(loaded.features, pd.DataFrame(store.features))
    assert loaded.provenance == store.provenance


def test_the_provenance_survives_the_parquet(store, tmp_path):
    provenance = load_store(tmp_path / "store.parquet").provenance
    assert provenance["gene_types"] == ["protein_coding"]
    assert provenance["reference_transcripts"] == 4
    assert provenance["gtf"].endswith("test.gtf")
    assert provenance["cau_version"]


def test_build_store_makes_the_parent_folder(gtf, tmp_path):
    build_store(gtf, tmp_path / "deep" / "down" / "store.parquet")
    assert (tmp_path / "deep" / "down" / "store.parquet").is_file()


def test_loading_a_store_that_is_not_there_says_how_to_build_one(tmp_path):
    with pytest.raises(StoreError, match="cau guide-annotation build-db"):
        load_store(tmp_path / "nowhere.parquet")


def test_loading_someone_elses_parquet_names_the_missing_columns(tmp_path):
    pd.DataFrame({"Chromosome": ["chr1"], "Start": [1]}).to_parquet(
        tmp_path / "other.parquet"
    )
    with pytest.raises(StoreError, match="is missing End, Feature"):
        load_store(tmp_path / "other.parquet")


# --- what a loaded store offers -------------------------------------------


def test_rows_selects_one_feature(store):
    assert (store.rows("exon")["Feature"] == "exon").all()
    assert len(store.rows("exon")) == 21


def test_rows_can_be_narrowed_to_the_reference_transcript(store):
    exons = store.rows("exon", reference_only=True)
    assert set(exons["enst"]) == set(REFERENCE.values())
    assert len(exons) == 3 + 3 + 3 + 1


def test_a_store_can_be_made_by_hand(store):
    assert Store(features=store.features, provenance={}).contigs == {"chr1", "chr2"}


# --- promoters from a BED -------------------------------------------------


def write_bed(path, lines):
    path.write_text("".join("\t".join(line) + "\n" for line in lines), encoding="utf-8")
    return path


BED6 = [
    ("chr1", "1000", "1001", "AAA", "0", "+"),
    ("chr1", "4499", "4500", "BBB", "0", "-"),
]


def test_a_six_column_bed_is_read_with_no_gene_id(tmp_path):
    frame = read_promoters(write_bed(tmp_path / "tss.bed", BED6))
    assert list(frame.columns) == [
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "gene_symbol",
        "ensg",
        "enst",
    ]
    assert frame["Start"].tolist() == [1000, 4499]
    assert frame["gene_symbol"].tolist() == ["AAA", "BBB"]
    assert frame["Strand"].tolist() == ["+", "-"]
    assert frame["ensg"].tolist() == ["", ""]
    assert frame["enst"].tolist() == ["", ""]


def test_a_leading_hash_line_is_a_header_and_not_a_site(tmp_path):
    header = ("#chr", "start", "end", "gene", "score", "strand")
    with_header = read_promoters(write_bed(tmp_path / "head.bed", [header, *BED6]))
    without = read_promoters(write_bed(tmp_path / "bare.bed", BED6))
    assert_frame_equal(with_header, without)


def test_a_seventh_column_is_the_gene_id_and_an_eighth_the_transcript(tmp_path):
    bed = write_bed(
        tmp_path / "tss8.bed",
        [
            (
                "chr1",
                "1000",
                "1001",
                "AAA",
                "0",
                "+",
                "ENSG00000000001.5",
                "ENST00000000011.2",
            )
        ],
    )
    frame = read_promoters(bed)
    assert frame["ensg"].tolist() == ["ENSG00000000001"]
    assert frame["enst"].tolist() == ["ENST00000000011"]
    assert list(frame.columns) == [
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "gene_symbol",
        "ensg",
        "enst",
    ]


def test_the_contig_stays_a_string(tmp_path):
    bed = write_bed(tmp_path / "plain.bed", [("1", "1000", "1001", "AAA", "0", "+")])
    assert read_promoters(bed)["Chromosome"].tolist() == ["1"]


def test_a_promoter_file_that_is_not_there_is_refused(tmp_path):
    with pytest.raises(PromoterError, match="does not exist"):
        read_promoters(tmp_path / "nowhere.bed")


@pytest.mark.parametrize("width", [1, 4, 5])
def test_too_few_columns_names_the_six_a_tss_bed_needs(tmp_path, width):
    bed = write_bed(tmp_path / "short.bed", [("chr1", "1000", "1001", "AAA")[:width]])
    with pytest.raises(PromoterError, match="a TSS BED needs at least 6"):
        read_promoters(bed)


# --- promoters from the store ---------------------------------------------


def test_the_store_gives_a_tss_for_every_transcript(store):
    frame = promoters_from_store(store)
    # EEE has a second transcript starting at 5500, and a guide at that
    # promoter is invisible to a one-site-per-gene file.
    assert sorted(frame[frame["ensg"] == "ENSG00000000005"]["Start"]) == [5000, 5500]
    assert len(frame) == len(store.rows("transcript"))
    assert (frame["End"] - frame["Start"]).unique().tolist() == [1]


def test_reference_only_gives_back_one_tss_per_gene(store):
    frame = promoters_from_store(store, reference_only=True)
    assert dict(zip(frame["ensg"], frame["Start"], strict=True)) == TSS
    assert (frame["End"] - frame["Start"]).tolist() == [1, 1, 1, 1]


def test_a_minus_strand_tss_is_the_last_base_of_the_transcript(store):
    frame = promoters_from_store(store)
    bbb = frame[frame["gene_symbol"] == "BBB"].iloc[0]
    transcript = store.rows("transcript", reference_only=True)
    b1 = transcript[transcript["enst"] == "ENST00000000021"].iloc[0]
    assert (b1["Strand"], b1["End"]) == ("-", 4500)
    assert (bbb["Start"], bbb["End"]) == (4499, 4500)


def test_promoters_from_the_store_carry_the_gene_columns(store):
    frame = promoters_from_store(store)
    assert list(frame.columns) == [
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "gene_symbol",
        "ensg",
        "enst",
    ]
    assert sorted(set(frame["gene_symbol"])) == ["AAA", "BBB", "CCC", "EEE"]
    assert frame["enst"].str.startswith("ENST").all()

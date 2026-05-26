import crispr_analysis_utils
import pandas as pd
import pysam
from crispr_analysis_utils.guide_qc import (
    _evaluate_alignment_layout,
    _mismatch_positions_from_md,
    _resolve_allowed_contigs,
    filter_guide_alignments,
)


def test_guides_to_fastq_from_top_level_namespace(tmp_path):
    guides = tmp_path / "guides.tsv"
    guides.write_text("guide_1\tACGT\nguide_2\tTGCA\n", encoding="utf-8")
    output = tmp_path / "guides.fastq"

    crispr_analysis_utils.guide_qc.guides_to_fastq(guides, output)

    text = output.read_text(encoding="utf-8")
    assert "@guide_1\nACGT\n+\nIIII\n" in text
    assert "@guide_2\nTGCA\n+\nIIII\n" in text


def test_guides_to_fastq_appends_pam(tmp_path):
    guides = tmp_path / "guides.tsv"
    guides.write_text("guide_1\tACGT\n", encoding="utf-8")
    output = tmp_path / "guides.fastq"

    crispr_analysis_utils.guide_qc.guides_to_fastq(guides, output, pam="NGG")

    text = output.read_text(encoding="utf-8")
    assert "@guide_1\nACGTNGG\n+\nIIIIIII\n" in text


def test_guides_to_fastq_adds_leading_g_only_if_missing(tmp_path):
    guides = tmp_path / "guides.tsv"
    guides.write_text("guide_1\tACGT\nguide_2\tGTTT\n", encoding="utf-8")
    output = tmp_path / "guides.fastq"

    crispr_analysis_utils.guide_qc.guides_to_fastq(
        guides, output, add_leading_g=True
    )

    text = output.read_text(encoding="utf-8")
    assert "@guide_1\nGACGT\n+\nIIIII\n" in text
    assert "@guide_2\nGTTT\n+\nIIII\n" in text


def test_guides_to_fastq_accepts_dataframe_with_column_overrides(tmp_path):
    guides = pd.DataFrame(
        {
            "id": ["guide_1", "guide_2"],
            "seq": ["ACGT", "TGCA"],
            "extra": [1, 2],
        }
    )
    output = tmp_path / "guides.fastq"

    crispr_analysis_utils.guide_qc.guides_to_fastq(
        guides,
        output,
        id_col="id",
        sequence_col="seq",
        pam="NGG",
    )

    text = output.read_text(encoding="utf-8")
    assert "@guide_1\nACGTNGG\n+\nIIIIIII\n" in text
    assert "@guide_2\nTGCANGG\n+\nIIIIIII\n" in text


def test_evaluate_alignment_layout_accepts_clean_match():
    # 1S20M3M -> leading G can be softclipped, spacer aligned, PAM aligned
    cigartuples = [(4, 1), (0, 20), (0, 3)]
    query_len = 24
    ok, reason = _evaluate_alignment_layout(
        cigartuples, query_len, set(), pam="NGG", allow_leading_g_softclip=True
    )
    assert ok
    assert reason == "valid"


def test_evaluate_alignment_layout_rejects_gg_pam_mismatch():
    cigartuples = [(0, 24)]
    # mismatch at second PAM base (index 22 for length 24)
    ok, reason = _evaluate_alignment_layout(
        cigartuples, 24, {22}, pam="NGG", allow_leading_g_softclip=True
    )
    assert not ok
    assert reason == "pam_gg_mismatch"


def test_evaluate_alignment_layout_rejects_spacer_indel():
    # insertion in spacer
    cigartuples = [(0, 10), (1, 1), (0, 13)]
    ok, reason = _evaluate_alignment_layout(
        cigartuples, 24, set(), pam="NGG", allow_leading_g_softclip=True
    )
    assert not ok
    assert reason == "spacer_insertion"


def test_mismatch_positions_from_md_parses_single_mismatch():
    # 10 matches, one mismatch, then 13 matches for 24bp read
    mm = _mismatch_positions_from_md("10A13", [(0, 24)])
    assert 10 in mm


def test_resolve_allowed_contigs_from_dataframe():
    df = pd.DataFrame({"chrom": ["chr1", "chr2"], "size": [1000, 2000]})
    allowed = _resolve_allowed_contigs(primary_contigs=None, chromsizes=df)
    assert allowed == {"chr1", "chr2"}


def test_filter_guide_alignments_uses_chromsizes_file(tmp_path):
    sam_path = tmp_path / "in.sam"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [
            {"SN": "chr1", "LN": 1000},
            {"SN": "chrUn", "LN": 1000},
        ],
    }
    with pysam.AlignmentFile(str(sam_path), "w", header=header) as out:
        a = pysam.AlignedSegment()
        a.query_name = "g1"
        a.query_sequence = "G" + "A" * 20 + "TGG"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigarstring = "24M"
        a.set_tag("MD", "24")
        a.set_tag("NM", 0)
        out.write(a)

        b = pysam.AlignedSegment()
        b.query_name = "g2"
        b.query_sequence = "G" + "A" * 20 + "TGG"
        b.flag = 0
        b.reference_id = 1
        b.reference_start = 100
        b.mapping_quality = 60
        b.cigarstring = "24M"
        b.set_tag("MD", "24")
        b.set_tag("NM", 0)
        out.write(b)

    chromsizes = tmp_path / "chrom.sizes"
    chromsizes.write_text("chr1\t1000\n", encoding="utf-8")

    summary = filter_guide_alignments(
        sam_path,
        tmp_path / "unique.sam",
        tmp_path / "multi.sam",
        chromsizes=chromsizes,
    )
    assert summary["valid_guides"] == 1
    assert summary["invalid_alignments"] == 1
    assert (tmp_path / "valid_alignments.bed").exists()
    assert (tmp_path / "discarded_alignments.tsv").exists()
    assert (tmp_path / "unmapped.tsv").exists()
    assert (tmp_path / "guide_alignment_log.tsv").exists()
    bed_line = (tmp_path / "valid_alignments.bed").read_text(encoding="utf-8").strip().splitlines()[0]
    assert len(bed_line.split("\t")) == 9


def test_filter_guide_alignments_guide_log_counts(tmp_path):
    sam_path = tmp_path / "in_counts.sam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(sam_path), "w", header=header) as out:
        # mapped valid
        a = pysam.AlignedSegment()
        a.query_name = "g1"
        a.query_sequence = "G" + "A" * 20 + "TGG"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigarstring = "24M"
        a.set_tag("MD", "24")
        a.set_tag("NM", 0)
        out.write(a)

        # mapped discarded (insertion in spacer)
        b = pysam.AlignedSegment()
        b.query_name = "g1"
        b.query_sequence = "A" * 21 + "TGG"
        b.flag = 0
        b.reference_id = 0
        b.reference_start = 200
        b.mapping_quality = 60
        b.cigarstring = "10M1I13M"
        b.set_tag("MD", "23")
        b.set_tag("NM", 1)
        out.write(b)

        # unmapped for another guide
        c = pysam.AlignedSegment()
        c.query_name = "g2"
        c.query_sequence = "G" + "A" * 20 + "TGG"
        c.flag = 4
        c.reference_id = -1
        c.reference_start = -1
        c.mapping_quality = 0
        c.cigarstring = None
        out.write(c)

    summary = filter_guide_alignments(sam_path, tmp_path / "unique.sam", tmp_path / "multi.sam")
    log_lines = (tmp_path / "guide_alignment_log.tsv").read_text(encoding="utf-8").strip().splitlines()
    assert log_lines[0] == "guide_id\tn_aligned\tn_valid\tn_discarded\tn_not_mapped"
    assert "g1\t2\t1\t1\t0" in log_lines
    assert "g2\t0\t0\t0\t1" in log_lines
    assert summary["guides_unique_valid"] == 1
    assert summary["guides_multi_valid"] == 0
    assert summary["guides_aligned_none_valid"] == 0
    assert summary["guides_unmapped"] == 1
    assert summary["guides_one_valid_plus_invalid"] == 1
    summary_lines = (tmp_path / "alignment_summary.tsv").read_text(encoding="utf-8").strip().splitlines()
    assert summary_lines[0] == "metric\tcount"
    assert "guides_one_valid_plus_invalid\t1" in summary_lines


def test_filter_guide_alignments_writes_alias_column(tmp_path):
    sam_path = tmp_path / "in_alias.sam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    seq = "G" + "A" * 20 + "TGG"
    with pysam.AlignmentFile(str(sam_path), "w", header=header) as out:
        a = pysam.AlignedSegment()
        a.query_name = seq
        a.query_sequence = seq
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigarstring = "24M"
        a.set_tag("MD", "24")
        a.set_tag("NM", 0)
        a.set_tag("AS", 24)
        out.write(a)

    filter_guide_alignments(
        sam_path,
        None,
        None,
        alias_by_guide_id={seq: "my_alias"},
        output_valid_bed=tmp_path / "valid_alignments.bed",
    )
    line = (tmp_path / "valid_alignments.bed").read_text(encoding="utf-8").strip()
    assert line.endswith("\tmy_alias")


def test_filter_guide_alignments_reuses_sequence_for_seq_star_records(tmp_path):
    sam_path = tmp_path / "in_seq_star.sam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr7", "LN": 200000000}]}
    seq = "CCNCTCTTCATTCATGTTTATTGC"
    with pysam.AlignmentFile(str(sam_path), "w", header=header) as out:
        # First alignment has sequence
        a = pysam.AlignedSegment()
        a.query_name = "guide_alias_1"
        a.query_sequence = seq
        a.flag = 16
        a.reference_id = 0
        a.reference_start = 45069710
        a.mapping_quality = 0
        a.cigarstring = "23M1S"
        a.set_tag("MD", "2A20")
        a.set_tag("NM", 1)
        a.set_tag("AS", 18)
        out.write(a)

        # Second alignment for same QNAME has SEQ=*
        b = pysam.AlignedSegment()
        b.query_name = "guide_alias_1"
        b.query_sequence = None
        b.flag = 256
        b.reference_id = 0
        b.reference_start = 153279377
        b.mapping_quality = 0
        b.cigarstring = "2S22M"
        b.set_tag("MD", "19G2")
        b.set_tag("NM", 1)
        b.set_tag("AS", 17)
        out.write(b)

    filter_guide_alignments(
        sam_path,
        None,
        None,
        output_discarded_tsv=tmp_path / "discarded_alignments.tsv",
    )
    discarded = (tmp_path / "discarded_alignments.tsv").read_text(encoding="utf-8")
    # Sequence should be used as guide_id/read_name even on SEQ=* line
    assert "CCNCTCTTCATTCATGTTTATTGC" in discarded


def test_filter_guide_alignments_reports_sequence_5prime_to_3prime(tmp_path):
    sam_path = tmp_path / "in_orientation.sam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr7", "LN": 200000000}]}
    seq_5to3 = "CCNCTCTTCATTCATGTTTATTGC"
    # What SAM stores for reverse-strand records is reverse-complemented SEQ.
    seq_stored_in_sam = "GCAATAAACATGAATGAAGAGNGG"
    with pysam.AlignmentFile(str(sam_path), "w", header=header) as out:
        a = pysam.AlignedSegment()
        a.query_name = "alias_orientation"
        a.query_sequence = seq_stored_in_sam
        a.flag = 16
        a.reference_id = 0
        a.reference_start = 45069710
        a.mapping_quality = 0
        a.cigarstring = "24M"
        a.set_tag("MD", "24")
        a.set_tag("NM", 0)
        a.set_tag("AS", 20)
        out.write(a)

    filter_guide_alignments(
        sam_path,
        None,
        None,
        output_valid_bed=tmp_path / "valid_alignments.bed",
    )
    bed = (tmp_path / "valid_alignments.bed").read_text(encoding="utf-8")
    assert seq_5to3 in bed

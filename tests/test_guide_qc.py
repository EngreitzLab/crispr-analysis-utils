import crispr_analysis_utils
import pandas as pd


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


def test_guides_to_fastq_always_adds_leading_g_when_enabled(tmp_path):
    guides = tmp_path / "guides.tsv"
    guides.write_text("guide_1\tACGT\nguide_2\tGTTT\n", encoding="utf-8")
    output = tmp_path / "guides.fastq"

    crispr_analysis_utils.guide_qc.guides_to_fastq(
        guides, output, add_leading_g=True
    )

    text = output.read_text(encoding="utf-8")
    assert "@guide_1\nGACGT\n+\nIIIII\n" in text
    assert "@guide_2\nGGTTT\n+\nIIIII\n" in text


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

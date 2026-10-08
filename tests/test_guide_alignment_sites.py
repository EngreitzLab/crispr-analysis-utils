import gzip
import shutil

import pandas as pd
import pytest
from guide_alignment_fixtures import (
    SPACER_A,
    SPACER_C,
    SPACER_G,
    Contig,
    hit,
    unmapped,
    write_fasta,
    write_sam,
)

from crispr_analysis_utils.guide_alignment import Site, read_guides, read_sites

SWAP = {"A": "C", "C": "G", "G": "T", "T": "A"}


def mutate(spacer, positions):
    """Change the bases at these 1-based positions."""
    bases = list(spacer)
    for position in positions:
        bases[position - 1] = SWAP[bases[position - 1]]
    return "".join(bases)


def library(*spacers, add_leading_g=False):
    rows = [(f"g{i + 1}", spacer) for i, spacer in enumerate(spacers)]
    return read_guides(pd.DataFrame(rows), add_leading_g=add_leading_g)


def run(tmp_path, reference, records, guides, **kwargs):
    """Write the reference and one SAM, then read the sites."""
    fasta = write_fasta(tmp_path / "ref.fa", reference)
    lengths = {name: len(sequence) for name, sequence in reference.items()}
    sam = write_sam(tmp_path / "sites.sam", lengths, records)
    return read_sites(sam, guides, fasta, **kwargs)


def test_plus_strand_site_without_g(tmp_path):
    contig = "A" * 10 + SPACER_A + "TGG" + "A" * 10

    result = run(
        tmp_path,
        {"chr1": contig},
        [hit("0:NGG", "chr1", 10, "+", 23, nm=1)],
        library(SPACER_A),
    )

    assert result.sites == {
        0: [Site("g1", "chr1", 10, 30, "+", "TGG", "NGG", 0, (), SPACER_A, None)]
    }
    assert result.n_nm_disagreements == 0


def test_minus_strand_site_counts_mismatches_from_the_guide_5prime_end(tmp_path):
    protospacer = mutate(SPACER_A, [1])  # guide orientation
    assert protospacer[0] == "C"
    contig = Contig()
    contig.add("A" * 10)
    assert contig.add_minus(protospacer + "AGG") == 10
    contig.add("A" * 10)
    # On the minus strand, the guide's 5' end is the protospacer's last
    # forward-strand base: [13, 33) holds it, so position 1 sits at 32.
    assert contig.sequence[10:13] == "CCT"
    assert contig.sequence[32] == "G"  # the complement of the mutated C

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [hit("0:NGG", "chr1", 10, "-", 23, nm=2)],
        library(SPACER_A),
    )

    assert result.sites == {
        0: [Site("g1", "chr1", 13, 33, "-", "AGG", "NGG", 1, (1,), protospacer, None)]
    }
    assert result.n_nm_disagreements == 0


def test_leading_g_on_both_strands(tmp_path):
    contig = Contig()
    contig.add("A" * 5)
    plus = contig.add("T" + SPACER_A + "CGG")
    contig.add("A" * 5)
    minus = contig.add_minus("C" + SPACER_A + "GGG")
    contig.add("A" * 5)
    assert (plus, minus) == (5, 34)
    assert contig.sequence[57] == "G"  # the minus-strand G sits after the protospacer

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [
            hit("0:NGG", "chr1", 5, "+", 24, nm=2),
            hit("0:NGG", "chr1", 34, "-", 24, nm=2, secondary=True),
        ],
        library(SPACER_A, add_leading_g=True),
    )

    # The G is not counted: its genomic base is reported instead.
    assert result.sites == {
        0: [
            Site("g1", "chr1", 6, 26, "+", "CGG", "NGG", 0, (), SPACER_A, "T"),
            Site("g1", "chr1", 37, 57, "-", "GGG", "NGG", 0, (), SPACER_A, "C"),
        ]
    }
    assert result.n_nm_disagreements == 0


def test_a_spacer_that_starts_with_g_counts_its_own_g(tmp_path):
    protospacer = mutate(SPACER_G, [1])
    contig = "A" * 5 + protospacer + "TGG" + "A" * 5
    guides = library(SPACER_G, add_leading_g=True)
    assert not guides[0].leading_g

    result = run(
        tmp_path, {"chr1": contig}, [hit("0:NGG", "chr1", 5, "+", 23, nm=2)], guides
    )

    (site,) = result.sites[0]
    assert (site.start, site.n_mismatches, site.mismatch_positions) == (5, 1, (1,))
    assert site.leading_g_base is None


def test_sites_running_off_a_contig_end_are_skipped(tmp_path):
    guides = library(SPACER_G, SPACER_C, add_leading_g=True)
    assert [guide.leading_g for guide in guides] == [False, True]
    g_cut = Contig()
    g_cut.add_minus(SPACER_C + "TGG")  # the minus-strand G would sit at 23
    contigs = {
        "pam_cut": "A" * 5 + SPACER_G + "TG",
        "pam_fits": "A" * 5 + SPACER_G + "TGG",
        "g_cut": g_cut.sequence,
        "g_first": "T" + SPACER_C + "AGG",
    }

    result = run(
        tmp_path,
        contigs,
        [
            hit("0:NGG", "pam_cut", 5, "+", 23),
            hit("0:NGG", "pam_fits", 5, "+", 23),
            hit("1:NGG", "g_cut", 0, "-", 24),
            hit("1:NGG", "g_first", 0, "+", 24),
        ],
        guides,
    )

    assert result.n_off_contig_end == 2
    assert result.sites == {
        0: [Site("g1", "pam_fits", 5, 25, "+", "TGG", "NGG", 0, (), SPACER_G, None)],
        1: [Site("g2", "g_first", 1, 21, "+", "AGG", "NGG", 0, (), SPACER_C, "T")],
    }


@pytest.mark.parametrize(
    ("alt_pams", "classes"),
    [
        (("NAG", "NGA"), ["NGG", "NAG", "NGA", "none", "NGG", "none"]),
        ((), ["NGG", "none", "none", "none", "NGG", "none"]),
        (("NRG",), ["NGG", "NRG", "none", "none", "NGG", "none"]),
    ],
)
def test_pam_classes_follow_the_genome(tmp_path, alt_pams, classes):
    genomic_pams = ["TGG", "TAG", "CGA", "TTT", "NGG", "TNG"]
    # The aligner counts the read's N, and each other PAM base that differs.
    nm = [1, 2, 2, 3, 1, 2]
    contig = Contig()
    positions = []
    for genomic_pam in genomic_pams:
        contig.add("A" * 5)
        positions.append(contig.add(SPACER_A + genomic_pam))

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [
            hit("0:NGG", "chr1", position, "+", 23, nm=n)
            for position, n in zip(positions, nm, strict=True)
        ],
        library(SPACER_A),
        alt_pams=alt_pams,
    )

    assert [site.pam for site in result.sites[0]] == genomic_pams
    assert [site.pam_class for site in result.sites[0]] == classes
    assert result.n_nm_disagreements == 0


def test_a_soft_masked_reference_is_uppercased(tmp_path):
    protospacer = mutate(SPACER_A, [10])
    contig = Contig()
    contig.add("a" * 5)
    plus = contig.add((protospacer + "tgg").lower())
    minus = contig.add_minus((protospacer + "agg").lower())

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [
            hit("0:NGG", "chr1", plus, "+", 23, nm=2),
            hit("0:NGG", "chr1", minus, "-", 23, nm=2),
        ],
        library(SPACER_A),
    )

    for site in result.sites[0]:
        assert site.genomic_protospacer == protospacer
        assert site.mismatch_positions == (10,)
        assert site.pam_class == "NGG"
    assert [site.pam for site in result.sites[0]] == ["TGG", "AGG"]
    assert result.n_nm_disagreements == 0


def test_a_genomic_n_in_the_protospacer_is_a_mismatch(tmp_path):
    protospacer = SPACER_A[:4] + "N" + SPACER_A[5:]
    contig = "A" * 5 + protospacer + "TGG"

    result = run(
        tmp_path,
        {"chr1": contig},
        [hit("0:NGG", "chr1", 5, "+", 23, nm=2)],
        library(SPACER_A),
    )

    (site,) = result.sites[0]
    assert (site.n_mismatches, site.mismatch_positions) == (1, (5,))
    assert site.genomic_protospacer == protospacer


@pytest.mark.parametrize(
    ("max_mismatches", "kept", "dropped"),
    [(3, [3], 1), (4, [3, 4], 0), (0, [], 2)],
)
def test_sites_over_max_mismatches_are_dropped(tmp_path, max_mismatches, kept, dropped):
    contig = Contig()
    three = contig.add(mutate(SPACER_A, [1, 2, 3]) + "TGG")
    four = contig.add_minus(mutate(SPACER_A, [1, 2, 3, 20]) + "TGG")

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [hit("0:NGG", "chr1", three, "+", 23), hit("0:NGG", "chr1", four, "-", 23)],
        library(SPACER_A),
        max_mismatches=max_mismatches,
    )

    assert [site.n_mismatches for site in result.sites.get(0, [])] == kept
    assert result.n_over_max_mismatches == dropped
    if 4 in kept:
        assert result.sites[0][1].mismatch_positions == (1, 2, 3, 20)


def test_sites_found_by_several_pam_reads_are_merged(tmp_path):
    contig = Contig()
    plus = contig.add(SPACER_A + "TGG")
    contig.add("A" * 5)
    minus = contig.add_minus(SPACER_A + "TAG")

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [
            hit("0:NGG", "chr1", plus, "+", 23, nm=1),
            hit("0:NAG", "chr1", plus, "+", 23, nm=2, secondary=True),
            hit("0:NGA", "chr1", plus, "+", 23, nm=2, secondary=True),
            hit("0:NGG", "chr1", minus, "-", 23, nm=2),
            hit("0:NAG", "chr1", minus, "-", 23, nm=1, secondary=True),
        ],
        library(SPACER_A),
    )

    assert [(s.start, s.strand, s.pam_class) for s in result.sites[0]] == [
        (0, "+", "NGG"),
        (31, "-", "NAG"),
    ]
    assert result.n_repeated == 3
    assert result.n_nm_disagreements == 0


def test_unmapped_reads_are_counted_and_skipped(tmp_path):
    contig = SPACER_A + "TGG"

    result = run(
        tmp_path,
        {"chr1": contig},
        [
            hit("0:NGG", "chr1", 0, "+", 23),
            unmapped("0:NAG"),
            unmapped("1:NGG"),
            unmapped("1:NAG"),
        ],
        library(SPACER_A, SPACER_C),
    )

    assert list(result.sites) == [0]
    assert result.counts() == {
        "n_records": 4,
        "n_unmapped": 3,
        "n_outside_contigs": 0,
        "n_off_contig_end": 0,
        "n_over_max_mismatches": 0,
        "n_repeated": 0,
        "n_nm_disagreements": 0,
        "n_sites": 1,
    }


def test_the_contig_filter(tmp_path):
    contigs = {"chr1": SPACER_A + "TGG", "chrUn": SPACER_A + "TGG"}
    records = [hit("0:NGG", "chr1", 0, "+", 23), hit("0:NGG", "chrUn", 0, "+", 23)]
    guides = library(SPACER_A)

    result = run(tmp_path, contigs, records, guides, contigs={"chr1"})

    assert [site.chr for site in result.sites[0]] == ["chr1"]
    assert result.n_outside_contigs == 1


@pytest.mark.parametrize(
    ("contigs", "error", "message"),
    [
        ({"chr1", "chrX"}, ValueError, "not in the reference: chrX"),
        (set(), ValueError, "at least one contig"),
        ("chr1", TypeError, "not a string"),
    ],
)
def test_the_contig_filter_rejects(tmp_path, contigs, error, message):
    with pytest.raises(error, match=message):
        run(
            tmp_path,
            {"chr1": SPACER_A + "TGG"},
            [hit("0:NGG", "chr1", 0, "+", 23)],
            library(SPACER_A),
            contigs=contigs,
        )


def test_sites_are_sorted_in_reference_order(tmp_path):
    contigs = {
        "chrB": SPACER_A + "TGG" + "A" * 5 + SPACER_A + "TGG",
        "chrA": SPACER_A + "TGG",
    }

    result = run(
        tmp_path,
        contigs,
        [
            hit("0:NGG", "chrA", 0, "+", 23),
            hit("0:NGG", "chrB", 28, "+", 23),
            hit("0:NGG", "chrB", 0, "+", 23),
        ],
        library(SPACER_A),
    )

    assert [(site.chr, site.start) for site in result.sites[0]] == [
        ("chrB", 0),
        ("chrB", 28),
        ("chrA", 0),
    ]


@pytest.mark.parametrize(
    ("sam_lengths", "message"),
    [
        ({"chr1": 99}, "lengths differ: chr1 \\(99 vs 23\\)"),
        ({"chr1": 23, "chr9": 50}, "only in the SAM: chr9"),
        ({}, "no @SQ header lines"),
    ],
)
def test_the_sam_header_must_match_the_reference(tmp_path, sam_lengths, message):
    fasta = write_fasta(tmp_path / "ref.fa", {"chr1": SPACER_A + "TGG"})
    sam = tmp_path / "sites.sam"
    if sam_lengths:
        write_sam(sam, sam_lengths, [unmapped("0:NGG")])
    else:
        sam.write_text("@HD\tVN:1.6\n0:NGG\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")

    with pytest.raises(ValueError, match=message):
        read_sites(sam, library(SPACER_A), fasta)


def test_a_reference_contig_missing_from_the_sam_is_an_error(tmp_path):
    fasta = write_fasta(tmp_path / "ref.fa", {"chr1": "ACGT", "chr2": "ACGT"})
    sam = write_sam(tmp_path / "sites.sam", {"chr1": 4}, [unmapped("0:NGG")])

    with pytest.raises(ValueError, match="only in the reference: chr2"):
        read_sites(sam, library(SPACER_A), fasta)


def test_nm_disagreements_are_counted(tmp_path):
    contig = Contig()
    first = contig.add(SPACER_A + "TGG")
    second = contig.add(mutate(SPACER_A, [7]) + "TGG")
    third = contig.add(SPACER_A + "AGG")

    result = run(
        tmp_path,
        {"chr1": contig.sequence},
        [
            hit("0:NGG", "chr1", first, "+", 23, nm=0),  # the N counts: 1
            hit("0:NGG", "chr1", second, "+", 23, nm=2),  # 1 + the N: 2
            hit("0:NGG", "chr1", third, "+", 23),  # no NM tag: not compared
        ],
        library(SPACER_A),
    )

    assert result.n_nm_disagreements == 1
    assert len(result.sites[0]) == 3


@pytest.mark.parametrize("cigar", ["20M3S", "3S20M", "10M1I12M", "10M1D13M"])
def test_a_gapped_or_clipped_alignment_is_an_error(tmp_path, cigar):
    with pytest.raises(ValueError, match="gapless, unclipped"):
        run(
            tmp_path,
            {"chr1": "A" * 5 + SPACER_A + "TGG"},
            [hit("0:NGG", "chr1", 5, "+", 23, cigar=cigar)],
            library(SPACER_A),
        )


def test_an_alignment_of_the_wrong_length_is_an_error(tmp_path):
    # The guide carries an added G, so its reads are 24 bases long.
    with pytest.raises(ValueError, match="all 24 bases"):
        run(
            tmp_path,
            {"chr1": "A" * 5 + SPACER_A + "TGG"},
            [hit("0:NGG", "chr1", 5, "+", 23)],
            library(SPACER_A, add_leading_g=True),
        )


@pytest.mark.parametrize(
    ("name", "message"),
    [
        ("g1", "not <row>:<pam>"),
        ("0:NNN", "not <row>:<pam>"),
        ("x:NGG", "does not carry a guide row"),
        ("5:NGG", "library has 2 guides"),
        ("1:NGG", "a duplicate of 'g1'"),
    ],
)
def test_read_names_must_match_the_library(tmp_path, name, message):
    with pytest.raises(ValueError, match=message):
        run(
            tmp_path,
            {"chr1": SPACER_A + "TGG"},
            [hit(name, "chr1", 0, "+", 23)],
            library(SPACER_A, SPACER_A),
        )


def test_several_sam_files_are_read_into_one_result(tmp_path):
    guides = library(SPACER_G, SPACER_C, add_leading_g=True)
    contig = Contig()
    first = contig.add(SPACER_G + "TGG")
    second = contig.add("G" + SPACER_C + "TGG")
    fasta = write_fasta(tmp_path / "ref.fa", {"chr1": contig.sequence})
    lengths = {"chr1": len(contig.sequence)}
    without_g = write_sam(
        tmp_path / "without_g.sam",
        lengths,
        [hit("0:NGG", "chr1", first, "+", 23), unmapped("0:NAG")],
    )
    with_g = write_sam(
        tmp_path / "with_g.sam", lengths, [hit("1:NGG", "chr1", second, "+", 24)]
    )

    result = read_sites([without_g, with_g], guides, fasta)

    assert [site.start for site in result.sites[0]] == [0]
    assert [site.leading_g_base for site in result.sites[1]] == ["G"]
    assert (result.n_records, result.n_unmapped) == (3, 1)


def test_bam_and_gzipped_sam_input(tmp_path):
    contig = "A" * 5 + SPACER_A + "TGG"
    fasta = write_fasta(tmp_path / "ref.fa", {"chr1": contig})
    records = [hit("0:NGG", "chr1", 5, "+", 23)]
    bam = write_sam(tmp_path / "sites.bam", {"chr1": len(contig)}, records)
    sam = write_sam(tmp_path / "sites.sam", {"chr1": len(contig)}, records)
    with open(sam, "rb") as source, gzip.open(tmp_path / "sites.sam.gz", "wb") as out:
        shutil.copyfileobj(source, out)

    for path in (bam, tmp_path / "sites.sam.gz"):
        result = read_sites(path, library(SPACER_A), fasta)
        assert [site.start for site in result.sites[0]] == [5]


def test_the_fasta_index_can_live_elsewhere(tmp_path):
    contig = "A" * 5 + SPACER_A + "TGG"
    fasta = write_fasta(tmp_path / "ref.fa", {"chr1": contig})
    sam = write_sam(
        tmp_path / "sites.sam",
        {"chr1": len(contig)},
        [hit("0:NGG", "chr1", 5, "+", 23)],
    )
    index = tmp_path / "outdir" / "ref.fa.fai"
    index.parent.mkdir()

    result = read_sites(sam, library(SPACER_A), fasta, fasta_index=index)

    assert [site.start for site in result.sites[0]] == [5]
    assert index.is_file()
    assert not (tmp_path / "ref.fa.fai").exists()


@pytest.mark.parametrize("max_mismatches", [-1, 1.5, True])
def test_max_mismatches_must_be_a_non_negative_integer(tmp_path, max_mismatches):
    with pytest.raises(ValueError, match="non-negative integer"):
        run(
            tmp_path,
            {"chr1": SPACER_A + "TGG"},
            [],
            library(SPACER_A),
            max_mismatches=max_mismatches,
        )

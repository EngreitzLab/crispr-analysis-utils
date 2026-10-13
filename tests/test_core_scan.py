"""The Rust scan finds exactly what a naive Python scan finds.

The oracle below is deliberately slow and obvious: every position, every
strand, compared base by base. It is the only thing standing between us and
trusting a search we cannot see inside.
"""

import random

import pytest

from crispr_analysis_utils import _core

IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def revcomp(sequence):
    return sequence[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def oracle(contigs, spacers, pam, upstream):
    """Every exact match of every spacer, by brute force."""
    found = []
    for name, sequence in contigs.items():
        text = sequence.upper()
        for index, spacer in enumerate(spacers):
            length = len(spacer)
            for strand, wanted in (("+", spacer), ("-", revcomp(spacer))):
                for start in range(len(text) - length + 1):
                    if text[start : start + length] != wanted:
                        continue
                    if strand == "+":
                        pam_bases = text[start + length : start + length + len(pam)]
                        up = text[max(0, start - upstream) : start]
                        if start - upstream < 0:
                            up = ""
                    else:
                        pam_bases = revcomp(text[max(0, start - len(pam)) : start])
                        if start - len(pam) < 0:
                            pam_bases = ""
                        up = revcomp(text[start + length : start + length + upstream])
                        if start + length + upstream > len(text):
                            up = ""
                    if len(pam_bases) != len(pam):
                        pam_bases = ""
                    matched = len(pam_bases) == len(pam) and all(
                        base in IUPAC[code]
                        for base, code in zip(pam_bases, pam, strict=True)
                    )
                    found.append((index, name, start, strand, pam_bases, up, matched))
    return sorted(found)


def write_fasta(path, contigs, width=60):
    with open(path, "w") as handle:
        for name, sequence in contigs.items():
            handle.write(f">{name} some description\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")


def run(tmp_path, contigs, spacers, pam="NGG", upstream=1):
    path = tmp_path / "genome.fa"
    write_fasta(path, contigs)
    got = sorted(_core.scan_fasta(str(path), spacers, pam, upstream))
    assert got == oracle(contigs, spacers, pam, upstream)
    return got


def random_sequence(rng, length, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(length))


def test_a_planted_site_is_found(tmp_path):
    rng = random.Random(0)
    spacer = "ACGTTGCAAGCTTCGATCGA"
    sequence = random_sequence(rng, 2000)
    sequence = sequence[:500] + spacer + "TGG" + sequence[523:]
    got = run(tmp_path, {"chrA": sequence}, [spacer])
    assert (0, "chrA", 500, "+", "TGG", sequence[499].upper(), True) in got


def test_a_site_on_the_minus_strand(tmp_path):
    rng = random.Random(1)
    spacer = "ACGTTGCAAGCTTCGATCGA"
    # Plant the reverse complement, with a CCA before it: in guide orientation
    # that is a TGG PAM 3' of the protospacer.
    planted = "CCA" + revcomp(spacer)
    sequence = random_sequence(rng, 2000)
    sequence = sequence[:800] + planted + sequence[800 + len(planted) :]
    got = run(tmp_path, {"chrA": sequence}, [spacer])
    assert any(hit[3] == "-" and hit[4] == "TGG" and hit[6] for hit in got)


@pytest.mark.parametrize("seed", range(6))
def test_random_genomes_agree_with_the_oracle(tmp_path, seed):
    rng = random.Random(seed)
    contigs = {
        f"chr{index}": random_sequence(rng, rng.randint(200, 900)) for index in range(3)
    }
    # Short spacers so matches actually happen by chance, which is the point.
    spacers = [random_sequence(rng, 6) for _ in range(5)]
    joined = "".join(contigs.values())
    spacers.append(joined[100:106])
    run(tmp_path, contigs, spacers, pam="NGG", upstream=1)


@pytest.mark.parametrize("pam", ["NGG", "NAG", "NRG", "TTTV", "N"])
def test_pam_patterns(tmp_path, pam):
    rng = random.Random(42)
    contigs = {"chrA": random_sequence(rng, 1500)}
    spacers = [contigs["chrA"][200:206], random_sequence(rng, 6)]
    run(tmp_path, contigs, spacers, pam=pam, upstream=2)


def test_sites_at_the_contig_ends(tmp_path):
    spacer = "ACGTAC"
    # One at the very start, one at the very end: their flanks run off.
    sequence = spacer + "TGG" + "AAAACCCCGGGGTTTT" + "CCA" + revcomp(spacer)
    run(tmp_path, {"chrEdge": sequence}, [spacer], upstream=1)


def test_n_runs_never_match(tmp_path):
    spacer = "ACGTAC"
    sequence = "N" * 50 + spacer + "TGGNNNN" + "N" * 50
    got = run(tmp_path, {"chrN": sequence}, [spacer])
    assert len(got) == 1


def test_lowercase_reference_matches(tmp_path):
    spacer = "ACGTTGCAAGCTTCGATCGA"
    sequence = "aaaa" + spacer.lower() + "tgg" + "cccc"
    got = run(tmp_path, {"chrSoft": sequence}, [spacer])
    assert got[0][4] == "TGG"


def test_a_palindromic_spacer_reads_on_both_strands(tmp_path):
    # Its own reverse complement, so the same bases are a guide either way --
    # with different flanks, so only one reading has a PAM.
    spacer = "ACGTACGT"
    assert revcomp(spacer) == spacer
    sequence = "TTTT" + spacer + "TGGAAAA"
    got = run(tmp_path, {"chrP": sequence}, [spacer])
    both = {hit[3]: hit for hit in got if hit[2] == 4}
    assert set(both) == {"+", "-"}
    assert both["+"][4] == "TGG" and both["+"][6]
    assert not both["-"][6]


def test_two_guides_sharing_a_spacer(tmp_path):
    spacer = "ACGTTGCAAGCTTCGATCGA"
    sequence = "TTTT" + spacer + "TGGAAAA"
    got = run(tmp_path, {"chrS": sequence}, [spacer, spacer])
    assert sorted(hit[0] for hit in got) == [0, 1]


def test_errors():
    with pytest.raises(ValueError, match="same length"):
        _core.scan_fasta("x", ["ACGT", "ACGTA"], "NGG", 1)
    with pytest.raises(ValueError, match="outside ACGT"):
        _core.scan_fasta("x", ["ACGN"], "NGG", 1)
    with pytest.raises(ValueError, match="not an IUPAC code"):
        _core.scan_fasta("x", ["ACGT"], "NGZ", 1)
    with pytest.raises(ValueError, match="must not be empty"):
        _core.scan_fasta("x", ["ACGT"], "", 1)
    with pytest.raises(OSError, match="Cannot read"):
        _core.scan_fasta("/nonexistent/genome.fa", ["ACGT"], "NGG", 1)


def test_contig_filter(tmp_path):
    spacer = "ACGTACGTACGTACGTACGT"
    contigs = {"keep": "TT" + spacer + "TGG", "drop": "AA" + spacer + "AGG"}
    path = tmp_path / "genome.fa"
    write_fasta(path, contigs)
    got = _core.scan_fasta(str(path), [spacer], "NGG", 1, ["keep"])
    assert {hit[1] for hit in got} == {"keep"}

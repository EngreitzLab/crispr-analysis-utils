"""The guide alignment pipeline, run on a small hand-planted reference.

There is nothing to mock: `run` scans a FASTA, so every test plants the sites
it wants and reads the tables back.
"""

import csv
import gzip
import json
import logging
import random

import pandas as pd
import pytest
from guide_alignment_fixtures import (
    SPACER_A,
    SPACER_C,
    SPACER_G,
    Contig,
    write_fasta,
)

from crispr_analysis_utils.guide_alignment import (
    OrientationError,
    read_guides,
    run,
    scan_guides,
)
from crispr_analysis_utils.guide_alignment.iupac import reverse_complement
from crispr_analysis_utils.guide_alignment.pipeline import METADATA_DIR, STEP

SPACER_NT = "TTAGCGCATAGCGTACGTAA"
GUIDES = [
    ("g_a", SPACER_A),
    ("g_c", SPACER_C),
    ("g_g", SPACER_G),
    ("g_nt", SPACER_NT),
    ("g_dup", SPACER_A),
]
OUTPUTS = (
    "guides_alignment_metadata.tsv",
    "guides_mapping.bed",
    "columns.yaml",
    "qc.tsv",
)
SITE_COLUMNS = (
    "chr",
    "start",
    "end",
    "strand",
    "cut_site",
    "pam_matches",
    "mapped_sequence",
)


def mutate(spacer, position):
    """Change the base at this 1-based position."""
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    return spacer[: position - 1] + swap[spacer[position - 1]] + spacer[position:]


def build_reference(path):
    """Two contigs: g_a once (+), g_c once (-), g_g once on each contig.

    g_c also has a 1-mismatch NAG site, which an exact scan never reports.
    Every site has TTT (or AAA) on its 5' side, so no flank of a
    reverse-complemented spacer reads as a PAM.
    """
    rng = random.Random(7)

    def filler(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    sites = {}
    chr1 = Contig()
    chr1.add(filler(40) + "TTT")
    sites["g_a"] = ("chr1", chr1.add(SPACER_A + "TGG"), "+")
    chr1.add(filler(40) + "TTT")
    # add_minus writes the PAM first, so the protospacer starts 3 bases later.
    sites["g_c"] = ("chr1", chr1.add_minus(SPACER_C + "AGG") + 3, "-")
    chr1.add("TTT" + filler(40) + "TTT")
    sites["g_c_nag"] = ("chr1", chr1.add(mutate(SPACER_C, 10) + "TAG"), "+")
    chr1.add(filler(40) + "TTT")
    sites["g_g_1"] = ("chr1", chr1.add(SPACER_G + "CGG"), "+")
    chr1.add(filler(40))
    chr2 = Contig()
    chr2.add(filler(30) + "TTT")
    sites["g_g_2"] = ("chr2", chr2.add(SPACER_G + "AGG"), "+")
    chr2.add(filler(30))
    path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(path, {"chr1": chr1.sequence, "chr2": chr2.sequence})
    return sites


def write_guides(path, guides=GUIDES):
    path.write_text("".join(f"{gid}\t{spacer}\n" for gid, spacer in guides))
    return path


@pytest.fixture
def setup(tmp_path):
    reference = tmp_path / "ref" / "genome.fa"
    sites = build_reference(reference)
    return {
        "reference": reference,
        "guides": write_guides(tmp_path / "guides.tsv"),
        "sites": sites,
        "outdir": tmp_path / "out",
    }


def go(setup, **kwargs):
    guides = kwargs.pop("guides", setup["guides"])
    outdir = kwargs.pop("outdir", setup["outdir"])
    return run(guides, setup["reference"], outdir, **kwargs)


def table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def by_alias(rows):
    """The rows grouped by alias: a guide with two sites has two rows."""
    grouped = {}
    for row in rows:
        grouped.setdefault(row["alias"], []).append(row)
    return grouped


def record(setup):
    """The most recent invocation record in the output folder."""
    records = sorted((setup["outdir"] / METADATA_DIR).glob(f"{STEP}_*.json"))
    return json.loads(records[-1].read_text())


def test_run_writes_every_output(setup):
    result = go(setup)
    out = setup["outdir"]

    # Results, metadata and logs each get their own folder, so nothing the run
    # produced is left sitting beside the inputs.
    for name in OUTPUTS:
        assert (out / "alignments" / name).is_file(), name
    assert (out / "logs" / "guide_alignment.log").is_file()
    records = list((out / METADATA_DIR).glob(f"{STEP}_*.json"))
    assert len(records) == 1
    assert sorted(p.name for p in out.iterdir()) == ["alignments", "logs", METADATA_DIR]

    guides = by_alias(table(out / "alignments" / "guides_alignment_metadata.tsv"))
    chrom, a_start, _ = setup["sites"]["g_a"]
    (g_a,) = guides["g_a"]
    assert g_a["site_class"] == "unique"
    assert (g_a["chr"], g_a["start"]) == (chrom, str(a_start))
    assert g_a["end"] == str(a_start + 20)
    assert g_a["strand"] == "+"
    assert g_a["pam_matches"] == "TRUE"
    assert g_a["cut_site"] == str(a_start + 17)
    assert g_a["guide_id"] == SPACER_A + "NGG"
    # The genomic PAM is the tail of mapped_sequence; the pattern ends guide_id.
    assert g_a["mapped_sequence"] == SPACER_A + "TGG"
    assert (g_a["n_matches"], g_a["n_cut_sites"]) == ("1", "1")
    assert result.summaries[0].site_pam == "TGG"
    assert result.summaries[0].has_perfect_site is True
    assert result.summaries[0].extra_5p is False
    _, c_start, _ = setup["sites"]["g_c"]
    (g_c,) = guides["g_c"]
    assert g_c["site_class"] == "unique"
    assert (g_c["start"], g_c["strand"]) == (str(c_start), "-")
    assert g_c["cut_site"] == str(c_start + 3)
    # Two sites, so two rows, with the guide's own columns repeated on each.
    assert {row["site_class"] for row in guides["g_g"]} == {"multi"}
    assert {row["n_cut_sites"] for row in guides["g_g"]} == {"2"}
    assert len(guides["g_g"]) == 2
    (g_nt,) = guides["g_nt"]
    assert g_nt["site_class"] == "no_site"
    assert [g_nt[name] for name in SITE_COLUMNS] == [""] * 7
    (g_dup,) = guides["g_dup"]
    assert g_dup["guide_id"] == g_a["guide_id"]
    assert g_dup["start"] == str(a_start)
    assert {g: {row["orientation"] for row in rows} for g, rows in guides.items()} == {
        "g_a": {"forward"},
        "g_c": {"forward"},
        "g_g": {"forward"},
        "g_nt": {"no_match"},
        "g_dup": {"forward"},
    }
    assert [s.site_class for s in result.summaries] == [
        "unique",
        "unique",
        "multi",
        "no_site",
        "unique",
    ]

    qc = {
        (r["section"], r["name"]): r["value"]
        for r in table(out / "alignments" / "qc.tsv")
    }
    assert qc[("guides", "n_entries")] == "5"
    assert qc[("guides", "n_guides")] == "4"
    assert qc[("class", "unique")] == "3"
    assert qc[("orientation", "forward")] == "4"
    assert qc[("scan", "n_sites")] == str(result.scan.n_sites)

    # The BED holds every site a guide could be cut at, named by its spacer.
    bed = [
        line.split("\t")
        for line in (out / "alignments" / "guides_mapping.bed").read_text().splitlines()
    ]
    assert [(row[0], row[1], row[3], row[5]) for row in bed] == [
        # One row per distinct guide: 'g_dup' is the same sequence as 'g_a' at
        # the same locus, so it must not appear twice.
        ("chr1", str(a_start), SPACER_A, "+"),
        ("chr1", str(c_start), SPACER_C, "-"),
        ("chr1", str(setup["sites"]["g_g_1"][1]), SPACER_G, "+"),
        ("chr2", str(setup["sites"]["g_g_2"][1]), SPACER_G, "+"),
    ]
    assert {int(row[2]) - int(row[1]) for row in bed} == {20}

    state = record(setup)
    assert state == result.metadata
    assert state["complete"] is True and state["error"] is None
    assert state["step"] == STEP
    assert state["counts"]["scan"] == result.scan.counts()
    assert state["classes"] == {"unique": 3, "multi": 1, "no_site": 1}
    assert state["orientation"] == {
        "forward": 4,
        "reverse_complement": 0,
        "no_pam": 0,
        "no_match": 1,
        "not_checked": 0,
    }
    assert state["parameters"]["pam"] == "NGG"
    assert state["parameters"]["orientation_check"] == "error"
    assert state["inputs"]["guides"] == str(setup["guides"])
    assert len(state["inputs"]["reference_md5"]) == 32
    assert state["orientation_fixed"] == 0


def test_the_scan_is_exact_so_a_mismatched_site_is_never_reported(setup):
    result = go(setup)

    _, nag_start, _ = setup["sites"]["g_c_nag"]
    rows = table(setup["outdir"] / "alignments" / "guides_alignment_metadata.tsv")
    assert str(nag_start) not in {row["start"] for row in rows}
    guides = by_alias(rows)
    assert guides["g_c"][0]["n_cut_sites"] == "1"
    assert guides["g_c"][0]["n_matches"] == "1"
    assert [s.min_mismatches for s in result.summaries if s.alias == "g_c"] == [0]
    assert {
        site.n_mismatches for sites in result.scan.sites.values() for site in sites
    } == {0}


def test_entries_sharing_a_sequence_are_scanned_once_and_written_per_entry(setup):
    result = go(setup)

    # The sequence is searched once, however many entries carry it: four
    # distinct spacers, four genomic matches between them.
    assert result.scan.n_matches == 4
    guides = table(setup["outdir"] / "alignments" / "guides_alignment_metadata.tsv")
    shared = [row for row in guides if row["guide_id"] == SPACER_A + "NGG"]
    # One row per entry per site: both entries get the one site of the sequence.
    assert sorted(row["alias"] for row in shared) == ["g_a", "g_dup"]
    assert {row["start"] for row in shared} == {str(setup["sites"]["g_a"][1])}
    assert result.guides[4].duplicate_of == "g_a"


def test_add_leading_g_records_the_base_under_the_g(setup):
    result = go(setup, add_leading_g=True)

    _, a_start, _ = setup["sites"]["g_a"]
    _, c_start, _ = setup["sites"]["g_c"]
    # TTT precedes g_a on +; TTT follows g_c's site on +, so A in guide orientation.
    bases = {
        (site.guide_id, site.start): site.genomic_5p_base
        for sites in result.scan.sites.values()
        for site in sites
    }
    assert bases[("G" + SPACER_A + "NGG", a_start)] == "T"
    assert bases[("G" + SPACER_C + "NGG", c_start)] == "A"
    guides = by_alias(
        table(setup["outdir"] / "alignments" / "guides_alignment_metadata.tsv")
    )
    assert result.guides[0].extra_5p is True
    assert guides["g_a"][0]["guide_id"] == "G" + SPACER_A + "NGG"
    # The G is dropped from mapped_sequence because the genome has a T there.
    assert guides["g_a"][0]["mapped_sequence"] == SPACER_A + "TGG"
    # The coordinates are the protospacer's: the added G is outside them.
    assert guides["g_a"][0]["start"] == str(a_start)


def test_contigs_keep_only_their_hits(setup):
    result = go(setup, contigs=["chr1"])

    guides = by_alias(
        table(setup["outdir"] / "alignments" / "guides_alignment_metadata.tsv")
    )
    (g_g,) = guides["g_g"]
    assert g_g["site_class"] == "unique"
    assert g_g["chr"] == "chr1"
    assert result.metadata["parameters"]["contigs"] == ["chr1"]


def test_an_unknown_contig_is_an_error_not_an_empty_scan(setup):
    # A typo here used to leave every guide no_site and the run marked
    # complete, which reads exactly like a library that maps nowhere.
    with pytest.raises(ValueError, match=r'has no contig named "chrUn"'):
        go(setup, contigs=["chr1", "chrUn"])


def test_a_pam_running_off_a_contig_end_is_counted_not_reported(tmp_path):
    reference = tmp_path / "ends.fa"
    write_fasta(
        reference,
        # One match ending at the last base on +, one starting at base 0 on -.
        {"c1": "AAAA" + SPACER_A, "c2": reverse_complement(SPACER_A) + "AAAA"},
    )
    result = run(
        write_guides(tmp_path / "g.tsv", [("g_a", SPACER_A)]), reference, "out"
    )

    counts = result.scan.counts()
    assert (counts["n_matches"], counts["n_off_contig_end"], counts["n_sites"]) == (
        2,
        2,
        0,
    )
    # Nothing was kept, so the guide gets its one row with no site in it.
    (row,) = table(result.outdir / "alignments" / "guides_alignment_metadata.tsv")
    assert [row[name] for name in SITE_COLUMNS] == [""] * 7
    assert (row["n_matches"], row["n_cut_sites"]) == ("0", "0")
    assert result.summaries[0].site_class == "no_site"
    assert result.summaries[0].orientation == "no_match"


def test_a_match_without_a_pam_is_kept_out_of_the_bed_but_written_to_the_tsv(tmp_path):
    reference = write_fasta(tmp_path / "nopam.fa", {"c1": "TTTT" + SPACER_A + "TTTTTT"})
    result = run(
        write_guides(tmp_path / "g.tsv", [("g_a", SPACER_A)]), reference, "out"
    )

    (site,) = result.scan.sites[0]
    assert (site.pam, site.pam_class) == ("TTT", "none")
    (row,) = table(result.outdir / "alignments" / "guides_alignment_metadata.tsv")
    assert (row["chr"], row["start"], row["pam_matches"]) == ("c1", "4", "FALSE")
    assert (row["n_matches"], row["n_cut_sites"]) == ("1", "0")
    assert row["mapped_sequence"] == SPACER_A + "TTT"
    # Nothing can be cut there, so the BED is empty.
    assert (result.outdir / "alignments" / "guides_mapping.bed").read_text() == ""
    assert result.scan.counts()["n_no_pam"] == 1
    assert result.scan.counts()["n_sites"] == 1
    assert result.summaries[0].site_class == "no_site"
    assert result.summaries[0].orientation == "no_pam"
    assert result.summaries[0].n_sites == (0,)


def test_a_spacer_that_matches_nowhere_is_no_match(setup):
    result = go(setup)

    (nt,) = [s for s in result.summaries if s.alias == "g_nt"]
    assert (nt.orientation, nt.site_class, nt.has_perfect_site) == (
        "no_match",
        "no_site",
        False,
    )
    assert nt.orientation_evidence is None


def test_a_reversed_entry_is_flagged_only_with_a_pam_bearing_reversed_site(tmp_path):
    contig = Contig()
    contig.add("TTT")
    contig.add(SPACER_A + "TGG")
    contig.add("TTT")
    # SPACER_C sits with TTT on both sides: neither reading finds a PAM.
    contig.add(SPACER_C + "TTT")
    reference = write_fasta(tmp_path / "ref.fa", {"chr1": contig.sequence})
    guides = read_guides(
        pd.DataFrame(
            [
                ("flagged", reverse_complement(SPACER_A)),
                ("no_pam", reverse_complement(SPACER_C)),
                ("absent", reverse_complement(SPACER_NT)),
            ]
        )
    )

    scan = scan_guides(guides, reference)

    # Only 'flagged' has a reversed reading with a PAM 3' of it.
    assert scan.verdicts == {
        0: "reverse_complement",
        1: "no_pam",
        2: "no_match",
    }


def test_orientation_evidence_is_ambiguous_without_a_base_5prime_of_the_spacer(
    tmp_path,
):
    reference = write_fasta(
        tmp_path / "ref.fa", {"chr1": "TTT" + SPACER_A + "TGG" + "TTT"}
    )
    bare = read_guides(pd.DataFrame([("bare", reverse_complement(SPACER_A))]))
    # A 21-base entry: its reversed reading is a different 20-mer, so the test
    # is independent of what happens to sit 5' of the match.
    longer = read_guides(pd.DataFrame([("longer", reverse_complement("G" + SPACER_A))]))

    assert bare[0].n_extra == 0
    assert longer[0].n_extra == 1
    assert scan_guides(bare, reference).evidence == {0: "ambiguous"}
    assert scan_guides(longer, reference).evidence == {0: "distinct"}
    assert scan_guides(longer, reference).verdicts == {0: "reverse_complement"}


def test_a_run_stops_only_when_reversed_entries_outnumber_forward_ones(setup, tmp_path):
    def reversed_rows(n):
        return [(gid, reverse_complement(spacer)) for gid, spacer in GUIDES[:n]] + [
            GUIDES[n]
        ]

    # Two reversed against one forward: the table reads 3' to 5' throughout.
    with pytest.raises(OrientationError, match=r"2 of 3 entries read 3' to 5'"):
        go(
            setup,
            guides=write_guides(tmp_path / "most.tsv", reversed_rows(2)),
            outdir=tmp_path / "most",
        )

    # One against one is not a majority, so the run carries on.
    tie = go(
        setup,
        guides=write_guides(tmp_path / "tie.tsv", reversed_rows(1)),
        outdir=tmp_path / "tie",
    )
    assert [s.orientation for s in tie.summaries] == ["reverse_complement", "forward"]
    assert (tmp_path / "tie" / "alignments" / "guides_alignment_metadata.tsv").is_file()


def test_a_reversed_library_under_error_leaves_a_failed_record(setup, tmp_path):
    reversed_guides = write_guides(
        tmp_path / "reversed.tsv",
        [(gid, reverse_complement(spacer)) for gid, spacer in GUIDES[:3]],
    )

    with pytest.raises(OrientationError, match="Pass --fix-orientation"):
        go(setup, guides=reversed_guides, orientation_check="error")

    out = setup["outdir"]
    # The tables are written after the check, so a stopped run leaves none.
    for name in OUTPUTS:
        assert not (out / "alignments" / name).exists(), name
    state = record(setup)
    assert state["complete"] is False
    assert state["error"].startswith("OrientationError")
    assert "read 3' to 5'" in (out / "logs" / "guide_alignment.log").read_text()


def test_a_reversed_library_is_only_logged_when_the_check_is_warn(setup, tmp_path):
    reversed_guides = write_guides(
        tmp_path / "reversed.tsv",
        [(gid, reverse_complement(spacer)) for gid, spacer in GUIDES[:3]],
    )

    result = go(setup, guides=reversed_guides, orientation_check="warn")

    # The warning goes to the run's own log, not to the terminal.
    log = (setup["outdir"] / "logs" / "guide_alignment.log").read_text()
    assert "3 of 3 entries read 3' to 5'" in log
    assert {s.orientation for s in result.summaries} == {"reverse_complement"}
    assert {s.orientation_evidence for s in result.summaries} == {"ambiguous"}
    assert (setup["outdir"] / "alignments" / "guides_alignment_metadata.tsv").is_file()
    assert result.metadata["orientation_fixed"] == 0


def test_orientation_check_off_reports_no_verdict_at_all(setup, tmp_path, caplog):
    reversed_guides = write_guides(
        tmp_path / "reversed.tsv",
        [(gid, reverse_complement(spacer)) for gid, spacer in GUIDES[:3]],
    )

    with caplog.at_level(logging.WARNING, logger="crispr_analysis_utils"):
        result = go(setup, guides=reversed_guides, orientation_check="off")

    # The scan holds both readings whatever happens, but "off" means off: the
    # verdict column says not_checked and a table
    # that is reversed throughout neither warns nor stops.
    assert {s.orientation for s in result.summaries} == {"not_checked"}
    assert {s.orientation_evidence for s in result.summaries} == {None}
    assert "read 3' to 5'" not in caplog.text
    assert result.metadata["complete"] is True
    assert result.metadata["orientation"]["not_checked"] == 3


def test_fix_orientation_reverse_complements_only_the_flagged_entries(setup, tmp_path):
    rows = [("g_a", reverse_complement(SPACER_A)), ("g_c", SPACER_C), ("g_g", SPACER_G)]

    result = go(
        setup,
        guides=write_guides(tmp_path / "one.tsv", rows),
        fix_orientation=True,
    )

    assert [guide.input_sequence for guide in result.guides] == [
        SPACER_A,
        SPACER_C,
        SPACER_G,
    ]
    assert [s.orientation for s in result.summaries] == ["forward"] * 3
    assert [s.site_class for s in result.summaries] == ["unique", "unique", "multi"]
    assert result.metadata["orientation_fixed"] == 1
    assert "rescan" in result.metadata["seconds"]
    assert result.metadata["counts"]["scan"] == result.scan.counts()


def test_fix_orientation_does_nothing_when_no_entry_is_flagged(setup):
    result = go(setup, fix_orientation=True)

    assert result.metadata["orientation_fixed"] == 0
    assert "rescan" not in result.metadata["seconds"]


def test_a_dataframe_library_is_accepted(setup):
    frame = pd.DataFrame(GUIDES, columns=["id", "spacer"])

    result = go(setup, guides=frame, id_col="id", spacer_col="spacer")

    assert [s.alias for s in result.summaries] == [gid for gid, _ in GUIDES]
    assert result.metadata["inputs"]["guides"] == "DataFrame"


def test_the_log_is_appended_across_runs(setup):
    go(setup)
    go(setup)

    log = (setup["outdir"] / "logs" / "guide_alignment.log").read_text()
    assert log.count("Guide alignment into") == 2
    assert log.count("Done in") == 2


def test_nothing_is_written_outside_outdir(setup, tmp_path):
    def snapshot():
        return {
            path.relative_to(tmp_path)
            for path in tmp_path.rglob("*")
            if not path.is_relative_to(setup["outdir"])
        }

    before = snapshot()
    go(setup)

    assert snapshot() == before
    # No .fai either: nothing reads the FASTA through pysam any more.
    assert not list(tmp_path.rglob("*.fai"))


def test_a_failure_inside_the_run_is_recorded(setup):
    with pytest.raises(ValueError, match="positive integer"):
        go(setup, spacer_length=0)

    state = record(setup)
    assert state["complete"] is False
    assert state["error"].startswith("ValueError")
    assert state["finished"] is not None
    log = (setup["outdir"] / "logs" / "guide_alignment.log").read_text()
    assert "Guide alignment failed" in log


@pytest.mark.parametrize(
    ("kwargs", "error"),
    [
        ({"orientation_check": "maybe"}, ValueError),
        ({"contigs": "chr1"}, TypeError),
        ({"pam": "NGZ"}, ValueError),
        ({"pam": ""}, ValueError),
    ],
)
def test_run_rejects_bad_parameters_before_it_writes_anything(setup, kwargs, error):
    with pytest.raises(error):
        go(setup, **kwargs)

    assert not setup["outdir"].exists()


def test_missing_and_compressed_references_are_reported(setup, tmp_path):
    squashed = tmp_path / "genome.fa.gz"
    squashed.write_bytes(gzip.compress(setup["reference"].read_bytes()))

    with pytest.raises(FileNotFoundError, match="The reference"):
        run(setup["guides"], tmp_path / "nope.fa", setup["outdir"])
    with pytest.raises(ValueError, match="gzip-compressed"):
        run(setup["guides"], squashed, setup["outdir"])
    assert not setup["outdir"].exists()


def test_a_5prime_base_at_a_contig_edge_does_not_crash(setup, tmp_path):
    # The core returns an empty flank where it would run off the contig. Taking
    # the last n_extra bases of "" gave "" rather than None, and the zip that
    # compares them is strict, so the whole run died on a site it had kept.
    contig = Contig()
    start = contig.add(SPACER_A + "TGG")
    contig.add("AAAACCCCGGGGTTTTAAAA")
    reference = write_fasta(tmp_path / "edge.fa", {"edge": contig.sequence})
    assert start == 0

    result = run(
        write_guides(tmp_path / "edge.tsv", [("g_a", SPACER_A)]),
        reference,
        tmp_path / "edge_out",
        add_leading_g=True,
        orientation_check="off",
    )

    (row,) = result.summaries
    assert row.site_class == "unique"
    assert result.scan.sites[0][0].genomic_5p_base is None
    assert result.scan.sites[0][0].mapped_sequence == SPACER_A + "TGG"

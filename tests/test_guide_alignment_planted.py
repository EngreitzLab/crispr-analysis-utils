"""The pipeline against a brute-force scan of a synthetic genome, with real GEM3.

`planted_genome` builds the genome, its guides and the truth: every site with
at most three spacer mismatches and a matching PAM, and every exact copy of a
bare spacer. The pipeline must find exactly those, with the same coordinates,
classes, cut sites and orientation verdicts.

Two kinds of site are left out of the comparison, each for a reason GEM cannot
get around:

- the sites of `planted_genome.GEM_MISSES`, whose protospacer holds a genomic
  N that the complete search cannot step through, or that sit in an N run
  gem-indexer strips;
- with an added 5' G, every site whose G position is outside the index, which
  `planted_genome.g_position_indexed` decides.

The run is marked ``gem`` and skipped without gem-mapper on the PATH; the
project's pixi environments have it.
"""

import shutil

import planted_genome as P
import pytest

from crispr_analysis_utils.guide_alignment import GemError, build_index, gem, run
from crispr_analysis_utils.guide_alignment.summary import fan_out_sites

pytestmark = [
    pytest.mark.gem,
    pytest.mark.skipif(
        shutil.which("gem-mapper") is None, reason="GEM3 is not on the PATH"
    ),
]

CRASHES = P.GEM_CRASHES_WITH_LEADING_G


def site_key(guide_id, chrom, start, strand):
    return (guide_id, chrom, start, strand)


def truth_site(site):
    """A truth site as the fields the pipeline reports for it."""
    return {
        "end": site.end,
        "n_mismatches": site.n_mismatches,
        "mismatch_positions": tuple(p + 1 for p in site.mismatch_positions),
        "genomic_protospacer": site.protospacer,
        "pam": site.pam,
        "pam_class": next(p for p in P.PAMS if p in site.pams_matched),
    }


def found_site(site):
    return {
        "end": site.end,
        "n_mismatches": site.n_mismatches,
        "mismatch_positions": site.mismatch_positions,
        "genomic_protospacer": site.genomic_protospacer,
        "pam": site.pam,
        "pam_class": site.pam_class,
    }


def expected_sites(genome, *, leading_g, guides=None):
    """The truth sites GEM can reach, by key.

    Sites of a guide whose read carries the added G need their G position
    inside the index. A guide whose spacer already starts with G never gets
    one, so its sites are always in reach.
    """
    names = {g.name for g in (guides if guides is not None else genome.guides)}
    with_g = {
        g.name
        for g in genome.guides
        if leading_g and not g.spacer.startswith("G") and g.name in names
    }
    expected = {}
    for site in P.scan_sites(genome.contigs, genome.guides):
        if site.guide not in names or site.guide in P.GEM_MISSES:
            continue
        if site.guide in with_g and not P.g_position_indexed(site, genome.contigs):
            continue
        expected[site_key(site.guide, site.chrom, site.start, site.strand)] = (
            truth_site(site)
        )
    return expected


def expected_classes(expected):
    """Each guide's class and unique site, from the sites GEM can reach."""
    perfect = {}
    any_site = set()
    for (guide, chrom, start, strand), fields in expected.items():
        if fields["pam_class"] != P.PRIMARY_PAM:
            continue
        any_site.add(guide)
        if fields["n_mismatches"] == 0:
            perfect.setdefault(guide, []).append((chrom, start, fields["end"], strand))
    classes, unique = {}, {}
    for guide in {g for g, *_ in expected} | any_site:
        hits = perfect.get(guide, [])
        if len(hits) == 1:
            classes[guide] = "unique"
            unique[guide] = hits[0]
        elif hits:
            classes[guide] = "multi"
        elif guide in any_site:
            classes[guide] = "off_target_only"
    return classes, unique


def expected_verdicts(genome, names):
    """Each guide's orientation verdict, from its exact hits and their flanks."""
    hits = {}
    for hit in P.scan_exact(genome.contigs, genome.guides):
        hits.setdefault(hit.guide, []).append(hit)
    reversed_pams = [P.revcomp(pam) for pam in P.PAMS]
    verdicts = {}
    for name in names:
        guide_hits = hits.get(name, [])
        if not guide_hits:
            verdicts[name] = "no_perfect_hit"
        elif any(
            P.pam_matches(pam, hit.flank3) for hit in guide_hits for pam in P.PAMS
        ):
            verdicts[name] = "ok"
        elif any(
            P.pam_matches(pam, hit.flank5)
            for hit in guide_hits
            for pam in reversed_pams
        ):
            verdicts[name] = "reversed"
        else:
            verdicts[name] = "no_pam"
    return verdicts


def found_sites(result):
    return {
        site_key(site.guide_id, site.chr, site.start, site.strand): found_site(site)
        for site in fan_out_sites(result.guides, result.sites.sites)
        if site.pam_class != "none"
    }


def difference(found, expected):
    """A message naming the missed and the extra sites."""
    missed = sorted(set(expected) - set(found))
    extra = sorted(set(found) - set(expected))
    differ = sorted(
        key for key in set(found) & set(expected) if found[key] != expected[key]
    )
    return (
        f"missed {len(missed)}: {missed[:10]}\n"
        f"extra {len(extra)}: {extra[:10]}\n"
        f"differ {len(differ)}: "
        + "; ".join(f"{key} {found[key]} != {expected[key]}" for key in differ[:5])
    )


def write_guides(path, guides):
    path.write_text("".join(f"{g.name}\t{g.spacer}\n" for g in guides))
    return path


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    """The genome, its index and its guide table, built once for this module."""
    folder = tmp_path_factory.mktemp("planted")
    genome = P.build_genome()
    fasta = folder / "genome.fa"
    P.write_fasta(genome, fasta)
    index = build_index(fasta, folder / "index" / "genome", threads=1)
    return {
        "genome": genome,
        "fasta": fasta,
        "index": index,
        "folder": folder,
        "no_crash": [g for g in genome.guides if g.name not in CRASHES],
    }


def align(planted, tmp_path, guides, **kwargs):
    table = write_guides(tmp_path / "guides.tsv", guides)
    return run(
        table,
        planted["fasta"],
        planted["index"],
        tmp_path / "out",
        threads=1,
        **kwargs,
    )


def test_every_site_is_found_without_the_added_g(planted, tmp_path):
    genome = planted["genome"]
    result = align(planted, tmp_path, genome.guides)
    expected = expected_sites(genome, leading_g=False)
    found = found_sites(result)
    assert found == expected, difference(found, expected)
    assert len(expected) > 350

    classes, unique = expected_classes(expected)
    summaries = {s.guide_id: s for s in result.summaries}
    assert {
        name: summary.alignment_class
        for name, summary in summaries.items()
        if summary.alignment_class != "no_site"
    } == classes
    for name, (chrom, start, end, strand) in unique.items():
        summary = summaries[name]
        assert (summary.guide_chr, summary.guide_start, summary.guide_end) == (
            chrom,
            start,
            end,
        )
        assert summary.strand == strand
        assert summary.targeting is True
        assert summary.cut_site == (end - 3 if strand == "+" else start + 3)
    assert summaries["multi60"].alignment_class == "multi"
    assert summaries["multi60"].n_ngg[0] == 60
    assert summaries["nt_00"].alignment_class in ("no_site", "off_target_only")

    bed = {
        line.split("\t")[3]: line.split("\t")
        for line in (result.outdir / "cut_sites.bed").read_text().splitlines()
    }
    assert set(bed) == set(unique)
    for name, (chrom, start, end, strand) in unique.items():
        cut = end - 3 if strand == "+" else start + 3
        assert bed[name][:3] == [chrom, str(cut), str(cut + 1)]


def test_orientation_verdicts_match_the_exact_hits(planted, tmp_path):
    genome = planted["genome"]
    result = align(planted, tmp_path, genome.guides)
    names = [g.name for g in genome.guides]
    assert {s.guide_id: s.orientation for s in result.summaries} == expected_verdicts(
        genome, names
    )
    verdicts = {s.guide_id: s.orientation for s in result.summaries}
    assert verdicts["revcomp_guide"] == "reversed"
    assert verdicts["nopam_plus"] == "no_pam"
    assert verdicts["uniq_plus"] == "ok"
    assert verdicts["nt_00"] == "no_perfect_hit"
    hits = {key: len(hits) for key, hits in result.orientation.hits.items()}
    assert sum(hits.values()) == len(P.scan_exact(genome.contigs, genome.guides))


def test_every_reachable_site_is_found_with_the_added_g(planted, tmp_path):
    genome = planted["genome"]
    guides = planted["no_crash"]
    result = align(planted, tmp_path, guides, add_leading_g=True)
    expected = expected_sites(genome, leading_g=True, guides=guides)
    found = found_sites(result)
    assert found == expected, difference(found, expected)

    # Both FASTQs were mapped: the four guides whose spacer starts with G keep
    # their reads in the one without the added G.
    passes = result.passes
    assert set(passes) == {"orientation", "sites", "sites_leading_g"}
    assert passes["sites"]["n_reads"] == 4 * len(P.PAMS)
    assert passes["sites_leading_g"]["n_reads"] == (len(guides) - 4) * len(P.PAMS)

    sites = {
        (s.guide_id, s.chr, s.start, s.strand): s
        for s in fan_out_sites(result.guides, result.sites.sites)
    }
    for (name, chrom, start, strand), site in sites.items():
        guide = next(g for g in guides if g.name == name)
        if guide.spacer.startswith("G"):
            assert site.leading_g_base is None
        else:
            position = start - 1 if strand == "+" else start + len(guide.spacer)
            base = genome.contigs[chrom][position].upper()
            assert site.leading_g_base == (base if strand == "+" else P.revcomp(base))


def test_a_leading_g_at_the_start_of_the_index_fails_the_run(planted, tmp_path):
    """GEM crashes on a read whose G would sit before the index's first base."""
    genome = planted["genome"]
    guides = [g for g in genome.guides if g.name in CRASHES]
    assert guides, "the planted genome must keep a guide at the first contig's start"
    with pytest.raises(GemError) as error:
        align(planted, tmp_path, guides, add_leading_g=True, orientation_check="off")
    message = str(error.value)
    assert "Signal raised" in message
    assert "sites_leading_g.fastq" in message
    assert f"guide {guides[0].name!r}, with the added G" in message
    assert "first contig" in message
    assert not (tmp_path / "out" / "alignments" / "sites_leading_g.sam.gz").exists()
    # Without the G the same guide maps, and its site is found.
    result = align(planted, tmp_path, guides, orientation_check="off")
    assert found_sites(result) == expected_sites(genome, leading_g=False, guides=guides)


@pytest.mark.parametrize("mode", ["fast", "sensitive"])
def test_the_comparison_fails_under_gems_own_mapping_modes(
    planted, tmp_path, monkeypatch, mode
):
    """The test discriminates: GEM's incomplete searches miss planted sites."""
    monkeypatch.setattr(
        gem,
        "MAPPER_OPTIONS",
        tuple(
            f"--mapping-mode={mode}" if option.startswith("--mapping-mode") else option
            for option in gem.MAPPER_OPTIONS
        ),
    )
    result = align(planted, tmp_path, planted["genome"].guides)
    found = found_sites(result)
    expected = expected_sites(planted["genome"], leading_g=False)
    assert set(found) < set(expected)

"""A synthetic genome with planted guide sites, and a brute-force scan of it.

The genome, its guides and the planted sites are fixed by the seed. The scan
is the truth the GEM integration test compares the pipeline with:

- `scan_sites`: every site, on both strands, whose spacer has at most ``k``
  mismatches and whose complete 3' PAM matches one of the PAM patterns.
- `scan_exact`: every exact copy of the bare spacer, with its flanks.

Coordinates are 0-based and half-open on the + strand, and a site's are the
protospacer's (no PAM, no leading G). Sequences, mismatch positions (0-based
here), PAMs and the leading-G base are in guide orientation, 5' to 3'. This
module shares no code with the package, so the scan is an independent truth.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np

SPACER_LENGTH = 20
MAX_MISMATCHES = 3
PRIMARY_PAM = "NGG"
ALT_PAMS = ("NAG", "NGA")
PAMS = (PRIMARY_PAM, *ALT_PAMS)
SEED = 20261007
# gem-indexer replaces runs of at least this many N by a single separator
# (--strip-unknown-bases-threshold, default 50).
STRIP_N_THRESHOLD = 50

# Sites GEM 3.6 cannot report with the tuned options: the sites of these guides
# whose protospacer holds a genomic N. A genomic N plus other mismatches is out of
# the complete search's reach at any -e/-E, and a stripped N run is not indexed.
GEM_MISSES = {
    "genomicN_3": "genomic N in the protospacer plus 2 mismatches",
    "near_plus_nrun": "protospacer starts inside a stripped N run",
}
# Guides whose read with a leading G makes gem-mapper crash ("Signal raised"): the G
# would sit before the first base of the first contig, the start of the index.
GEM_CRASHES_WITH_LEADING_G = {"edge_first_plus0"}

# Pattern base -> genomic bases it accepts. A genomic N passes only at a pattern N.
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
    "N": "ACGTN",
}
_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn"
)
_CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[_b + 32] = _i


def revcomp(seq: str) -> str:
    """Reverse complement, IUPAC-aware, case-preserving."""
    return seq.translate(_COMPLEMENT)[::-1]


def pam_matches(pattern: str, genomic: str) -> bool:
    """True when ``genomic`` (uppercase, guide orientation) matches ``pattern``."""
    if len(genomic) != len(pattern):
        return False
    return all(b in IUPAC[p] for p, b in zip(pattern, genomic, strict=True))


def _encode(seq: str) -> np.ndarray:
    return _CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]


@dataclasses.dataclass(frozen=True)
class Guide:
    name: str
    category: str
    spacer: str


@dataclasses.dataclass(frozen=True)
class Planted:
    """A site written into the genome on purpose (the scanner is still the truth)."""

    guide: str
    chrom: str
    start: int
    end: int
    strand: str
    n_mismatches: int
    pam: str  # genomic, guide orientation; shorter than the pattern at a contig end
    g_base: str  # genomic base at the leading-G position; '' when outside the contig
    note: str


@dataclasses.dataclass(frozen=True)
class Site:
    guide: str
    chrom: str
    start: int
    end: int
    strand: str
    n_mismatches: int
    mismatch_positions: tuple[int, ...]
    protospacer: str
    pam: str
    pams_matched: tuple[str, ...]
    g_base: str  # '' when the leading-G position is outside the contig

    @property
    def key(self) -> tuple[str, str, int, str]:
        return (self.guide, self.chrom, self.start, self.strand)

    @property
    def pam_class(self) -> str:
        if PRIMARY_PAM in self.pams_matched:
            return "primary"
        return "alt" if self.pams_matched else "none"


@dataclasses.dataclass(frozen=True)
class ExactHit:
    guide: str
    chrom: str
    start: int
    end: int
    strand: str
    flank5: str  # guide orientation, truncated at a contig end
    flank3: str

    @property
    def key(self) -> tuple[str, str, int, str]:
        return (self.guide, self.chrom, self.start, self.strand)


@dataclasses.dataclass
class Genome:
    contigs: dict[str, str]  # FASTA order; soft-masked bases are lowercase
    descriptions: dict[str, str]
    guides: list[Guide]
    planted: list[Planted]
    pass1_planted: list[tuple[str, str, int, int, str, str]]


def _windows(
    contigs: dict[str, str], spacer: str, max_mismatches: int
) -> Iterable[tuple[str, int, str, str, int, int]]:
    """Yield (chrom, length, strand, strand_seq, j, start) for spacer windows <= k mm."""
    s = _encode(spacer.upper())
    width = len(s)
    for chrom, seq in contigs.items():
        plus = seq.upper()
        length = len(plus)
        n = length - width + 1
        if n <= 0:
            continue
        for strand, strand_seq in (("+", plus), ("-", revcomp(plus))):
            a = _encode(strand_seq)
            mm = np.zeros(n, dtype=np.int16)
            for i, code in enumerate(s):
                mm += a[i : i + n] != code
            for j in np.flatnonzero(mm <= max_mismatches).tolist():
                start = j if strand == "+" else length - j - width
                yield chrom, length, strand, strand_seq, j, start


def scan_sites(
    contigs: dict[str, str],
    guides: Sequence[Guide],
    pams: Sequence[str] = PAMS,
    max_mismatches: int = MAX_MISMATCHES,
) -> list[Site]:
    """Every site with <= ``max_mismatches`` spacer mismatches and a complete PAM.

    A genomic N (or any non-ACGT base) counts as a mismatch in the spacer, and
    fails the PAM except at a pattern N.
    """
    sites = []
    for guide in guides:
        spacer = guide.spacer.upper()
        width = len(spacer)
        for chrom, _length, strand, sseq, j, start in _windows(
            contigs, spacer, max_mismatches
        ):
            matched = tuple(
                p for p in pams if pam_matches(p, sseq[j + width : j + width + len(p)])
            )
            if not matched:
                continue
            proto = sseq[j : j + width]
            positions = tuple(i for i in range(width) if proto[i] != spacer[i])
            sites.append(
                Site(
                    guide=guide.name,
                    chrom=chrom,
                    start=start,
                    end=start + width,
                    strand=strand,
                    n_mismatches=len(positions),
                    mismatch_positions=positions,
                    protospacer=proto,
                    pam=sseq[j + width : j + width + max(len(p) for p in pams)],
                    pams_matched=matched,
                    g_base=sseq[j - 1] if j >= 1 else "",
                )
            )
    return sites


def scan_exact(
    contigs: dict[str, str], guides: Sequence[Guide], flank: int = 3
) -> list[ExactHit]:
    """Every exact occurrence of each bare spacer, on both strands, with flanks."""
    hits = []
    for guide in guides:
        width = len(guide.spacer)
        for chrom, _length, strand, sseq, j, start in _windows(
            contigs, guide.spacer, 0
        ):
            hits.append(
                ExactHit(
                    guide=guide.name,
                    chrom=chrom,
                    start=start,
                    end=start + width,
                    strand=strand,
                    flank5=sseq[max(0, j - flank) : j],
                    flank3=sseq[j + width : j + width + flank],
                )
            )
    return hits


class _Builder:
    """Mutable genome under construction (uppercase until the soft-mask pass)."""

    def __init__(self, seed: int):
        self.rng = np.random.default_rng(seed)
        self.seqs: dict[str, bytearray] = {}
        self.descriptions: dict[str, str] = {}
        self.guides: list[Guide] = []
        self.planted: list[Planted] = []
        self.pass1_planted: list[tuple[str, str, int, int, str, str]] = []
        self.lowercase: list[tuple[str, int, int]] = []
        self.reserved: dict[str, list[tuple[int, int]]] = {}
        self._slots: list[tuple[str, int]] = []

    # -- sequences -------------------------------------------------------------
    def bases(self, n: int, alphabet: str = "ACGT") -> str:
        letters = np.frombuffer(alphabet.encode(), dtype=np.uint8)
        return letters[self.rng.integers(0, len(alphabet), size=n)].tobytes().decode()

    def add_contig(self, name: str, length: int, description: str = "") -> None:
        self.seqs[name] = bytearray(self.bases(length).encode())
        self.descriptions[name] = description

    def put(self, chrom: str, pos: int, text: str) -> None:
        assert 0 <= pos and pos + len(text) <= len(self.seqs[chrom]), (chrom, pos)
        self.seqs[chrom][pos : pos + len(text)] = text.encode()

    def reserve(self, chrom: str, start: int, end: int) -> None:
        self.reserved.setdefault(chrom, []).append((start, end))

    def make_slots(self, chroms: Sequence[str], width: int = 64) -> None:
        slots = []
        for chrom in chroms:
            for pos in range(width, len(self.seqs[chrom]) - 2 * width, width):
                if any(
                    pos < e and s < pos + width for s, e in self.reserved.get(chrom, [])
                ):
                    continue
                slots.append((chrom, pos))
        order = self.rng.permutation(len(slots))
        self._slots = [slots[i] for i in order]

    def slot(self) -> tuple[str, int]:
        chrom, pos = self._slots.pop()
        return chrom, pos + 20

    # -- guides and sites --------------------------------------------------------
    def guide(
        self, name: str, category: str, spacer: str | None = None, first: str = "H"
    ) -> Guide:
        """``first``: 'G' forces a 5' G, 'H' forces a non-G 5' base."""
        if spacer is None:
            spacer = self.bases(SPACER_LENGTH)
            if first == "G":
                spacer = "G" + spacer[1:]
            elif spacer[0] == "G":
                spacer = self.bases(1, "ACT") + spacer[1:]
        assert name not in {g.name for g in self.guides}, name
        g = Guide(name, category, spacer)
        self.guides.append(g)
        return g

    def mutate(self, spacer: str, positions: Iterable[int]) -> str:
        s = list(spacer)
        for p in positions:
            s[p] = self.bases(1, "ACGT".replace(spacer[p], ""))
        return "".join(s)

    def positions(self, n: int) -> list[int]:
        return sorted(self.rng.choice(SPACER_LENGTH, size=n, replace=False).tolist())

    def pam(self, pattern: str | None) -> str:
        """A concrete PAM for ``pattern``; ``None`` gives a 3-mer matching no PAM."""
        if pattern is None:
            while True:
                flank = self.bases(3)
                if not any(pam_matches(p, flank) for p in PAMS):
                    return flank
        return "".join(self.bases(1, IUPAC[p].replace("N", "")) for p in pattern)

    def g_base(self, match: bool) -> str:
        return "G" if match else self.bases(1, "ACT")

    def plant(
        self,
        guide: Guide,
        chrom: str,
        start: int,
        strand: str,
        *,
        n_mm: int = 0,
        positions: Sequence[int] | None = None,
        pam: str | None = PRIMARY_PAM,
        pam_seq: str | None = None,
        g_match: bool = True,
        g_base: str | None = None,
        protospacer: str | None = None,
        note: str = "",
    ) -> Planted:
        """Write [G] + protospacer + PAM (guide orientation) with the protospacer at start."""
        width = len(guide.spacer)
        if protospacer is None:
            if positions is None:
                positions = self.positions(n_mm)
            protospacer = self.mutate(guide.spacer, positions)
        if pam_seq is None:
            pam_seq = self.pam(pam)
        if g_base is None:
            g_base = self.g_base(g_match)
        length = len(self.seqs[chrom])
        if strand == "+":
            g_room, pam_room = start, length - (start + width)
        else:
            g_room, pam_room = length - (start + width), start
        left = g_base if g_room >= 1 else ""
        right = pam_seq[: max(0, pam_room)]
        oriented = left + protospacer + right
        if strand == "+":
            self.put(chrom, start - len(left), oriented)
        else:
            self.put(chrom, start - len(right), revcomp(oriented))
        n = sum(a != b for a, b in zip(protospacer, guide.spacer, strict=True))
        site = Planted(
            guide.name, chrom, start, start + width, strand, n, right, left, note
        )
        self.planted.append(site)
        return site

    def plant_random(
        self, guide: Guide, strand: str | None = None, **kwargs
    ) -> Planted:
        chrom, start = self.slot()
        if strand is None:
            strand = "+" if self.rng.random() < 0.5 else "-"
        return self.plant(guide, chrom, start, strand, **kwargs)

    def finish(self) -> Genome:
        for chrom, start, end in self.lowercase:
            self.seqs[chrom][start:end] = self.seqs[chrom][start:end].lower()
        contigs = {name: seq.decode() for name, seq in self.seqs.items()}
        return Genome(
            contigs, self.descriptions, self.guides, self.planted, self.pass1_planted
        )


def build_genome(seed: int = SEED) -> Genome:
    """Build the synthetic genome, its guides and the planted sites."""
    b = _Builder(seed)
    # FASTA order matters: GEM concatenates the contigs in this order.
    b.add_contig("edgeFirst", 200)
    b.add_contig("chr1", 150_000)
    b.add_contig("chr2", 80_000, "AC:CM000664.2 LN:80000 rl:Chromosome AS:synthetic")
    b.add_contig("chr3", 40_000)
    b.add_contig("chr4", 20_000)
    b.add_contig("tiny", 10)
    b.add_contig("HLA-A*01:01:01:01", 400)
    b.add_contig("bndX1", 300)
    b.add_contig("bndX2", 300)
    b.add_contig("bndY1", 300)
    b.add_contig("bndY2", 300)
    b.add_contig("edgeA", 200)
    b.add_contig("edgeB", 200)
    b.add_contig("edgeC", 200)
    b.add_contig("edgeD", 200)
    b.add_contig("edgeLast", 200)

    # Special regions: a long N run (GEM strips runs >= 50 by default), a short one,
    # a soft-masked stretch, and room for the low-complexity probes.
    for run in (20_000, 22_000, 24_000):
        b.put("chr3", run, "N" * 200)
    b.put("chr3", 30_000, "N" * 10)
    b.reserve("chr3", 19_900, 24_300)
    b.reserve("chr3", 29_900, 30_100)
    b.reserve("chr4", 4_900, 9_100)
    b.lowercase.append(("chr4", 5_000, 9_000))
    b.reserve("chr2", 39_900, 40_200)
    b.reserve("chr2", 49_900, 50_200)
    b.reserve("chr2", 60_000, 60_100)
    b.make_slots(["chr1", "chr2", "chr3", "chr4"])

    # Unique perfect NGG sites.
    b.plant_random(b.guide("uniq_plus", "unique"), "+", g_match=True)
    b.plant_random(b.guide("uniq_minus", "unique"), "-", g_match=False)
    b.plant_random(b.guide("uniq_plus_startG", "unique", first="G"), "+")
    b.plant_random(
        b.guide("uniq_minus_startG", "unique", first="G"), "-", g_match=False
    )
    b.plant(b.guide("chr2_described", "unique"), "chr2", 60_020, "+")
    b.plant(b.guide("hla_name", "unique"), "HLA-A*01:01:01:01", 100, "-")

    # Perfect copies: 2, 8 (> GEM's default -M 5), 60.
    g = b.guide("multi2", "multi")
    b.plant_random(g, "+")
    b.plant_random(g, "-", g_match=False)
    g = b.guide("multi8", "multi")
    for _ in range(8):
        b.plant_random(g, g_match=bool(b.rng.random() < 0.5))
    g = b.guide("multi60", "multi_high")
    for _ in range(60):
        b.plant_random(g, g_match=bool(b.rng.random() < 0.5))

    # Perfect plus 1/2/3-mismatch NGG off-targets, both strands.
    g = b.guide("off0123_plus", "offtarget")
    b.plant_random(g, "+", g_match=True)
    b.plant_random(g, "+", n_mm=1)
    b.plant_random(g, "-", n_mm=2, g_match=False)
    b.plant_random(g, "+", n_mm=3, g_match=False)
    g = b.guide("off0123_minus", "offtarget")
    b.plant_random(g, "-", g_match=False)
    b.plant_random(g, "-", n_mm=1, g_match=False)
    b.plant_random(g, "+", n_mm=2)
    b.plant_random(g, "-", n_mm=3, g_match=False)

    # Strata probes for -E against -s.
    g = b.guide("strata_13", "strata")
    b.plant_random(g, "+", n_mm=1, g_match=True)
    b.plant_random(g, "-", n_mm=3, g_match=False)
    b.plant_random(b.guide("strata_3_plus", "strata"), "+", n_mm=3, g_match=False)
    b.plant_random(b.guide("strata_3_minus", "strata"), "-", n_mm=3, g_match=False)
    g = b.guide("strata_03", "strata")
    b.plant_random(g, "+", g_match=True)
    b.plant_random(g, "-", n_mm=3, g_match=False)

    # Three mismatches at adversarial positions, beside a perfect site; the
    # 3-mismatch site has a non-G at the G position: the worst-case read.
    placements = {
        "5p": (0, 3, 6),
        "seed": (13, 16, 19),
        "adjacent": (9, 10, 11),
        "ends": (0, 10, 19),
        "spread": (2, 10, 17),
    }
    for label, pos in placements.items():
        for strand in "+-":
            name = f"mm3_{label}_{'plus' if strand == '+' else 'minus'}"
            g = b.guide(name, "mm3_worst")
            b.plant_random(g, g_match=True)
            b.plant_random(
                g, strand, positions=pos, g_match=False, note="worst-case read"
            )

    # Alternative PAMs.
    b.plant_random(b.guide("nag_only", "alt_pam"), "+", pam="NAG")
    b.plant_random(b.guide("nga_only", "alt_pam"), "-", pam="NGA", g_match=False)
    b.plant_random(b.guide("nag_2mm", "alt_pam"), "-", n_mm=2, pam="NAG")
    g = b.guide("alt_mixed", "alt_pam")
    b.plant_random(g, "+")
    b.plant_random(g, "-", n_mm=3, pam="NAG", g_match=False)
    b.plant_random(g, "+", n_mm=3, pam="NGA", g_match=False)

    # Perfect protospacers with no PAM.
    b.plant_random(b.guide("nopam_plus", "no_pam"), "+", pam=None)
    b.plant_random(b.guide("nopam_minus", "no_pam"), "-", pam=None)

    # A guide given as the reverse complement of a genomic NGG protospacer: in guide
    # orientation it sits on the - strand with CCN 5' of it and no PAM 3' of it.
    g = b.guide("revcomp_guide", "reversed")
    chrom, start = b.slot()
    flank3 = b.pam(None)
    b.put(chrom, start - 3, revcomp(flank3) + revcomp(g.spacer) + b.pam("NGG"))
    b.pass1_planted.append(
        (g.name, chrom, start, start + SPACER_LENGTH, "-", "reversed")
    )

    # Contig edges: no room / room for the leading G, the PAM ending at the last
    # base, and a PAM cut short by the contig end (not a site).
    n = 200
    b.plant(b.guide("edge_first_plus0", "edge"), "edgeFirst", 0, "+", note="no G room")
    b.plant(
        b.guide("edge_first_minus_end", "edge"),
        "edgeFirst",
        n - 20,
        "-",
        note="no G room",
    )
    b.plant(
        b.guide("edge_plus1", "edge"), "edgeA", 1, "+", g_match=False, note="G room"
    )
    b.plant(
        b.guide("edge_minus_room", "edge"),
        "edgeA",
        n - 21,
        "-",
        g_match=False,
        note="G room",
    )
    b.plant(
        b.guide("edge_plus_pamlast", "edge"),
        "edgeB",
        n - 23,
        "+",
        note="PAM at the end",
    )
    b.plant(
        b.guide("edge_minus_pamfirst", "edge"), "edgeB", 3, "-", note="PAM at the start"
    )
    b.plant(
        b.guide("edge_plus_pamcut", "edge_pamcut"), "edgeC", n - 22, "+", note="PAM cut"
    )
    b.plant(
        b.guide("edge_minus_pamcut", "edge_pamcut"), "edgeC", 1, "-", note="PAM cut"
    )
    b.plant(b.guide("edge_mid_plus0", "edge"), "edgeD", 0, "+", note="no G room")
    b.plant(
        b.guide("edge_mid_minus_end", "edge"), "edgeD", n - 20, "-", note="no G room"
    )
    b.plant(b.guide("edge_last_plus0", "edge"), "edgeLast", 0, "+", note="no G room")
    b.plant(
        b.guide("edge_last_minus_end", "edge"),
        "edgeLast",
        n - 20,
        "-",
        note="no G room",
    )

    # Around the 200-nt N runs at chr3:[20000, 20200), [22000, 22200), [24000, 24200).
    b.plant(b.guide("nrun_plus_after_noroom", "n_run"), "chr3", 20_200, "+", g_base="N")
    b.plant(
        b.guide("nrun_minus_before_noroom", "n_run"), "chr3", 19_980, "-", g_base="N"
    )
    b.plant(
        b.guide("nrun_plus_after_room", "n_run"), "chr3", 22_201, "+", g_match=False
    )
    b.plant(b.guide("nrun_plus_before", "n_run"), "chr3", 21_977, "+")
    b.plant(b.guide("nrun_minus_after", "n_run"), "chr3", 24_203, "-", g_match=False)
    b.plant(
        b.guide("nrun_minus_before_room", "n_run"), "chr3", 23_979, "-", g_match=False
    )

    # Genomic N inside the spacer (a mismatch) and at the PAM's N (still a PAM).
    g = b.guide("genomicN_1", "genomic_n")
    b.plant_random(g, "+", protospacer=g.spacer[:10] + "N" + g.spacer[11:])
    g = b.guide("genomicN_3", "genomic_n")
    proto = b.mutate(g.spacer, (12, 17))
    b.plant_random(g, "-", protospacer=proto[:5] + "N" + proto[6:], g_match=False)
    b.plant_random(b.guide("genomicN_pam", "genomic_n"), "+", pam_seq="NGG")

    # Soft-masked (lowercase) chr4:[5000, 9000).
    b.plant(b.guide("soft_perfect", "softmask"), "chr4", 6_000, "+")
    b.plant(b.guide("soft_2mm_minus", "softmask"), "chr4", 7_000, "-", n_mm=2)
    b.plant(b.guide("soft_boundary", "softmask"), "chr4", 8_990, "+", g_match=False)

    # Repeat families: 1-3 mismatches, mixed PAMs, both strands.
    def family(guide: Guide, copies: int) -> None:
        for _ in range(copies):
            r = b.rng.random()
            pam = "NGG" if r < 0.7 else "NAG" if r < 0.8 else "NGA" if r < 0.9 else None
            b.plant_random(
                guide,
                n_mm=int(b.rng.integers(1, 4)),
                pam=pam,
                g_match=bool(b.rng.random() < 0.5),
            )

    g = b.guide("family_offonly", "family")
    family(g, 120)
    for _ in range(20):  # 4 mismatches: never a site, but inside GEM's budget with a G
        b.plant_random(g, n_mm=4, g_match=True, note="4-mismatch decoy")
    g = b.guide("family_withperfect", "family")
    b.plant_random(g, "+")
    family(g, 110)

    # Low complexity: overlapping sites 2 nt and 1 nt apart.
    g = b.guide("tandem_GA", "low_complexity", spacer="GA" * 10)
    b.put("chr2", 40_020, "GA" * 30 + "TGG")
    g = b.guide("polyG", "low_complexity", spacer="G" * 20)
    b.put("chr2", 50_020, "A" + "G" * 40 + "TCA")

    # A read that spans two contigs (the index separates them by one separator):
    # within budget with no indel (X) or with one deletion (Y). Not sites.
    g = b.guide("bnd_hamming", "boundary")
    b.put("bndX1", 300 - 11, g.spacer[:11])
    b.put("bndX2", 0, g.spacer[12:] + "AGG")
    g = b.guide("bnd_gap", "boundary")
    b.put("bndY1", 300 - 11, g.spacer[:11])
    b.put("bndY2", 0, g.spacer[11:] + "AGG")

    for i in range(10):
        b.guide(f"nt_{i:02d}", "non_targeting")

    # Near-sites: the guide minus its 5'-most base, at a contig start (+), at a contig
    # end (-) and right after a stripped N run (+). The read overhangs by one base, so
    # these are not sites, but GEM reports them with out-of-contig or wrong positions.
    # Kept off the first contig, where an overhang can crash gem-mapper.
    g = b.guide("near_plus_start", "near_edge")
    b.put("chr2", 0, g.spacer[1:] + "TGG")
    g = b.guide("near_minus_end", "near_edge")
    b.put("chr3", 40_000 - 22, revcomp(g.spacer[1:] + "TGG"))
    b.put("chr3", 23_000, "N" * 100)
    g = b.guide("near_plus_nrun", "near_edge")
    b.put("chr3", 23_100, g.spacer[1:] + "TGG")

    return b.finish()


def stripped_runs(
    seq: str, threshold: int = STRIP_N_THRESHOLD
) -> list[tuple[int, int]]:
    """Runs of uppercase N at least ``threshold`` long, as [start, end) intervals.

    gem-indexer strips these; it keeps runs of lowercase n.
    """
    runs, start = [], None
    for i, base in enumerate(seq + "$"):
        if base == "N":
            start = i if start is None else start
        elif start is not None:
            if i - start >= threshold:
                runs.append((start, i))
            start = None
    return runs


def g_position_indexed(site: Site, contigs: dict[str, str]) -> bool:
    """True when the leading-G position exists and is not in a stripped N run.

    A read carrying the leading G aligns only where this holds.
    """
    if not site.g_base:
        return False
    pos = site.start - 1 if site.strand == "+" else site.end
    return not any(s <= pos < e for s, e in stripped_runs(contigs[site.chrom]))


def write_fasta(genome: Genome, path: Path, width: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in genome.contigs.items():
            desc = genome.descriptions.get(name, "")
            fh.write(f">{name}{' ' + desc if desc else ''}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")

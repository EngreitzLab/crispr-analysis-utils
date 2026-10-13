"""Per-guide classes, and the guide alignment's output tables.

A guide's class counts the matches whose genomic PAM matches: ``unique`` with
exactly one, ``multi`` with more, ``no_site`` with none. A match whose PAM does
not is still written, with ``pam_matches`` FALSE, and counted nowhere. An entry
that duplicates an earlier one's sequence gets that entry's sites.

Coordinates are 0-based and half-open on the reference's forward strand, and
the protospacer excludes the PAM and any base 5' of it. The cut site is a
boundary ``b``: the cut lies between ``b - 1`` and ``b``, 3 bases 5' of the PAM
(``end - 3`` on ``+``, ``start + 3`` on ``-``). In the tables, a value that
does not apply is an empty field, and booleans are ``TRUE`` or ``FALSE``.
"""

from __future__ import annotations

import os
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, TextIO

from .iupac import validate_pams
from .library import Guide
from .orientation import NOT_CHECKED, VERDICTS
from .sites import NO_PAM_CLASS, Site

if TYPE_CHECKING:
    pass

UNIQUE = "unique"
"""Exactly one perfect site with a matching PAM."""
MULTI = "multi"
"""More than one perfect site with a matching PAM."""
NO_SITE = "no_site"
"""No site with a matching PAM within `max_mismatches`."""
CLASSES = (UNIQUE, MULTI, NO_SITE)
"""Every site class, in report order."""

CUT_OFFSET = 3
"""Bases between the cut and the PAM."""


_QC_SECTIONS = ("guides", "class", "orientation")


@dataclass(frozen=True)
class GuideSummary:
    """One guide's class, site counts and, when unique, its site.

    Attributes
    ----------
    guide_id
        The sequence that was aligned: the input sequence and the PAM.
    alias
        The name the guide table gives this entry.
    spacer
        The spacer, 5' to 3'.
    mapped_sequence
        The guide's bases that pair with the genome at the unique site, as in
        `Site.mapped_sequence`; None unless the class is ``unique``.
    site_pam
        The genomic PAM of the unique site, in guide orientation; None unless
        the class is ``unique``.
    has_perfect_site
        Whether the guide has at least one perfect site with a matching PAM.
    guide_chr
        The unique site's contig; None unless unique.
    guide_start
        The unique site's protospacer start; None unless unique.
    guide_end
        The unique site's protospacer end; None unless unique.
    strand
        The unique site's strand; None unless unique.
    cut_site
        The unique site's cut boundary (see the module description); None
        unless unique.
    site_class
        ``unique``, ``multi`` or ``no_site``.
    orientation
        The orientation verdict, or ``not_checked``.
    orientation_evidence
        For a ``reverse_complement`` guide, whether that reading was a
        distinct sequence or only the spacer's reverse complement; None
        otherwise.
    extra_5p
        Whether the guide's site reads carry an added 5' G.
    min_mismatches
        The fewest spacer mismatches among the sites with a matching PAM;
        None without any.
    n_matches
        Every exact match of the spacer, whether or not a PAM sits beside it.
        More than one means the spacer occurs more than once in the genome,
        which `site_class` does not say: that counts cuttable sites only.
    n_sites
        Sites with a matching PAM by spacer mismatches: ``n_sites[m]`` sites
        have ``m``, for ``m`` from 0 to `max_mismatches`.
    input_sequence
        The sequence as given in the guide table, uppercased.
    """

    guide_id: str
    alias: str
    spacer: str
    mapped_sequence: str | None
    site_pam: str | None
    has_perfect_site: bool
    guide_chr: str | None
    guide_start: int | None
    guide_end: int | None
    strand: str | None
    cut_site: int | None
    site_class: str
    orientation: str
    orientation_evidence: str | None
    extra_5p: bool
    min_mismatches: int | None
    n_matches: int
    n_sites: tuple[int, ...]
    input_sequence: str


def summarize(
    guides: Sequence[Guide],
    sites: Mapping[int, Sequence[Site]],
    orientation: Mapping[int, str] | None = None,
    *,
    evidence: Mapping[int, str | None] | None = None,
    pam: str = "NGG",
    max_mismatches: int = 3,
) -> list[GuideSummary]:
    """Classify every guide from its sites.

    Parameters
    ----------
    guides
        The library, as returned by `read_guides`.
    sites
        The kept sites by guide row, as in `SiteResult.sites`.
    orientation
        The orientation verdicts by guide row, as in
        `ScanResult.verdicts`; None when the check did not run.
    pam
        The PAM pattern the sites were classified with.
    max_mismatches
        The mismatch limit the sites were kept with.

    Returns
    -------
    list of GuideSummary
        One per guide, in library order, duplicates included.

    Raises
    ------
    ValueError
        If the PAM patterns fail `validate_pams`, `max_mismatches` is
        negative, a site has more mismatches than `max_mismatches` or a PAM
        class that is neither one of the patterns nor ``none``, a guide is a
        duplicate of an alias not in the library, or `orientation` lacks a
        guide's verdict.
    """
    patterns = validate_pams([pam])
    _check_max_mismatches(max_mismatches)
    primary = patterns[0]
    summaries = []
    for guide, source in zip(guides, _source_rows(guides), strict=True):
        n_matches = 0
        n_sites = [0] * (max_mismatches + 1)
        perfect = []
        for site in sites.get(source, ()):
            n_matches += 1
            if site.n_mismatches > max_mismatches:
                raise ValueError(
                    f"A site of guide {site.guide_id!r} has {site.n_mismatches} "
                    f"mismatches, more than max_mismatches={max_mismatches}."
                )
            if site.pam_class == primary:
                n_sites[site.n_mismatches] += 1
                if site.n_mismatches == 0:
                    perfect.append(site)
            elif site.pam_class != NO_PAM_CLASS:
                raise ValueError(
                    f"A site of guide {site.guide_id!r} has PAM class "
                    f"{site.pam_class!r}, which is not one of the PAM patterns "
                    f"{', '.join(patterns)}."
                )

        if orientation is None:
            verdict = NOT_CHECKED
        elif source in orientation:
            verdict = orientation[source]
        else:
            raise ValueError(f"No orientation verdict for guide {guide.alias!r}.")

        if len(perfect) == 1:
            site_class = UNIQUE
        elif perfect:
            site_class = MULTI
        else:
            site_class = NO_SITE
        unique = perfect[0] if site_class == UNIQUE else None

        summaries.append(
            GuideSummary(
                guide_id=guide.guide_id,
                alias=guide.alias,
                spacer=guide.spacer,
                mapped_sequence=unique.mapped_sequence if unique else None,
                site_pam=unique.pam if unique else None,
                has_perfect_site=bool(perfect),
                guide_chr=unique.chr if unique else None,
                guide_start=unique.start if unique else None,
                guide_end=unique.end if unique else None,
                strand=unique.strand if unique else None,
                cut_site=_cut_site(unique) if unique else None,
                site_class=site_class,
                orientation=verdict,
                orientation_evidence=(
                    None if evidence is None else evidence.get(source)
                ),
                extra_5p=guide.extra_5p,
                min_mismatches=next(
                    (count for count, n in enumerate(n_sites) if n), None
                ),
                n_matches=n_matches,
                n_sites=tuple(n_sites),
                input_sequence=guide.input_sequence,
            )
        )
    return summaries


def fan_out_sites(
    guides: Sequence[Guide], sites: Mapping[int, Sequence[Site]]
) -> list[Site]:
    """Every guide's sites, in library order, duplicates included.

    Parameters
    ----------
    guides
        The library, as returned by `read_guides`.
    sites
        The kept sites by guide row, as in `SiteResult.sites`.

    Returns
    -------
    list of Site
        For each distinct `Guide.guide_id` in turn, the sites of its row, or
        of the entry it duplicates, with `Site.guide_id` set to it. Entries
        that share a sequence are one guide, so their sites are written once.
    """
    fanned = []
    seen = set()
    for guide, source in zip(guides, _source_rows(guides), strict=True):
        if guide.guide_id in seen:
            continue
        seen.add(guide.guide_id)
        fanned.extend(
            site
            if site.guide_id == guide.guide_id
            else replace(site, guide_id=guide.guide_id)
            for site in sites.get(source, ())
        )
    return fanned


GROUPS = (
    ("from the guide table, before anything is aligned", 3),
    ("about the whole guide, repeated on each of its rows", 5),
    ("about this row's match; empty when the guide matches nowhere", 7),
)
"""The column groups, in order, and how many columns each has."""

COLUMNS = {
    # from the guide table
    "guide_id": (
        "The sequence that was aligned, ending in the PAM pattern."
        " Identifies the guide."
    ),
    "alias": "Name for the guide from the input table.",
    "spacer": (
        "The protospacer: guide_id without the PAM pattern or any base 5' of it."
    ),
    # about the whole guide
    "site_class": "unique, multi or no_site: how many places this guide can cut.",
    "n_matches": (
        "How many places the protospacer matches, which is how many rows"
        " this guide has below."
    ),
    "n_cut_sites": "How many of those have a PAM, so Cas9 could cut there.",
    "orientation": (
        "Whether the sequence reads 5' to 3': forward, reverse_complement,"
        " no_pam, no_match, or not_checked."
    ),
    "orientation_evidence": (
        "Only for reverse_complement. distinct: the reversed reading was a"
        " different sequence. ambiguous: it was just the reverse complement,"
        " which a chance CCN upstream looks identical to."
    ),
    # about this row's match
    "chr": "Contig of this match.",
    "start": "First base of the protospacer, 0-based.",
    "end": "One past its last base, so the protospacer is [start, end).",
    "strand": "Strand the protospacer reads 5' to 3' on.",
    "cut_site": (
        "Where Cas9 cuts: between cut_site - 1 and cut_site, 3 bases from the PAM."
    ),
    "pam_matches": (
        "TRUE when the genomic PAM fits the pattern, so this match is cuttable."
    ),
    "mapped_sequence": (
        "What the genome holds here, to compare against guide_id: the"
        " protospacer and the real PAM, keeping the 5' base only where the"
        " genome has that same base."
    ),
}
"""Every column of guides_alignment_metadata.tsv, and what it holds."""

BED_COLUMNS = {
    "chr": "Contig of the match.",
    "start": "First base of the protospacer, 0-based.",
    "end": "One past its last base.",
    "name": "The protospacer, as spacer in the metadata table.",
    "score": "Always 0; BED6 wants the field.",
    "strand": "Strand the protospacer reads 5' to 3' on.",
}
"""The six fields of guides_mapping.bed, which has no header of its own."""


def guide_columns() -> list[str]:
    """The columns of guides_alignment_metadata.tsv, in order.

    One row per library entry per site, so the guide's own columns repeat on
    every site it has, and an entry with no site gets one row whose site
    columns are empty.

    Returns
    -------
    list of str
        The column names: the guide, then what the whole guide is, then the
        one site the row is about.
    """
    return list(COLUMNS)


def write_guides_tsv(
    path: str | os.PathLike[str],
    guides: Sequence[Guide],
    summaries: Sequence[GuideSummary],
    sites: Mapping[int, Sequence[Site]],
) -> None:
    """Write guides_alignment_metadata.tsv: every entry, a row per site.

    Parameters
    ----------
    path
        The file to write.
    guides
        The library, as returned by `read_guides`.
    summaries
        The guides, as returned by `summarize`, in the same order.
    sites
        The sites by guide row, as in `ScanResult.sites`.

    Raises
    ------
    ValueError
        If `summaries` is not one per guide.

    Notes
    -----
    `write_guides_bed` writes the same sites as a BED.

    ``n_matches`` is not ``site_class``: it counts every place the protospacer
    occurs, while ``multi`` counts only the matches with a PAM beside them, so
    a guide can match many places and still be cuttable at one.
    """
    if len(summaries) != len(guides):
        raise ValueError(f"Got {len(summaries)} summaries for {len(guides)} guides.")
    with _open_table(path, guide_columns()) as handle:
        for summary, source in zip(summaries, _source_rows(guides), strict=True):
            head = [
                summary.guide_id,
                summary.alias,
                summary.spacer,
                summary.site_class,
                summary.n_matches,
                summary.n_sites[0],
                summary.orientation,
                summary.orientation_evidence,
            ]
            found = sites.get(source, ())
            if not found:
                _write_row(handle, [*head, *[None] * 7])
            for site in found:
                _write_row(
                    handle,
                    [
                        *head,
                        site.chr,
                        site.start,
                        site.end,
                        site.strand,
                        _cut_site(site),
                        site.pam_class != NO_PAM_CLASS,
                        site.mapped_sequence,
                    ],
                )


def write_guides_bed(
    path: str | os.PathLike[str],
    guides: Sequence[Guide],
    sites: Mapping[int, Sequence[Site]],
) -> None:
    """Write guides_mapping.bed: every site a guide could be cut at.

    Parameters
    ----------
    path
        The file to write.
    guides
        The library, as returned by `read_guides`.
    sites
        The sites by guide row, as in `ScanResult.sites`.

    Notes
    -----
    BED6 without a header, fields as in `BED_COLUMNS`, one row per site whose
    genomic PAM matches, so a spacer that cuts in two places gets a row for
    each. The name field is the protospacer, which is the same on every row
    of a guide that matches more than once. Sorted by contig, then position.
    A match with no PAM
    is not somewhere anything can be cut, so it is left out;
    guides_alignment_metadata.tsv has them all.
    """
    # One row per distinct guide, not per entry: two entries of one sequence
    # would otherwise give byte-identical rows and bedtools would count twice.
    seen: set[str] = set()
    rows = []
    for guide, source in zip(guides, _source_rows(guides), strict=True):
        if guide.guide_id in seen:
            continue
        seen.add(guide.guide_id)
        rows.extend(
            (site.chr, site.start, site.end, guide.spacer, site.strand)
            for site in sites.get(source, ())
            if site.pam_class != NO_PAM_CLASS
        )
    rows.sort()
    with _open_table(path, None) as handle:
        for chrom, start, end, name, strand in rows:
            _write_row(handle, [chrom, start, end, name, 0, strand])


def write_columns_yaml(path: str | os.PathLike[str]) -> None:
    """Write columns.yaml: what every column of the outputs holds.

    Parameters
    ----------
    path
        The file to write.

    Notes
    -----
    Generated from `COLUMNS` and `BED_COLUMNS`, so it cannot drift from the
    tables themselves.
    """
    lines = ["# What each column of this folder's outputs holds.", ""]
    lines.append("guides_alignment_metadata.tsv:")
    columns = iter(COLUMNS.items())
    for heading, count in GROUPS:
        lines.append(f"  # {heading}")
        for _ in range(count):
            column, description = next(columns)
            lines.append(f'  {column}: "{_one_line(description)}"')
    lines.append("")
    lines.append("guides_mapping.bed:")
    for column, description in BED_COLUMNS.items():
        lines.append(f'  {column}: "{_one_line(description)}"')
    lines.append("")
    _write_text(path, "\n".join(lines))


def write_qc_tsv(
    path: str | os.PathLike[str],
    summaries: Sequence[GuideSummary],
    *,
    counts: Mapping[str, Mapping[str, int]] | None = None,
) -> None:
    """Write qc.tsv: the counts worth checking on a run.

    Parameters
    ----------
    path
        The file to write.
    summaries
        The guides, as returned by `summarize`.
    counts
        Further counts by section name, such as
        ``{"site_search": SiteResult.counts()}``.

    Raises
    ------
    ValueError
        If a summary has an unknown class or verdict, or a `counts` section
        reuses a built-in section name.

    Notes
    -----
    Three columns, ``section``, ``name`` and ``value``. Section ``guides``
    (``n_entries``, ``n_guides``, ``n_shared_spacers``); section ``class``,
    one row per site class; section ``orientation``, one row per verdict,
    ``not_checked`` included; then each `counts` section. The class rows and
    the orientation rows each sum to ``n_entries``. The run's parameters are
    not here: they belong to the step's metadata file.
    """
    classes = Counter(s.site_class for s in summaries)
    verdicts = Counter(s.orientation for s in summaries)
    unknown = sorted(set(classes) - set(CLASSES)) + sorted(
        set(verdicts) - set(VERDICTS)
    )
    if unknown:
        raise ValueError(f"Unknown classes or verdicts: {', '.join(unknown)}.")
    clashes = sorted(set(counts or {}) & set(_QC_SECTIONS))
    if clashes:
        raise ValueError(f"counts reuses built-in section names: {', '.join(clashes)}.")

    n_spacers = len({s.spacer for s in summaries})
    rows = [
        ("guides", "n_entries", len(summaries)),
        ("guides", "n_guides", len({s.guide_id for s in summaries})),
        ("guides", "n_shared_spacers", len(summaries) - n_spacers),
        *(("class", name, classes[name]) for name in CLASSES),
        *(("orientation", name, verdicts[name]) for name in VERDICTS),
    ]
    for section, section_counts in (counts or {}).items():
        rows.extend((section, name, value) for name, value in section_counts.items())
    with _open_table(path, ("section", "name", "value")) as handle:
        for row in rows:
            _write_row(handle, row)


def _source_rows(guides: Sequence[Guide]) -> list[int]:
    """For each guide, the row whose alignments it uses: its own, or the original's."""
    row_of = {guide.alias: row for row, guide in enumerate(guides)}
    rows = []
    for row, guide in enumerate(guides):
        if guide.duplicate_of is None:
            rows.append(row)
        elif guide.duplicate_of in row_of:
            rows.append(row_of[guide.duplicate_of])
        else:
            raise ValueError(
                f"Guide {guide.alias!r} is a duplicate of "
                f"{guide.duplicate_of!r}, which is not in the library."
            )
    return rows


def _cut_site(site: Site) -> int:
    """The cut boundary of a site, `CUT_OFFSET` bases 5' of its PAM."""
    return site.end - CUT_OFFSET if site.strand == "+" else site.start + CUT_OFFSET


def _open_table(path: str | os.PathLike[str], columns: Sequence[str] | None) -> TextIO:
    """Open a tab-separated file for writing, with a header unless BED."""
    handle = open(path, "w", encoding="utf-8", newline="\n")
    if columns is not None:
        handle.write("\t".join(columns) + "\n")
    return handle


def _write_row(handle: TextIO, values: Sequence[object]) -> None:
    """Write one tab-separated row (see `_cell`)."""
    handle.write("\t".join(_cell(value) for value in values) + "\n")


def _cell(value: object) -> str:
    """A table value: empty for None, TRUE or FALSE, collections comma-separated."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, (set, frozenset)):
        return ",".join(sorted(str(item) for item in value))
    if isinstance(value, (list, tuple)):
        return ",".join(str(item) for item in value)
    return str(value)


def _check_max_mismatches(value: object) -> None:
    """A mismatch limit must be a non-negative integer."""
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("max_mismatches must be a non-negative integer.")


def _write_text(path: str | os.PathLike[str], text: str) -> None:
    """Write text with the newline the other writers use."""
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _one_line(text: str) -> str:
    """A description as one line, safe inside a double-quoted YAML scalar."""
    return " ".join(text.split()).replace('"', "'")

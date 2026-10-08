"""Per-guide classes, and the guide alignment's output tables.

A guide's class counts only the sites whose genomic PAM matches: ``unique``
with exactly one perfect site, ``multi`` with more than one, ``mismatch_only``
with none but at least one site with 1 to `max_mismatches` mismatches, and
``no_site`` otherwise. A site whose genomic PAM matches no pattern is written
to sites.tsv with ``pam_class`` ``none`` and counted nowhere. Alternative PAMs
are off by default; given any, their sites are counted in their own columns
but never decide the class. A guide that duplicates an earlier entry's spacer
gets that entry's sites.

Coordinates are 0-based and half-open on the reference's forward strand, and
the protospacer excludes the PAM and the added G. The cut site is a boundary
``b``: the cut lies between bases ``b - 1`` and ``b``, 3 bases 5' of the PAM
(``end - 3`` on ``+``, ``start + 3`` on ``-``). In the tables, a value that
does not apply is an empty field, and booleans are ``TRUE`` or ``FALSE``.
"""

from __future__ import annotations

import os
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import TextIO

from ._alignment import check_max_mismatches
from .iupac import validate_pams
from .library import Guide
from .orientation import NOT_CHECKED, VERDICTS, OrientationResult
from .sites import NO_PAM_CLASS, Site

UNIQUE = "unique"
"""Exactly one perfect site with a matching PAM."""
MULTI = "multi"
"""More than one perfect site with a matching PAM."""
MISMATCH_ONLY = "mismatch_only"
"""No perfect site with a matching PAM, but at least one with mismatches."""
NO_SITE = "no_site"
"""No site with a matching PAM within `max_mismatches`."""
CLASSES = (UNIQUE, MULTI, MISMATCH_ONLY, NO_SITE)
"""Every site class, in report order."""

CUT_OFFSET = 3
"""Bases between the cut and the PAM."""

SITE_COLUMNS = (
    "guide_id",
    "chr",
    "start",
    "end",
    "strand",
    "pam",
    "pam_class",
    "n_mismatches",
    "mismatch_positions",
    "mapped_sequence",
    "genomic_protospacer",
    "added_g_base",
)
"""The columns of sites.tsv, in order: the fields of `Site`."""

ORIENTATION_COLUMNS = (
    "guide_id",
    "alias",
    "chr",
    "start",
    "end",
    "strand",
    "flank_5p",
    "flank_3p",
    "verdict",
)
"""The columns of orientation.tsv, in order."""

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
        ``unique``, ``multi``, ``mismatch_only`` or ``no_site``.
    orientation
        The orientation verdict, or ``not_checked``.
    g_added
        Whether the guide's site reads carry an added 5' G.
    min_mismatches
        The fewest spacer mismatches among the sites with a matching PAM;
        None without any.
    n_alignments
        The alignments the aligner reported for this guide in the site passes,
        before any were verified: a repeat shows a large number here and few
        sites. None when the counts were not given.
    n_sites
        Sites with a matching PAM by spacer mismatches: ``n_sites[m]`` sites
        have ``m``, for ``m`` from 0 to `max_mismatches`.
    n_alt_pam
        Alternative-PAM sites, counted the same way; all zero without any.
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
    g_added: bool
    min_mismatches: int | None
    n_alignments: int | None
    n_sites: tuple[int, ...]
    n_alt_pam: tuple[int, ...]
    input_sequence: str


def summarize(
    guides: Sequence[Guide],
    sites: Mapping[int, Sequence[Site]],
    orientation: Mapping[int, str] | None = None,
    *,
    pam: str = "NGG",
    alt_pams: Iterable[str] = (),
    max_mismatches: int = 3,
    n_alignments: Mapping[int, int] | None = None,
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
        `OrientationResult.verdicts`; None when the check did not run.
    pam
        The PAM pattern the sites were classified with.
    alt_pams
        The alternative PAM patterns the sites were classified with.
    max_mismatches
        The mismatch limit the sites were kept with.
    n_alignments
        The aligner's alignment count by guide row, for the `n_alignments`
        column. None leaves the column empty.

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
    patterns = validate_pams([pam, *alt_pams])
    check_max_mismatches(max_mismatches)
    primary, alternatives = patterns[0], frozenset(patterns[1:])
    summaries = []
    for guide, source in zip(guides, _source_rows(guides), strict=True):
        n_sites = [0] * (max_mismatches + 1)
        n_alt_pam = [0] * (max_mismatches + 1)
        perfect = []
        for site in sites.get(source, ()):
            if site.n_mismatches > max_mismatches:
                raise ValueError(
                    f"A site of guide {site.guide_id!r} has {site.n_mismatches} "
                    f"mismatches, more than max_mismatches={max_mismatches}."
                )
            if site.pam_class == primary:
                n_sites[site.n_mismatches] += 1
                if site.n_mismatches == 0:
                    perfect.append(site)
            elif site.pam_class in alternatives:
                n_alt_pam[site.n_mismatches] += 1
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
        elif any(n_sites):
            site_class = MISMATCH_ONLY
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
                g_added=guide.g_added,
                n_alignments=(
                    None if n_alignments is None else n_alignments.get(source, 0)
                ),
                min_mismatches=next(
                    (count for count, n in enumerate(n_sites) if n), None
                ),
                n_sites=tuple(n_sites),
                n_alt_pam=tuple(n_alt_pam),
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


def guide_columns(max_mismatches: int, alt_pams: Sequence[str] = ()) -> list[str]:
    """The columns of guides.tsv, in order.

    Parameters
    ----------
    max_mismatches
        The mismatch limit, which sets the count columns: ``n_sites_0mm`` to
        ``n_sites_<max_mismatches>mm``.
    alt_pams
        The alternative PAMs the run used. Given any, the same count columns
        follow for them, named ``n_alt_pam_<m>mm``; without them, which is the
        default, there are no alternative-PAM columns.

    Returns
    -------
    list of str
        The column names. ``spacer`` to ``strand`` follow IGVF guide-metadata
        names, except that ``guide_id`` here is the aligned sequence and the
        table's own name for the guide is ``alias``; the rest are this
        tool's.
    """
    check_max_mismatches(max_mismatches)
    counts = range(max_mismatches + 1)
    return [
        "guide_id",
        "alias",
        "spacer",
        "mapped_sequence",
        "site_pam",
        "has_perfect_site",
        "guide_chr",
        "guide_start",
        "guide_end",
        "strand",
        "cut_site",
        "site_class",
        "orientation",
        "g_added",
        "n_alignments",
        "min_mismatches",
        *(f"n_sites_{m}mm" for m in counts),
        *(f"n_alt_pam_{m}mm" for m in counts if alt_pams),
        "input_sequence",
    ]


def write_guides_tsv(
    path: str | os.PathLike[str],
    summaries: Sequence[GuideSummary],
    *,
    max_mismatches: int,
    alt_pams: Sequence[str] = (),
) -> None:
    """Write guides.tsv: one row per guide, with the columns of `guide_columns`.

    Parameters
    ----------
    path
        The file to write.
    summaries
        The guides, as returned by `summarize`.
    max_mismatches
        The mismatch limit `summarize` ran with.
    alt_pams
        The alternative PAMs `summarize` ran with, which add their count
        columns.

    Raises
    ------
    ValueError
        If a summary's count tuples do not have ``max_mismatches + 1`` entries.
    """
    columns = guide_columns(max_mismatches, alt_pams)
    with _open_table(path, columns) as handle:
        for summary in summaries:
            if len(summary.n_sites) != max_mismatches + 1:
                raise ValueError(
                    f"The counts of guide {summary.alias!r} do not match "
                    f"max_mismatches={max_mismatches}."
                )
            _write_row(
                handle,
                [
                    summary.guide_id,
                    summary.alias,
                    summary.spacer,
                    summary.mapped_sequence,
                    summary.site_pam,
                    summary.has_perfect_site,
                    summary.guide_chr,
                    summary.guide_start,
                    summary.guide_end,
                    summary.strand,
                    summary.cut_site,
                    summary.site_class,
                    summary.orientation,
                    summary.g_added,
                    summary.n_alignments,
                    summary.min_mismatches,
                    *summary.n_sites,
                    *(summary.n_alt_pam if alt_pams else ()),
                    summary.input_sequence,
                ],
            )


def write_sites_tsv(
    path: str | os.PathLike[str],
    guides: Sequence[Guide],
    sites: Mapping[int, Sequence[Site]],
) -> None:
    """Write sites.tsv: every kept site of every guide, duplicates included.

    Parameters
    ----------
    path
        The file to write.
    guides
        The library, as returned by `read_guides`.
    sites
        The kept sites by guide row, as in `SiteResult.sites`.

    Notes
    -----
    The columns are `SITE_COLUMNS`, the fields of `Site`. Rows follow
    `fan_out_sites`; ``mismatch_positions`` is comma-separated.
    """
    with _open_table(path, SITE_COLUMNS) as handle:
        for site in fan_out_sites(guides, sites):
            _write_row(handle, [getattr(site, column) for column in SITE_COLUMNS])


def write_orientation_tsv(
    path: str | os.PathLike[str],
    guides: Sequence[Guide],
    orientation: OrientationResult,
) -> None:
    """Write orientation.tsv: each guide's perfect hits, flanks and verdict.

    Parameters
    ----------
    path
        The file to write.
    guides
        The library, as returned by `read_guides`.
    orientation
        The result of `check_orientation`.

    Raises
    ------
    ValueError
        If `orientation` lacks a guide's verdict.

    Notes
    -----
    The columns are `ORIENTATION_COLUMNS`. Each guide, duplicates included,
    has one row per perfect hit, or a single row with empty hit fields when it
    has none; ``verdict`` repeats on every row.
    """
    with _open_table(path, ORIENTATION_COLUMNS) as handle:
        for guide, source in zip(guides, _source_rows(guides), strict=True):
            if source not in orientation.verdicts:
                raise ValueError(
                    f"No orientation verdict for guide {guide.guide_id!r}."
                )
            verdict = orientation.verdicts[source]
            hits = orientation.hits.get(source, ())
            if not hits:
                _write_row(handle, [guide.guide_id, guide.alias, *[None] * 6, verdict])
            for hit in hits:
                _write_row(
                    handle,
                    [
                        guide.guide_id,
                        guide.alias,
                        hit.chr,
                        hit.start,
                        hit.end,
                        hit.strand,
                        hit.flank_5p,
                        hit.flank_3p,
                        verdict,
                    ],
                )


def write_guides_bed(
    path: str | os.PathLike[str], summaries: Sequence[GuideSummary]
) -> None:
    """Write guides.bed: the protospacer of every unique guide.

    Parameters
    ----------
    path
        The file to write.
    summaries
        The guides, as returned by `summarize`.

    Notes
    -----
    BED6 without a header: contig, start, end, guide id, score 0, strand,
    sorted by contig name, then start, end and guide id.
    """
    _write_bed(
        path,
        (
            (s.guide_chr, s.guide_start, s.guide_end, s.alias, s.strand)
            for s in summaries
            if s.site_class == UNIQUE
        ),
    )


def write_cut_sites_bed(
    path: str | os.PathLike[str], summaries: Sequence[GuideSummary]
) -> None:
    """Write cut_sites.bed: the cut site of every unique guide.

    Parameters
    ----------
    path
        The file to write.
    summaries
        The guides, as returned by `summarize`.

    Notes
    -----
    BED6 without a header, one ``[b, b + 1)`` row per guide for its cut
    boundary ``b``: contig, ``b``, ``b + 1``, guide id, score 0, strand,
    sorted by contig name, then position and guide id.
    """
    _write_bed(
        path,
        (
            (s.guide_chr, s.cut_site, s.cut_site + 1, s.alias, s.strand)
            for s in summaries
            if s.site_class == UNIQUE
        ),
    )


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


def _open_table(path: str | os.PathLike[str], columns: Sequence[str]) -> TextIO:
    """Open a tab-separated file for writing and write its header line."""
    handle = open(path, "w", encoding="utf-8", newline="\n")
    handle.write("\t".join(columns) + "\n")
    return handle


def _write_bed(
    path: str | os.PathLike[str], rows: Iterable[tuple[str, int, int, str, str]]
) -> None:
    """Write sorted BED6 rows with a score of 0."""
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for chrom, start, end, name, strand in sorted(rows):
            handle.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")


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

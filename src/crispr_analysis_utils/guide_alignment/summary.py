"""Per-guide classes, and the guide alignment's output tables.

A guide's class is decided on its primary-PAM sites only (the ``ngg`` columns
count those, whatever the primary PAM is): ``unique`` with exactly one perfect
site, ``multi`` with more than one, ``off_target_only`` with none but at least
one site with 1 to `max_mismatches` mismatches, and ``no_site`` otherwise.
Alternative-PAM sites are counted but never decide the class, and sites whose
genomic PAM matches no pattern are not counted at all. A guide that duplicates
an earlier guide's spacer gets that guide's sites.

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
"""Exactly one perfect primary-PAM site."""
MULTI = "multi"
"""More than one perfect primary-PAM site."""
OFF_TARGET_ONLY = "off_target_only"
"""No perfect primary-PAM site, but at least one with mismatches."""
NO_SITE = "no_site"
"""No primary-PAM site within `max_mismatches`."""
CLASSES = (UNIQUE, MULTI, OFF_TARGET_ONLY, NO_SITE)
"""Every alignment class, in report order."""

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
    "genomic_protospacer",
    "leading_g_base",
)
"""The columns of sites.tsv, in order: the fields of `Site`."""

ORIENTATION_COLUMNS = (
    "guide_id",
    "chr",
    "start",
    "end",
    "strand",
    "flank_5p",
    "flank_3p",
    "verdict",
)
"""The columns of orientation.tsv, in order."""

_SUMMARY_SECTIONS = ("guides", "class", "orientation", "parameter")


@dataclass(frozen=True)
class GuideSummary:
    """One guide's class, site counts and, when unique, its site.

    Attributes
    ----------
    guide_id
        The guide's id.
    spacer
        The spacer, 5' to 3'.
    pam
        The genomic PAM of the unique site, in guide orientation; None unless
        the class is ``unique``.
    targeting
        Whether the guide has at least one perfect primary-PAM site.
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
    alignment_class
        ``unique``, ``multi``, ``off_target_only`` or ``no_site``.
    orientation
        The orientation verdict, or ``not_checked``.
    leading_g
        Whether the guide's site reads carry an added 5' G.
    min_ngg_mismatches
        The fewest spacer mismatches among primary-PAM sites; None without
        any.
    n_ngg
        Primary-PAM sites by spacer mismatches: ``n_ngg[m]`` sites have ``m``,
        for ``m`` from 0 to `max_mismatches`.
    n_alt_pam
        Alternative-PAM sites, counted the same way.
    duplicate_of
        The id of the first guide with the same spacer, or None.
    input_sequence
        The sequence as given in the guide table, uppercased.
    """

    guide_id: str
    spacer: str
    pam: str | None
    targeting: bool
    guide_chr: str | None
    guide_start: int | None
    guide_end: int | None
    strand: str | None
    cut_site: int | None
    alignment_class: str
    orientation: str
    leading_g: bool
    min_ngg_mismatches: int | None
    n_ngg: tuple[int, ...]
    n_alt_pam: tuple[int, ...]
    duplicate_of: str | None
    input_sequence: str


def summarize(
    guides: Sequence[Guide],
    sites: Mapping[int, Sequence[Site]],
    orientation: Mapping[int, str] | None = None,
    *,
    pam: str = "NGG",
    alt_pams: Iterable[str] = ("NAG", "NGA"),
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
        `OrientationResult.verdicts`; None when the check did not run.
    pam
        The primary PAM pattern the sites were classified with.
    alt_pams
        The alternative PAM patterns the sites were classified with.
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
        duplicate of an id not in the library, or `orientation` lacks a
        guide's verdict.
    """
    patterns = validate_pams([pam, *alt_pams])
    check_max_mismatches(max_mismatches)
    primary, alternatives = patterns[0], frozenset(patterns[1:])
    summaries = []
    for guide, source in zip(guides, _source_rows(guides), strict=True):
        n_ngg = [0] * (max_mismatches + 1)
        n_alt_pam = [0] * (max_mismatches + 1)
        perfect = []
        for site in sites.get(source, ()):
            if site.n_mismatches > max_mismatches:
                raise ValueError(
                    f"A site of guide {site.guide_id!r} has {site.n_mismatches} "
                    f"mismatches, more than max_mismatches={max_mismatches}."
                )
            if site.pam_class == primary:
                n_ngg[site.n_mismatches] += 1
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
            raise ValueError(f"No orientation verdict for guide {guide.guide_id!r}.")

        if len(perfect) == 1:
            alignment_class = UNIQUE
        elif perfect:
            alignment_class = MULTI
        elif any(n_ngg):
            alignment_class = OFF_TARGET_ONLY
        else:
            alignment_class = NO_SITE
        unique = perfect[0] if alignment_class == UNIQUE else None

        summaries.append(
            GuideSummary(
                guide_id=guide.guide_id,
                spacer=guide.spacer,
                pam=unique.pam if unique else None,
                targeting=bool(perfect),
                guide_chr=unique.chr if unique else None,
                guide_start=unique.start if unique else None,
                guide_end=unique.end if unique else None,
                strand=unique.strand if unique else None,
                cut_site=_cut_site(unique) if unique else None,
                alignment_class=alignment_class,
                orientation=verdict,
                leading_g=guide.leading_g,
                min_ngg_mismatches=next(
                    (count for count, n in enumerate(n_ngg) if n), None
                ),
                n_ngg=tuple(n_ngg),
                n_alt_pam=tuple(n_alt_pam),
                duplicate_of=guide.duplicate_of,
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
        For each guide in turn, the sites of its row, or of the guide it
        duplicates, with `Site.guide_id` set to the guide's own id.
    """
    fanned = []
    for guide, source in zip(guides, _source_rows(guides), strict=True):
        fanned.extend(
            site
            if site.guide_id == guide.guide_id
            else replace(site, guide_id=guide.guide_id)
            for site in sites.get(source, ())
        )
    return fanned


def guide_columns(max_mismatches: int) -> list[str]:
    """The columns of guides.tsv, in order.

    Parameters
    ----------
    max_mismatches
        The mismatch limit, which sets the count columns:
        ``n_ngg_0mm`` to ``n_ngg_<max_mismatches>mm``, then the same for
        ``n_alt_pam``.

    Returns
    -------
    list of str
        The column names. The first eight (``guide_id`` to ``strand``) follow
        IGVF guide-metadata names; the rest are this tool's.
    """
    check_max_mismatches(max_mismatches)
    counts = range(max_mismatches + 1)
    return [
        "guide_id",
        "spacer",
        "pam",
        "targeting",
        "guide_chr",
        "guide_start",
        "guide_end",
        "strand",
        "cut_site",
        "alignment_class",
        "orientation",
        "leading_g",
        "min_ngg_mismatches",
        *(f"n_ngg_{m}mm" for m in counts),
        *(f"n_alt_pam_{m}mm" for m in counts),
        "duplicate_of",
        "input_sequence",
    ]


def write_guides_tsv(
    path: str | os.PathLike[str],
    summaries: Sequence[GuideSummary],
    *,
    max_mismatches: int,
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

    Raises
    ------
    ValueError
        If a summary's count tuples do not have ``max_mismatches + 1`` entries.
    """
    columns = guide_columns(max_mismatches)
    with _open_table(path, columns) as handle:
        for summary in summaries:
            n_counts = max_mismatches + 1
            if len(summary.n_ngg) != n_counts or len(summary.n_alt_pam) != n_counts:
                raise ValueError(
                    f"The counts of guide {summary.guide_id!r} do not match "
                    f"max_mismatches={max_mismatches}."
                )
            _write_row(
                handle,
                [
                    summary.guide_id,
                    summary.spacer,
                    summary.pam,
                    summary.targeting,
                    summary.guide_chr,
                    summary.guide_start,
                    summary.guide_end,
                    summary.strand,
                    summary.cut_site,
                    summary.alignment_class,
                    summary.orientation,
                    summary.leading_g,
                    summary.min_ngg_mismatches,
                    *summary.n_ngg,
                    *summary.n_alt_pam,
                    summary.duplicate_of,
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
                _write_row(handle, [guide.guide_id, *[None] * 6, verdict])
            for hit in hits:
                _write_row(
                    handle,
                    [
                        guide.guide_id,
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
            (s.guide_chr, s.guide_start, s.guide_end, s.guide_id, s.strand)
            for s in summaries
            if s.alignment_class == UNIQUE
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
            (s.guide_chr, s.cut_site, s.cut_site + 1, s.guide_id, s.strand)
            for s in summaries
            if s.alignment_class == UNIQUE
        ),
    )


def write_summary_tsv(
    path: str | os.PathLike[str],
    summaries: Sequence[GuideSummary],
    *,
    counts: Mapping[str, Mapping[str, int]] | None = None,
    parameters: Mapping[str, object] | None = None,
) -> None:
    """Write summary.tsv: library, class and orientation counts, and the run.

    Parameters
    ----------
    path
        The file to write.
    summaries
        The guides, as returned by `summarize`.
    counts
        Further counts by section name, such as
        ``{"site_search": SiteResult.counts()}``.
    parameters
        The run's parameters by name.

    Raises
    ------
    ValueError
        If a summary has an unknown class or verdict, or a `counts` section
        reuses a built-in section name.

    Notes
    -----
    Three columns, ``section``, ``name`` and ``value``, in this order:
    section ``guides`` (``n_guides``, ``n_distinct_spacers``,
    ``n_duplicates``); section ``class``, one row per alignment class; section
    ``orientation``, one row per verdict, ``not_checked`` included; then each
    `counts` section; then section ``parameter``. The class rows and the
    orientation rows each sum to ``n_guides``. A list parameter is written
    comma-separated.
    """
    classes = Counter(s.alignment_class for s in summaries)
    verdicts = Counter(s.orientation for s in summaries)
    unknown = sorted(set(classes) - set(CLASSES)) + sorted(
        set(verdicts) - set(VERDICTS)
    )
    if unknown:
        raise ValueError(f"Unknown classes or verdicts: {', '.join(unknown)}.")
    clashes = sorted(set(counts or {}) & set(_SUMMARY_SECTIONS))
    if clashes:
        raise ValueError(f"counts reuses built-in section names: {', '.join(clashes)}.")

    n_distinct = sum(s.duplicate_of is None for s in summaries)
    rows = [
        ("guides", "n_guides", len(summaries)),
        ("guides", "n_distinct_spacers", n_distinct),
        ("guides", "n_duplicates", len(summaries) - n_distinct),
        *(("class", name, classes[name]) for name in CLASSES),
        *(("orientation", name, verdicts[name]) for name in VERDICTS),
    ]
    for section, section_counts in (counts or {}).items():
        rows.extend((section, name, value) for name, value in section_counts.items())
    rows.extend(
        ("parameter", name, value) for name, value in (parameters or {}).items()
    )
    with _open_table(path, ("section", "name", "value")) as handle:
        for row in rows:
            _write_row(handle, row)


def _source_rows(guides: Sequence[Guide]) -> list[int]:
    """For each guide, the row whose alignments it uses: its own, or the original's."""
    row_of = {guide.guide_id: row for row, guide in enumerate(guides)}
    rows = []
    for row, guide in enumerate(guides):
        if guide.duplicate_of is None:
            rows.append(row)
        elif guide.duplicate_of in row_of:
            rows.append(row_of[guide.duplicate_of])
        else:
            raise ValueError(
                f"Guide {guide.guide_id!r} is a duplicate of "
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

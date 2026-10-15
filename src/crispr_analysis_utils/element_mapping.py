"""Map guides to the elements a screen targets.

A site is `targeting` when the protospacer overlaps an element and
`safe_targeting` when it overlaps none. A guide that matches nowhere is
`non_targeting` and has no site.

How safe a safe-targeting site really is depends on what the screen reads
out, so that is left as `distance_to_tss` rather than decided by a threshold
here. The column is filled for every site, since it is worth having for a
targeting one too.

The class is per site, not per guide, so a guide whose sites disagree reports
each of them.
"""

from __future__ import annotations

import logging
import os
from collections import Counter
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from ._files import file_md5, write_json
from ._version import __version__

logger = logging.getLogger(__name__)

TARGETING = "targeting"
NON_TARGETING = "non_targeting"
SAFE_TARGETING = "safe_targeting"
CLASSES = (TARGETING, SAFE_TARGETING, NON_TARGETING)
"""Every class, in the order the counts report them."""

SITE_COLUMNS = ("guide_id", "chr", "start", "end", "strand")
"""What `map_to_elements` needs of every site."""

COLUMNS = {
    "guide_id": (
        "The sequence that was aligned. With the four site columns, the key"
        " back to guide_genome_annotation.tsv, which holds the rest."
    ),
    "chr": "Contig of this site. Empty when the guide matches nowhere.",
    "start": "First base of the protospacer, 0-based.",
    "end": "One past its last base.",
    "strand": "Strand the protospacer reads 5' to 3' on.",
    "element_id": (
        "The element this site overlaps, from the BED's fourth column when it"
        " has one, else chr:start-end. Empty unless targeting."
    ),
    "element_category": "The BED's fifth column, when it has one.",
    "target_class": (
        "targeting: the protospacer overlaps this element. safe_targeting:"
        " it maps but overlaps no element. non_targeting: the guide matches"
        " nowhere in the genome."
    ),
    "distance_to_tss": (
        "Bases from the cut site to the nearest protein-coding TSS, signed in"
        " that gene's orientation. How far a safe-targeting site has to be to"
        " be safe depends on the screen, so threshold this yourself rather"
        " than trusting the class alone. From guide_gene_annotation.tsv, so"
        " empty without it."
    ),
}
"""Every column of the mapping table, and what it holds."""

TABLE = "guide_element_mapping.tsv"
"""The one results table."""


class ElementError(RuntimeError):
    """The elements cannot be compared with the guides at all."""


def read_elements(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read the elements BED.

    Parameters
    ----------
    path
        A BED of the regions the screen targets. Three columns at least; a
        fourth is the element's id, a fifth its category.

    Returns
    -------
    pandas.DataFrame
        ``Chromosome``, ``Start``, ``End``, ``element_id`` and
        ``element_category``.

    Raises
    ------
    ElementError
        The file is missing, empty, or has fewer than three columns.
    """
    path = Path(path)
    if not path.is_file():
        raise ElementError(f"The elements file {path} does not exist.")

    with open(path, encoding="utf-8") as handle:
        first = handle.readline()
    skip = 1 if first.startswith(("#", "track", "browser")) else 0
    try:
        frame = pd.read_csv(path, sep="\t", header=None, skiprows=skip, dtype={0: str})
    except pd.errors.EmptyDataError as error:
        raise ElementError(f"{path} holds no element.") from error
    if frame.shape[1] < 3:
        raise ElementError(
            f"{path} has {frame.shape[1]} columns; a BED needs at least three."
        )
    if frame.empty:
        raise ElementError(f"{path} holds no element.")

    names = ["Chromosome", "Start", "End", "element_id", "element_category"]
    frame = frame.iloc[:, : len(names)]
    frame.columns = names[: frame.shape[1]]
    if "element_id" not in frame:
        # A BED with a name keeps it, so an id written elsewhere keeps its
        # spelling. Otherwise 1-based, which is how a locus is read aloud.
        frame["element_id"] = (
            frame["Chromosome"].astype(str)
            + ":"
            + (frame["Start"] + 1).astype(str)
            + "-"
            + frame["End"].astype(str)
        )
    if "element_category" not in frame:
        frame["element_category"] = ""
    return frame[names].fillna({"element_category": ""})


def map_to_elements(
    sites: pd.DataFrame,
    elements: pd.DataFrame,
    guides: pd.Series | None = None,
    distances: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Map every site to the elements it overlaps, and classify it.

    Parameters
    ----------
    sites
        The mapped sites, with `SITE_COLUMNS`.
    elements
        The elements, as `read_elements` returns them.
    guides
        Every guide in the library, so the ones that map nowhere are
        classified too. Defaults to the guides the sites name.
    distances
        ``guide_gene_annotation.tsv``, for its ``distance_to_tss``. Without
        it that column is empty; the classes do not depend on it.

    Returns
    -------
    pandas.DataFrame
        One row per site per element it overlaps, one row for a site that
        overlaps none, and one row for every guide that matches nowhere.

    Raises
    ------
    ElementError
        No site contig appears among the elements.
    """
    import pyranges1 as pr

    missing = set(SITE_COLUMNS) - set(sites.columns)
    if missing:
        raise ElementError(f"The sites are missing {', '.join(sorted(missing))}.")

    sites = sites[list(SITE_COLUMNS)].drop_duplicates().reset_index(drop=True)
    if len(sites):
        _check_contigs(sites, elements)

    hits = _overlaps(sites, elements, pr)
    targeting = set(hits["_row"]) if len(hits) else set()

    plain = sites[~sites.index.isin(targeting)].copy()
    plain["element_id"] = ""
    plain["element_category"] = ""
    plain["target_class"] = SAFE_TARGETING

    rows = [plain]
    if len(hits):
        rows.insert(
            0,
            sites.loc[hits["_row"].to_numpy()]
            .reset_index(drop=True)
            .assign(
                element_id=hits["element_id"].to_numpy(),
                element_category=hits["element_category"].to_numpy(),
                target_class=TARGETING,
            ),
        )

    every = (
        pd.Index(guides).unique()
        if guides is not None
        else pd.Index(sites["guide_id"]).unique()
    )
    unmapped = every.difference(pd.Index(sites["guide_id"]).unique())
    if len(unmapped):
        rows.append(
            pd.DataFrame(
                {
                    "guide_id": np.asarray(unmapped),
                    "chr": "",
                    "start": pd.NA,
                    "end": pd.NA,
                    "strand": "",
                    "element_id": "",
                    "element_category": "",
                    "target_class": NON_TARGETING,
                }
            )
        )

    frame = pd.concat(rows, ignore_index=True)
    order = {name: i for i, name in enumerate(every)}
    frame["_order"] = frame["guide_id"].map(order)
    frame = frame.sort_values(["_order", "chr", "start", "element_id"])
    for column in ("start", "end"):
        frame[column] = frame[column].astype("Float64").astype("Int64")
    frame = _add_distance(frame, distances)
    return frame.reset_index(drop=True).reindex(columns=list(COLUMNS))


def counts(mapping: pd.DataFrame) -> dict:
    """Count the sites in each class, and the elements that got a guide.

    Parameters
    ----------
    mapping
        What `map_to_elements` returned.

    Returns
    -------
    dict
        The class counts over sites, the number of guides in each class, and
        how many elements have at least one guide.
    """
    key = ["guide_id", "chr", "start", "end", "strand"]
    per_site = mapping.drop_duplicates(key)
    sites = Counter(per_site["target_class"])
    # A guide that matches nowhere has no site, so it is not counted as one.
    mapped = per_site[per_site["chr"].notna() & (per_site["chr"] != "")]
    guides = {
        name: int(mapping.loc[mapping["target_class"] == name, "guide_id"].nunique())
        for name in CLASSES
    }
    return {
        "rows": len(mapping),
        "guides": int(mapping["guide_id"].nunique()),
        "sites": int(len(mapped)),
        "classes": {name: int(sites.get(name, 0)) for name in CLASSES},
        "guides_by_class": guides,
        "elements_with_a_guide": int(
            mapping.loc[mapping["element_id"] != "", "element_id"].nunique()
        ),
    }


def write_mapping_tsv(path: str | os.PathLike[str], mapping: pd.DataFrame) -> None:
    """Write the guide-to-element mapping.

    Coordinates and distances are whole numbers, so they are written without
    the trailing .0 that a merge against an absent row leaves behind.
    """
    frame = mapping.reindex(columns=list(COLUMNS))
    for column in ("start", "end", "distance_to_tss"):
        frame[column] = frame[column].astype("Float64").astype("Int64")
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def write_qc_tsv(
    path: str | os.PathLike[str],
    mapping: pd.DataFrame,
    elements: pd.DataFrame,
    sites: pd.DataFrame | None = None,
) -> None:
    """Write qc.tsv: how the library covers the elements.

    Parameters
    ----------
    path
        The file to write.
    mapping
        What `map_to_elements` returned.
    elements
        The elements, as `read_elements` returns them, so the ones no guide
        reaches are counted rather than missing.
    sites
        The mapped sites, used to say how far an element with no guide is
        from the nearest one.

    Notes
    -----
    Three columns, ``section``, ``name`` and ``value``. Guides per element
    counts distinct guides, and the density beside it is per 100 bases,
    because a wide element collects more guides without being better covered.
    """
    found = counts(mapping)
    every = elements.drop_duplicates("element_id")
    per = (
        mapping.loc[mapping["target_class"] == TARGETING]
        .groupby("element_id")["guide_id"]
        .nunique()
        .reindex(every["element_id"])
        .fillna(0)
        .astype(int)
    )
    width = (every["End"] - every["Start"]).to_numpy()
    density = 100 * per.to_numpy() / np.where(width > 0, width, 1)
    covered = per > 0
    middle = float(np.median(density[covered])) if covered.any() else 0.0

    rows: list[tuple[str, str, object]] = [("classes", "guides", found["guides"])]
    rows += [("classes", f"sites_{k}", v) for k, v in found["classes"].items()]
    rows += [
        ("elements", "total", len(every)),
        ("elements", "with_a_guide", int(covered.sum())),
        ("elements", "without_a_guide", int((~covered).sum())),
        ("elements", "below_half_median_density", int((density < middle / 2).sum())),
        (
            "guides_per_element",
            "median",
            _round(np.median(per[covered]) if covered.any() else 0),
        ),
        ("guides_per_element", "min", int(per[covered].min()) if covered.any() else 0),
        ("guides_per_element", "max", int(per.max()) if len(per) else 0),
        ("guides_per_element", "median_per_100bp", _round(middle)),
    ]
    for name, group in every.groupby("element_category", dropna=False):
        if not str(name):
            continue
        got = per.reindex(group["element_id"])
        rows += [
            ("category", f"{name}_elements", len(group)),
            ("category", f"{name}_without_a_guide", int((got == 0).sum())),
            (
                "category",
                f"{name}_median_guides",
                _round(np.median(got[got > 0]) if (got > 0).any() else 0),
            ),
        ]
    rows += _untiled_rows(every[~covered.to_numpy()], sites)

    lines = ["section\tname\tvalue"]
    lines += [f"{section}\t{name}\t{value}" for section, name, value in rows]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _untiled_rows(
    untiled: pd.DataFrame, sites: pd.DataFrame | None
) -> list[tuple[str, str, object]]:
    """How far the elements with no guide are from the nearest mapped one.

    A few hundred bases means the guides missed the boundary; tens of
    kilobases means the library never targeted the element at all.
    """
    if sites is None or untiled.empty or not len(sites):
        return []
    import pyranges1 as pr

    near = pd.DataFrame(
        pr.PyRanges(
            untiled[["Chromosome", "Start", "End"]].reset_index(drop=True)
        ).nearest_ranges(
            pr.PyRanges(
                {
                    "Chromosome": sites["chr"].to_numpy(),
                    "Start": sites["start"].to_numpy().astype("int64"),
                    "End": sites["end"].to_numpy().astype("int64"),
                }
            ),
            strand_behavior="ignore",
            dist_col="_d",
            ties="first",
        )
    )
    if near.empty:
        return []
    far = near["_d"].to_numpy()
    return [
        ("untiled_elements", "median_bases_to_nearest_guide", int(np.median(far))),
        ("untiled_elements", "within_1kb_of_a_guide", int((far < 1_000).sum())),
        ("untiled_elements", "beyond_10kb_of_any_guide", int((far >= 10_000).sum())),
    ]


def _round(value: float) -> float | int:
    """A number for the QC table, without a pointless decimal."""
    value = float(value)
    return int(value) if value == int(value) else round(value, 1)


def write_columns_yaml(path: str | os.PathLike[str]) -> None:
    """Write guide_element_mapping_columns.yaml, generated from `COLUMNS` so it cannot drift."""
    lines = ["# What each column of this folder's output holds.", "", f"{TABLE}:"]
    lines += [f'  {name}: "{_one_line(text)}"' for name, text in COLUMNS.items()]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _overlaps(sites: pd.DataFrame, elements: pd.DataFrame, pr) -> pd.DataFrame:
    """Every site and element whose protospacer and region overlap."""
    if not len(sites):
        return pd.DataFrame(columns=["_row", "element_id", "element_category"])
    ranges = pr.PyRanges(
        {
            "Chromosome": sites["chr"].to_numpy(),
            "Start": sites["start"].to_numpy().astype("int64"),
            "End": sites["end"].to_numpy().astype("int64"),
            "_row": np.arange(len(sites)),
        }
    )
    regions = pr.PyRanges(
        elements[["Chromosome", "Start", "End", "element_id", "element_category"]]
        .reset_index(drop=True)
        .astype({"Start": "int64", "End": "int64"})
    )
    hits = pd.DataFrame(
        ranges.join_overlaps(
            regions, join_type="inner", strand_behavior="ignore", suffix="_element"
        )
    )
    if hits.empty:
        return pd.DataFrame(columns=["_row", "element_id", "element_category"])
    return (
        hits[["_row", "element_id", "element_category"]]
        .drop_duplicates()
        .sort_values(["_row", "element_id"])
        .reset_index(drop=True)
    )


def _add_distance(frame: pd.DataFrame, distances: pd.DataFrame | None) -> pd.DataFrame:
    """Carry each site's distance to the nearest protein-coding TSS across."""
    if distances is None or "distance_to_tss" not in getattr(distances, "columns", ()):
        return frame.assign(
            distance_to_tss=pd.Series(pd.NA, index=frame.index, dtype="Int64")
        )
    key = list(SITE_COLUMNS)
    per_site = (
        distances.drop_duplicates(key)[[*key, "distance_to_tss"]]
        .dropna(subset=["chr"])
        .copy()
    )
    for column in ("start", "end"):
        per_site[column] = per_site[column].astype("Float64").astype("Int64")
    merged = frame.merge(per_site, on=key, how="left")
    merged["distance_to_tss"] = (
        merged["distance_to_tss"].astype("Float64").astype("Int64")
    )
    merged.index = frame.index
    return merged


def _check_contigs(sites: pd.DataFrame, elements: pd.DataFrame) -> None:
    """Fail when nothing could ever overlap, rather than call every guide safe."""
    shared = set(sites["chr"].unique()) & set(elements["Chromosome"].unique())
    if not shared:
        raise ElementError(
            "No contig of the sites appears among the elements. The sites use "
            f"{_sample(sites['chr'])} and the elements use "
            f"{_sample(elements['Chromosome'])}. Do they name contigs the same way?"
        )


def _sample(values: pd.Series, n: int = 3) -> str:
    """A few distinct values, for an error message."""
    names = list(dict.fromkeys(values.astype(str)))[:n]
    return ", ".join(names) + ("..." if len(names) == n else "")


def _one_line(text: str) -> str:
    """One line, with any quote escaped, so the YAML stays valid."""
    return " ".join(text.split()).replace('"', '\\"')


# --- the pipeline -------------------------------------------------------


RESULTS_DIR = "elements"
"""Where the step's results go, so nothing of its own sits beside the inputs."""
METADATA_DIR = "run_metadata"
STEP = "element_mapping"
LOG_FILE = "element_mapping.log"
_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"


@dataclass
class RunResult:
    """What `run` produced.

    Attributes
    ----------
    outdir
        The output folder.
    mapping
        The mapping table: one row per site per element it overlaps.
    metadata
        The record written under ``run_metadata/``.
    """

    outdir: Path
    mapping: pd.DataFrame
    metadata: dict = field(default_factory=dict)


def run(
    alignment: str | os.PathLike[str] | pd.DataFrame,
    elements: str | os.PathLike[str],
    outdir: str | os.PathLike[str],
    *,
    annotation: str | os.PathLike[str] | pd.DataFrame | None = None,
) -> RunResult:
    """Map a library to the elements a screen targets, and write the table.

    Parameters
    ----------
    alignment
        ``guide_genome_annotation.tsv`` from `guide_alignment.run`, or the
        same table as a DataFrame. Every guide in it is classified, including
        the ones that match nowhere.
    elements
        A BED of the regions the screen targets.
    outdir
        Where the table goes. Created if it does not exist; nothing is ever
        written outside it.
    annotation
        ``guide_gene_annotation.tsv`` from `guide_annotation.run`, which
        fills the ``distance_to_tss`` column.

    Returns
    -------
    RunResult
        The mapping and the record.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)
    (outdir / METADATA_DIR).mkdir(exist_ok=True)
    (outdir / RESULTS_DIR).mkdir(exist_ok=True)
    record = outdir / METADATA_DIR / f"{STEP}_{_stamp()}.json"

    state: dict = {
        "step": STEP,
        "complete": False,
        "error": None,
        "started": _now(),
        "finished": None,
        "cau_version": __version__,
        "inputs": {},
        "parameters": {},
        "counts": {},
    }
    handler = logging.FileHandler(outdir / "logs" / LOG_FILE, encoding="utf-8")
    handler.setFormatter(logging.Formatter(_LOG_FORMAT))
    package = logging.getLogger(__name__.rsplit(".", 1)[0])
    package.addHandler(handler)
    if package.level == logging.NOTSET:
        package.setLevel(logging.INFO)
    # The run's log belongs in the file, not on the terminal. pyranges calls
    # logging.basicConfig on import, which puts a handler on the root logger.
    propagate = package.propagate
    package.propagate = False
    try:
        library, sites = _read_alignment(alignment)
        state["inputs"]["alignment"] = (
            "DataFrame" if isinstance(alignment, pd.DataFrame) else str(alignment)
        )
        if not isinstance(alignment, pd.DataFrame):
            state["inputs"]["alignment_md5"] = file_md5(alignment)

        regions = read_elements(elements)
        state["inputs"]["elements"] = str(elements)
        state["inputs"]["elements_md5"] = file_md5(elements)
        logger.info(
            "%d elements, %d guides, %d mapped sites",
            len(regions),
            len(library),
            len(sites),
        )

        distances = None
        if annotation is not None:
            distances = (
                annotation
                if isinstance(annotation, pd.DataFrame)
                else pd.read_csv(
                    annotation, sep="\t", dtype={"chr": str}, low_memory=False
                )
            )
            state["inputs"]["annotation"] = (
                "DataFrame" if isinstance(annotation, pd.DataFrame) else str(annotation)
            )
        else:
            logger.warning("No annotation given, so distance_to_tss will be empty.")

        mapping = map_to_elements(
            sites,
            regions,
            guides=library,
            distances=distances,
        )
        write_mapping_tsv(outdir / RESULTS_DIR / TABLE, mapping)
        write_columns_yaml(outdir / RESULTS_DIR / "guide_element_mapping_columns.yaml")
        write_qc_tsv(outdir / RESULTS_DIR / "qc.tsv", mapping, regions, sites)

        state["counts"] = counts(mapping)
        state["complete"] = True
        state["finished"] = _now()
        write_json(record, state)
        logger.info("Element mapping complete: %s", state["counts"]["classes"])
        return RunResult(outdir=outdir, mapping=mapping, metadata=state)
    except Exception as error:
        state["error"] = f"{type(error).__name__}: {error}"
        state["finished"] = _now()
        write_json(record, state)
        logger.error("Element mapping failed: %s", state["error"])
        raise
    finally:
        package.removeHandler(handler)
        package.propagate = propagate
        handler.close()


def _read_alignment(
    alignment: str | os.PathLike[str] | pd.DataFrame,
) -> tuple[pd.Series, pd.DataFrame]:
    """Every guide in the library, and the sites that have a PAM."""
    if isinstance(alignment, pd.DataFrame):
        frame = alignment.copy()
    else:
        path = Path(alignment)
        if not path.is_file():
            raise FileNotFoundError(f"The alignment table {path} does not exist.")
        frame = pd.read_csv(path, sep="\t", dtype={"chr": str}, low_memory=False)
    if "guide_id" not in frame:
        raise ElementError("The alignment table has no guide_id column.")

    library = frame["guide_id"].drop_duplicates()
    sites = frame
    if "pam_matches" in sites:
        sites = sites[sites["pam_matches"].map(_truthy)]
    sites = sites[sites["chr"].notna() & (sites["chr"].astype(str) != "")]
    sites = sites.copy()
    for column in ("start", "end"):
        sites[column] = sites[column].astype("int64")
    return library, sites[list(SITE_COLUMNS)].reset_index(drop=True)


def _truthy(value: object) -> bool:
    """TRUE, true and True all mean the same; an empty cell does not."""
    if isinstance(value, str):
        return value.strip().upper() == "TRUE"
    if value is None or pd.isna(value):
        return False
    return bool(value)


def _stamp() -> str:
    """A UTC timestamp for a file name."""
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def _now() -> str:
    """The time now, in UTC, as ISO 8601."""
    return datetime.now(UTC).isoformat(timespec="seconds")

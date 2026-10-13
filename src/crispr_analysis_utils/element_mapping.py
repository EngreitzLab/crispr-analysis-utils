"""Map guides to the elements a screen targets.

A site is `targeting` when the protospacer overlaps an element, and
`safe_targeting` when it overlaps none and lies further than the safe
distance from any protein-coding TSS. A site that overlaps no element but
sits nearer than that is `unassigned`: it maps, but to none of the elements
given. A guide that matches nowhere is `non_targeting` and has no site.

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
UNASSIGNED = "unassigned"
CLASSES = (TARGETING, SAFE_TARGETING, NON_TARGETING, UNASSIGNED)
"""Every class, in the order the counts report them."""

SAFE_DISTANCE = 1_000_000
"""How far from any protein-coding TSS a site must be to count as safe."""

SITE_COLUMNS = ("guide_id", "chr", "start", "end", "strand")
"""What `map_to_elements` needs of every site."""

COLUMNS = {
    "guide_id": (
        "The sequence that was aligned. With the four site columns, the key"
        " back to guides_alignment_metadata.tsv, which holds the rest."
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
        " it overlaps none and is further than the safe distance from any"
        " protein-coding TSS. unassigned: it overlaps none but is nearer than"
        " that. non_targeting: the guide matches nowhere in the genome."
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
    *,
    safe_distance: int = SAFE_DISTANCE,
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
        ``guides_gene_annotation.tsv``, for its ``distance_to_tss``. Without
        it nothing can be called safe, and a site that overlaps no element is
        `unassigned`.
    safe_distance
        How far from a protein-coding TSS a site must be to count as safe.

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
    if safe_distance < 0:
        raise ValueError(f"safe_distance must not be negative, got {safe_distance}.")

    sites = sites[list(SITE_COLUMNS)].drop_duplicates().reset_index(drop=True)
    if len(sites):
        _check_contigs(sites, elements)

    hits = _overlaps(sites, elements, pr)
    targeting = set(hits["_row"]) if len(hits) else set()

    far = _far_from_a_tss(sites, distances, safe_distance)
    plain = sites[~sites.index.isin(targeting)].copy()
    plain["element_id"] = ""
    plain["element_category"] = ""
    plain["target_class"] = np.where(plain.index.isin(far), SAFE_TARGETING, UNASSIGNED)

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
    guides = {
        name: int(mapping.loc[mapping["target_class"] == name, "guide_id"].nunique())
        for name in CLASSES
    }
    return {
        "rows": len(mapping),
        "guides": int(mapping["guide_id"].nunique()),
        "sites": int(len(per_site)),
        "classes": {name: int(sites.get(name, 0)) for name in CLASSES},
        "guides_by_class": guides,
        "elements_with_a_guide": int(
            mapping.loc[mapping["element_id"] != "", "element_id"].nunique()
        ),
    }


def write_mapping_tsv(path: str | os.PathLike[str], mapping: pd.DataFrame) -> None:
    """Write the guide-to-element mapping."""
    mapping.reindex(columns=list(COLUMNS)).to_csv(
        path, sep="\t", index=False, na_rep=""
    )


def write_columns_yaml(path: str | os.PathLike[str]) -> None:
    """Write columns.yaml, generated from `COLUMNS` so it cannot drift."""
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


def _far_from_a_tss(
    sites: pd.DataFrame, distances: pd.DataFrame | None, safe_distance: int
) -> set[int]:
    """The site rows further than `safe_distance` from any protein-coding TSS."""
    if distances is None or "distance_to_tss" not in getattr(distances, "columns", ()):
        return set()
    key = list(SITE_COLUMNS)
    near = (
        distances.assign(_d=distances["distance_to_tss"].abs())
        .groupby(key, dropna=False)["_d"]
        .min()
    )
    joined = sites.merge(
        near.rename("_d").reset_index(), on=key, how="left", validate="one_to_one"
    )
    return set(joined.index[(joined["_d"] >= safe_distance).fillna(False)])


def _check_contigs(sites: pd.DataFrame, elements: pd.DataFrame) -> None:
    """Fail when nothing could ever overlap, rather than call everything unassigned."""
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
    safe_distance: int = SAFE_DISTANCE,
) -> RunResult:
    """Map a library to the elements a screen targets, and write the table.

    Parameters
    ----------
    alignment
        ``guides_alignment_metadata.tsv`` from `guide_alignment.run`, or the
        same table as a DataFrame. Every guide in it is classified, including
        the ones that match nowhere.
    elements
        A BED of the regions the screen targets.
    outdir
        Where the table goes. Created if it does not exist; nothing is ever
        written outside it.
    annotation
        ``guides_gene_annotation.tsv`` from `guide_annotation.run`. Without
        it no site can be called safe.
    safe_distance
        How far from a protein-coding TSS a site must be to count as safe.

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
        "parameters": {"safe_distance": safe_distance},
        "counts": {},
    }
    handler = logging.FileHandler(outdir / "logs" / LOG_FILE, encoding="utf-8")
    handler.setFormatter(logging.Formatter(_LOG_FORMAT))
    package = logging.getLogger(__name__.rsplit(".", 1)[0])
    package.addHandler(handler)
    if package.level == logging.NOTSET:
        package.setLevel(logging.INFO)
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
            logger.warning(
                "No annotation given, so no site can be called %s.", SAFE_TARGETING
            )

        mapping = map_to_elements(
            sites,
            regions,
            guides=library,
            distances=distances,
            safe_distance=safe_distance,
        )
        write_mapping_tsv(outdir / RESULTS_DIR / TABLE, mapping)
        write_columns_yaml(outdir / RESULTS_DIR / "columns.yaml")

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
    if value is None or (isinstance(value, float) and value != value):
        return False
    if isinstance(value, str):
        return value.strip().upper() == "TRUE"
    return bool(value)


def _stamp() -> str:
    """A UTC timestamp for a file name."""
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def _now() -> str:
    """The time now, in UTC, as ISO 8601."""
    return datetime.now(UTC).isoformat(timespec="seconds")

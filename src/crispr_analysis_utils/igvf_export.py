"""Write the guide design in the IGVF per-guide metadata format.

The three steps before this one between them know everything the format
asks: alignment gives each guide its site and its PAM, annotation gives the
genes, and element mapping gives the element and the guide's role.

What a guide targets is decided from the genome, not from the design's own
label for the element. The exception is a validated enhancer, since that an
element was shown to regulate a gene is experimental knowledge no annotation
can supply.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from ._files import file_md5, write_json
from ._version import __version__

logger = logging.getLogger(__name__)

COLUMNS = {
    "guide_id": "Unique name of the guide.",
    "spacer": "The targeting sequence of the gRNA, read 5' to 3', in DNA bases.",
    "oligo": (
        "The sequence the library ordered, without the PAM: the spacer plus"
        " any base 5' of it, such as a G added for U6 transcription."
    ),
    "n_oligos": (
        "How many library entries share this spacer. Above one, they are one"
        " feature in the counts and the reads cannot be told apart, so the"
        " counts are a sum."
    ),
    "targeting": (
        "TRUE when the guide targets a locus relevant to the experiment."
        " FALSE for a safe-targeting or non-targeting control."
    ),
    "type": (
        "The guide's role: targeting, safe-targeting, non-targeting or"
        " positive control."
    ),
    "guide_chr": "Chromosome of the guide's target sequence.",
    "guide_start": "Start of the guide's target sequence, 0-based.",
    "guide_end": "End of it, half-open, and without the PAM.",
    "strand": "Strand of the genomic sequence being targeted.",
    "pam": (
        "The protospacer adjacent motif as the genome holds it, 5' to 3' on"
        " the targeted strand, immediately 3' of the spacer."
    ),
    "genomic_element": (
        "What the targeted element is: promoter, enhancer, distal element,"
        " exon or intron. Decided from the genome."
    ),
    "intended_target_name": (
        "Names the target: the gene's ENSEMBL id for a promoter, the exon's"
        " ENSEMBL id for an exon, and the coordinates otherwise."
    ),
    "intended_target_chr": "Chromosome of the target element.",
    "intended_target_start": "Start of the target element, 0-based.",
    "intended_target_end": "End of it, half-open.",
    "putative_target_genes": (
        "Every gene in ENSEMBL associated with the element, comma-separated."
        " A measured pair can name several, as a bidirectional promoter does;"
        " an inferred one names the nearest TSS's gene alone."
    ),
    "reporter": "Sequence of the reporter element. Unused here.",
    "imperfect": "Sequence of the matching perfect guide. Unused here.",
}
"""The columns of the IGVF per-guide metadata table, in order."""

PROMOTER = "promoter"
ENHANCER = "enhancer"
DISTAL = "distal element"
EXON = "exon"
INTRON = "intron"
ELEMENTS = (PROMOTER, ENHANCER, DISTAL, EXON, INTRON)
"""The `genomic_element` values this writes. The format allows more."""

TARGETING = "targeting"
SAFE = "safe-targeting"
NON = "non-targeting"
POSITIVE = "positive control"
TYPES = (TARGETING, SAFE, NON, POSITIVE)
"""The `type` values this writes."""

PROMOTER_WINDOW = 2000
"""Within this of a TSS, an element is a promoter.

The format defines `distal element` as lying more than 2 kb from any
annotated TSS, so the same threshold decides the other side.
"""

VALIDATED = "e2g_validated"
"""The one element category taken from the design rather than the genome."""

POSITIVE_CATEGORIES = ("promoter", VALIDATED)
"""The design's element kinds meant as positive controls."""

TABLE = "guide_metadata.tsv"
RESULTS_DIR = "igvf"
METADATA_DIR = "run_metadata"
STEP = "igvf_export"
LOG_FILE = "igvf_export.log"
_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"


class ExportError(RuntimeError):
    """The inputs cannot produce a valid table."""


@dataclass
class RunResult:
    """What `run` produced.

    Attributes
    ----------
    outdir
        The output folder.
    table
        The IGVF table, one row per guide.
    elements
        How each element was classified.
    metadata
        The record written under ``run_metadata/``.
    """

    outdir: Path
    table: pd.DataFrame
    elements: pd.DataFrame
    metadata: dict = field(default_factory=dict)


def classify_elements(
    elements: pd.DataFrame,
    store,
    promoters: pd.DataFrame,
    *,
    promoter_window: int = PROMOTER_WINDOW,
    genes: Mapping[str, list[str]] | None = None,
) -> pd.DataFrame:
    """Decide what each element is, and what it targets.

    Parameters
    ----------
    elements
        The elements, as `element_mapping.read_elements` returns them.
    store
        The annotation store, as `guide_annotation.load_store` returns it.
    promoters
        Transcription start sites, as `guide_annotation.read_promoters`
        returns them.
    promoter_window
        Within this many bases of a TSS, an element is a promoter.
    genes
        Element id to the ENSEMBL ids known to be regulated by it, for the
        elements where that was measured rather than inferred.

    Returns
    -------
    pandas.DataFrame
        One row per element: ``genomic_element``, ``intended_target_name``,
        the target's coordinates and ``putative_target_genes``.

    Notes
    -----
    Everything but `VALIDATED` is decided from the genome. The design's own
    label for an element is not evidence about what the element is.
    """
    import pyranges1 as pr

    if elements.empty:
        raise ExportError("There are no elements to classify.")

    frame = elements.drop_duplicates("element_id").reset_index(drop=True)
    ranges = pr.PyRanges(frame[["Chromosome", "Start", "End", "element_id"]])

    in_gene = _overlapping(ranges, store.rows("gene"), pr)
    in_exon = _overlapping(ranges, store.rows("exon"), pr)
    near, nearest_gene = _nearest_tss(ranges, promoters, pr)

    frame["genomic_element"] = [
        _element_kind(row, in_gene, in_exon, near, promoter_window)
        for row in frame.itertuples()
    ]
    frame["nearest_ensg"] = frame["element_id"].map(nearest_gene)

    exon_ids = _exon_ids(ranges, store, pr)
    introns = _intron_spans(ranges, store, pr)

    known = {} if genes is None else dict(genes)
    frame["putative_target_genes"] = [
        ",".join(_genes_for(row, known)) for row in frame.itertuples()
    ]
    frame["intended_target_name"] = [
        _target_name(row, known, exon_ids, introns) for row in frame.itertuples()
    ]
    frame["intended_target_chr"] = frame["Chromosome"]
    frame["intended_target_start"] = frame["Start"]
    frame["intended_target_end"] = frame["End"]
    return frame.drop(columns=["nearest_ensg"])


def _element_kind(row, in_gene: set, in_exon: set, near: Mapping, window: int) -> str:
    """One element's class, from the genome."""
    if getattr(row, "element_category", "") == VALIDATED:
        return ENHANCER
    if near.get(row.element_id, np.inf) <= window:
        return PROMOTER
    if row.element_id in in_exon:
        return EXON
    if row.element_id in in_gene:
        return INTRON
    return DISTAL


def _genes_for(row, known: Mapping[str, list[str]]) -> list[str]:
    """Every gene associated with an element, measured ones first."""
    measured = known.get(row.element_id, [])
    if measured:
        return list(dict.fromkeys(measured))
    if row.genomic_element == PROMOTER and isinstance(row.nearest_ensg, str):
        return [row.nearest_ensg]
    return []


def _target_name(row, known: Mapping, exon_ids: Mapping, introns: Mapping) -> str:
    """What `intended_target_name` holds, which the element's class decides."""
    if row.genomic_element == PROMOTER:
        genes = _genes_for(row, known)
        return genes[0] if genes else ""
    if row.genomic_element == EXON:
        return exon_ids.get(row.element_id, "")
    if row.genomic_element == INTRON:
        span = introns.get(row.element_id)
        if span:
            return f"{row.Chromosome}:{span[0]}-{span[1]}"
    return f"{row.Chromosome}:{row.Start}-{row.End}"


def _overlapping(ranges, rows: pd.DataFrame, pr) -> set[str]:
    """The element ids overlapping any of `rows`."""
    if rows.empty:
        return set()
    hits = ranges.join_overlaps(
        pr.PyRanges(rows[["Chromosome", "Start", "End"]].reset_index(drop=True)),
        join_type="inner",
        strand_behavior="ignore",
        suffix="_x",
    )
    frame = pd.DataFrame(hits)
    return set(frame["element_id"]) if len(frame) else set()


def _nearest_tss(ranges, promoters: pd.DataFrame, pr) -> tuple[dict, dict]:
    """Each element's distance to the nearest TSS, and that TSS's gene."""
    columns = ["Chromosome", "Start", "End", "ensg"]
    near = pd.DataFrame(
        ranges.nearest_ranges(
            pr.PyRanges(promoters[columns].reset_index(drop=True)),
            strand_behavior="ignore",
            dist_col="_d",
            ties="first",
            suffix="_tss",
        )
    )
    if near.empty:
        return {}, {}
    near = near.drop_duplicates("element_id")
    gene = near["ensg_tss"] if "ensg_tss" in near else near["ensg"]
    return (
        dict(zip(near["element_id"], near["_d"], strict=True)),
        dict(zip(near["element_id"], gene, strict=True)),
    )


def _exon_ids(ranges, store, pr) -> dict[str, str]:
    """The ENSEMBL exon id for each element overlapping an exon.

    The reference transcript's exon is preferred, so a rebuild chooses the
    same one; otherwise the lowest id, which is just as stable.
    """
    rows = store.rows("exon")
    if rows.empty or "exon_id" not in rows:
        return {}
    columns = ["Chromosome", "Start", "End", "exon_id", "is_reference"]
    hits = pd.DataFrame(
        ranges.join_overlaps(
            pr.PyRanges(rows[columns].reset_index(drop=True)),
            join_type="inner",
            strand_behavior="ignore",
            suffix="_x",
        )
    )
    if hits.empty:
        return {}
    hits = hits.dropna(subset=["exon_id"]).sort_values(
        ["element_id", "is_reference", "exon_id"], ascending=[True, False, True]
    )
    first = hits.drop_duplicates("element_id")
    return dict(zip(first["element_id"], first["exon_id"], strict=True))


def _intron_spans(ranges, store, pr) -> dict[str, tuple[int, int]]:
    """The intron each element falls in, as (start, end).

    Introns are the gaps between a reference transcript's exons, which the
    store holds rather than the GTF.
    """
    exons = store.rows("exon", reference_only=True)
    if exons.empty:
        return {}
    exons = exons.sort_values(["enst", "Start"])
    starts = exons["End"].to_numpy()[:-1]
    ends = exons["Start"].to_numpy()[1:]
    same = exons["enst"].to_numpy()[:-1] == exons["enst"].to_numpy()[1:]
    keep = same & (ends > starts)
    introns = pd.DataFrame(
        {
            "Chromosome": exons["Chromosome"].to_numpy()[:-1][keep],
            "Start": starts[keep],
            "End": ends[keep],
        }
    )
    if introns.empty:
        return {}
    hits = pd.DataFrame(
        ranges.join_overlaps(
            pr.PyRanges(introns.reset_index(drop=True)),
            join_type="inner",
            strand_behavior="ignore",
            suffix="_i",
        )
    )
    if hits.empty:
        return {}
    hits = hits.sort_values(["element_id", "Start_i"]).drop_duplicates("element_id")
    return {
        row.element_id: (int(row.Start_i), int(row.End_i)) for row in hits.itertuples()
    }


# --- the table ----------------------------------------------------------


def build_table(
    alignment: pd.DataFrame,
    mapping: pd.DataFrame,
    elements: pd.DataFrame,
    *,
    positive_categories: tuple[str, ...] = POSITIVE_CATEGORIES,
) -> pd.DataFrame:
    """Assemble the IGVF table, one row per guide.

    Parameters
    ----------
    alignment
        ``guide_genome_annotation.tsv``.
    mapping
        ``guide_element_mapping.tsv``.
    elements
        What `classify_elements` returned.
    positive_categories
        The design's element categories meant as positive controls. A guide
        is only called one when the genome agrees the element is a promoter,
        or the element is a validated enhancer.

    Returns
    -------
    pandas.DataFrame
        The columns of `COLUMNS`, one row per guide, in library order.

    Raises
    ------
    ExportError
        A guide has more than one spacer, which the format forbids.
    """
    guides = alignment.drop_duplicates("guide_id").set_index("guide_id")
    if guides["alias"].duplicated().any():
        raise ExportError("Two guides share an alias, so guide_id would not be unique.")

    rows = _chosen_rows(mapping, elements)
    missing = guides.index.difference(rows.index)
    if len(missing):
        raise ExportError(
            f"{len(missing)} guides are in the alignment but not in the element"
            f" mapping, starting with {list(missing[:3])}. The two tables"
            " describe different libraries."
        )
    chosen = rows["element_id"]
    info = elements.set_index("element_id")

    out = pd.DataFrame(index=guides.index)
    out["guide_id"] = guides["alias"]
    out["spacer"] = guides["spacer"]
    # Carried so the submitted table says which guides share a feature: their
    # counts are a sum and are not comparable with a guide ordered once.
    for column in ("oligo", "n_oligos"):
        # An alignment table written before these existed has neither, and
        # the column is then written empty rather than the table refused.
        out[column] = guides[column] if column in guides.columns else pd.NA

    element = chosen.reindex(out.index)
    kind = element.map(info["genomic_element"])
    category = (
        element.map(info.set_index(info.index)["element_category"])
        if "element_category" in info
        else pd.Series(dtype=object)
    )
    klass = rows["target_class"].reindex(out.index)

    positive = category.isin(positive_categories) & (
        (category == VALIDATED) | (kind == PROMOTER)
    )
    out["targeting"] = (klass == "targeting").map({True: "TRUE", False: "FALSE"})
    out["type"] = np.select(
        [klass == "safe_targeting", klass == "non_targeting", positive],
        [SAFE, NON, POSITIVE],
        default=TARGETING,
    )

    mapped = klass != "non_targeting"
    for column, source in (
        ("guide_chr", "chr"),
        ("guide_start", "start"),
        ("guide_end", "end"),
        ("strand", "strand"),
    ):
        out[column] = rows[source].reindex(out.index).where(mapped)
    out["pam"] = (
        guides["mapped_sequence"]
        .str[-3:]
        .where(mapped & guides["mapped_sequence"].notna())
    )

    for column, source in (
        ("genomic_element", "genomic_element"),
        ("intended_target_name", "intended_target_name"),
        ("intended_target_chr", "intended_target_chr"),
        ("intended_target_start", "intended_target_start"),
        ("intended_target_end", "intended_target_end"),
        ("putative_target_genes", "putative_target_genes"),
    ):
        out[column] = element.map(info[source])
    out["reporter"] = ""
    out["imperfect"] = ""

    for column in (
        "guide_start",
        "guide_end",
        "intended_target_start",
        "intended_target_end",
    ):
        out[column] = out[column].astype("Float64").astype("Int64")
    return out.reset_index(drop=True)[list(COLUMNS)]


def _chosen_rows(mapping: pd.DataFrame, elements: pd.DataFrame) -> pd.DataFrame:
    """One row of the mapping per guide, the one its output row describes.

    The format allows one row per guide, so a guide inside two overlapping
    elements has to name one: the element sharing the most bases with the
    protospacer wins, and its id breaks a tie, so the choice is the same on
    every run. The class and the coordinates come from that same row, or the
    row would name one site's element beside another site's coordinates.
    """
    named = mapping[mapping["element_id"].notna() & (mapping["element_id"] != "")]
    best = named.head(0)
    if len(named):
        bounds = elements.drop_duplicates("element_id").set_index("element_id")
        joined = named.join(bounds[["Start", "End"]], on="element_id")
        overlap = np.minimum(joined["End"], joined["end"]) - np.maximum(
            joined["Start"], joined["start"]
        )
        best = (
            joined.assign(_overlap=overlap)
            .sort_values(
                ["guide_id", "_overlap", "element_id"], ascending=[True, False, True]
            )
            .drop_duplicates("guide_id")
            .drop(columns=["Start", "End", "_overlap"])
        )
    rest = mapping[~mapping["guide_id"].isin(set(best["guide_id"]))].drop_duplicates(
        "guide_id"
    )
    return pd.concat([best, rest], ignore_index=True).set_index("guide_id")


def write_table(path: str | os.PathLike[str], table: pd.DataFrame) -> None:
    """Write the IGVF table, with empty cells where a field does not apply."""
    table.reindex(columns=list(COLUMNS)).to_csv(path, sep="\t", index=False, na_rep="")


def write_columns_yaml(path: str | os.PathLike[str]) -> None:
    """Write guide_metadata_columns.yaml, generated from `COLUMNS` so it cannot drift."""
    lines = ["# What each column of this folder's output holds.", "", f"{TABLE}:"]
    lines += [f'  {name}: "{" ".join(text.split())}"' for name, text in COLUMNS.items()]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def counts(table: pd.DataFrame, elements: pd.DataFrame) -> dict:
    """What the table holds, by type and by element class."""
    return {
        "guides": int(len(table)),
        "type": {name: int((table["type"] == name).sum()) for name in TYPES},
        "genomic_element": {
            name: int((table["genomic_element"] == name).sum()) for name in ELEMENTS
        },
        "elements": int(len(elements)),
        "with_a_gene": int((elements["putative_target_genes"] != "").sum()),
        "missing_pam": int(((table["type"] != NON) & table["pam"].isna()).sum()),
    }


# --- the pipeline -------------------------------------------------------


def run(
    alignment: str | os.PathLike[str],
    mapping: str | os.PathLike[str],
    elements: str | os.PathLike[str],
    store: str | os.PathLike[str],
    outdir: str | os.PathLike[str],
    *,
    promoters: str | os.PathLike[str] | None = None,
    genes: str | os.PathLike[str] | None = None,
    promoter_window: int = PROMOTER_WINDOW,
) -> RunResult:
    """Write the guide design in the IGVF per-guide metadata format.

    Parameters
    ----------
    alignment
        ``guide_genome_annotation.tsv``.
    mapping
        ``guide_element_mapping.tsv``.
    elements
        The elements BED the mapping used.
    store
        The annotation store, which supplies the genes, exons and exon ids.
    outdir
        Where the table goes. Nothing is written outside it.
    promoters
        A TSS BED. Without one, the store's reference transcript starts.
    genes
        A TSV of measured element-to-gene pairs, with an element id column
        and an ENSEMBL id column. Used for the elements whose gene was
        measured rather than inferred.
    promoter_window
        Within this many bases of a TSS, an element is a promoter.

    Returns
    -------
    RunResult
        The table, the element classification and the record.
    """
    from .element_mapping import read_elements
    from .guide_annotation import load_store, promoters_from_store, read_promoters

    outdir = Path(outdir)
    for name in (RESULTS_DIR, METADATA_DIR, "logs"):
        (outdir / name).mkdir(parents=True, exist_ok=True)
    record = outdir / METADATA_DIR / f"{STEP}_{_stamp()}.json"

    state: dict = {
        "step": STEP,
        "complete": False,
        "error": None,
        "started": _now(),
        "cau_version": __version__,
        "inputs": {
            "alignment": str(alignment),
            "mapping": str(mapping),
            "elements": str(elements),
            "store": str(store),
            "genes": None if genes is None else str(genes),
        },
        "parameters": {"promoter_window": promoter_window},
        "counts": {},
    }
    handler = logging.FileHandler(outdir / "logs" / LOG_FILE, encoding="utf-8")
    handler.setFormatter(logging.Formatter(_LOG_FORMAT))
    package = logging.getLogger(__name__.rsplit(".", 1)[0])
    package.addHandler(handler)
    if package.level == logging.NOTSET:
        package.setLevel(logging.INFO)
    propagate = package.propagate
    package.propagate = False
    try:
        loaded = load_store(store)
        tss = read_promoters(promoters) if promoters else promoters_from_store(loaded)
        regions = read_elements(elements)
        align = pd.read_csv(alignment, sep="\t", dtype={"chr": str}, low_memory=False)
        maps = pd.read_csv(mapping, sep="\t", dtype={"chr": str}, low_memory=False)
        state["inputs"]["alignment_md5"] = file_md5(alignment)
        state["inputs"]["elements_md5"] = file_md5(elements)

        info = classify_elements(
            regions,
            loaded,
            tss,
            promoter_window=promoter_window,
            genes=_read_genes(genes),
        )
        table = build_table(align, maps, info)

        write_table(outdir / RESULTS_DIR / TABLE, table)
        write_columns_yaml(outdir / RESULTS_DIR / "guide_metadata_columns.yaml")
        state["counts"] = counts(table, info)
        state["complete"] = True
        state["finished"] = _now()
        write_json(record, state)
        logger.info("IGVF export complete: %s", state["counts"]["type"])
        return RunResult(outdir=outdir, table=table, elements=info, metadata=state)
    except Exception as error:
        state["error"] = f"{type(error).__name__}: {error}"
        state["finished"] = _now()
        write_json(record, state)
        logger.error("IGVF export failed: %s", state["error"])
        raise
    finally:
        package.removeHandler(handler)
        package.propagate = propagate
        handler.close()


def _read_genes(path: str | os.PathLike[str] | None) -> dict[str, list[str]] | None:
    """Measured element-to-gene pairs, as element id to ENSEMBL ids.

    The first column is taken as the element and the first column holding
    ENSEMBL gene ids as the gene, so a pairs table written for another tool
    can be used as it is.
    """
    if path is None:
        return None
    frame = pd.read_csv(path, sep="\t" if str(path).endswith((".tsv", ".txt")) else ",")
    if frame.empty:
        return {}
    ensembl = [
        column
        for column in frame.columns
        if frame[column].astype(str).str.startswith("ENSG").any()
    ]
    if not ensembl:
        raise ExportError(f"{path} has no column of ENSEMBL gene ids.")
    element, gene = frame.columns[0], ensembl[0]
    pairs: dict[str, list[str]] = {}
    for key, rows in frame.groupby(element)[gene]:
        pairs[str(key)] = list(dict.fromkeys(rows.astype(str)))
    return pairs


def _stamp() -> str:
    """A UTC timestamp for a file name."""
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def _now() -> str:
    """The time now, in UTC, as ISO 8601."""
    return datetime.now(UTC).isoformat(timespec="seconds")

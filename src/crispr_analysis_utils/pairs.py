"""The guide-gene pairs a screen tests: discovery and positive control.

Built from the IGVF per-guide table alone, so the pairs name exactly the
elements and genes the submission does.

- **positive_control**: a guide typed ``positive control`` with each gene in
  its ``putative_target_genes``. Those genes are the ones the design measured
  the element against, so a gene outside the readout panel is kept and
  flagged rather than dropped.
- **discovery**: a targeting guide with each readout gene whose TSS lies
  within ``max_distance`` of the element's midpoint. Where that element and
  gene are a positive control's, the pair is labelled positive_control
  instead, so a promoter still pairs with its other neighbours as discovery.

Control guides get no pairs; the inference calibrates on them.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from ._files import file_md5, write_json
from ._version import __version__

logger = logging.getLogger(__name__)

DISCOVERY = "discovery"
POSITIVE_CONTROL = "positive_control"
PAIR_TYPES = (DISCOVERY, POSITIVE_CONTROL)

TARGET_COLUMNS = (
    "intended_target_name",
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end",
)
"""The element's identity, as the per-guide table and the pipeline name it."""

COLUMNS = {
    "guide_id": "The guide, as guide_metadata.tsv names it.",
    "gene_id": "ENSEMBL gene id, without its version.",
    "gene_symbol": "The gene's symbol in the annotation store.",
    "intended_target_name": "The element, as guide_metadata.tsv names it.",
    "intended_target_chr": "Chromosome of the element.",
    "intended_target_start": "Start of the element, 0-based.",
    "intended_target_end": "End of the element, half-open.",
    "pair_type": (
        "positive_control: the gene the design measured this element against."
        " discovery: a readout gene whose TSS is within the distance limit."
    ),
    "distance_to_tss": (
        "Bases from the element's midpoint to the gene's TSS, signed in the"
        " gene's orientation: negative upstream, positive downstream."
    ),
    "in_panel": "TRUE when the gene is in the readout panel given to --genes.",
}
"""Every column of the pairs table, and what it holds."""

TABLE = "pairs_to_test.csv"
"""Comma-separated, because the CRISPR pipeline reads its pairs with read_csv."""

RESULTS_DIR = "pairs"
METADATA_DIR = "run_metadata"
STEP = "pairs"
LOG_FILE = "pairs.log"
_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"

DEFAULT_MAX_DISTANCE = 1_000_000


class PairsError(RuntimeError):
    """The inputs cannot be turned into pairs."""


@dataclass
class RunResult:
    """What `run` produced."""

    outdir: Path
    pairs: pd.DataFrame
    metadata: dict = field(default_factory=dict)


def read_panel(path: str | os.PathLike[str]) -> pd.Index:
    """The readout genes, as unversioned ENSEMBL ids.

    Parameters
    ----------
    path
        A CSV or TSV. The first column holding ENSEMBL gene ids is used.

    Returns
    -------
    pandas.Index
        Unique ids, in file order.

    Raises
    ------
    PairsError
        No column holds ENSEMBL gene ids.
    """
    sep = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    frame = pd.read_csv(path, sep=sep, dtype=str)
    for column in frame.columns:
        values = frame[column].dropna()
        if len(values) and values.str.startswith("ENSG").all():
            return pd.Index(_unversioned(values).unique())
    raise PairsError(f"{path} has no column of ENSEMBL gene ids.")


def gene_tss(store) -> pd.DataFrame:
    """One TSS per gene, from the transcript the annotation chooses.

    Parameters
    ----------
    store
        A loaded `guide_annotation.db.Store`.

    Returns
    -------
    pandas.DataFrame
        Indexed by ``gene_id``: ``gene_symbol``, ``chr``, ``tss`` (0-based) and
        ``strand``. The reference transcript is MANE Select, else Ensembl
        canonical, else the longest CDS; a gene without one falls back to its
        own row.
    """
    transcripts = store.rows("transcript", reference_only=True, primary=False)
    genes = store.rows("gene", primary=False)
    rows = pd.concat(
        [transcripts, genes[~genes["ensg"].isin(transcripts["ensg"])]],
        ignore_index=True,
    ).drop_duplicates("ensg")
    plus = rows["Strand"] == "+"
    return pd.DataFrame(
        {
            "gene_symbol": rows["gene_symbol"].to_numpy(),
            "chr": rows["Chromosome"].astype(str).to_numpy(),
            "tss": np.where(plus, rows["Start"], rows["End"] - 1).astype("int64"),
            "strand": rows["Strand"].astype(str).to_numpy(),
        },
        index=pd.Index(rows["ensg"].to_numpy(), name="gene_id"),
    )


def build_pairs(
    guides: pd.DataFrame,
    panel: pd.Index,
    tss: pd.DataFrame,
    *,
    max_distance: int = DEFAULT_MAX_DISTANCE,
) -> pd.DataFrame:
    """Every positive-control and discovery pair, one row per guide and gene.

    Parameters
    ----------
    guides
        The IGVF per-guide table: ``guide_id``, ``targeting``, ``type``, the
        four `TARGET_COLUMNS` and ``putative_target_genes``.
    panel
        The readout genes, unversioned.
    tss
        `gene_tss` output.
    max_distance
        The largest distance, in bases, from an element's midpoint to a
        readout gene's TSS for a discovery pair.

    Returns
    -------
    pandas.DataFrame
        The columns of `COLUMNS`, positive controls first.

    Raises
    ------
    PairsError
        A required column is missing, or a panel gene or a positive-control
        gene is not in the annotation.
    """
    missing = {
        "guide_id",
        "targeting",
        "type",
        "putative_target_genes",
        *TARGET_COLUMNS,
    } - set(guides)
    if missing:
        raise PairsError(f"The guide table lacks {', '.join(sorted(missing))}.")
    absent = panel.difference(tss.index)
    if len(absent):
        raise PairsError(
            f"{len(absent)} panel genes are not in the annotation, e.g. {', '.join(absent[:3])}."
            " The panel and the annotation must be the same GENCODE release."
        )

    frame = guides.copy()
    frame["targeting"] = frame["targeting"].map(_truthy)
    for column in ("intended_target_start", "intended_target_end"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce").astype("Int64")
    targeting = frame[frame["targeting"] & frame["intended_target_start"].notna()]

    positive = _positive_pairs(targeting, tss)
    discovery = _discovery_pairs(targeting, panel, tss, max_distance)
    key = [
        "intended_target_chr",
        "intended_target_start",
        "intended_target_end",
        "gene_id",
    ]
    taken = pd.MultiIndex.from_frame(positive[key])
    discovery.loc[pd.MultiIndex.from_frame(discovery[key]).isin(taken), "pair_type"] = (
        POSITIVE_CONTROL
    )

    pairs = pd.concat([positive, discovery], ignore_index=True)
    pairs = pairs.drop_duplicates(["guide_id", "gene_id"])
    pairs["in_panel"] = pairs["gene_id"].isin(panel)
    return pairs[list(COLUMNS)].reset_index(drop=True)


def _positive_pairs(targeting: pd.DataFrame, tss: pd.DataFrame) -> pd.DataFrame:
    """Each positive-control guide with each gene the design names for it."""
    rows = targeting[targeting["type"] == "positive control"].copy()
    rows["gene_id"] = rows["putative_target_genes"].fillna("").str.split(",")
    rows = rows.explode("gene_id")
    rows["gene_id"] = _unversioned(rows["gene_id"].str.strip())
    unnamed = rows[rows["gene_id"] == ""]
    if len(unnamed):
        logger.warning(
            "%d positive-control guides name no gene and give no pair, e.g. %s",
            unnamed["guide_id"].nunique(),
            ", ".join(unnamed["guide_id"].head(3)),
        )
    rows = rows[rows["gene_id"] != ""]
    unknown = sorted(set(rows["gene_id"]) - set(tss.index))
    if unknown:
        raise PairsError(
            f"{len(unknown)} positive-control genes are not in the annotation, e.g. {', '.join(unknown[:3])}."
        )
    return _with_gene(rows, tss, POSITIVE_CONTROL)


def _discovery_pairs(
    targeting: pd.DataFrame, panel: pd.Index, tss: pd.DataFrame, max_distance: int
) -> pd.DataFrame:
    """Each targeting guide with each readout gene whose TSS is near its element."""
    elements = targeting[list(TARGET_COLUMNS)].drop_duplicates()
    readout = tss.loc[panel].rename_axis("gene_id").reset_index()
    near = elements.merge(readout, left_on="intended_target_chr", right_on="chr")
    midpoint = (near["intended_target_start"] + near["intended_target_end"]) // 2
    near = near[(midpoint - near["tss"]).abs() <= max_distance]
    rows = targeting.merge(
        near[[*TARGET_COLUMNS, "gene_id"]], on=list(TARGET_COLUMNS), how="inner"
    )
    return _with_gene(rows, tss, DISCOVERY)


def _with_gene(rows: pd.DataFrame, tss: pd.DataFrame, pair_type: str) -> pd.DataFrame:
    """Add the gene's name and the signed distance from the element to its TSS."""
    gene = tss.loc[rows["gene_id"].to_numpy()]
    midpoint = (
        (rows["intended_target_start"] + rows["intended_target_end"]) // 2
    ).to_numpy(dtype="int64")
    offset = midpoint - gene["tss"].to_numpy()
    out = rows[["guide_id", "gene_id", *TARGET_COLUMNS]].reset_index(drop=True)
    out["gene_symbol"] = gene["gene_symbol"].to_numpy()
    out["pair_type"] = pair_type
    out["distance_to_tss"] = np.where(gene["strand"].to_numpy() == "-", -offset, offset)
    return out


def counts(pairs: pd.DataFrame, panel: pd.Index) -> dict:
    """What the run record and qc.tsv report."""
    found: dict = {"pairs": len(pairs)}
    for pair_type in PAIR_TYPES:
        part = pairs[pairs["pair_type"] == pair_type]
        found[pair_type] = {
            "pairs": len(part),
            "guides": int(part["guide_id"].nunique()),
            "elements": int(part[list(TARGET_COLUMNS)].drop_duplicates().shape[0]),
            "genes": int(part["gene_id"].nunique()),
            "genes_outside_panel": int(
                part.loc[~part["in_panel"], "gene_id"].nunique()
            ),
        }
    found["panel_genes"] = len(panel)
    found["panel_genes_without_pairs"] = int((~panel.isin(pairs["gene_id"])).sum())
    return found


def write_qc_tsv(path: str | os.PathLike[str], found: dict) -> None:
    """Counts, one per row, in the section/name/value shape the other steps use."""
    rows = [("pairs", "total", found["pairs"])]
    for pair_type in PAIR_TYPES:
        rows += [(pair_type, name, value) for name, value in found[pair_type].items()]
    rows += [
        ("panel", "genes", found["panel_genes"]),
        ("panel", "genes_without_pairs", found["panel_genes_without_pairs"]),
    ]
    pd.DataFrame(rows, columns=["section", "name", "value"]).to_csv(
        path, sep="\t", index=False
    )


def write_columns_yaml(path: str | os.PathLike[str]) -> None:
    """Write pairs_to_test_columns.yaml, generated from `COLUMNS` so it cannot drift."""
    lines = ["# What each column of this folder's output holds.", "", f"{TABLE}:"]
    lines += [f'  {name}: "{text}"' for name, text in COLUMNS.items()]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(
    guide_metadata: str | os.PathLike[str],
    panel: str | os.PathLike[str],
    store: str | os.PathLike[str],
    outdir: str | os.PathLike[str],
    *,
    max_distance: int = DEFAULT_MAX_DISTANCE,
) -> RunResult:
    """Build the pairs and write them.

    Parameters
    ----------
    guide_metadata
        ``guide_metadata.tsv`` from `igvf_export.run`.
    panel
        The readout genes: a CSV or TSV with a column of ENSEMBL gene ids.
    store
        The annotation store the guide table was built with.
    outdir
        Where the table goes. Nothing is written outside it.
    max_distance
        The largest distance from an element's midpoint to a readout gene's
        TSS for a discovery pair.

    Returns
    -------
    RunResult
        The pairs and the run record.
    """
    from .guide_annotation.db import load_store

    outdir = Path(outdir)
    for sub in ("logs", METADATA_DIR, RESULTS_DIR):
        (outdir / sub).mkdir(parents=True, exist_ok=True)
    record = outdir / METADATA_DIR / f"{STEP}_{_stamp()}.json"
    state: dict = {
        "step": STEP,
        "complete": False,
        "error": None,
        "started": _now(),
        "finished": None,
        "cau_version": __version__,
        "inputs": {},
        "parameters": {"max_distance": max_distance},
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
        for name, path in (
            ("guide_metadata", guide_metadata),
            ("panel", panel),
            ("store", store),
        ):
            if not Path(path).is_file():
                raise FileNotFoundError(
                    f"The {name.replace('_', ' ')} {path} does not exist."
                )
            state["inputs"][name] = str(path)
            state["inputs"][f"{name}_md5"] = file_md5(path)

        guides = pd.read_csv(guide_metadata, sep="\t", dtype=str, keep_default_na=False)
        genes = read_panel(panel)
        tss = gene_tss(load_store(store))
        logger.info(
            "%d guides, %d panel genes, %d annotated genes",
            len(guides),
            len(genes),
            len(tss),
        )

        pairs = build_pairs(guides, genes, tss, max_distance=max_distance)
        write_table(outdir / RESULTS_DIR / TABLE, pairs)
        write_columns_yaml(outdir / RESULTS_DIR / "pairs_to_test_columns.yaml")
        found = counts(pairs, genes)
        write_qc_tsv(outdir / RESULTS_DIR / "qc.tsv", found)
        unpaired = genes[~genes.isin(pairs["gene_id"])]
        tss.loc[unpaired, ["gene_symbol", "chr", "tss", "strand"]].reset_index().to_csv(
            outdir / RESULTS_DIR / "panel_genes_without_pairs.tsv",
            sep="\t",
            index=False,
        )

        state["counts"] = found
        state["complete"] = True
        state["finished"] = _now()
        write_json(record, state)
        logger.info(
            "Pairs complete: %d positive control, %d discovery",
            found[POSITIVE_CONTROL]["pairs"],
            found[DISCOVERY]["pairs"],
        )
        return RunResult(outdir=outdir, pairs=pairs, metadata=state)
    except Exception as error:
        state["error"] = f"{type(error).__name__}: {error}"
        state["finished"] = _now()
        write_json(record, state)
        logger.error("Pairs failed: %s", state["error"])
        raise
    finally:
        package.removeHandler(handler)
        package.propagate = propagate
        handler.close()


def write_table(path: str | os.PathLike[str], pairs: pd.DataFrame) -> None:
    """Write the pairs, with booleans as TRUE and FALSE like the other tables."""
    out = pairs.copy()
    out["in_panel"] = out["in_panel"].map({True: "TRUE", False: "FALSE"})
    out.to_csv(path, index=False)


def _unversioned(ids: pd.Series) -> pd.Series:
    """ENSG00000186092.7 -> ENSG00000186092; a _PAR_Y suffix stays."""
    return ids.astype(str).str.replace(r"\.\d+$", "", regex=True)


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

"""The guide alignment: read a guide table, scan a reference, write the tables.

One pass of the compiled core finds every exact site of every guide, on both
strands, and reads the PAM 3' of each and the base 5' of it. There is no
aligner, no index and no error budget: a guide is designed, not sequenced, so
its own site is an exact match, and everything the tables report follows from
that.

`outdir` gets three folders and nothing else: ``alignments/`` with the
results (``guide_genome_annotation.tsv``, ``guide_mapping.bed``, or
``guide_mapping.bed.gz`` and its ``.tbi`` under ``bgzip``,
``guide_genome_annotation_columns.yaml``, ``qc.tsv``), ``run_metadata/`` with
this invocation's record, and ``logs/``.
"""

from __future__ import annotations

import logging
import os
import time
from collections import Counter
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .._files import (
    check_plain_fasta,
    file_md5,
    run_record,
    timed,
    write_json,
)
from .iupac import reverse_complement, validate_pams
from .library import Guide, read_guides
from .orientation import FORWARD, REVERSE_COMPLEMENT, VERDICTS
from .scan import ScanResult, scan_guides
from .summary import (
    CLASSES,
    GuideSummary,
    summarize,
    write_columns_yaml,
    write_guides_bed,
    write_guides_tsv,
    write_qc_tsv,
)

logger = logging.getLogger(__name__)

RESULTS_DIR = "alignments"
"""Where the step's results go, so nothing of its own sits beside the inputs."""
STEP = "guide_alignment"
"""The step its run records and log are named for."""

ORIENTATION_CHECKS = ("warn", "error", "off")
"""The values of `run`'s ``orientation_check``."""


class OrientationError(RuntimeError):
    """The library reads 3' to 5' and `orientation_check` is ``error``."""


@dataclass
class RunResult:
    """What `run` produced.

    Attributes
    ----------
    outdir
        The output folder.
    guides
        The library, as returned by `read_guides`.
    summaries
        One per guide, in library order, as returned by `summarize`.
    scan
        The genome scan: each guide's sites, its orientation, and the counts.
    metadata
        The record written under ``run_metadata/``.
    """

    outdir: Path
    guides: list[Guide]
    summaries: list[GuideSummary]
    scan: ScanResult
    metadata: dict = field(default_factory=dict)


def run(
    guides: str | os.PathLike[str] | pd.DataFrame,
    reference: str | os.PathLike[str],
    outdir: str | os.PathLike[str],
    *,
    pam: str = "NGG",
    spacer_length: int = 20,
    tolerate_5p: int = 0,
    add_leading_g: bool = False,
    orientation_check: str = "error",
    fix_orientation: bool = False,
    contigs: Collection[str] | None = None,
    pam_only: bool = False,
    bgzip: bool = False,
    id_col: int | str = 0,
    spacer_col: int | str = 1,
    sep: str = "\t",
    header: bool = False,
) -> RunResult:
    """Find every guide's genomic site and write the tables.

    Parameters
    ----------
    guides
        The guide table, as `read_guides` takes it.
    reference
        The reference FASTA, uncompressed.
    outdir
        Where the tables go. Created if it does not exist; nothing is ever
        written outside it.
    pam
        The PAM, in IUPAC codes, 3' of the protospacer.
    spacer_length
        See `read_guides`.
    tolerate_5p
        See `scan_guides`. The protospacer and its coordinates are unchanged;
        only the exact search is shorter.
    add_leading_g
        See `read_guides`.
    orientation_check
        ``error``: stop when entries reading 3' to 5' outnumber those reading
        5' to 3', which is a table given reverse-complemented throughout. A
        handful is reported either way, since one guide can look reversed for
        want of nothing more than a PAM that differs between assemblies.
        ``warn``: only report. ``off``: skip the check.
    fix_orientation
        Align the reverse complement of the entries that read 3' to 5',
        instead of leaving them as given.
    contigs
        If given, read only these contigs.
    pam_only
        Write only the sites whose genomic PAM matches, which is every site
        a nuclease could act at. A short spacer matches far more places than
        it can cut, and the rest are what makes the table unusable.
    bgzip
        Write ``guide_mapping.bed.gz`` and its tabix index instead of a
        plain BED, so a browser or `tabix` can seek into it.
    id_col, spacer_col, sep, header
        See `read_guides`.

    Returns
    -------
    RunResult
        The library, the summaries, the scan and this invocation's metadata.

    Raises
    ------
    FileNotFoundError
        If `reference` does not exist.
    ValueError
        If `reference` is compressed, `pam` is not IUPAC codes,
        `orientation_check` is not one of `ORIENTATION_CHECKS`, or the guide
        table fails `read_guides`.
    OrientationError
        If the library reads 3' to 5' and `orientation_check` is ``error``.
    """
    if isinstance(contigs, str):
        raise TypeError("contigs must be a collection of contig names, not a string.")
    if orientation_check not in ORIENTATION_CHECKS:
        raise ValueError(
            f"orientation_check must be one of {', '.join(ORIENTATION_CHECKS)}, "
            f"not {orientation_check!r}."
        )
    patterns = validate_pams([pam])
    reference = check_plain_fasta(reference, "reference")

    outdir = Path(outdir)
    parameters = {
        "pam": patterns[0],
        "spacer_length": spacer_length,
        "tolerate_5p": tolerate_5p,
        "add_leading_g": add_leading_g,
        "orientation_check": orientation_check,
        "fix_orientation": fix_orientation,
        "contigs": None if contigs is None else list(contigs),
        "pam_only": pam_only,
        "bgzip": bgzip,
        "id_col": id_col,
        "spacer_col": spacer_col,
        "sep": sep,
        "header": header,
    }
    started = time.perf_counter()
    with run_record(
        outdir, STEP, RESULTS_DIR, parameters, logger, "Guide alignment"
    ) as (state, record):
        state.update(classes={}, orientation={}, orientation_fixed=0)
        logger.info("Guide alignment into %s, with %s", outdir, parameters)
        return _run(guides, reference, outdir, parameters, state, record, started)


def _run(
    guides: str | os.PathLike[str] | pd.DataFrame,
    reference: Path,
    outdir: Path,
    parameters: dict,
    state: dict,
    record: Path,
    started: float,
) -> RunResult:
    """The run itself, inside the log handler and the metadata record."""
    with timed(state, "guides"):
        source = guides if isinstance(guides, pd.DataFrame) else os.fspath(guides)
        library = read_guides(
            guides,
            id_col=parameters["id_col"],
            spacer_col=parameters["spacer_col"],
            sep=parameters["sep"],
            header=parameters["header"],
            spacer_length=parameters["spacer_length"],
            add_leading_g=parameters["add_leading_g"],
            pam=parameters["pam"],
        )
    state["inputs"] = {
        "guides": "DataFrame" if isinstance(guides, pd.DataFrame) else source,
        "reference": str(reference),
    }
    logger.info(
        "%d entries, %d distinct guides, %d with a base 5' of the spacer",
        len(library),
        len({guide.guide_id for guide in library}),
        sum(guide.extra_5p for guide in library),
    )

    with timed(state, "checksums"):
        state["inputs"]["reference_md5"] = file_md5(reference)
    write_json(record, state)

    with timed(state, "scan"):
        scan = scan_guides(
            library,
            reference,
            pam=parameters["pam"],
            contigs=parameters["contigs"],
            pam_only=parameters["pam_only"],
            tolerate_5p=parameters["tolerate_5p"],
        )
    state["counts"]["scan"] = scan.counts()
    write_json(record, state)

    library, scan = _check_orientation(library, scan, reference, parameters, state)

    with timed(state, "outputs"):
        summaries = summarize(
            library,
            scan.sites,
            scan.verdicts or None,
            evidence=scan.evidence,
            pam=parameters["pam"],
            max_mismatches=parameters["tolerate_5p"],
            tolerate_5p=parameters["tolerate_5p"],
            matches=scan.matches,
        )
        tables = outdir / RESULTS_DIR
        write_guides_tsv(
            tables / "guide_genome_annotation.tsv", library, summaries, scan.sites
        )
        write_guides_bed(
            tables / "guide_mapping.bed",
            library,
            scan.sites,
            bgzip=parameters["bgzip"],
        )
        write_columns_yaml(tables / "guide_genome_annotation_columns.yaml")
        write_qc_tsv(tables / "qc.tsv", summaries, counts=state["counts"])

    state["classes"] = _ordered((summary.site_class for summary in summaries), CLASSES)
    state["orientation"] = _ordered(
        (summary.orientation for summary in summaries), VERDICTS
    )
    logger.info(
        "Done in %.1f s: classes %s; outputs in %s",
        time.perf_counter() - started,
        state["classes"],
        outdir,
    )
    return RunResult(
        outdir=outdir,
        guides=list(library),
        summaries=summaries,
        scan=scan,
        metadata=state,
    )


def _check_orientation(
    library: list[Guide],
    scan: ScanResult,
    reference: Path,
    parameters: dict,
    state: dict,
) -> tuple[list[Guide], ScanResult]:
    """Act on the orientation verdicts, flipping entries when asked to."""
    check = parameters["orientation_check"]
    if check == "off":
        # The scan holds both readings whatever happens, but "off" has to mean
        # off: every entry's verdict is not_checked.
        logger.info("Orientation check: off")
        scan.verdicts = {}
        scan.evidence = {}
        return library, scan

    counts = Counter(scan.verdicts.values())
    backwards = [
        row for row, verdict in scan.verdicts.items() if verdict == REVERSE_COMPLEMENT
    ]
    if not backwards:
        return library, scan

    aliases = ", ".join(repr(library[row].alias) for row in backwards[:5])
    if len(backwards) > 5:
        aliases += ", ..."
    message = (
        f"{len(backwards)} of {len(scan.verdicts)} entries read 3' to 5': their "
        f"reverse complement has a site with a PAM and the sequence as given "
        f"does not ({aliases})."
    )
    # Asked to fix them, fix them: the refusal below names this flag, so it has
    # to be acted on before we refuse, or the advice would be impossible to take.
    if not parameters["fix_orientation"]:
        # One guide can look reversed because its PAM differs between
        # assemblies, so only a table reversed throughout stops the run.
        if counts[REVERSE_COMPLEMENT] > counts[FORWARD] and check == "error":
            raise OrientationError(
                f"{message} Pass --fix-orientation to align their reverse "
                "complement instead, or --orientation-check warn to carry on."
            )
        logger.warning(message)
        return library, scan

    logger.warning("%s Reverse-complementing those entries.", message)
    with timed(state, "rescan"):
        library = _reverse_complemented(library, backwards, parameters)
        scan = scan_guides(
            library,
            reference,
            pam=parameters["pam"],
            contigs=parameters["contigs"],
            pam_only=parameters["pam_only"],
            tolerate_5p=parameters["tolerate_5p"],
        )
    state["orientation_fixed"] = len(backwards)
    state["counts"]["scan"] = scan.counts()
    return library, scan


def _reverse_complemented(
    library: Sequence[Guide], rows: Sequence[int], parameters: Mapping[str, object]
) -> list[Guide]:
    """The library with only `rows` read again from their reverse complement."""
    flip = set(rows)
    table = pd.DataFrame(
        {
            "alias": [guide.alias for guide in library],
            "sequence": [
                reverse_complement(guide.input_sequence)
                if row in flip
                else guide.input_sequence
                for row, guide in enumerate(library)
            ],
        }
    )
    return read_guides(
        table,
        id_col="alias",
        spacer_col="sequence",
        spacer_length=parameters["spacer_length"],
        add_leading_g=parameters["add_leading_g"],
        pam=parameters["pam"],
    )


def _ordered(values: Iterable[str], order: Sequence[str]) -> dict[str, int]:
    counts = Counter(values)
    return {name: counts[name] for name in order}

"""The guide alignment pipeline: both GEM passes, the site checks and the outputs.

`run` reads the guide table, checks the library's orientation from the exact
hits of the bare spacers (pass 1), searches each guide's sites with one read
per PAM (pass 2), verifies every hit against the reference, classifies the
guides and writes every output into one folder. A GEM pass is reused from an
earlier run into the same folder only when ``run.json`` records the same
reads, options, index and GEM version for it, and its alignments are
unchanged.
"""

from __future__ import annotations

import json
import logging
import os
import time
from collections import Counter
from collections.abc import Collection, Iterable, Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .._version import __version__
from ._alignment import check_contigs, check_max_mismatches, open_reference
from ._files import check_plain_fasta, file_md5, write_json
from .gem import (
    GemError,
    gem_version,
    index_file,
    map_reads,
    mapper_options,
    read_sidecar,
    sidecar_path,
)
from .iupac import count_ambiguous, validate_pams
from .library import Guide, read_guides, write_site_reads, write_spacer_reads
from .orientation import (
    OK,
    REVERSED,
    VERDICTS,
    OrientationResult,
    check_orientation,
    library_reversed,
)
from .sites import SiteResult, read_sites
from .summary import (
    CLASSES,
    GuideSummary,
    summarize,
    write_cut_sites_bed,
    write_guides_bed,
    write_guides_tsv,
    write_orientation_tsv,
    write_sites_tsv,
    write_summary_tsv,
)

logger = logging.getLogger(__name__)

ORIENTATION_CHECKS = ("warn", "error", "off")
"""The values of `run`'s ``orientation_check``."""

PASS_ORIENTATION = "orientation"
"""Pass 1: the bare spacers, matched exactly."""
PASS_SITES = "sites"
"""Pass 2, the reads without an added G."""
PASS_SITES_LEADING_G = "sites_leading_g"
"""Pass 2, the reads with an added 5' G."""

RUN_JSON = "run.json"
LOG_FILE = "guide_alignment.log"
_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"


class OrientationError(RuntimeError):
    """Reversed guides outnumber ok ones, and the orientation check is ``error``."""


@dataclass(frozen=True)
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
    sites
        The site search: the kept sites by guide row, and its counts.
    orientation
        The orientation check, or None when it was off.
    passes
        The record of each GEM pass that ran or was reused, by pass name, as
        in ``run.json``.
    """

    outdir: Path
    guides: list[Guide]
    summaries: list[GuideSummary]
    sites: SiteResult
    orientation: OrientationResult | None
    passes: dict[str, dict]


def run(
    guides: str | os.PathLike[str] | pd.DataFrame,
    reference: str | os.PathLike[str],
    index: str | os.PathLike[str],
    outdir: str | os.PathLike[str],
    *,
    max_mismatches: int = 3,
    pam: str = "NGG",
    alt_pams: Iterable[str] = ("NAG", "NGA"),
    spacer_length: int = 20,
    add_leading_g: bool = False,
    orientation_check: str = "warn",
    contigs: Collection[str] | None = None,
    threads: int = 8,
    id_col: int | str = 0,
    spacer_col: int | str = 1,
    sep: str = "\t",
    header: bool = False,
) -> RunResult:
    """Align a guide library to a reference and write every output to `outdir`.

    Parameters
    ----------
    guides
        The guide table: a delimited text file or a DataFrame, read by
        `read_guides` with `id_col`, `spacer_col`, `sep`, `header`,
        `spacer_length` and `add_leading_g`.
    reference
        The reference FASTA the index was built from, uncompressed. Its
        ``.fai`` is used when there is one; otherwise one is built in
        `outdir`.
    index
        The GEM index, ``<prefix>.gem`` or ``<prefix>``. When it has the
        sidecar `build_index` writes, the reference's MD5 checksum must match
        the one recorded there.
    outdir
        The output folder, created if missing. Nothing is written elsewhere.
    max_mismatches
        The spacer mismatches a site may have.
    pam
        The primary PAM, in IUPAC codes, 3' of the protospacer.
    alt_pams
        The alternative PAMs, in IUPAC codes; empty for none. Every PAM must
        have the same length.
    spacer_length
        See `read_guides`.
    add_leading_g
        See `read_guides`.
    orientation_check
        ``warn``: run pass 1 and log a warning when ``reversed`` guides
        outnumber ``ok`` ones. ``error``: stop after pass 1 in that case.
        ``off``: skip pass 1; every guide is ``not_checked``. Guides are
        counted per distinct spacer.
    contigs
        Keep only the hits on these contigs; None keeps every contig.
    threads
        gem-mapper threads.
    id_col
        See `read_guides`.
    spacer_col
        See `read_guides`.
    sep
        See `read_guides`.
    header
        See `read_guides`.

    Returns
    -------
    RunResult
        The library, the guide summaries, the sites, the orientation check
        and the GEM passes.

    Raises
    ------
    FileNotFoundError
        If the guide table, the reference or the index does not exist.
    ValueError
        If a parameter or the guide table is invalid, the reference is
        gzip-compressed or does not match the index's sidecar, `contigs`
        names a contig the reference lacks, or an alignment does not fit the
        reference (see `read_sites`).
    OrientationError
        If `orientation_check` is ``error`` and ``reversed`` guides outnumber
        ``ok`` ones.
    GemError
        If GEM is missing, fails or crashes; the message names the guide when
        gem-mapper names the read it crashed on.

    Notes
    -----
    `outdir` receives ``guides.tsv``, ``sites.tsv``, ``orientation.tsv``
    (unless the check is off), ``guides.bed``, ``cut_sites.bed``,
    ``summary.tsv`` and ``run.json``; the reads and the alignments of each
    pass go in ``alignments/``, and the GEM logs and ``guide_alignment.log``
    (appended to) in ``logs/``. ``run.json`` records the inputs and their
    checksums, the parameters, the GEM version, each pass's command, counts
    and wall time, the alignments per guide and the step times; it is
    rewritten after each pass, and says whether the run completed.
    """
    if isinstance(alt_pams, str):
        raise TypeError("alt_pams must be a list of PAM patterns, not a string.")
    if isinstance(contigs, str):
        raise TypeError("contigs must be a collection of contig names, not a string.")
    patterns = validate_pams([pam, *alt_pams])
    check_max_mismatches(max_mismatches)
    if orientation_check not in ORIENTATION_CHECKS:
        raise ValueError(
            f"orientation_check must be one of {', '.join(ORIENTATION_CHECKS)}, "
            f"not {orientation_check!r}."
        )
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer.")
    contig_list = None if contigs is None else list(contigs)
    parameters = {
        "max_mismatches": max_mismatches,
        "pam": patterns[0],
        "alt_pams": list(patterns[1:]),
        "spacer_length": spacer_length,
        "add_leading_g": add_leading_g,
        "orientation_check": orientation_check,
        "contigs": contig_list,
        "threads": threads,
        "id_col": id_col,
        "spacer_col": spacer_col,
        "sep": sep,
        "header": header,
    }

    outdir = Path(outdir)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    (outdir / "alignments").mkdir(exist_ok=True)
    with _log_to(outdir / "logs" / LOG_FILE):
        return _Run(outdir, parameters, patterns).execute(guides, reference, index)


class _Run:
    """One run of the pipeline, and its ``run.json`` record."""

    def __init__(
        self, outdir: Path, parameters: dict, patterns: tuple[str, ...]
    ) -> None:
        self.outdir = outdir
        self.alignments = outdir / "alignments"
        self.logs = outdir / "logs"
        self.parameters = parameters
        self.patterns = patterns
        self.run_json = outdir / RUN_JSON
        self.previous = _read_previous(self.run_json)
        self.previous_passes = _mapping(self.previous.get("passes"))
        self.attempted: set[str] = set()
        self.started = time.perf_counter()
        self.state: dict = {
            "complete": False,
            "error": None,
            "started": _now(),
            "finished": None,
            "cau_version": __version__,
            "gem_version": None,
            "inputs": {},
            "parameters": parameters,
            "passes": {},
            "counts": {},
            "classes": {},
            "orientation": {},
            "seconds": {},
            "alignments_per_guide": {},
        }

    def save(self) -> None:
        """Write ``run.json``.

        Until the run completes, the earlier run's records of the passes this
        run has not reached are kept, marked ``carried_over``, so a failed run
        does not lose alignments a later run could reuse.
        """
        state = self.state
        if not state["complete"]:
            carried = {
                name: {**record, "carried_over": True}
                for name, record in self.previous_passes.items()
                if name not in self.attempted and isinstance(record, dict)
            }
            if carried:
                state = {**state, "passes": {**carried, **state["passes"]}}
        write_json(self.run_json, state)

    @contextmanager
    def timed(self, step: str) -> Iterator[None]:
        started = time.perf_counter()
        try:
            yield
        finally:
            self.state["seconds"][step] = round(time.perf_counter() - started, 3)

    def execute(
        self,
        guides: str | os.PathLike[str] | pd.DataFrame,
        reference: str | os.PathLike[str],
        index: str | os.PathLike[str],
    ) -> RunResult:
        logger.info("Guide alignment into %s, with %s", self.outdir, self.parameters)
        try:
            result = self._execute(guides, reference, index)
        except BaseException as error:
            self.state["error"] = f"{type(error).__name__}: {error}"
            self.state["finished"] = _now()
            self.state["seconds"]["total"] = round(
                time.perf_counter() - self.started, 3
            )
            self.save()
            logger.error("Guide alignment failed: %s", self.state["error"])
            raise
        return result

    def _execute(
        self,
        guides: str | os.PathLike[str] | pd.DataFrame,
        reference: str | os.PathLike[str],
        index: str | os.PathLike[str],
    ) -> RunResult:
        parameters = self.parameters
        with self.timed("guides"):
            library = read_guides(
                guides,
                id_col=parameters["id_col"],
                spacer_col=parameters["spacer_col"],
                sep=parameters["sep"],
                header=parameters["header"],
                spacer_length=parameters["spacer_length"],
                add_leading_g=parameters["add_leading_g"],
            )
        n_distinct = sum(guide.duplicate_of is None for guide in library)
        n_leading_g = sum(guide.leading_g for guide in library)
        logger.info(
            "%d guides, %d distinct spacers, %d with an added G",
            len(library),
            n_distinct,
            n_leading_g,
        )

        reference = check_plain_fasta(reference, "reference").resolve()
        index = index_file(index).resolve()
        if not index.is_file():
            raise FileNotFoundError(f"The GEM index {index} does not exist.")
        version = gem_version()
        self.state["gem_version"] = version
        with self.timed("checksums"):
            reference_md5 = file_md5(reference)
            guides_input = (
                {"path": None, "md5": None}
                if isinstance(guides, pd.DataFrame)
                else {
                    "path": str(Path(guides).resolve()),
                    "md5": file_md5(guides),
                }
            )
        sidecar = read_sidecar(index)
        _check_sidecar(sidecar, index, reference, reference_md5, version)
        fasta_index, fasta_index_input = self._fasta_index(reference, reference_md5)
        index_stat = index.stat()
        index_identity = {
            "path": str(index),
            "size": index_stat.st_size,
            "mtime_ns": index_stat.st_mtime_ns,
            "reference_md5": sidecar.get("reference_md5") if sidecar else None,
        }
        self.state["inputs"] = {
            "guides": guides_input,
            "reference": {
                "path": str(reference),
                "size": reference.stat().st_size,
                "md5": reference_md5,
            },
            "index": index_identity,
            "fasta_index": fasta_index_input,
        }
        with open_reference(reference, fasta_index) as fasta:
            check_contigs(parameters["contigs"], fasta)
        self.save()

        context = {
            "library": library,
            "index": index,
            "index_identity": index_identity,
            "version": version,
        }
        orientation = self._orientation(context, reference, fasta_index)
        sites = self._sites(context, reference, fasta_index)

        with self.timed("outputs"):
            summaries = summarize(
                library,
                sites.sites,
                None if orientation is None else orientation.verdicts,
                pam=self.patterns[0],
                alt_pams=self.patterns[1:],
                max_mismatches=parameters["max_mismatches"],
            )
            self._write_outputs(
                library,
                summaries,
                sites,
                orientation,
                {
                    "guides": guides_input["path"] or "DataFrame",
                    "reference": str(reference),
                    "index": str(index),
                    "gem_version": version,
                    "cau_version": __version__,
                },
            )

        self.state["classes"] = _ordered_counts(
            (s.alignment_class for s in summaries), CLASSES
        )
        self.state["orientation"] = _ordered_counts(
            (s.orientation for s in summaries), VERDICTS
        )
        self.state["alignments_per_guide"] = self._alignments_per_guide(library)
        self.state["complete"] = True
        self.state["finished"] = _now()
        self.state["seconds"]["total"] = round(time.perf_counter() - self.started, 3)
        self.save()
        logger.info(
            "Done in %.1f s: classes %s; outputs in %s",
            self.state["seconds"]["total"],
            self.state["classes"],
            self.outdir,
        )
        return RunResult(
            outdir=self.outdir,
            guides=library,
            summaries=summaries,
            sites=sites,
            orientation=orientation,
            passes=dict(self.state["passes"]),
        )

    def _fasta_index(
        self, reference: Path, reference_md5: str
    ) -> tuple[Path | None, dict]:
        """The FASTA index to use: the reference's own, else one in `outdir`."""
        own = Path(f"{reference}.fai")
        if own.is_file():
            return None, {"path": str(own), "built_in_outdir": False}
        path = (self.outdir / own.name).resolve()
        earlier = _mapping(_mapping(self.previous.get("inputs")).get("fasta_index"))
        built_for_this = (
            earlier.get("path") == str(path)
            and earlier.get("reference_md5") == reference_md5
        )
        if path.exists() and not built_for_this:
            path.unlink()
        return path, {
            "path": str(path),
            "built_in_outdir": True,
            "reference_md5": reference_md5,
        }

    def _orientation(
        self, context: dict, reference: Path, fasta_index: Path | None
    ) -> OrientationResult | None:
        """Pass 1 and the library-level orientation check."""
        fastq = self.alignments / f"{PASS_ORIENTATION}.fastq"
        sam = self.alignments / f"{PASS_ORIENTATION}.sam.gz"
        table = self.outdir / "orientation.tsv"
        check = self.parameters["orientation_check"]
        if check == "off":
            self.attempted.add(PASS_ORIENTATION)
            for stale in (fastq, sam, table):
                stale.unlink(missing_ok=True)
            logger.info("Orientation check: off")
            return None

        with self.timed(PASS_ORIENTATION):
            n_reads = write_spacer_reads(context["library"], fastq)
            self._gem_pass(
                PASS_ORIENTATION,
                context,
                fastq=fastq,
                n_reads=n_reads,
                sam=sam,
                max_mismatches=0,
                n_ambiguous=0,
                leading_g=False,
            )
        with self.timed("orientation_check"):
            orientation = check_orientation(
                sam,
                context["library"],
                reference,
                pam=self.patterns[0],
                alt_pams=self.patterns[1:],
                contigs=self.parameters["contigs"],
                fasta_index=fasta_index,
            )
            write_orientation_tsv(table, context["library"], orientation)
        self.state["counts"]["orientation_check"] = {
            "n_out_of_range": self.state["passes"][PASS_ORIENTATION]["counts"][
                "n_out_of_range"
            ],
            **orientation.counts(),
        }
        self.save()

        verdicts = Counter(orientation.verdicts.values())
        if library_reversed(orientation.verdicts.values()):
            message = (
                f"Reversed guides outnumber ok ones ({verdicts[REVERSED]} against "
                f"{verdicts[OK]}, counted per distinct spacer): the guide table may "
                "give the spacers 3' to 5' or reverse-complemented. Each guide's "
                f"perfect hits and their flanks are in {table}."
            )
            if check == "error":
                raise OrientationError(message)
            logger.warning(message)
        return orientation

    def _sites(
        self, context: dict, reference: Path, fasta_index: Path | None
    ) -> SiteResult:
        """Pass 2, with and without the added G, then the site checks."""
        library = context["library"]
        n_ambiguous = max(count_ambiguous(pattern) for pattern in self.patterns)
        sams = []
        for name, leading_g in ((PASS_SITES, False), (PASS_SITES_LEADING_G, True)):
            fastq = self.alignments / f"{name}.fastq"
            sam = self.alignments / f"{name}.sam.gz"
            with self.timed(name):
                n_reads = write_site_reads(
                    library, fastq, pams=self.patterns, leading_g=leading_g
                )
                if not n_reads:
                    self.attempted.add(name)
                    fastq.unlink(missing_ok=True)
                    sam.unlink(missing_ok=True)
                    continue
                self._gem_pass(
                    name,
                    context,
                    fastq=fastq,
                    n_reads=n_reads,
                    sam=sam,
                    max_mismatches=self.parameters["max_mismatches"],
                    n_ambiguous=n_ambiguous,
                    leading_g=leading_g,
                )
            sams.append(sam)

        with self.timed("site_check"):
            sites = read_sites(
                sams,
                library,
                reference,
                pam=self.patterns[0],
                alt_pams=self.patterns[1:],
                max_mismatches=self.parameters["max_mismatches"],
                contigs=self.parameters["contigs"],
                fasta_index=fasta_index,
            )
        passes = self.state["passes"]
        self.state["counts"]["site_search"] = {
            "n_out_of_range": sum(
                passes[name]["counts"]["n_out_of_range"]
                for name in (PASS_SITES, PASS_SITES_LEADING_G)
                if name in passes
            ),
            **sites.counts(),
        }
        self.save()
        return sites

    def _gem_pass(
        self,
        name: str,
        context: dict,
        *,
        fastq: Path,
        n_reads: int,
        sam: Path,
        max_mismatches: int,
        n_ambiguous: int,
        leading_g: bool,
    ) -> None:
        """Map one FASTQ, or reuse an earlier run's alignments of it."""
        key = {
            "fastq_md5": file_md5(fastq),
            "options": mapper_options(
                max_mismatches=max_mismatches,
                n_ambiguous=n_ambiguous,
                leading_g=leading_g,
            ),
            "index": context["index_identity"],
            "gem_version": context["version"],
        }
        self.attempted.add(name)
        earlier = _mapping(self.previous_passes.get(name))
        if (
            earlier.get("key") == key
            and sam.is_file()
            and file_md5(sam) == earlier.get("sam_md5")
        ):
            logger.info(
                "%s: reusing %s; its reads, options, index and GEM version match %s",
                name,
                sam,
                RUN_JSON,
            )
            earlier.pop("carried_over", None)
            self.state["passes"][name] = {**earlier, "reused": True}
            self.save()
            return

        sam.unlink(missing_ok=True)
        log = self.logs / f"gem-mapper.{name}.log"
        threads = self.parameters["threads"]
        logger.info(
            "%s: mapping %d reads with gem-mapper (%d threads); log in %s",
            name,
            n_reads,
            threads,
            log,
        )
        try:
            result = map_reads(
                context["index"],
                fastq,
                sam,
                log=log,
                max_mismatches=max_mismatches,
                n_ambiguous=n_ambiguous,
                leading_g=leading_g,
                threads=threads,
            )
        except GemError as error:
            raise _naming_the_guide(error, context["library"]) from error
        self.state["passes"][name] = {
            "fastq": self._relative(fastq),
            "n_reads": n_reads,
            "sam": self._relative(sam),
            "sam_md5": file_md5(sam),
            "log": self._relative(log),
            "key": key,
            "argv": list(result.argv),
            "reused": False,
            "seconds": round(result.seconds, 3),
            "counts": result.counts(),
            "mapped_per_row": _per_row(result.mapped_per_read, len(context["library"])),
        }
        self.save()

    def _write_outputs(
        self,
        library: Sequence[Guide],
        summaries: Sequence[GuideSummary],
        sites: SiteResult,
        orientation: OrientationResult | None,
        inputs: dict,
    ) -> None:
        max_mismatches = self.parameters["max_mismatches"]
        write_guides_tsv(
            self.outdir / "guides.tsv", summaries, max_mismatches=max_mismatches
        )
        write_sites_tsv(self.outdir / "sites.tsv", library, sites.sites)
        write_guides_bed(self.outdir / "guides.bed", summaries)
        write_cut_sites_bed(self.outdir / "cut_sites.bed", summaries)
        write_summary_tsv(
            self.outdir / "summary.tsv",
            summaries,
            counts=self.state["counts"],
            parameters={
                name: _printable(value)
                for name, value in {**self.parameters, **inputs}.items()
            },
        )

    def _alignments_per_guide(self, library: Sequence[Guide]) -> dict[str, dict]:
        """The kept alignment records of each distinct guide, by pass."""
        passes = self.state["passes"]
        orientation = passes.get(PASS_ORIENTATION, {}).get("mapped_per_row")
        site_rows = [
            passes[name]["mapped_per_row"]
            for name in (PASS_SITES, PASS_SITES_LEADING_G)
            if name in passes
        ]
        per_guide = {}
        for row, guide in enumerate(library):
            if guide.duplicate_of is not None:
                continue
            entry = {}
            if orientation is not None:
                entry["orientation"] = _at(orientation, row)
            entry["sites"] = sum(_at(rows, row) for rows in site_rows)
            per_guide[guide.guide_id] = entry
        return per_guide

    def _relative(self, path: Path) -> str:
        return os.path.relpath(path, self.outdir)


def _check_sidecar(
    sidecar: dict | None,
    index: Path,
    reference: Path,
    reference_md5: str,
    version: str,
) -> None:
    """Check the reference against the index's sidecar, when there is one."""
    if sidecar is None:
        logger.warning(
            "The index %s has no sidecar (%s), so the reference cannot be checked "
            "against it; only the contig names and lengths are. Indexes built with "
            "`cau guide-alignment index` have one.",
            index,
            sidecar_path(index).name,
        )
        return
    expected = sidecar.get("reference_md5")
    if expected != reference_md5:
        raise ValueError(
            f"The index {index} was built from a reference with MD5 {expected} "
            f"({sidecar.get('reference')}), but {reference} has MD5 "
            f"{reference_md5}. Pass the reference the index was built from, or "
            "rebuild the index."
        )
    built_with = sidecar.get("gem_version")
    if built_with and built_with != version:
        logger.warning(
            "The index was built with GEM %s; gem-mapper is %s.", built_with, version
        )


def _naming_the_guide(error: GemError, library: Sequence[Guide]) -> GemError:
    """The same error, naming the guide of the read gem-mapper crashed on."""
    read_name = error.read_name
    row_text = (read_name or "").partition(":")[0]
    if not (row_text.isascii() and row_text.isdigit()):
        return error
    row = int(row_text)
    if row >= len(library):
        return error
    guide = library[row]
    described = f"read {read_name!r}"
    detail = f"guide {guide.guide_id!r}"
    if guide.leading_g:
        detail += ", with the added G"
    message = str(error)
    if described in message:
        message = message.replace(described, f"{described} ({detail})", 1)
    else:
        message += f" Read {read_name!r} is {detail}."
    return GemError(message, log=error.log, read_name=read_name)


def _per_row(mapped_per_read: dict[str, int], n_guides: int) -> list[int]:
    """Mapped records by guide row, from read names ``<row>`` or ``<row>:<pam>``."""
    rows = [0] * n_guides
    for name, count in mapped_per_read.items():
        row_text = name.partition(":")[0]
        if row_text.isascii() and row_text.isdigit() and int(row_text) < n_guides:
            rows[int(row_text)] += count
    return rows


def _at(rows: Sequence[int], row: int) -> int:
    return rows[row] if row < len(rows) else 0


def _ordered_counts(values: Iterable[str], order: Sequence[str]) -> dict[str, int]:
    counts = Counter(values)
    return {name: counts[name] for name in order}


def _printable(value: object) -> object:
    """A parameter value for a TSV cell: tabs and line breaks escaped."""
    if isinstance(value, str):
        return value.replace("\t", "\\t").replace("\n", "\\n").replace("\r", "\\r")
    if isinstance(value, list):
        return [_printable(item) for item in value]
    return value


def _read_previous(path: Path) -> dict:
    """An earlier run's ``run.json``, or an empty dict."""
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        logger.warning("Ignoring %s, which cannot be read: %s", path, error)
        return {}
    return data if isinstance(data, dict) else {}


def _mapping(value: object) -> dict:
    """A copy of `value` if it is a dict (as read from JSON), else an empty dict."""
    return dict(value) if isinstance(value, dict) else {}


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


@contextmanager
def _log_to(path: Path) -> Iterator[None]:
    """Append this package's log records, from INFO up, to `path`."""
    package = logging.getLogger(__package__)
    handler = logging.FileHandler(path, mode="a", encoding="utf-8")
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(_LOG_FORMAT))
    level = package.level
    if package.getEffectiveLevel() > logging.INFO:
        package.setLevel(logging.INFO)
    package.addHandler(handler)
    try:
        yield
    finally:
        package.removeHandler(handler)
        handler.close()
        package.setLevel(level)

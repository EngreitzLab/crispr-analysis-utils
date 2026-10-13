"""Guide alignment: find every guide's genomic site.

A guide is designed, not sequenced, so its site is an exact match of the
spacer. The library is read and validated (`read_guides`), `scan_guides` walks
the reference once holding both readings of every guide, and `summarize`
classifies each one for the output tables (the ``write_*`` functions). `run`
chains it into one folder of outputs. The IUPAC helpers live in `iupac`, and
the scan itself is the package's Rust core.
"""

from . import iupac, pipeline, scan
from .library import Guide, read_guides
from .orientation import EVIDENCE, VERDICTS
from .pipeline import OrientationError, RunResult, run
from .scan import ScanResult, scan_guides
from .sites import Site
from .summary import (
    GuideSummary,
    fan_out_sites,
    guide_columns,
    summarize,
    write_columns_yaml,
    write_guides_bed,
    write_guides_tsv,
    write_qc_tsv,
)

__all__ = [
    "EVIDENCE",
    "Guide",
    "GuideSummary",
    "OrientationError",
    "RunResult",
    "ScanResult",
    "Site",
    "VERDICTS",
    "fan_out_sites",
    "guide_columns",
    "iupac",
    "pipeline",
    "read_guides",
    "run",
    "scan",
    "scan_guides",
    "summarize",
    "write_columns_yaml",
    "write_guides_bed",
    "write_guides_tsv",
    "write_qc_tsv",
]

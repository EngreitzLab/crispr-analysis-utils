"""Guide alignment: map a CRISPR guide library to its genomic sites.

The library is read and validated (`read_guides`), then written as reads for
the aligner: bare spacers for the orientation check (`write_spacer_reads`) and
one ``[G] + spacer + PAM`` read per guide and PAM for the site search
(`write_site_reads`). GEM3 maps them (`build_index`, `map_reads`, in `gem`).
The alignments come back through `check_orientation` and `read_sites`, which
verify every hit against the reference, and `summarize` classifies each guide
for the output tables (the ``write_*`` functions). `run` chains it all into
one folder of outputs. The IUPAC helpers live in `iupac`.
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
    write_cut_sites_bed,
    write_guides_bed,
    write_guides_tsv,
    write_orientation_tsv,
    write_qc_tsv,
    write_sites_tsv,
)

__all__ = [
    "Guide",
    "GuideSummary",
    "OrientationError",
    "RunResult",
    "EVIDENCE",
    "EVIDENCE",
    "ScanResult",
    "VERDICTS",
    "VERDICTS",
    "Site",
    "fan_out_sites",
    "guide_columns",
    "iupac",
    "pipeline",
    "read_guides",
    "run",
    "scan",
    "scan_guides",
    "summarize",
    "write_cut_sites_bed",
    "write_guides_bed",
    "write_guides_tsv",
    "write_orientation_tsv",
    "write_sites_tsv",
    "write_qc_tsv",
]

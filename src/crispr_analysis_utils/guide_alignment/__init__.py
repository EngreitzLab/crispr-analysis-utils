"""Guide alignment: map a CRISPR guide library to its genomic sites.

The library is read and validated (`read_guides`), then written as reads for
the aligner: bare spacers for the orientation check (`write_spacer_reads`) and
one ``[G] + spacer + PAM`` read per guide and PAM for the site search
(`write_site_reads`). The alignments come back through `check_orientation`
and `read_sites`, which verify every hit against the reference, and
`summarize` classifies each guide for the output tables (the ``write_*``
functions). The IUPAC helpers live in `iupac`.
"""

from . import iupac
from .library import Guide, read_guides, write_site_reads, write_spacer_reads
from .orientation import (
    OrientationResult,
    PerfectHit,
    check_orientation,
    library_reversed,
)
from .sites import Site, SiteResult, read_sites
from .summary import (
    GuideSummary,
    fan_out_sites,
    guide_columns,
    summarize,
    write_cut_sites_bed,
    write_guides_bed,
    write_guides_tsv,
    write_orientation_tsv,
    write_sites_tsv,
    write_summary_tsv,
)

__all__ = [
    "Guide",
    "GuideSummary",
    "OrientationResult",
    "PerfectHit",
    "Site",
    "SiteResult",
    "check_orientation",
    "fan_out_sites",
    "guide_columns",
    "iupac",
    "library_reversed",
    "read_guides",
    "read_sites",
    "summarize",
    "write_cut_sites_bed",
    "write_guides_bed",
    "write_guides_tsv",
    "write_orientation_tsv",
    "write_site_reads",
    "write_sites_tsv",
    "write_spacer_reads",
    "write_summary_tsv",
]

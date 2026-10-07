"""Guide alignment: map a CRISPR guide library to its genomic sites.

The library is read and validated (`read_guides`), then written as reads for
the aligner: bare spacers for the orientation check (`write_spacer_reads`) and
one ``[G] + spacer + PAM`` read per guide and PAM for the site search
(`write_site_reads`). The alignments come back through `check_orientation`
and `read_sites`, which verify every hit against the reference. The IUPAC
helpers live in `iupac`.
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

__all__ = [
    "Guide",
    "OrientationResult",
    "PerfectHit",
    "Site",
    "SiteResult",
    "check_orientation",
    "iupac",
    "library_reversed",
    "read_guides",
    "read_sites",
    "write_site_reads",
    "write_spacer_reads",
]

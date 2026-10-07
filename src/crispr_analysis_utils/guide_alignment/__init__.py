"""Guide alignment: map a CRISPR guide library to its genomic sites.

The library is read and validated (`read_guides`), then written as reads for
the aligner: bare spacers for the orientation check (`write_spacer_reads`) and
one ``[G] + spacer + PAM`` read per guide and PAM for the site search
(`write_site_reads`). The IUPAC helpers live in `iupac`.
"""

from . import iupac
from .library import Guide, read_guides, write_site_reads, write_spacer_reads

__all__ = [
    "Guide",
    "iupac",
    "read_guides",
    "write_site_reads",
    "write_spacer_reads",
]

"""Annotate aligned guides against genes.

`guide_alignment` answers where in the genome a guide is. This answers where
in a gene: exon, intron or intergenic, whether it falls in a promoter, which
gene and transcript, and how far it is from the nearest protein-coding TSS.
"""

from __future__ import annotations

from .annotate import AnnotationError, annotate_sites
from .db import Store, StoreError, build_store, load_store
from .pipeline import RunResult, run
from .promoters import PromoterError, promoters_from_store, read_promoters

__all__ = [
    "AnnotationError",
    "PromoterError",
    "RunResult",
    "Store",
    "StoreError",
    "annotate_sites",
    "build_store",
    "load_store",
    "promoters_from_store",
    "read_promoters",
    "run",
]

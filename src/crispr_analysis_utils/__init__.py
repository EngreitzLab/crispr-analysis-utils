"""Reusable helpers for CRISPR screen analysis."""

from importlib.metadata import PackageNotFoundError, version

from . import gem_mapper, guide_qc, utils
from .normalization import counts_per_million

try:
    __version__ = version("crispr-analysis-utils")
except PackageNotFoundError:  # a source tree that was never installed
    __version__ = "0+unknown"

__all__ = ["__version__", "counts_per_million", "gem_mapper", "guide_qc", "utils"]

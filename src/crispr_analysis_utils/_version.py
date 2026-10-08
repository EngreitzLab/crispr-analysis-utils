"""The installed package's version, apart so that subpackages can import it."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("crispr-analysis-utils")
except PackageNotFoundError:  # a source tree that was never installed
    __version__ = "0+unknown"

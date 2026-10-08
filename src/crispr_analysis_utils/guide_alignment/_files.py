"""File helpers shared by the GEM wrappers and the pipeline."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

_GZIP_MAGIC = b"\x1f\x8b"


def file_md5(path: str | os.PathLike[str]) -> str:
    """The MD5 checksum of a file, as 32 hexadecimal digits."""
    with open(path, "rb") as handle:
        digest = hashlib.file_digest(handle, lambda: hashlib.md5(usedforsecurity=False))
    return digest.hexdigest()


def is_gzip(path: str | os.PathLike[str]) -> bool:
    """Whether a file starts with the gzip magic bytes."""
    with open(path, "rb") as handle:
        return handle.read(2) == _GZIP_MAGIC


def check_plain_fasta(path: str | os.PathLike[str], what: str) -> Path:
    """The path of an existing, uncompressed FASTA, or an error naming `what`."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"The {what} {path} does not exist.")
    if is_gzip(path):
        raise ValueError(
            f"The {what} {path} is gzip-compressed; GEM and the site checks need "
            "an uncompressed FASTA."
        )
    return path


def write_json(path: str | os.PathLike[str], data: object) -> None:
    """Write JSON through a temporary file and a rename, so no reader sees half."""
    path = Path(path)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)

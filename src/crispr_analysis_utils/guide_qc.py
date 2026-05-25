"""Guide-level QC helpers."""

from __future__ import annotations

from pathlib import Path
import gzip

import pandas as pd


def guides_to_fastq(
    guides_input: str | Path | pd.DataFrame,
    output_path: str | Path,
    *,
    quality_char: str = "I",
    pam: str = "",
    add_leading_g: bool = False,
    id_col: int | str = 0,
    sequence_col: int | str = 1,
) -> Path:
    """Write guide sequences into a synthetic FASTQ file.

    Parameters
    ----------
    guides_input
        Either:
        - Path to a tab-separated file
        - pandas DataFrame
    output_path
        FASTQ output path. Use `.gz` suffix for gzip-compressed output.
    quality_char
        Single character repeated to sequence length for quality scores.
        Defaults to `I` (Q = 40 in Phred+33 encoding).
    pam
        Optional PAM sequence appended to each guide sequence.
    add_leading_g
        If `True`, always prepend `G` to each guide sequence.
    id_col
        Guide id column in `guides_input` (integer position or column name).
    sequence_col
        Guide sequence column in `guides_input` (integer position or column
        name).

    Returns
    -------
    pathlib.Path
        Path to the written FASTQ file.
    """
    if isinstance(guides_input, pd.DataFrame):
        guides = guides_input
    else:
        guides = pd.read_csv(guides_input, sep="\t", header=None)

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    open_fn = gzip.open if path.suffix == ".gz" else open
    with open_fn(path, "wt", encoding="utf-8") as handle:
        for _, row in guides.iterrows():
            if isinstance(id_col, int):
                record_id = str(row.iloc[id_col])
            else:
                record_id = str(row[id_col])

            if isinstance(sequence_col, int):
                sequence = str(row.iloc[sequence_col]).strip().upper()
            else:
                sequence = str(row[sequence_col]).strip().upper()
            if add_leading_g:
                sequence = f"G{sequence}"
            sequence = f"{sequence}{pam}"
            quality = quality_char * len(sequence)
            handle.write(f"@{record_id}\n{sequence}\n+\n{quality}\n")

    return path

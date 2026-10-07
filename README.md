# CRISPR Analysis Utils

Documentation: <https://engreitzlab.github.io/crispr-analysis-utils/>

Utilities and reusable scripts for CRISPR screen analysis, preprocessing, plotting,
and downstream workflows.

The Python package is `crispr_analysis_utils`, conventionally imported as `cau`.

The goal is that collaborators can clone the repository, install the functions
locally, and read rendered examples/API documentation on GitHub Pages.

## Install

From a local clone:

```bash
git clone https://github.com/EngreitzLab/crispr-analysis-utils.git
cd crispr-analysis-utils
python -m pip install -e .
```

For SAM/BAM guide filtering, install the `alignment` extra (pysam):

```bash
python -m pip install -e ".[alignment]"
```

## Use

```python
import pandas as pd
from crispr_analysis_utils import counts_per_million

counts = pd.DataFrame(
    {"sample_a": [100, 300], "sample_b": [50, 50]},
    index=["guide_1", "guide_2"],
)

cpm = counts_per_million(counts)
```

## Documentation

Functions are documented from their docstrings with MkDocs and mkdocstrings.

Build the documentation site locally with:

```bash
uv run --only-group docs mkdocs build --strict --site-dir site
```

GitHub Pages is configured through `.github/workflows/docs.yml`. In the GitHub
repository settings, set Pages to deploy from GitHub Actions.

## Repository layout

```text
src/crispr_analysis_utils/  Python package
docs/                      MkDocs source pages
tests/                     Tests
.github/workflows/         GitHub Actions workflows
```

## Development

Install pre-commit hooks:

```bash
python -m pip install pre-commit
pre-commit install
pre-commit run --all-files
```

Run the tests:

```bash
uv sync
uv run pytest
```

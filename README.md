# CRISPR Analysis Utils

Documentation: https://engreitzlab.github.io/crispr-analysis-utils/

Utilities and reusable scripts for CRISPR screen analysis, preprocessing, plotting,
and downstream workflows.

This repository is set up as a small two-language package:

- Python package: `crispr_analysis_utils`
- R package: `crisprAnalysisUtils`

The goal is that collaborators can clone the repository, install the functions
locally, and read rendered examples/API documentation on GitHub Pages.

## Install

From a local clone:

```bash
git clone https://github.com/EngreitzLab/crispr-analysis-utils.git
cd crispr-analysis-utils
python -m pip install -e .
```

```r
install.packages("remotes")
remotes::install_github("EngreitzLab/crispr-analysis-utils")
```

For Python plotting helpers, install the optional plotting dependencies:

```bash
python -m pip install -e ".[plots]"
```

## Use

Python:

```python
import pandas as pd
from crispr_analysis_utils import counts_per_million

counts = pd.DataFrame(
    {"sample_a": [100, 300], "sample_b": [50, 50]},
    index=["guide_1", "guide_2"],
)

cpm = counts_per_million(counts)
```

R:

```r
library(crisprAnalysisUtils)

counts <- matrix(
  c(100, 300, 50, 50),
  nrow = 2,
  dimnames = list(c("guide_1", "guide_2"), c("sample_a", "sample_b"))
)

cpm <- counts_per_million(counts)
```

## Documentation

Python functions are documented from docstrings with MkDocs and mkdocstrings.
R functions are documented from roxygen comments with pkgdown.

Build the combined documentation site locally with:

```bash
python -m pip install -e ".[docs]"
mkdocs build --strict --site-dir site
Rscript -e 'roxygen2::roxygenise(); pkgdown::build_site(new_process = FALSE, override = list(destination = "site/r"))'
```

GitHub Pages is configured through `.github/workflows/docs.yml`. In the GitHub
repository settings, set Pages to deploy from GitHub Actions.

## Repository layout

```text
src/crispr_analysis_utils/  Python package
R/                         R package functions
docs/                      MkDocs source pages
tests/                     Python and R tests
.github/workflows/         GitHub Actions workflows
```

## Development

Install pre-commit hooks:

```bash
python -m pip install pre-commit
pre-commit install
pre-commit run --all-files
```

Run Python tests:

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

Run R tests:

```r
testthat::test_local(load_package = "source")
```

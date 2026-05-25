# Getting Started

## Install Python helpers

```bash
git clone https://github.com/EngreitzLab/crispr-analysis-utils.git
cd crispr-analysis-utils
python -m pip install -e .
```

```python
from crispr_analysis_utils import counts_per_million
```

## Install R helpers

```r
install.packages("remotes")
remotes::install_github("EngreitzLab/crispr-analysis-utils")
```

```r
library(crisprAnalysisUtils)
```

## Add a New Utility

For Python:

1. Add the function under `src/crispr_analysis_utils/`.
2. Export it from `src/crispr_analysis_utils/__init__.py` if it is public.
3. Add a NumPy-style docstring with parameters, returns, and examples.
4. Add a focused test under `tests/`.

For R:

1. Add the function under `R/`.
2. Add roxygen comments with `@param`, `@return`, `@examples`, and `@export`.
3. Run `roxygen2::roxygenise()` when documentation changes.
4. Add a focused test under `tests/testthat/`.

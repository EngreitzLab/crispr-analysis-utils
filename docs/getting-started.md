# Getting Started

## Install

```bash
git clone https://github.com/EngreitzLab/crispr-analysis-utils.git
cd crispr-analysis-utils
python -m pip install -e .
```

```python
from crispr_analysis_utils import counts_per_million
```

## Add a New Utility

1. Add the function under `src/crispr_analysis_utils/`.
2. Export it from `src/crispr_analysis_utils/__init__.py` if it is public.
3. Add a NumPy-style docstring with parameters, returns, and examples.
4. Add a focused test under `tests/`.

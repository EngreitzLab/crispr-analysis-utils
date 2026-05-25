# Contributing

Keep reusable analysis helpers small, documented, and tested.

## Documentation Style

Python functions should use NumPy-style docstrings:

```python
def my_function(x: float) -> float:
    """Short summary.

    Parameters
    ----------
    x
        Description of the input.

    Returns
    -------
    float
        Description of the output.
    """
```

R functions should use roxygen comments:

```r
#' Short summary
#'
#' @param x Description of the input.
#' @return Description of the output.
#' @export
my_function <- function(x) {
  x
}
```

## Checks

Run these before opening a pull request:

```bash
python -m pip install -e ".[dev]"
python -m pytest
Rscript -e 'testthat::test_local(load_package = "source")'
```

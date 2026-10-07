# Contributing

Keep reusable analysis helpers small, documented, and tested.

## Documentation Style

Functions use NumPy-style docstrings:

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

## Checks

Run these before opening a pull request:

```bash
uv sync
uv run pytest
```

# Contributing

Keep analysis helpers small, documented and tested. `CLAUDE.md` in the
repository holds the full set of conventions and the reasons behind them; this
page is the working summary.

## Set up

```bash
git clone https://github.com/EngreitzLab/crispr-analysis-utils.git
cd crispr-analysis-utils
pixi install -e dev
pixi run -e dev install-hooks  # the commit and push hooks
```

With uv instead: `uv sync`, then `uv run pre-commit install`.

## What the hooks and CI enforce

- 4-space indentation and no tab characters. ruff formats Python, including
  the Python blocks in Markdown; markdownlint checks Markdown; JSON is
  formatted at 4 spaces. YAML stays at 2.
- `ruff check` passes, with the rules set in `pyproject.toml`.
- No file over 512 KB, except the generated `uv.lock` and `pixi.lock`. Data
  never goes in git; build test fixtures inside the tests.
- No commit on `main`: work on a branch and open a pull request. `main`
  accepts only pull requests whose CI checks pass.
- `uv.lock` and `pixi.lock` match `pyproject.toml`: after changing a
  dependency, run `uv lock` and `pixi lock` and commit both.

## Commit messages

Start with an UPPERCASE verb: `ADD`, `FIX`, `UPDATE`, `REWRITE`, `RELEASE` or
`REMOVE`, then say what changed, for example
`FIX counts_per_million accepting a NaN pseudocount`.

## Adding a module

1. Write the module under `src/crispr_analysis_utils/`, with NumPy-style
   docstrings.
2. Re-export its public functions from `src/crispr_analysis_utils/__init__.py`
   if they belong at the top level.
3. Add its skill, `src/crispr_analysis_utils/skills/cau-<module>/SKILL.md`,
   and, only if using the module is a long multi-step job, an agent in
   `src/crispr_analysis_utils/agents/` that preloads that skill. See
   [Agents and skills](agents-and-skills.md).
4. Map the new skill or agent to the module in `COVERAGE`, in
   `tests/test_skills.py`.
5. Add an API reference page under `docs/python/reference/`, and its entry in
   the `mkdocs.yml` nav and on the [Agents and skills](agents-and-skills.md)
   page.
6. Add focused tests under `tests/`. Tests must not write outside `tmp_path`.
7. If the module has a command, add a `cau` subcommand module under
   `src/crispr_analysis_utils/cli/` and register it in `COMMANDS`. Import
   optional dependencies inside its `run()`.

## Documentation style

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

A docstring says what the function does, takes, returns and raises. Design
rationale and measurements belong in the documentation pages.

## Checks

Run these before opening a pull request:

```bash
pixi run -e dev test
pixi run -e dev lint
pixi run -e docs docs-build
```

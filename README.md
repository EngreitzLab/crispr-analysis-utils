# crispr-analysis-utils

Python utilities and command-line tools for analysing CRISPR perturbation
screens.

Documentation: <https://engreitzlab.github.io/crispr-analysis-utils/>

<!-- --8<-- [start:install] -->
## Install

The package is not on PyPI yet. Install it from GitHub:

```bash
pip install "crispr-analysis-utils @ git+https://github.com/EngreitzLab/crispr-analysis-utils"
```

Filtering guide alignments (SAM/BAM) needs pysam, from the `alignment` extra:

```bash
pip install "crispr-analysis-utils[alignment] @ git+https://github.com/EngreitzLab/crispr-analysis-utils"
```

To work on the package, or to get the GEM3 aligner as well (linux-64 only),
clone the repository and use [pixi](https://pixi.sh) or
[uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/EngreitzLab/crispr-analysis-utils.git
cd crispr-analysis-utils
pixi install  # conda + PyPI environment; includes gem3-mapper on linux-64
uv sync       # or a uv virtual environment, with the development tools
```
<!-- --8<-- [end:install] -->

<!-- --8<-- [start:usage] -->
## Use

```python
import pandas as pd

import crispr_analysis_utils as cau

counts = pd.DataFrame(
    {"sample_a": [100, 300], "sample_b": [50, 50]},
    index=["guide_1", "guide_2"],
)
cpm = cau.counts_per_million(counts)
```

| Module | What it does |
| --- | --- |
| `cau.counts_per_million` | Counts-per-million normalization of count matrices |
| `cau.guide_qc` | Guide FASTQs for alignment, and QC filtering of guide alignments |
| `cau.gem_mapper` | GEM3 genome indexing and guide mapping |

The `cau` command line tool comes with the package:

```bash
cau --help
```
<!-- --8<-- [end:usage] -->

<!-- --8<-- [start:skills] -->
## Claude Code skills and agents

Every analysis module comes with a [Claude Code](https://claude.com/claude-code)
skill that teaches Claude to use it: what to ask you before running anything,
the real defaults, and how to read the outputs. Modules whose use is a long
multi-step job can also get an agent. Install them one of two ways, not both
(with both, every skill is listed twice):

- **With the package**, so the skills match the version you have installed:

    ```bash
    cau install-skills                   # into ~/.claude (or $CLAUDE_CONFIG_DIR)
    cau install-skills --target .claude  # into the current project only
    cau install-skills --list            # what is bundled
    ```

    After upgrading the package, run it again with `--force`.

- **As a Claude Code plugin**, which follows the `main` branch:

    ```text
    /plugin marketplace add EngreitzLab/crispr-analysis-utils
    /plugin install crispr-analysis-utils@crispr-analysis-utils
    ```

    Update it with `/plugin marketplace update crispr-analysis-utils`, or add
    the marketplace as `EngreitzLab/crispr-analysis-utils#<tag>` to follow a
    tag instead of `main`.
<!-- --8<-- [end:skills] -->

## Development

```bash
pixi run -e dev install-hooks  # once per clone: the commit and push hooks
pixi run -e dev test
pixi run -e dev lint           # every hook on every tracked file
pixi run -e docs docs          # the documentation site, served locally
```

With uv instead: `uv sync`, `uv run pytest`, `uv run pre-commit install`.

`main` changes only through pull requests whose CI checks pass. The
conventions (4-space indentation, `ruff check`, no file over 512 KB, UPPERCASE
commit verbs) are in the
[contributing guide](https://engreitzlab.github.io/crispr-analysis-utils/contributing/)
and in `CLAUDE.md`.

## License

MIT. See `LICENSE`.

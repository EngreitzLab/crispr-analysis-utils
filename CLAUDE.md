# crispr-analysis-utils

Python utilities and the `cau` command-line tool for analysing CRISPR
perturbation screens. The distribution is `crispr-analysis-utils`; the import
is `import crispr_analysis_utils as cau`.

**Python only.** The repository began as a Python and R pair. The R package
held one function, `counts_per_million`, which was ported (with R's input
validation) and then removed. Don't reintroduce R: port it.

**Agent-ready by construction.** Every analysis module ships a Claude Code
skill inside the package, and the repository is also a Claude Code plugin
marketplace serving the same files. "Skills and agents" below has the rule,
and `tests/test_skills.py` enforces it.

`README.md` is the user-facing reference (install, usage, the two skill
install routes). This file is the contributor-facing complement: don't
duplicate the README here, point at it.

## Layout

src layout: the package is importable only once installed, and editable is
fine (`uv sync` and `pixi install` both do it).

- `src/crispr_analysis_utils/` -- the package, **and the Claude Code plugin
  root**: `.claude-plugin/marketplace.json` points its one plugin here.
    - `normalization.py`  -- `counts_per_million`.
    - `guide_qc.py`       -- `guides_to_fastq` and `filter_guide_alignments`
                             (pysam, imported inside the function).
    - `gem_mapper.py`     -- GEM3 indexing and mapping, through `run_shell_cmd`.
    - `utils.py`          -- `run_shell_cmd`. Infrastructure: no skill.
    - `cli/`              -- the `cau` console script, one module per
                             subcommand (`install_skills.py`). Infrastructure:
                             no skill.
    - `skills/<name>/`    -- one skill per analysis module: `SKILL.md`, plus
                             `references/` for detail. Data, not packages.
    - `agents/<name>.md`  -- subagents. The directory does not exist until the
                             first agent.
- `tests/`              -- pytest. `conftest.py` runs every test inside its
                           own `tmp_path`.
- `docs/`               -- the MkDocs-Material site, published to GitHub Pages
                           by `.github/workflows/docs.yml`.
- `scripts/`            -- the guide-alignment pipeline script and its sbatch
                           template. Not shipped in the wheel or the sdist.
- `.claude-plugin/marketplace.json` -- makes the repository a Claude Code
  plugin marketplace.

## Commands

```bash
uv sync                          # .venv with the dev group, editable install
uv run pytest

pixi install -e dev              # conda + PyPI; gem3-mapper on linux-64
pixi run -e dev test
pixi run -e dev lint             # every pre-commit hook on every tracked file
pixi run -e dev install-hooks    # once per clone: commit and push hooks
pixi run -e docs docs            # serve the site locally
pixi run -e docs docs-build      # mkdocs build --strict, as CI runs it

uv build && uvx twine check --strict dist/*
uv lock && pixi lock             # after any dependency change; commit both

cau install-skills --list
claude plugin validate --strict . && claude plugin validate --strict src/crispr_analysis_utils
```

Python 3.11+ (`requires-python`). CI tests 3.11 to 3.14.

## Conventions

- **Indentation: 4 spaces, no tabs.** What enforces it depends on the file
  type: ruff format for Python and for the Python blocks inside Markdown (ruff
  0.16 formats Markdown by default, and `*.md` is deliberately not excluded);
  markdownlint MD007 for nested Markdown lists, which MkDocs needs at 4;
  `pretty-format-json` for JSON. TOML and shell follow `.editorconfig` by
  convention only. YAML is the exception, at 2. The `forbid-tabs` hook rejects
  a tab anywhere except tab-separated data fixtures (`.tsv`, `.bed`, `.sam`,
  `.vcf`, `.fastq`).
- **The ruff config lives in `pyproject.toml`.** Its rule selection is
  explicit: ruff 0.16 replaced its 59-rule default set with 413 rules, and
  `select` replaces the default. `E501` is off because the formatter owns line
  length. ruff runs from pre-commit, whose `rev` is the only ruff pin.
- **No file over 512 KB.** `check-added-large-files --maxkb=512 --enforce-all`
  runs at commit and again at push. Without `--enforce-all` it would check only
  newly added files, and CI's `--all-files` run would check nothing. CI also
  checks every blob a pull request adds, so a large file added and then
  deleted inside the PR cannot reach history. `uv.lock` and `pixi.lock` are
  exempt because they are generated. Data never goes in git (`.gitignore`
  covers the usual formats); build test fixtures inside the tests.
- **`main` is PR-only**, by branch protection (Settings, Branches) that
  requires the nine CI checks, admins included. Work on a branch. The
  `no-commit-to-branch` hook catches a local commit to `main`; it always runs,
  so the pixi `lint` task and CI set `SKIP=no-commit-to-branch`. The required
  checks are listed by job name (`test-pip (3.11)` and so on): renaming a job,
  or changing the Python matrix, means updating them, or every pull request
  waits on a check that never reports.
- **Commit messages: UPPERCASE verb prefix** -- `ADD`, `FIX`, `UPDATE`,
  `REWRITE`, `RELEASE`, `REMOVE`.
- **Docs accuracy is a hard rule.** Every concrete detail (defaults, column
  meanings, versions, flags) must be confirmable from source. If you can't
  verify it, omit it. The guide-QC output tables in the docs drifted from the
  code from their first version until the Python-only setup, because the code
  changes that followed left the docs alone.
- **Docstrings: NumPy style** (`docstring_style: numpy` in `mkdocs.yml`). A
  docstring says what a function does, takes, returns and raises. Design
  rationale and measurements belong in the docs pages; leave at most a
  one-line pointer in the source.
- **The docs include, they don't copy.** `docs/` pages pull `README.md`
  sections through `pymdownx.snippets` markers (`<!-- --8<-- [start:name] -->`),
  and `docs/python/guide-qc.md` pulls its output tables from the
  `cau-guide-qc` skill's `references/outputs.md`. Don't delete a marker. Links
  inside a marked section must be absolute URLs: a relative one resolves
  differently on GitHub and in the site, and `mkdocs build --strict` fails.
- **Tests never write outside `tmp_path`.** `tests/conftest.py` changes into
  it because `filter_guide_alignments` writes "auto" outputs to the working
  directory. CI fails if a test leaves files in the checkout.
- **Plots must be colorblind-safe**: Okabe-Ito for categories, `cividis` for
  continuous scales, never `jet` or `rainbow`. Nothing plots yet; this applies
  to whatever does first.

## Skills and agents

- **Every analysis module has a skill**, `skills/cau-<module>/SKILL.md`. An
  analysis module is any module of the package except those starting with `_`
  and the infrastructure listed in `tests/test_skills.py`
  (`INFRASTRUCTURE = ("cli", "utils")`).
- **An agent is optional.** Add one when using a module means a long
  multi-step job whose intermediate output would flood the conversation (run a
  pipeline, read its logs, summarize). It gets its own name, such as
  `cau-<module>-runner`, and preloads the module's skill through its
  `skills:` frontmatter instead of repeating it.
- **`COVERAGE` in `tests/test_skills.py` maps** each skill and agent to the
  modules it covers. Merging two modules' skills, or splitting one, happens
  there.
- **Keep skill frontmatter to the six Agent Skills fields** (`name`,
  `description`, `license`, `compatibility`, `metadata`, `allowed-tools`):
  claude.ai uploads reject anything else. `name` equals the directory;
  `description` fits in 1024 characters and says when to use the skill;
  SKILL.md stays under 500 lines, with detail in `references/`. Claude
  Code-only behaviour belongs in an agent.
- **A skill quotes the real defaults** and tells Claude what to ask the user
  before running anything (cherimoya's "ask, don't guess"). Update a skill in
  the same commit as the API it describes.
- **Adding a module**: follow the checklist in `docs/contributing.md`.
- **Two install routes, one set of files.** `cau install-skills` copies them
  from the installed package, so they match that version. The plugin serves
  `src/crispr_analysis_utils/` from `main`. Users should pick one, or every
  skill is listed twice.

## Gotchas

- **The package directory is the plugin root.** Nothing in
  `src/crispr_analysis_utils/` may use a plugin component name (`commands`,
  `hooks`, `bin`, `workflows`, `themes`, `monitors`, `output-styles`,
  `settings.json`, `.mcp.json`, `.lsp.json`), or Claude Code loads it as part
  of the plugin. A test enforces it. `workflows` is the likely collision, so
  name such a subpackage something else.
- **The plugin pins no version, on purpose.** Plugin users get each commit on
  `main` (`/plugin marketplace update`), so a skill fix needs no release. The
  marketplace entry is the plugin manifest; there is no `plugin.json`. Adding a
  `version` would freeze plugin users until the next bump.
- **Plugin agents ignore** `hooks`, `mcpServers`, `permissionMode` and
  `initialPrompt` in their frontmatter; those take effect only for agents
  installed by `cau install-skills`. Before relying on a plugin agent's
  `skills:` preload, check that a bare `cau-x` resolves to the plugin's
  namespaced `crispr-analysis-utils:cau-x`.
- **`cau install-skills` reads frontmatter without PyYAML**, which is not a
  runtime dependency. `tests/test_skills.py` checks the reader against PyYAML
  on every bundled file. Write frontmatter in the shapes it handles:
  `key: value`, folded `key: >-` blocks, and `[a, b]` or `- item` lists.
- **Hatchling ships every file under `src/crispr_analysis_utils/` that the
  root `.gitignore` does not exclude, tracked or not.** It ignores
  `.git/info/exclude`, nested `.gitignore` files and global excludes, so keep
  scratch files out of the package directory. CI builds from a clean checkout
  and checks that the wheel holds no `tests/`, `scripts/` or stray scripts, and
  that it holds every skill file; a local `uv build` has no such guard.
- **gem3-mapper is linux-64 only here.** bioconda builds it for linux-64 and
  osx-64, not Apple silicon, so the pixi environments include it only on
  linux-64. The unit tests mock the binaries; CI's `test` job checks that they
  run.
- **pysam is optional** (the `alignment` extra) and imported inside
  `filter_guide_alignments`. Keep it out of module-level imports. The dev group
  installs it for the tests.
- **Tool versions are pinned in one place each.** uv: `[tool.uv]
  required-version` (setup-uv reads it in CI). pixi: `requires-pixi`, plus
  `PIXI_VERSION` in `.github/workflows/ci.yml`, where setup-pixi needs the
  leading `v`. ruff, markdownlint and the hygiene hooks: the revs in
  `.pre-commit-config.yaml`. Bump deliberately.
- **`uv.lock` is committed, but `test-pip` installs unlocked** against the
  `pyproject.toml` ranges, so an upstream release that breaks the package shows
  up in CI instead of hiding behind pins. The `lockfiles` job keeps both
  lockfiles current.
- **The docs stay on MkDocs-Material**, which is in maintenance mode (its
  9.7.7 release notes schedule end of life for 2026-11-05). Material caps
  `mkdocs<2`, and that cap matters: MkDocs 2.0, in development, is a rewrite
  without plugins. The successor, Zensical, is still alpha; moving to it is a
  separate change.
- **`guide_qc` behaviour that trips readers up** (may change in the
  guide-alignment redesign; the `cau-guide-qc` skill states it for users too):
    - per-guide counts are keyed by read name (QNAME), not by sequence;
    - inside `filter_guide_alignments`, `guide_id` is the read's sequence
      recovered 5' to 3' (reverse-complemented for reverse-strand records);
    - a secondary record with `SEQ='*'` reuses the sequence seen earlier for the
      same QNAME;
    - BED column 4 is the protospacer (the read minus the PAM and any
      soft-clipped leading base) since 81c6db4;
    - on the reverse strand the PAM is at the start of the read as SAM stores
      it, and the rules are mirrored except at the spacer's two ends: a deletion
      between spacer and PAM passes on `+` but fails on `-`, and one right after
      the read's first base fails on `+` but passes on `-`;
    - the "auto" outputs go next to the unique SAM, else the multi SAM, else the
      working directory, and the SAM outputs are written only when both are
      given;
    - `alias_by_guide_id` is looked up by that full read sequence, PAM included,
      not by the spacer.
- **`test_filter_guide_alignments_reports_sequence_5prime_to_3prime` is a
  strict xfail.** It has been stale since f3ba75e, and its input stores the
  PAM-first orientation as SAM SEQ. Fix or replace it in the guide-alignment
  redesign; being strict, it fails loudly if it starts passing.

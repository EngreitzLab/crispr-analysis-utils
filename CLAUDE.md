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
fine (`pixi install` does it, building the Rust core as it goes).

The build backend is **maturin**, not hatchling, because the package ships a
compiled extension. `Cargo.toml` at the root declares the crate, its sources
are in `rust/`, and `[tool.maturin]` keeps `python-source = "src"` so the
Python layout is unchanged. The wheel is **abi3** (`abi3-py311`), so one wheel
per platform serves every supported Python rather than one per version.

- `src/crispr_analysis_utils/` -- the package, **and the Claude Code plugin
  root**: `.claude-plugin/marketplace.json` points its one plugin here.
    - `accession.py`      -- `guide_accession` encodes a guide (`G...NGG`,
                             SpCas9) rather than hashing it, so two guides
                             never collide and `guide_sequence` decodes one
                             back. It validates and never repairs. Changing
                             the scheme renames every guide.
    - `normalization.py`  -- `counts_per_million`.
    - `guide_alignment/`  -- maps a guide library to its genomic sites with
                             an exact scan: `library` (the guide table),
                             `scan` (the genome pass, over `_core`), `sites`
                             (the `Site` record), `orientation` (the
                             verdicts), `summary` (classes, tables),
                             `pipeline` (`run`), `iupac`. One analysis module,
                             one skill.
    - `utils.py`          -- `run_shell_cmd`. Infrastructure: no skill.
    - `_core.abi3.so`     -- the compiled Rust core, built in place by an
                             editable install and gitignored. Its sources are
                             in `rust/src/`.
    - `_version.py`       -- `__version__`. Import it from here: `__init__`
                             imports `guide_alignment` first, so `from ..
                             import __version__` inside a subpackage is a
                             circular import.
    - `cli/`              -- the `cau` console script, one module per
                             subcommand (`guide_alignment.py`,
                             `install_skills.py`). Infrastructure: no skill.
    - `skills/<name>/`    -- one skill per analysis module: `SKILL.md`, plus
                             `references/` for detail. Data, not packages.
    - `agents/`           -- subagents, each preloading the skill of the
                             module it runs. Empty: no module's use is long
                             enough to need one.
- `Cargo.toml`, `rust/` -- the Rust crate behind `_core`, and its sources.
                           `Cargo.lock` is committed, like the other two
                           lockfiles.
- `tests/`              -- pytest. `conftest.py` runs every test inside its
                           own `tmp_path`. `test_core_scan.py` checks the Rust
                           scan against a naive Python scan of every position
                           on every strand; that oracle is what makes the scan
                           trustworthy.
- `docs/`               -- the MkDocs-Material site, published to GitHub Pages
                           by `.github/workflows/docs.yml`.
- `.claude-plugin/marketplace.json` -- makes the repository a Claude Code
  plugin marketplace.

## Commands

```bash
pixi install -e dev              # conda + PyPI: the package, its Rust core and the dev tools
pixi run -e dev test
pixi run -e dev lint             # every pre-commit hook on every tracked file
pixi run -e dev install-hooks    # once per clone: commit and push hooks
pixi run -e dev cau guide-alignment run --help
pixi run -e docs docs            # serve the site locally
pixi run -e docs docs-build      # mkdocs build --strict, as CI runs it

uv lock && pixi lock             # after any dependency change; commit both
uv build && uvx twine check --strict dist/*

pixi run -e dev cau install-skills --list
claude plugin validate --strict . && claude plugin validate --strict src/crispr_analysis_utils
```

pixi is the default, and its `dev` environment is the one environment for
everything: run the tests, the linters and `cau` there, not in a uv `.venv`.
Every environment is native on all three platforms, since the genome scan is
the package's own Rust core. uv locks `uv.lock`, builds, and runs the CI jobs
that need no conda: the `test-pip` matrix, `build` and `docs-build`.

Python 3.12+ (`requires-python`). CI tests 3.12 to 3.14; the pixi environments
use 3.13. The floor is 3.12 because `pyranges1` has never shipped for 3.11.
The Rust core stays `abi3-py311`: that wheel serves every supported Python,
and narrowing it would buy nothing.

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
  verify it, omit it. The old guide-QC output tables drifted from the code
  from their first version until the Python-only setup, because the code
  changes that followed left the docs alone.
- **Docstrings: NumPy style** (`docstring_style: numpy` in `mkdocs.yml`). A
  docstring says what a function does, takes, returns and raises. Design
  rationale and measurements belong in the docs pages; leave at most a
  one-line pointer in the source.
- **The docs include, they don't copy.** `docs/` pages pull `README.md`
  sections through `pymdownx.snippets` markers (`<!-- --8<-- [start:name] -->`),
  and `docs/guide-alignment.md` pulls its output tables and its migration
  note from the `cau-guide-alignment` skill's `references/outputs.md` and
  `references/migration.md`. Don't delete a marker. Links
  inside a marked section must be absolute URLs: a relative one resolves
  differently on GitHub and in the site, and `mkdocs build --strict` fails.
- **Tests never write outside `tmp_path`.** `tests/conftest.py` changes into
  it, so a test that writes to a relative path stays inside it. CI fails if a
  test leaves files in the checkout.
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
- **maturin ships every file under `src/crispr_analysis_utils/` that the root
  `.gitignore` does not exclude, tracked or not** (checked, not assumed). So
  keep scratch files out of the package directory. CI builds from a clean
  checkout and checks that the wheel holds no `tests/`, `scripts/` or stray
  scripts, and that it holds every skill file; a local `uv build` has no such
  guard.
- **pysam is a core dependency, imported inside functions.** The modules that
  read SAM, BAM and FASTA import it where they open a file, never at module
  level, so `import crispr_analysis_utils` and `cau --help` do not load it.
  Keep it that way.
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
- **`guide_alignment` behaviour that trips readers up** (the
  `cau-guide-alignment` skill states the user-facing part):
    - a guide is designed, not sequenced, so its site is an **exact** match of
      the spacer. There is no aligner, no index and no error budget. The three
      questions are where the spacer matches, what sits 3' of it (the PAM) and
      what sits 5' of it;
    - `guide_id` is the sequence that is aligned plus the PAM, and `alias` is
      the table's name for the entry. Two entries of one sequence are one
      guide; `duplicate_of` only records that, and is not an output column;
    - `g_added` comes from the length (`n_extra = len(aligned) - spacer_length`),
      never from whether the first base is G. The old rule dropped the real
      added G of a `GG` 21-mer and then called it a duplicate;
    - `mapped_sequence` is the guide's own bases that pair with the genome:
      the spacer, keeping a base 5' of it only where the genome holds the
      same base. The coordinates stay on the declared protospacer, so a G that
      pairs by chance never shifts a start;
    - coordinates are 0-based half-open protospacers, without the PAM or the
      added base, and `cut_site` is a column, the boundary `b` (`end - 3` on
      `+`, `start + 3` on `-`). The IGVF-named columns keep these
      definitions, which differ from the lab's DC_TAP tables (the guide page
      has the comparison);
    - classes count sites whose genomic PAM matches. A site whose PAM does
      not stays in the table with `pam_matches` `FALSE` and counts nowhere,
      so a spacer matching twice with one PAM is `unique`;
    - `run` writes only inside `outdir`, and only three folders: the results
      in `alignments/`, this invocation's record in `run_metadata/` and the
      log in `logs/`. Nothing of the run's own sits beside an input;
    - `alignments/` holds four files: `guides_alignment_metadata.tsv`, one
      row per library entry per site and the only one with the whole picture;
      `guides_mapping.bed`, its PAM-bearing rows for a browser, named by
      spacer; `columns.yaml`, generated from the same definitions as the
      table, so it always describes the folder it sits in; and `qc.tsv`,
      counts;
    - orientation is decided from both readings, held in the same scan: an
      entry is `reverse_complement` only when its reversed reading has a
      PAM-bearing site and the sequence as given does not. A `CCN` 5' of a
      match is not evidence; it falls there by chance one time in sixteen.
      Only a table reversed throughout stops a run.

"""The bundled Claude Code skills and agents, and the plugin that serves them.

Every analysis module must be covered by a skill. An agent is optional, for a
module whose use is a long multi-step job, and it preloads that module's skill
instead of repeating it. CLAUDE.md has the full rule.
"""

import json
import pkgutil
import re
from pathlib import Path

import pytest
import yaml

import crispr_analysis_utils
from crispr_analysis_utils.cli.install_skills import read_frontmatter

PACKAGE = Path(crispr_analysis_utils.__file__).parent
SKILLS = PACKAGE / "skills"
AGENTS = PACKAGE / "agents"
REPO = Path(__file__).resolve().parent.parent

# Which analysis modules each skill or agent covers.
COVERAGE = {
    "cau-normalization": ("normalization",),
    "cau-guide-alignment": ("guide_alignment",),
    "cau-guide-alignment-runner": ("guide_alignment",),
}

# Modules that are plumbing rather than analysis, so need no skill.
INFRASTRUCTURE = ("cli", "utils")

# The package directory is also the Claude Code plugin root, where these names
# would be read as plugin components.
PLUGIN_DEFAULT_NAMES = (
    "commands",
    "hooks",
    "bin",
    "workflows",
    "themes",
    "monitors",
    "output-styles",
    "settings.json",
    ".mcp.json",
    ".lsp.json",
)

NAME_PATTERN = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
MD_MENTION = re.compile(r"`([A-Za-z0-9_./-]*\.md)`")

SKILL_DIRS = sorted(
    p for p in SKILLS.iterdir() if p.is_dir() and p.name != "__pycache__"
)
AGENT_FILES = sorted(AGENTS.glob("*.md")) if AGENTS.is_dir() else []


def frontmatter(path):
    text = path.read_text(encoding="utf-8")
    assert text.startswith("---\n"), f"{path} must open with a frontmatter block"
    return yaml.safe_load(text.split("---\n", 2)[1])


def analysis_modules():
    return sorted(
        info.name
        for info in pkgutil.iter_modules([str(PACKAGE)])
        if not info.name.startswith("_") and info.name not in INFRASTRUCTURE
    )


def test_every_analysis_module_is_covered_by_a_skill():
    skills = {path.name for path in SKILL_DIRS}
    covered = {
        module
        for name, modules in COVERAGE.items()
        if name in skills
        for module in modules
    }
    missing = sorted(set(analysis_modules()) - covered)
    assert not missing, f"analysis modules without a skill: {missing}"


def test_every_skill_and_agent_is_in_the_coverage_map():
    bundled = [path.name for path in SKILL_DIRS] + [path.stem for path in AGENT_FILES]
    assert sorted(bundled) == sorted(COVERAGE), (
        "COVERAGE must list exactly the bundled skills and agents"
    )
    modules = set(analysis_modules())
    for name, covered in COVERAGE.items():
        stale = sorted(set(covered) - modules)
        assert not stale, f"{name} covers modules that do not exist: {stale}"


def test_skill_and_agent_names_are_unique():
    names = [path.name for path in SKILL_DIRS] + [path.stem for path in AGENT_FILES]
    assert len(names) == len(set(names))


def test_every_package_directory_is_a_regular_package():
    """A directory without __init__.py would escape the coverage rule."""
    for path in PACKAGE.rglob("*"):
        if not path.is_dir() or "__pycache__" in path.parts:
            continue
        relative = path.relative_to(PACKAGE)
        if relative.parts[0] in ("skills", "agents"):
            assert not (path / "__init__.py").exists(), (
                f"{relative} is data, not a Python package"
            )
        else:
            assert (path / "__init__.py").is_file(), f"{relative} has no __init__.py"


def test_package_root_uses_no_plugin_default_names():
    clashes = sorted(
        p.name for p in PACKAGE.iterdir() if p.name in PLUGIN_DEFAULT_NAMES
    )
    assert not clashes, f"plugin component names in the package root: {clashes}"


@pytest.mark.parametrize("skill_dir", SKILL_DIRS, ids=lambda p: p.name)
def test_skill_frontmatter(skill_dir):
    fields = frontmatter(skill_dir / "SKILL.md")
    name = fields.get("name")
    assert name == skill_dir.name, "name must equal the directory name"
    assert NAME_PATTERN.match(name) and len(name) <= 64
    assert name.startswith("cau-")
    description = fields.get("description")
    assert description, "a skill needs a description"
    assert len(description) <= 1024


@pytest.mark.parametrize("skill_dir", SKILL_DIRS, ids=lambda p: p.name)
def test_skill_body_stays_short(skill_dir):
    lines = (skill_dir / "SKILL.md").read_text(encoding="utf-8").count("\n")
    assert lines <= 500, "move detail into references/"


@pytest.mark.parametrize("skill_dir", SKILL_DIRS, ids=lambda p: p.name)
def test_skill_mentions_resolve(skill_dir):
    """Every `*.md` a skill names is SKILL.md or a references/ path that exists."""
    documents = [skill_dir / "SKILL.md", *sorted(skill_dir.glob("references/*.md"))]
    for document in documents:
        for target in MD_MENTION.findall(document.read_text(encoding="utf-8")):
            assert target == "SKILL.md" or target.startswith("references/"), (
                f"{document.name} names `{target}`: write the skill-relative path"
            )
            assert (skill_dir / target).is_file(), (
                f"{document.name} points at `{target}`, which does not exist"
            )


@pytest.mark.parametrize("skill_dir", SKILL_DIRS, ids=lambda p: p.name)
def test_every_reference_is_linked_from_its_skill(skill_dir):
    router = (skill_dir / "SKILL.md").read_text(encoding="utf-8")
    for reference in sorted(skill_dir.glob("references/*.md")):
        assert f"`references/{reference.name}`" in router, (
            f"references/{reference.name} is not linked from SKILL.md"
        )


def test_agents():
    skills = {path.name for path in SKILL_DIRS}
    for path in AGENT_FILES:
        fields = frontmatter(path)
        assert fields.get("name") == path.stem, (
            f"{path.name}: name must equal the file stem"
        )
        assert path.stem.startswith("cau-")
        assert fields.get("description"), f"{path.name} needs a description"
        preloads = fields.get("skills") or []
        if isinstance(preloads, str):
            preloads = [s.strip() for s in preloads.split(",")]
        assert preloads, f"{path.name} must preload its module's skill"
        unknown = sorted(set(preloads) - skills)
        assert not unknown, (
            f"{path.name} preloads skills that are not bundled: {unknown}"
        )
        modules = set(COVERAGE.get(path.stem, ()))
        preloaded_modules = {m for skill in preloads for m in COVERAGE.get(skill, ())}
        assert modules & preloaded_modules, (
            f"{path.name} must preload a skill covering one of its modules"
        )


@pytest.mark.parametrize(
    "path",
    [d / "SKILL.md" for d in SKILL_DIRS] + AGENT_FILES,
    ids=lambda p: p.parent.name if p.name == "SKILL.md" else p.stem,
)
def test_installer_reads_frontmatter_like_yaml(path):
    """cau install-skills parses frontmatter without PyYAML; keep them in step."""
    expected = frontmatter(path)
    parsed = read_frontmatter(path)
    for key in ("name", "description", "skills"):
        if key in expected:
            assert parsed[key] == expected[key], f"{path}: {key} differs"


def test_marketplace_serves_the_package_directory():
    manifest = REPO / ".claude-plugin" / "marketplace.json"
    if not manifest.is_file():
        pytest.skip("not a source checkout (the sdist does not ship the marketplace)")
    marketplace = json.loads(manifest.read_text(encoding="utf-8"))
    assert marketplace["name"] == "crispr-analysis-utils"
    assert marketplace["owner"]["name"]
    (plugin,) = marketplace["plugins"]
    assert plugin["name"] == "crispr-analysis-utils"
    assert (REPO / plugin["source"]).resolve() == (
        REPO / "src" / "crispr_analysis_utils"
    )
    # No version: plugin users follow main, by commit, and need no release to update.
    assert "version" not in plugin
    assert not (REPO / plugin["source"] / ".claude-plugin").exists(), (
        "the marketplace entry is the plugin manifest"
    )

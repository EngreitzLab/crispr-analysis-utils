"""``cau install-skills``: install the bundled Claude Code skills and agents.

Skills are copied (or linked) to ``<target>/skills/<name>/`` and agents to
``<target>/agents/<name>.md``, where ``<target>`` is the Claude Code
configuration directory: ``$CLAUDE_CONFIG_DIR``, else ``~/.claude``.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from dataclasses import dataclass
from importlib import metadata, resources
from pathlib import Path

from . import CommandError

_COPY_IGNORE = shutil.ignore_patterns("__pycache__", ".ipynb_checkpoints", ".DS_Store")


@dataclass(frozen=True)
class Bundled:
    """A skill directory or an agent file shipped inside the package."""

    name: str
    kind: str  # "skill" or "agent"
    description: str
    source: Path
    preloads: tuple[str, ...] = ()  # the skills an agent loads at startup


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "names",
        nargs="*",
        metavar="NAME",
        help="Skills or agents to install (default: all of them). An agent also "
        "brings the skills it preloads.",
    )
    parser.add_argument(
        "--target",
        type=Path,
        help="Claude Code configuration directory to install into (default: "
        "$CLAUDE_CONFIG_DIR, else ~/.claude). Use .claude to install into the "
        "current project.",
    )
    parser.add_argument(
        "--symlink",
        action="store_true",
        help="Link to the installed package instead of copying. Follows edits to "
        "an editable install, and breaks when the environment is rebuilt.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace skills or agents of the same name already installed.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_only",
        help="List the bundled skills and agents, and exit.",
    )


def run(args: argparse.Namespace) -> int:
    bundled = bundled_items()
    if args.list_only:
        for item in bundled.values():
            print(f"{item.name}\t{item.kind}\t{item.description}")
        return 0

    selected = _select(bundled, args.names)
    target = _resolve_target(args.target)
    plan = [(item, _destination(target, item)) for item in selected]
    _check_destinations(plan, force=args.force)
    if args.symlink and not _is_editable_install():
        print(
            "cau: warning: crispr-analysis-utils is not an editable install, so the "
            "links point into this environment and break when it is rebuilt or "
            "removed. Copying (the default) avoids that.",
            file=sys.stderr,
        )

    for item, destination in plan:
        try:
            _install(item, destination, symlink=args.symlink)
        except OSError as error:
            raise CommandError(
                f"could not install {item.name} to {destination}: {error}"
            ) from error
        print(f"Installed {item.kind} {item.name}: {destination}")
    print(
        "New Claude Code sessions pick these up. Install either this way or as the "
        "Claude Code plugin, not both: with both, every skill is listed twice."
    )
    return 0


def bundled_items() -> dict[str, Bundled]:
    """The skills and agents shipped in the package, keyed by name."""
    root = _bundled_root()
    items: dict[str, Bundled] = {}
    skills_dir = root / "skills"
    if skills_dir.is_dir():
        for path in sorted(skills_dir.iterdir()):
            if (path / "SKILL.md").is_file():
                fields = read_frontmatter(path / "SKILL.md")
                items[path.name] = Bundled(
                    name=path.name,
                    kind="skill",
                    description=str(fields.get("description", "")),
                    source=path,
                )
    agents_dir = root / "agents"  # absent until the first agent exists
    if agents_dir.is_dir():
        for path in sorted(agents_dir.glob("*.md")):
            fields = read_frontmatter(path)
            preloads = fields.get("skills", [])
            if isinstance(preloads, str):
                preloads = [s.strip() for s in preloads.split(",") if s.strip()]
            items[path.stem] = Bundled(
                name=path.stem,
                kind="agent",
                description=str(fields.get("description", "")),
                source=path,
                preloads=tuple(preloads),
            )
    return items


def read_frontmatter(path: Path) -> dict[str, str | list[str]]:
    """Read the YAML frontmatter of a skill or agent file.

    Handles the shapes the bundled files use: ``key: value``, a folded
    ``key: >-`` block, and lists written ``key: [a, b]`` or as ``- item``
    lines. ``tests/test_skills.py`` checks it against PyYAML on every bundled
    file, so PyYAML is not a runtime dependency.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "---":
        raise CommandError(f"{path} does not start with a frontmatter block")
    try:
        end = next(i for i in range(1, len(lines)) if lines[i].strip() == "---")
    except StopIteration:
        raise CommandError(f"{path} has an unterminated frontmatter block") from None

    fields: dict[str, str | list[str]] = {}
    key = None
    for line in lines[1:end]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("- ") and isinstance(fields.get(key), list):
            fields[key].append(_unquote(stripped[2:]))
            continue
        if line[0] in " \t" and key is not None:
            value = fields[key]
            if not isinstance(value, list):
                fields[key] = f"{value} {stripped}".strip()
            continue
        key, _, raw = line.partition(":")
        key, raw = key.strip(), raw.strip()
        if raw in (">", ">-", "|", "|-"):
            fields[key] = ""
        elif raw == "":
            fields[key] = []
        elif raw.startswith("[") and raw.endswith("]"):
            fields[key] = [_unquote(v) for v in raw[1:-1].split(",") if v.strip()]
        else:
            fields[key] = _unquote(raw)
    return fields


def _unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
        return value[1:-1]
    return value


def _bundled_root() -> Path:
    root = resources.files("crispr_analysis_utils")
    if not isinstance(root, Path):
        raise CommandError(
            "crispr-analysis-utils is not installed as plain files, so its skills "
            "cannot be copied; reinstall it with pip, uv or pixi."
        )
    return root


def _select(bundled: dict[str, Bundled], names: list[str]) -> list[Bundled]:
    if not names:
        return list(bundled.values())
    unknown = [name for name in names if name not in bundled]
    if unknown:
        raise CommandError(
            f"unknown skill or agent: {', '.join(unknown)}. "
            f"Available: {', '.join(bundled)}."
        )
    selected: dict[str, Bundled] = {}
    for name in names:
        item = bundled[name]
        selected[name] = item
        for skill in item.preloads:
            if skill not in bundled:
                raise CommandError(
                    f"agent {name} preloads {skill}, which is not bundled"
                )
            selected.setdefault(skill, bundled[skill])
    return list(selected.values())


def _resolve_target(target: Path | None) -> Path:
    if target is None:
        target = Path(os.environ.get("CLAUDE_CONFIG_DIR") or "~/.claude")
    target = target.expanduser()
    if target.exists() and not target.is_dir():
        raise CommandError(f"--target {target} is not a directory")
    if target.name in ("skills", "agents"):
        raise CommandError(
            f"--target is the Claude Code configuration directory, not its "
            f"{target.name}/ folder: did you mean --target {target.parent}?"
        )
    return target


def _destination(target: Path, item: Bundled) -> Path:
    if item.kind == "skill":
        return target / "skills" / item.name
    return target / "agents" / f"{item.name}.md"


def _check_destinations(plan: list[tuple[Bundled, Path]], *, force: bool) -> None:
    """Refuse the whole install before anything is written."""
    existing = []
    for item, destination in plan:
        # Compare file identity, not path strings: on a case-insensitive
        # filesystem (macOS) or through a linked parent, a different-looking path
        # can still be the bundled copy, which --force would delete. An earlier
        # --symlink install is a link to the bundled copy and may be replaced.
        if (
            os.path.lexists(destination)
            and not destination.is_symlink()
            and destination.samefile(item.source)
        ):
            raise CommandError(
                f"{destination} is the bundled copy itself; choose another --target."
            )
        if os.path.lexists(destination):
            existing.append(str(destination))
    if existing and not force:
        raise CommandError(
            "already installed: "
            + ", ".join(existing)
            + ". Re-run with --force to replace them."
        )


def _install(item: Bundled, destination: Path, *, symlink: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink() or destination.is_file():
        destination.unlink()  # a link is removed, never followed
    elif destination.exists():
        shutil.rmtree(destination)
    if symlink:
        destination.symlink_to(item.source, target_is_directory=item.kind == "skill")
    elif item.kind == "skill":
        shutil.copytree(item.source, destination, ignore=_COPY_IGNORE)
    else:
        shutil.copy2(item.source, destination)


def _is_editable_install() -> bool:
    """Whether the package was installed in editable mode (PEP 610)."""
    try:
        direct_url = metadata.distribution("crispr-analysis-utils").read_text(
            "direct_url.json"
        )
    except metadata.PackageNotFoundError:
        return False
    if not direct_url:
        return False
    return bool(json.loads(direct_url).get("dir_info", {}).get("editable", False))

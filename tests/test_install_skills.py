"""`cau install-skills`: copy or link the bundled skills and agents.

Most tests install a small stand-in bundle (two skills, one agent) so they do
not depend on which skills the package ships; a few use the real bundle.
"""

import pytest

from crispr_analysis_utils import cli
from crispr_analysis_utils.cli import install_skills

SKILL = """---
name: {name}
description: >-
  Test skill {name}.
---

# {name}
"""

AGENT = """---
name: cau-alpha-runner
description: Test agent.
skills: [cau-alpha]
---

Run alpha.
"""


def run_cau(*args):
    """Run `cau install-skills ...`; return the exit status."""
    try:
        return cli.main(["install-skills", *args])
    except SystemExit as exit_info:
        return exit_info.code


@pytest.fixture
def bundle(tmp_path, monkeypatch):
    """A stand-in package root with skills cau-alpha and cau-beta, agent cau-alpha-runner."""
    root = tmp_path / "package"
    for name in ("cau-alpha", "cau-beta"):
        skill = root / "skills" / name
        (skill / "references").mkdir(parents=True)
        (skill / "SKILL.md").write_text(SKILL.format(name=name))
        (skill / "references" / "notes.md").write_text("notes\n")
        (skill / "__pycache__").mkdir()
        (skill / "__pycache__" / "junk.pyc").write_bytes(b"\0")
        (skill / ".ipynb_checkpoints").mkdir()
    (root / "agents").mkdir()
    (root / "agents" / "cau-alpha-runner.md").write_text(AGENT)
    monkeypatch.setattr(install_skills, "_bundled_root", lambda: root)
    return root


@pytest.fixture
def target(tmp_path):
    return tmp_path / "claude"


def test_installs_every_skill_and_agent(bundle, target):
    assert run_cau("--target", str(target)) == 0
    for name in ("cau-alpha", "cau-beta"):
        assert (target / "skills" / name / "SKILL.md").is_file()
        assert (target / "skills" / name / "references" / "notes.md").is_file()
    assert (target / "agents" / "cau-alpha-runner.md").is_file()


def test_copy_leaves_out_caches_and_autosaves(bundle, target):
    run_cau("--target", str(target))
    installed = target / "skills" / "cau-alpha"
    assert not (installed / "__pycache__").exists()
    assert not (installed / ".ipynb_checkpoints").exists()


def test_installs_the_real_bundle(target):
    assert run_cau("--target", str(target)) == 0
    bundled = install_skills.bundled_items()
    assert bundled, "the package ships no skills"
    for item in bundled.values():
        if item.kind == "skill":
            assert (target / "skills" / item.name / "SKILL.md").is_file()


def test_subset_by_name(bundle, target):
    assert run_cau("cau-beta", "--target", str(target)) == 0
    assert (target / "skills" / "cau-beta").is_dir()
    assert not (target / "skills" / "cau-alpha").exists()
    assert not (target / "agents").exists()


def test_an_agent_brings_the_skills_it_preloads(bundle, target):
    assert run_cau("cau-alpha-runner", "--target", str(target)) == 0
    assert (target / "agents" / "cau-alpha-runner.md").is_file()
    assert (target / "skills" / "cau-alpha" / "SKILL.md").is_file()
    assert not (target / "skills" / "cau-beta").exists()


def test_no_agents_directory_means_no_agents(bundle, target):
    (bundle / "agents" / "cau-alpha-runner.md").unlink()
    (bundle / "agents").rmdir()
    assert run_cau("--target", str(target)) == 0
    assert (target / "skills" / "cau-alpha").is_dir()
    assert not (target / "agents").exists()


def test_unknown_name_is_an_error(bundle, target, capsys):
    assert run_cau("cau-gamma", "--target", str(target)) == 1
    err = capsys.readouterr().err
    assert "unknown skill or agent: cau-gamma" in err
    assert "cau-alpha" in err


def test_list(bundle, capsys):
    assert run_cau("--list") == 0
    lines = capsys.readouterr().out.splitlines()
    assert "cau-alpha\tskill\tTest skill cau-alpha." in lines
    assert "cau-alpha-runner\tagent\tTest agent." in lines


def test_collision_without_force_writes_nothing(bundle, target):
    existing = target / "skills" / "cau-beta"
    existing.mkdir(parents=True)
    assert run_cau("--target", str(target)) == 1
    assert not (target / "skills" / "cau-alpha").exists()
    assert not (target / "agents").exists()


def test_force_replaces_an_existing_install(bundle, target):
    run_cau("--target", str(target))
    stray = target / "skills" / "cau-alpha" / "stray.txt"
    stray.write_text("x")
    assert run_cau("--target", str(target), "--force") == 0
    assert not stray.exists()
    assert (target / "skills" / "cau-alpha" / "SKILL.md").is_file()


def test_symlink_install_points_at_the_bundle(bundle, target):
    assert run_cau("--target", str(target), "--symlink") == 0
    link = target / "skills" / "cau-alpha"
    assert link.is_symlink()
    assert link.resolve() == (bundle / "skills" / "cau-alpha").resolve()
    assert (target / "agents" / "cau-alpha-runner.md").is_symlink()


def test_symlink_install_can_be_redone_with_force(bundle, target):
    run_cau("--target", str(target), "--symlink")
    assert run_cau("--target", str(target), "--symlink", "--force") == 0
    assert (target / "skills" / "cau-alpha").is_symlink()


def test_copy_over_a_symlink_install_leaves_the_bundle_intact(bundle, target):
    run_cau("--target", str(target), "--symlink")
    assert run_cau("--target", str(target), "--force") == 0
    installed = target / "skills" / "cau-alpha"
    assert installed.is_dir() and not installed.is_symlink()
    assert (bundle / "skills" / "cau-alpha" / "SKILL.md").is_file()


def test_refuses_to_install_over_the_bundle_itself(bundle):
    assert run_cau("--target", str(bundle), "--force") == 1
    assert (bundle / "skills" / "cau-alpha" / "SKILL.md").is_file()
    assert (bundle / "agents" / "cau-alpha-runner.md").is_file()


@pytest.mark.parametrize("folder", ["skills", "agents"])
def test_target_must_be_the_configuration_directory(bundle, target, folder, capsys):
    assert run_cau("--target", str(target / folder)) == 1
    assert f"did you mean --target {target}" in capsys.readouterr().err


def test_default_target_honours_claude_config_dir(bundle, tmp_path, monkeypatch):
    config = tmp_path / "config"
    monkeypatch.setenv("CLAUDE_CONFIG_DIR", str(config))
    assert run_cau() == 0
    assert (config / "skills" / "cau-alpha" / "SKILL.md").is_file()


def test_symlink_warns_when_not_an_editable_install(
    bundle, target, monkeypatch, capsys
):
    monkeypatch.setattr(install_skills, "_is_editable_install", lambda: False)
    assert run_cau("--target", str(target), "--symlink") == 0
    assert "not an editable install" in capsys.readouterr().err


def test_refuses_a_case_variant_of_the_bundle(bundle):
    """On a case-insensitive filesystem (the macOS default) PACKAGE is package."""
    variant = bundle.with_name(bundle.name.upper())
    if not variant.exists():
        pytest.skip("case-sensitive filesystem")
    assert run_cau("--target", str(variant), "--force") == 1
    assert (bundle / "skills" / "cau-alpha" / "SKILL.md").is_file()


def test_refuses_a_target_whose_skills_folder_links_to_the_bundle(bundle, target):
    target.mkdir()
    (target / "skills").symlink_to(bundle / "skills", target_is_directory=True)
    assert run_cau("--target", str(target), "--force") == 1
    assert (bundle / "skills" / "cau-alpha" / "SKILL.md").is_file()


def test_a_dangling_link_counts_as_installed(bundle, target, tmp_path):
    """What --symlink leaves behind once its environment is rebuilt."""
    (target / "skills").mkdir(parents=True)
    link = target / "skills" / "cau-alpha"
    link.symlink_to(tmp_path / "gone", target_is_directory=True)
    assert run_cau("--target", str(target)) == 1
    assert link.is_symlink() and not link.exists()
    assert not (target / "skills" / "cau-beta").exists()
    assert not (target / "agents").exists()


def test_target_that_is_a_file_is_an_error(bundle, tmp_path, capsys):
    not_a_dir = tmp_path / "file"
    not_a_dir.write_text("x")
    assert run_cau("--target", str(not_a_dir)) == 1
    assert "is not a directory" in capsys.readouterr().err


def test_filesystem_errors_are_reported_without_traceback(bundle, target, capsys):
    agents = target / "agents"
    agents.mkdir(parents=True)
    agents.chmod(0o500)
    try:
        assert run_cau("--target", str(target)) == 1
    finally:
        agents.chmod(0o700)
    assert "cau: error: could not install cau-alpha-runner" in capsys.readouterr().err


def test_reads_a_block_list_written_at_column_zero(bundle, target):
    agent = bundle / "agents" / "cau-alpha-runner.md"
    agent.write_text(AGENT.replace("skills: [cau-alpha]", "skills:\n- cau-alpha"))
    assert install_skills.bundled_items()["cau-alpha-runner"].preloads == ("cau-alpha",)
    assert run_cau("cau-alpha-runner", "--target", str(target)) == 0
    assert (target / "skills" / "cau-alpha" / "SKILL.md").is_file()

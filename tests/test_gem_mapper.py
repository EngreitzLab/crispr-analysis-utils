import crispr_analysis_utils.gem_mapper as gem_mapper


def test_build_gem_index_command(monkeypatch):
    captured = {}

    def fake_run_shell_cmd(cmd):
        captured["cmd"] = cmd
        return ""

    monkeypatch.setattr(gem_mapper, "run_shell_cmd", fake_run_shell_cmd)
    gem_mapper.build_gem_index("ref.fa", "idx/ref", threads=8)

    assert captured["cmd"] == "gem-indexer -i ref.fa -o idx/ref -t 8"


def test_map_guides_with_gem_command(monkeypatch):
    captured = {}

    def fake_run_shell_cmd(cmd):
        captured["cmd"] = cmd
        return ""

    monkeypatch.setattr(gem_mapper, "run_shell_cmd", fake_run_shell_cmd)
    gem_mapper.map_guides_with_gem(
        "idx/ref.gem",
        "input.fastq",
        "out.sam",
        mapping_mode="sensitive",
        threads=8,
        sam_compact=False,
    )

    assert captured["cmd"] == (
        "gem-mapper -I idx/ref.gem -i input.fastq -o out.sam "
        "--sam-compact=false --mapping-mode sensitive --threads 8 "
        "> out.log 2>&1"
    )


def test_map_guides_with_gem_command_custom_log_path(monkeypatch):
    captured = {}

    def fake_run_shell_cmd(cmd):
        captured["cmd"] = cmd
        return ""

    monkeypatch.setattr(gem_mapper, "run_shell_cmd", fake_run_shell_cmd)
    gem_mapper.map_guides_with_gem(
        "idx/ref.gem",
        "input.fastq",
        "out.sam",
        log_path="logs/map.log",
    )

    assert captured["cmd"].endswith("> logs/map.log 2>&1")

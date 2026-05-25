import crispr_analysis_utils as cau


def test_run_shell_cmd_returns_stdout():
    stdout = cau.utils.run_shell_cmd("printf 'hello'")
    assert stdout == "hello"

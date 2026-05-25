"""General utility helpers."""

from __future__ import annotations

import logging
import os
import signal
import subprocess
import time


def run_shell_cmd(cmd: str) -> str:
    """Run a shell command through bash with pipefail enabled.

    Parameters
    ----------
    cmd
        Shell command string.

    Returns
    -------
    str
        Command stdout with trailing newlines stripped.
    """
    proc = subprocess.Popen(
        ["/bin/bash", "-o", "pipefail"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid,
    )
    pid = proc.pid
    pgid = os.getpgid(pid)
    logging.info("run_shell_cmd: PID=%s, PGID=%s, CMD=%s", pid, pgid, cmd)

    started = time.perf_counter()
    stdout, stderr = proc.communicate(cmd)
    elapsed_s = time.perf_counter() - started
    rc = proc.returncode

    err_str = (
        f"PID={pid}, PGID={pgid}, RC={rc}, ELAPSED_S={elapsed_s:.3f}\n"
        f"STDERR={stderr.strip()}\n"
        f"STDOUT={stdout.strip()}"
    )
    if rc:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except Exception:
            pass
        raise RuntimeError(err_str)

    logging.info(err_str)
    return stdout.strip("\n")

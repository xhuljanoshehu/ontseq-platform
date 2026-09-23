from __future__ import annotations

import os
import shlex
import shutil
import signal
import sys
import time
from pathlib import Path

import pytest

from ontseq_platform.execution import SubprocessRunner, ToolExecutionError


def _linux_pid_is_live(pid: int) -> bool:
    """Return False for absent or zombie Linux processes."""
    try:
        fields = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8").split()
    except FileNotFoundError:
        return False
    return len(fields) > 2 and fields[2] != "Z"


def _wait_until_not_live(pid: int, *, timeout_seconds: float = 2.0) -> bool:
    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        if not _linux_pid_is_live(pid):
            return True
        time.sleep(0.02)
    return not _linux_pid_is_live(pid)


def _descendant_script(pid_path: Path, *, emit_stdout: bool = False) -> str:
    prefix = "printf 'partial-output\\n'; " if emit_stdout else ""
    return (
        prefix
        + "sleep 30 & child=$!; "
        + f"printf '%s' \"$child\" > {shlex.quote(str(pid_path))}; "
        + 'wait "$child"'
    )


def _kill_if_still_live(pid: int) -> None:
    if _linux_pid_is_live(pid):
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-group cleanup contract")
def test_timeout_kills_descendant_process_tree(tmp_path: Path) -> None:
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "descendant.pid"

    with pytest.raises(ToolExecutionError, match="Command timed out after 1 seconds"):
        SubprocessRunner().run(
            [shell, "-c", _descendant_script(pid_path)], timeout_seconds=1
        )

    assert pid_path.is_file()
    pid = int(pid_path.read_text(encoding="utf-8"))
    try:
        assert _wait_until_not_live(pid), "timed-out tool left a live descendant process"
    finally:
        _kill_if_still_live(pid)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-group cleanup contract")
def test_streaming_timeout_kills_descendants_and_removes_staged_output(tmp_path: Path) -> None:
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "stream-descendant.pid"
    output = tmp_path / "stream.bin"

    with pytest.raises(ToolExecutionError, match="Command timed out after 1 seconds"):
        SubprocessRunner().run_to_file(
            [shell, "-c", _descendant_script(pid_path, emit_stdout=True)],
            output,
            timeout_seconds=1,
        )

    assert pid_path.is_file()
    pid = int(pid_path.read_text(encoding="utf-8"))
    try:
        assert _wait_until_not_live(pid), "streaming timeout left a live descendant process"
    finally:
        _kill_if_still_live(pid)
    assert not output.exists()
    assert list(tmp_path.glob(f".{output.name}.*.tmp")) == []


def test_captured_command_result_semantics_are_preserved() -> None:
    result = SubprocessRunner().run(
        [
            sys.executable,
            "-c",
            "import sys; print('ok'); print('diagnostic', file=sys.stderr); raise SystemExit(7)",
        ],
        timeout_seconds=30,
    )

    assert result.returncode == 7
    assert result.stdout == "ok\n"
    assert result.stderr == "diagnostic\n"

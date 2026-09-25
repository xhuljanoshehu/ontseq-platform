from __future__ import annotations

import contextlib
import os
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

from ontseq_platform import execution
from ontseq_platform.execution import SubprocessRunner, ToolExecutionError


def _linux_pid_is_live(pid: int) -> bool:
    """Return False for absent or zombie Linux processes.

    A process that exits between opening and reading ``/proc/<pid>/stat`` makes the read
    fail with ESRCH (``ProcessLookupError``) rather than ``FileNotFoundError``; both mean
    the process is gone.
    """
    try:
        fields = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8").split()
    except (FileNotFoundError, ProcessLookupError):
        return False
    return len(fields) > 2 and fields[2] != "Z"


def _linux_state(pid: int) -> str | None:
    """Return the ``/proc/<pid>/stat`` state letter, or None once the process is gone."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
    except (FileNotFoundError, ProcessLookupError):
        return None
    return stat.rpartition(")")[2].split()[0]


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
        with contextlib.suppress(ProcessLookupError):
            os.kill(pid, signal.SIGKILL)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
def test_timeout_kills_descendant_process_tree(tmp_path: Path) -> None:
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "descendant.pid"

    with pytest.raises(ToolExecutionError, match="Command timed out after 1 seconds"):
        SubprocessRunner().run([shell, "-c", _descendant_script(pid_path)], timeout_seconds=1)

    assert pid_path.is_file()
    pid = int(pid_path.read_text(encoding="utf-8"))
    try:
        assert _wait_until_not_live(pid), "timed-out tool left a live descendant process"
    finally:
        _kill_if_still_live(pid)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
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
            "-S",
            "-c",
            "import sys; print('ok'); print('diagnostic', file=sys.stderr); raise SystemExit(7)",
        ],
        timeout_seconds=30,
    )

    assert result.returncode == 7
    assert result.stdout == "ok\n"
    assert result.stderr == "diagnostic\n"


def _interrupt_once_descendant_starts(monkeypatch: pytest.MonkeyPatch, pid_path: Path) -> None:
    """Make the wait raise ``KeyboardInterrupt`` once the tool's descendant is running.

    This is what Ctrl+C does to the parent. The terminal's SIGINT reaches the tool's
    process group as well, but a tool may ignore it (a non-interactive shell's background
    jobs do), so the runner itself must stop the tree.
    """

    def interrupted(self: subprocess.Popen[str], *args: object, **kwargs: object) -> None:
        deadline = time.monotonic() + 10
        while not (pid_path.is_file() and pid_path.read_text(encoding="utf-8")):
            if time.monotonic() > deadline:
                raise AssertionError("the synthetic tool never started its descendant")
            time.sleep(0.02)
        raise KeyboardInterrupt

    monkeypatch.setattr(subprocess.Popen, "communicate", interrupted)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
@pytest.mark.parametrize("streaming", [False, True])
def test_interrupt_kills_descendant_process_tree(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, streaming: bool
) -> None:
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "interrupted-descendant.pid"
    output = tmp_path / "interrupted.bin"
    _interrupt_once_descendant_starts(monkeypatch, pid_path)
    argv = [shell, "-c", _descendant_script(pid_path, emit_stdout=streaming)]

    with pytest.raises(KeyboardInterrupt):
        if streaming:
            SubprocessRunner().run_to_file(argv, output, timeout_seconds=60)
        else:
            SubprocessRunner().run(argv, timeout_seconds=60)

    pid = int(pid_path.read_text(encoding="utf-8"))
    try:
        assert _wait_until_not_live(pid), "an interrupted run left a live descendant process"
    finally:
        _kill_if_still_live(pid)
    assert not output.exists()
    assert list(tmp_path.glob(f".{output.name}.*.tmp")) == []


def _record_started_processes(monkeypatch: pytest.MonkeyPatch) -> list[subprocess.Popen[str]]:
    started: list[subprocess.Popen[str]] = []
    real_popen = subprocess.Popen

    def recording_popen(*args: object, **kwargs: object) -> subprocess.Popen[str]:
        process: subprocess.Popen[str] = real_popen(*args, **kwargs)  # type: ignore[call-overload]
        started.append(process)
        return process

    monkeypatch.setattr(subprocess, "Popen", recording_popen)
    return started


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
@pytest.mark.parametrize("streaming", [False, True])
def test_aborted_runs_close_their_pipes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, streaming: bool
) -> None:
    """A retained traceback must not keep the aborted process's pipe descriptors open."""
    started = _record_started_processes(monkeypatch)
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "pipes-descendant.pid"
    argv = [shell, "-c", _descendant_script(pid_path, emit_stdout=streaming)]

    with pytest.raises(ToolExecutionError, match="timed out"):
        if streaming:
            SubprocessRunner().run_to_file(argv, tmp_path / "pipes.bin", timeout_seconds=1)
        else:
            SubprocessRunner().run(argv, timeout_seconds=1)

    _kill_if_still_live(int(pid_path.read_text(encoding="utf-8")))
    (process,) = started
    pipes = [stream for stream in (process.stdout, process.stderr) if stream is not None]
    assert pipes
    assert all(stream.closed for stream in pipes)


_SRC = Path(__file__).resolve().parents[1] / "src"
_PIPELINE = (
    "import sys\n"
    "sys.path.insert(0, sys.argv[1])\n"
    "from ontseq_platform.execution import SubprocessRunner\n"
    "SubprocessRunner().run(['sh', '-c', sys.argv[2]], timeout_seconds=60)\n"
)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-group contract")
@pytest.mark.parametrize("signum", [signal.SIGTERM, signal.SIGHUP])
def test_signal_to_the_pipeline_process_group_reaches_the_tool(
    tmp_path: Path, signum: signal.Signals
) -> None:
    """A service manager or WSL teardown signals the pipeline's group; the tool must go too.

    The pipeline runs in its own session here only to isolate the test runner.
    """
    pid_path = tmp_path / "group-descendant.pid"
    pipeline = subprocess.Popen(
        [sys.executable, "-c", _PIPELINE, str(_SRC), _descendant_script(pid_path)],
        start_new_session=True,
    )
    try:
        deadline = time.monotonic() + 20
        while not (pid_path.is_file() and pid_path.read_text(encoding="utf-8")):
            assert time.monotonic() < deadline, "the synthetic tool never started"
            time.sleep(0.02)
        pid = int(pid_path.read_text(encoding="utf-8"))
        os.killpg(pipeline.pid, signum)
        pipeline.wait(timeout=10)
        try:
            assert _wait_until_not_live(pid), "a group signal left the tool's process alive"
        finally:
            _kill_if_still_live(pid)
    finally:
        with contextlib.suppress(ProcessLookupError):
            os.killpg(pipeline.pid, signal.SIGKILL)
        pipeline.wait(timeout=10)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree contract")
def test_timeout_cleanup_leaves_unrelated_processes_alone(tmp_path: Path) -> None:
    bystander = subprocess.Popen(["sleep", "30"])
    try:
        shell = shutil.which("sh")
        assert shell is not None
        pid_path = tmp_path / "bystander-descendant.pid"
        with pytest.raises(ToolExecutionError, match="timed out"):
            SubprocessRunner().run([shell, "-c", _descendant_script(pid_path)], timeout_seconds=1)
        _kill_if_still_live(int(pid_path.read_text(encoding="utf-8")))
        assert _linux_pid_is_live(bystander.pid), "timeout cleanup killed an unrelated process"
    finally:
        bystander.kill()
        bystander.wait()


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
def test_a_reused_pid_found_by_the_scan_is_resumed_not_killed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A scanned PID that now belongs to an unrelated process is resumed and left alone."""
    real_scan = execution._linux_descendants
    bystander = subprocess.Popen(["sleep", "30"])
    try:
        monkeypatch.setattr(
            execution, "_linux_descendants", lambda root: [*real_scan(root), bystander.pid]
        )
        shell = shutil.which("sh")
        assert shell is not None
        pid_path = tmp_path / "reused-pid-descendant.pid"
        with pytest.raises(ToolExecutionError, match="timed out"):
            SubprocessRunner().run([shell, "-c", _descendant_script(pid_path)], timeout_seconds=1)
        pid = int(pid_path.read_text(encoding="utf-8"))
        try:
            assert _wait_until_not_live(pid), "the tool's own descendant was not killed"
        finally:
            _kill_if_still_live(pid)
        assert _linux_state(bystander.pid) not in (None, "Z", "T", "t"), (
            "timeout cleanup stopped or killed a process outside the tool's tree"
        )
    finally:
        bystander.kill()
        bystander.wait()


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
def test_interrupted_cleanup_still_kills_what_it_stopped(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A second Ctrl+C during the descendant scan must not leave stopped processes behind."""
    real_scan = execution._linux_descendants
    scans = 0

    def scan_then_interrupt(root: int) -> list[int]:
        nonlocal scans
        scans += 1
        if scans > 1:
            raise KeyboardInterrupt
        return real_scan(root)

    monkeypatch.setattr(execution, "_linux_descendants", scan_then_interrupt)
    started = _record_started_processes(monkeypatch)
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "rescan-descendant.pid"

    with pytest.raises(KeyboardInterrupt):
        SubprocessRunner().run([shell, "-c", _descendant_script(pid_path)], timeout_seconds=1)

    (process,) = started
    pid = int(pid_path.read_text(encoding="utf-8"))
    try:
        assert scans >= 2
        assert _wait_until_not_live(pid), "an interrupted cleanup left a stopped descendant"
        assert process.returncode == -signal.SIGKILL, "the stopped tool was not killed and reaped"
    finally:
        _kill_if_still_live(pid)
        if process.returncode is None:
            process.kill()
            process.wait()


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
def test_interrupt_right_after_a_descendant_is_stopped_still_kills_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The window between stopping a descendant and classifying it must not leak it."""
    started = _record_started_processes(monkeypatch)
    real_signal = execution._signal
    interrupted = False

    def stop_then_interrupt(pid: int, signum: int) -> None:
        nonlocal interrupted
        real_signal(pid, signum)
        if signum == signal.SIGSTOP and pid != started[0].pid and not interrupted:
            interrupted = True
            raise KeyboardInterrupt

    monkeypatch.setattr(execution, "_signal", stop_then_interrupt)
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "stop-window-descendant.pid"

    with pytest.raises(KeyboardInterrupt):
        SubprocessRunner().run([shell, "-c", _descendant_script(pid_path)], timeout_seconds=1)

    (process,) = started
    pid = int(pid_path.read_text(encoding="utf-8"))
    try:
        assert interrupted
        assert _wait_until_not_live(pid), "a descendant stopped just before an interrupt leaked"
        assert process.returncode == -signal.SIGKILL
    finally:
        _kill_if_still_live(pid)
        if process.returncode is None:
            process.kill()
            process.wait()


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
@pytest.mark.parametrize("receiver", ["runner-thread", "other-thread"])
def test_ctrl_c_during_cleanup_is_deferred_until_the_tree_is_dead(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, receiver: str
) -> None:
    """A further Ctrl+C must not cut the cleanup short; it is raised once the tree is dead.

    Python raises ``KeyboardInterrupt`` in the main thread even when the kernel delivered
    SIGINT to another thread, so both deliveries are exercised.
    """
    idle = threading.Event()
    other = threading.Thread(target=idle.wait, daemon=True)
    other.start()
    target = threading.get_ident() if receiver == "runner-thread" else other.ident
    assert target is not None
    real_scan = execution._linux_descendants
    scans = 0

    def scan_after_ctrl_c(root: int) -> list[int]:
        nonlocal scans
        scans += 1
        if scans == 1:
            signal.pthread_kill(target, signal.SIGINT)
            time.sleep(0.2)  # let the receiving thread trip Python's signal flag
        return real_scan(root)

    monkeypatch.setattr(execution, "_linux_descendants", scan_after_ctrl_c)
    started = _record_started_processes(monkeypatch)
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "deferred-ctrl-c-descendant.pid"

    try:
        with pytest.raises(KeyboardInterrupt):
            SubprocessRunner().run([shell, "-c", _descendant_script(pid_path)], timeout_seconds=1)
    finally:
        idle.set()
        other.join()

    (process,) = started
    pid = int(pid_path.read_text(encoding="utf-8"))
    try:
        assert scans >= 2, "the Ctrl+C cut the descendant scan short"
        assert _wait_until_not_live(pid), "a Ctrl+C during cleanup left a descendant behind"
        assert process.returncode == -signal.SIGKILL
    finally:
        _kill_if_still_live(pid)
        if process.returncode is None:
            process.kill()
            process.wait()


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
def test_reaping_the_killed_tool_stays_interruptible(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Only the signalling defers Ctrl+C, not the wait for a tool stuck in I/O."""
    handler_before = signal.getsignal(signal.SIGINT)
    handlers: list[object] = []
    real_wait = subprocess.Popen.wait

    def recording_wait(self: subprocess.Popen[str], timeout: float | None = None) -> int:
        handlers.append(signal.getsignal(signal.SIGINT))
        return real_wait(self, timeout)

    monkeypatch.setattr(subprocess.Popen, "wait", recording_wait)
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "reaped-descendant.pid"

    with pytest.raises(ToolExecutionError, match="timed out"):
        SubprocessRunner().run([shell, "-c", _descendant_script(pid_path)], timeout_seconds=1)

    _kill_if_still_live(int(pid_path.read_text(encoding="utf-8")))
    assert handlers
    assert all(handler is handler_before for handler in handlers)
    assert signal.getsignal(signal.SIGINT) is handler_before


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
def test_pipes_are_closed_even_if_reaping_is_interrupted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ctrl+C during the interruptible wait must not leak the aborted process's pipes."""
    popen = subprocess.Popen
    real_wait = popen.wait

    def interrupted_wait(self: subprocess.Popen[str], timeout: float | None = None) -> int:
        raise KeyboardInterrupt

    monkeypatch.setattr(popen, "wait", interrupted_wait)
    started = _record_started_processes(monkeypatch)
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "interrupted-reap-descendant.pid"

    with pytest.raises(KeyboardInterrupt):
        SubprocessRunner().run([shell, "-c", _descendant_script(pid_path)], timeout_seconds=1)

    (process,) = started
    try:
        pipes = [stream for stream in (process.stdout, process.stderr) if stream is not None]
        assert pipes
        assert all(stream.closed for stream in pipes)
    finally:
        _kill_if_still_live(int(pid_path.read_text(encoding="utf-8")))
        real_wait(process)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
def test_timeout_cleanup_works_off_the_main_thread(tmp_path: Path) -> None:
    """Worker threads cannot swap signal handlers; their cleanup must still kill the tree."""
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "worker-descendant.pid"
    argv = [shell, "-c", _descendant_script(pid_path)]
    errors: list[BaseException] = []

    def run_in_worker() -> None:
        try:
            SubprocessRunner().run(argv, timeout_seconds=1)
        except BaseException as exc:
            errors.append(exc)

    worker = threading.Thread(target=run_in_worker)
    worker.start()
    worker.join(timeout=30)
    assert not worker.is_alive()
    pid = int(pid_path.read_text(encoding="utf-8"))
    try:
        assert len(errors) == 1
        assert isinstance(errors[0], ToolExecutionError), errors[0]
        assert _wait_until_not_live(pid), "a worker-thread timeout left a descendant alive"
    finally:
        _kill_if_still_live(pid)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux/WSL process-tree cleanup contract")
def test_ctrl_c_before_the_deferral_takes_effect_still_kills_the_tree(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A Ctrl+C that lands before the deferring handler is installed must not leak the tree."""
    real_getsignal = signal.getsignal
    lookups = 0

    def getsignal_interrupted_once(signalnum: int) -> object:
        nonlocal lookups
        lookups += 1
        if lookups == 1:
            raise KeyboardInterrupt
        return real_getsignal(signalnum)

    monkeypatch.setattr(signal, "getsignal", getsignal_interrupted_once)
    started = _record_started_processes(monkeypatch)
    shell = shutil.which("sh")
    assert shell is not None
    pid_path = tmp_path / "early-ctrl-c-descendant.pid"

    with pytest.raises(KeyboardInterrupt):
        SubprocessRunner().run([shell, "-c", _descendant_script(pid_path)], timeout_seconds=1)

    (process,) = started
    pid = int(pid_path.read_text(encoding="utf-8"))
    try:
        assert lookups >= 2, "the cleanup was not retried under the deferral"
        assert _wait_until_not_live(pid), "an early Ctrl+C left a descendant behind"
        assert process.returncode == -signal.SIGKILL
    finally:
        _kill_if_still_live(pid)
        if process.returncode is None:
            process.kill()
            process.wait()


def test_cleanup_never_signals_an_already_reaped_tool(monkeypatch: pytest.MonkeyPatch) -> None:
    """Once reaped, the tool's PID may belong to another process; it must not be frozen."""
    frozen_roots: list[int] = []
    monkeypatch.setattr(execution, "_kill_descendants", frozen_roots.append)
    real_communicate = subprocess.Popen.communicate

    def finish_then_interrupt(
        self: subprocess.Popen[str], input: str | None = None, timeout: float | None = None
    ) -> tuple[str, str]:
        real_communicate(self, input, timeout)
        raise KeyboardInterrupt

    monkeypatch.setattr(subprocess.Popen, "communicate", finish_then_interrupt)

    with pytest.raises(KeyboardInterrupt):
        SubprocessRunner().run([sys.executable, "-S", "-c", "pass"], timeout_seconds=30)

    assert frozen_roots == []

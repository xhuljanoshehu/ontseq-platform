from __future__ import annotations

import contextlib
import os
import signal
import subprocess
import tempfile
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol


@dataclass(frozen=True)
class CommandResult:
    argv: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


class CommandRunner(Protocol):
    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult: ...


class StreamingCommandRunner(CommandRunner, Protocol):
    """A runner that can send a command's stdout straight to a file.

    Some tools emit binary output on stdout rather than accepting an output path. Reading
    that through the text-mode capture used elsewhere would corrupt it, so those adapters
    depend on this narrower protocol instead.
    """

    def run_to_file(
        self, argv: Sequence[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult: ...


class ToolExecutionError(RuntimeError):
    """Raised when a required local executable cannot be started or times out."""


def _normalize(argv: Sequence[str]) -> tuple[str, ...]:
    normalized = tuple(str(item) for item in argv)
    if not normalized or any(not item or "\x00" in item for item in normalized):
        raise ValueError("Command arguments must be non-empty and cannot contain NUL bytes")
    return normalized


_PROC = Path("/proc")
#: Rounds of "find the tool's descendants, stop them" before the frozen tree is killed.
#: A stopped process cannot fork, so the tree stops growing after a few rounds.
_FREEZE_ROUNDS = 8


def _linux_parent_pid(pid: int) -> int | None:
    """Return the parent of a live Linux process, or None if it is gone or a zombie."""
    try:
        stat = (_PROC / str(pid) / "stat").read_text(encoding="utf-8", errors="replace")
    except OSError:  # ENOENT, or ESRCH when the process exits mid-read
        return None
    # The command name is parenthesised and may itself contain spaces or parentheses.
    fields = stat.rpartition(")")[2].split()
    if len(fields) < 2 or fields[0] == "Z":
        return None
    try:
        return int(fields[1])
    except ValueError:
        return None


def _linux_descendants(root: int) -> list[int]:
    """Live descendants of ``root`` found through ``/proc``, parents before children."""
    try:
        entries = [int(entry.name) for entry in _PROC.iterdir() if entry.name.isdigit()]
    except OSError:
        return []
    children: dict[int, list[int]] = {}
    for pid in entries:
        parent = _linux_parent_pid(pid)
        if parent is not None:
            children.setdefault(parent, []).append(pid)
    ordered: list[int] = []
    pending = [root]
    while pending:
        for child in children.get(pending.pop(0), ()):
            ordered.append(child)
            pending.append(child)
    return ordered


def _signal(pid: int, signum: int) -> None:
    with contextlib.suppress(ProcessLookupError, PermissionError):
        os.kill(pid, signum)


def _kill_descendants(root: int) -> None:
    """Freeze and kill the processes the tool started, without touching anything else.

    The tool stays in the pipeline's process group, so a signal to that group (a service
    manager, ``kill -- -PGID``, the WSL session teardown behind the Desktop app, a terminal
    hangup) still reaches the tool and its descendants, as it did with ``subprocess.run``.
    Its descendants are therefore found through ``/proc`` instead of a dedicated process
    group. The root is stopped first so it cannot fork more children, each descendant is
    stopped as it is found, and one whose parent is no longer in the frozen tree (a reused
    PID) is resumed and left alone. Only the frozen tree is killed. A descendant that
    outlived its own parent has already been re-parented away and is not found. Where
    ``/proc`` is unavailable this finds nothing and only the direct child is killed, as
    before.
    """
    _signal(root, signal.SIGSTOP)
    members = {root}
    # Every candidate is recorded before it is stopped (``None``: not checked yet), so an
    # exception at any point still leaves no stopped process behind: the ``finally`` block
    # checks what is unchecked, then kills the members children first.
    candidates: dict[int, bool | None] = {}
    try:
        for _ in range(_FREEZE_ROUNDS):
            found = [pid for pid in _linux_descendants(root) if pid not in members]
            if not found:
                break
            for pid in found:
                candidates[pid] = None
                candidates[pid] = _freeze_if_member(pid, members)
    finally:
        for pid, member in reversed(candidates.items()):
            if member is None:
                member = _freeze_if_member(pid, members)
            if member:
                _signal(pid, signal.SIGKILL)


def _freeze_if_member(pid: int, members: set[int]) -> bool:
    """Stop ``pid``; keep it stopped only if its parent belongs to the frozen tree."""
    _signal(pid, signal.SIGSTOP)
    if _linux_parent_pid(pid) in members:
        members.add(pid)
        return True
    _signal(pid, signal.SIGCONT)
    return False


@contextlib.contextmanager
def _sigint_held() -> Iterator[None]:
    """Hold SIGINT for the calling thread; a Ctrl+C meanwhile is raised afterwards.

    A cleanup that a further Ctrl+C interrupts half-way could leave stopped processes
    behind, so that Ctrl+C stays pending until every stopped process has been killed or
    resumed and then raises ``KeyboardInterrupt`` as usual. Only the calling thread is
    masked: in a multi-threaded process the kernel may deliver the signal to another
    thread, which is why ``_kill_descendants`` also records every process before stopping
    it.
    """
    if not hasattr(signal, "pthread_sigmask"):  # native Windows
        yield
        return
    previous = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGINT})
    try:
        yield
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous)


def _terminate_process_tree(process: subprocess.Popen[bytes] | subprocess.Popen[str]) -> None:
    """Kill the tool and every process it started, then reap it and release its pipes.

    Called on a timeout and on any other abort while waiting (``KeyboardInterrupt``,
    ``SystemExit``). ``subprocess.run`` killed only the direct child, so a timed-out
    multiprocessing tool could leave CPU/RAM-consuming workers behind. Native Windows keeps
    direct-child termination until a separately tested job-object contract exists.
    """
    try:
        # Only the signalling is held off from Ctrl+C. The wait below stays interruptible
        # in case the killed tool is stuck in uninterruptible I/O.
        with _sigint_held():
            try:
                if os.name == "posix":
                    _kill_descendants(process.pid)
            finally:
                # On POSIX the tool is already stopped here, so it must be killed even if
                # the descendant scan was itself interrupted.
                process.kill()
    finally:
        process.wait()
        # The Popen object can outlive this call inside a retained traceback, so release
        # its pipes now rather than at garbage collection. Closing (not draining) cannot
        # block on a descendant that still holds the write end.
        for stream in (process.stdin, process.stdout, process.stderr):
            if stream is not None:
                stream.close()


class SubprocessRunner:
    """Execute an argument vector locally without a shell."""

    def run_to_file(
        self, argv: Sequence[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult:
        """Run a command, streaming its stdout into ``output_path`` in binary.

        The file is written under a temporary name in the destination directory and moved
        into place only after a zero exit code, so a failed or interrupted run never
        leaves a plausible-looking but truncated artifact behind.
        """
        normalized = _normalize(argv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, staged_name = tempfile.mkstemp(
            dir=output_path.parent, prefix=f".{output_path.name}.", suffix=".tmp"
        )
        staged = Path(staged_name)
        try:
            with os.fdopen(descriptor, "wb") as handle:
                process = subprocess.Popen(
                    normalized,
                    stdout=handle,
                    stderr=subprocess.PIPE,
                    stdin=subprocess.DEVNULL,
                )
                try:
                    _stdout, stderr_bytes = process.communicate(timeout=timeout_seconds)
                except subprocess.TimeoutExpired as exc:
                    _terminate_process_tree(process)
                    raise ToolExecutionError(
                        f"Command timed out after {timeout_seconds} seconds: {normalized[0]}"
                    ) from exc
                except BaseException:
                    _terminate_process_tree(process)
                    raise
                handle.flush()
                os.fsync(handle.fileno())
        except FileNotFoundError as exc:
            staged.unlink(missing_ok=True)
            raise ToolExecutionError(f"Required executable not found: {normalized[0]}") from exc
        except ToolExecutionError:
            staged.unlink(missing_ok=True)
            raise
        except OSError as exc:
            staged.unlink(missing_ok=True)
            raise ToolExecutionError(f"Could not execute {normalized[0]}: {exc}") from exc
        except BaseException:
            staged.unlink(missing_ok=True)
            raise

        stderr = stderr_bytes.decode("utf-8", "replace") if stderr_bytes else ""
        if process.returncode != 0:
            staged.unlink(missing_ok=True)
        else:
            os.replace(staged, output_path)
        return CommandResult(
            argv=normalized, returncode=process.returncode, stdout="", stderr=stderr
        )

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        """Run a command and capture its output as text.

        Decoding is pinned to UTF-8 with replacement rather than left to ``text=True``,
        which would use the platform locale. Under the ``C``/``POSIX`` locale that
        resolves to ASCII, and a single non-ASCII byte anywhere in a tool banner — an
        em dash in a version string, a micro sign in a progress line — would raise
        ``UnicodeDecodeError`` from inside the runner. That would make the outcome of a
        stage depend on the operator's environment, which is the one thing a
        reproducible pipeline cannot allow. Replacement matches ``run_to_file``, and
        losing a character from a diagnostic line is preferable to losing the run.
        """
        normalized = _normalize(argv)
        try:
            process = subprocess.Popen(
                normalized,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdin=subprocess.DEVNULL,
                encoding="utf-8",
                errors="replace",
            )
            try:
                stdout, stderr = process.communicate(timeout=timeout_seconds)
            except subprocess.TimeoutExpired as exc:
                _terminate_process_tree(process)
                raise ToolExecutionError(
                    f"Command timed out after {timeout_seconds} seconds: {normalized[0]}"
                ) from exc
            except BaseException:
                _terminate_process_tree(process)
                raise
        except FileNotFoundError as exc:
            raise ToolExecutionError(f"Required executable not found: {normalized[0]}") from exc
        except ToolExecutionError:
            raise
        except OSError as exc:
            raise ToolExecutionError(f"Could not execute {normalized[0]}: {exc}") from exc
        return CommandResult(
            argv=normalized,
            returncode=process.returncode,
            stdout=stdout,
            stderr=stderr,
        )

from __future__ import annotations

import os
import signal
import subprocess
import tempfile
from collections.abc import Sequence
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


def _terminate_process_tree(process: subprocess.Popen[bytes] | subprocess.Popen[str]) -> None:
    """Terminate the owned process tree and reap the direct child.

    Called on a timeout and on any other abort while waiting (``KeyboardInterrupt``,
    ``SystemExit``). POSIX tools run in a dedicated session, so the terminal's SIGINT no
    longer reaches them; signalling that process group here is what stops the tool and its
    multiprocessing descendants instead of leaving them orphaned. Native Windows keeps the
    previous direct-child termination semantics until a separately tested
    job-object/process-tree contract is introduced.
    """
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except OSError:
            process.kill()
    else:
        process.kill()
    process.wait()


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
                    start_new_session=os.name == "posix",
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
                start_new_session=os.name == "posix",
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

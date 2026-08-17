from __future__ import annotations

import os
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
                completed = subprocess.run(
                    normalized,
                    check=False,
                    stdout=handle,
                    stderr=subprocess.PIPE,
                    stdin=subprocess.DEVNULL,
                    timeout=timeout_seconds,
                )
                handle.flush()
                os.fsync(handle.fileno())
        except FileNotFoundError as exc:
            staged.unlink(missing_ok=True)
            raise ToolExecutionError(f"Required executable not found: {normalized[0]}") from exc
        except subprocess.TimeoutExpired as exc:
            staged.unlink(missing_ok=True)
            raise ToolExecutionError(
                f"Command timed out after {timeout_seconds} seconds: {normalized[0]}"
            ) from exc
        except OSError as exc:
            staged.unlink(missing_ok=True)
            raise ToolExecutionError(f"Could not execute {normalized[0]}: {exc}") from exc

        stderr = completed.stderr.decode("utf-8", "replace") if completed.stderr else ""
        if completed.returncode != 0:
            staged.unlink(missing_ok=True)
        else:
            os.replace(staged, output_path)
        return CommandResult(
            argv=normalized, returncode=completed.returncode, stdout="", stderr=stderr
        )

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = _normalize(argv)
        try:
            completed = subprocess.run(
                normalized,
                check=False,
                capture_output=True,
                stdin=subprocess.DEVNULL,
                text=True,
                timeout=timeout_seconds,
            )
        except FileNotFoundError as exc:
            raise ToolExecutionError(f"Required executable not found: {normalized[0]}") from exc
        except subprocess.TimeoutExpired as exc:
            raise ToolExecutionError(
                f"Command timed out after {timeout_seconds} seconds: {normalized[0]}"
            ) from exc
        except OSError as exc:
            raise ToolExecutionError(f"Could not execute {normalized[0]}: {exc}") from exc
        return CommandResult(
            argv=normalized,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
        )

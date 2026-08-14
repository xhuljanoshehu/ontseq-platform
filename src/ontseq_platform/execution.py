from __future__ import annotations

import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol


@dataclass(frozen=True)
class CommandResult:
    argv: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


class CommandRunner(Protocol):
    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult: ...


class ToolExecutionError(RuntimeError):
    """Raised when a required local executable cannot be started or times out."""


class SubprocessRunner:
    """Execute an argument vector locally without a shell."""

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(str(item) for item in argv)
        if not normalized or any(not item or "\x00" in item for item in normalized):
            raise ValueError("Command arguments must be non-empty and cannot contain NUL bytes")
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

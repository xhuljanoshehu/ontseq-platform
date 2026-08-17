"""An exclusive lock over one run envelope.

Nothing prevented two ``ontseq run`` invocations from working on the same envelope. Under a
person typing commands that is a theoretical problem. Under a watch folder it is not: a
watcher that fires twice — a duplicate filesystem event, a restart mid-run, two operators
pointing at the same drop directory — puts two processes into the same directory, where
they race on the run report and on each other's resume decisions.

The failure is quiet, which is what makes it worth preventing rather than documenting.
Atomic writes stop an artifact from being *truncated*, and content-addressed resume stops a
stale artifact from being *accepted*. Neither notices that a second process is rewriting the
same run report from a different set of stage records; the loser simply disappears from the
history with nothing recording that it existed.

What this does and does not promise
-----------------------------------

**Acquisition is atomic** via ``O_CREAT | O_EXCL``, so two processes cannot both believe
they hold the lock, even if they arrive in the same microsecond.

**A crashed run does not block the next one forever.** The lock records the host and PID
that took it. When the holder is on this machine and that process is gone, the lock is
reclaimed and the reclaim is reported — a run interrupted by a power cut must not need
manual cleanup before it can be resumed.

**A lock held from another host is never reclaimed automatically.** On shared storage this
process cannot tell a crashed remote run from a running one, and guessing wrong means two
live runs in one envelope — the exact thing the lock exists to prevent. It fails closed and
names the file to remove.

**PID reuse fails closed.** If the recorded PID has been recycled by an unrelated process,
the lock looks held and the run refuses. Refusing a run that could have proceeded is a
delay; proceeding on a run that should have refused is a corrupted envelope.

This module is dependency-free so it can be tested without pydantic.
"""

from __future__ import annotations

import contextlib
import json
import os
import socket
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

#: Sits at the envelope root rather than under ``work/``, which callers may clean.
LOCK_FILENAME = ".ontseq-run.lock"


@dataclass(frozen=True)
class LockHolder:
    """Who took the lock, recorded so a human can decide what to do about it."""

    pid: int
    hostname: str
    acquired_at: str
    run_id: str
    sample_id: str
    pipeline_version: str

    def describe(self) -> str:
        return (
            f"pid {self.pid} on {self.hostname}, "
            f"run {self.run_id}/{self.sample_id}, "
            f"acquired {self.acquired_at}"
        )


class LockError(RuntimeError):
    """Raised when the run envelope cannot be locked."""


class RunAlreadyRunning(LockError):
    """Raised when another process holds the lock.

    Carries the holder so a caller can report *who* rather than only *that* it is held.
    """

    def __init__(self, holder: LockHolder | None, lock_path: Path, detail: str) -> None:
        self.holder = holder
        self.lock_path = lock_path
        super().__init__(detail)


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _process_is_alive(pid: int) -> bool:
    """Return whether a PID exists on this machine.

    ``PermissionError`` means the process exists and belongs to someone else, which counts
    as alive. Any other OS error is treated as alive too: an unreadable answer must not be
    read as permission to proceed.
    """
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return True
    return True


def read_holder(lock_path: Path) -> LockHolder | None:
    """Read the lock file, returning ``None`` when it is absent, empty or unreadable.

    A lock file written by a process that died between ``O_EXCL`` and its first write is
    empty. Treating that as "unknown holder" rather than as an error is what lets the stale
    path reclaim it instead of demanding manual cleanup for a zero-byte file.
    """
    try:
        payload = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    try:
        return LockHolder(
            pid=int(payload["pid"]),
            hostname=str(payload["hostname"]),
            acquired_at=str(payload["acquired_at"]),
            run_id=str(payload["run_id"]),
            sample_id=str(payload["sample_id"]),
            pipeline_version=str(payload["pipeline_version"]),
        )
    except (KeyError, TypeError, ValueError):
        return None


def holder_is_running(holder: LockHolder) -> bool | None:
    """Return whether the lock holder still exists, or ``None`` when that is unknowable.

    ``None`` is not a shrug — it is the honest answer for a lock taken on another host,
    where this process has no way to look. A caller reporting run state must show that
    difference rather than collapse it into "probably fine".
    """
    if holder.hostname != socket.gethostname():
        return None
    return _process_is_alive(holder.pid)


def _stale_reason(holder: LockHolder | None, lock_path: Path) -> str | None:
    """Return why an existing lock may be reclaimed, or ``None`` when it may not."""
    if holder is None:
        return f"{lock_path.name} was unreadable or incomplete, so no live holder is recorded"
    if holder.hostname != socket.gethostname():
        return None
    if _process_is_alive(holder.pid):
        return None
    return f"the process that held it is gone ({holder.describe()})"


def _write_lock(descriptor: int, holder: LockHolder) -> None:
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        json.dump(asdict(holder), handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


@contextlib.contextmanager
def run_lock(
    envelope_root: Path,
    *,
    run_id: str,
    sample_id: str,
    pipeline_version: str,
    reclaim_stale: bool = True,
) -> Iterator[list[str]]:
    """Hold an exclusive lock on one run envelope for the duration of the block.

    Yields a list of warnings — empty on a clean acquisition, one entry when a stale lock
    was reclaimed. The caller is expected to carry those into the run report: silently
    stepping over another run's lock file is exactly the kind of event that must survive
    into the record rather than scroll past in a terminal.

    Raises :class:`RunAlreadyRunning` when the envelope is in use.
    """
    envelope_root.mkdir(parents=True, exist_ok=True)
    lock_path = envelope_root / LOCK_FILENAME
    holder = LockHolder(
        pid=os.getpid(),
        hostname=socket.gethostname(),
        acquired_at=_now(),
        run_id=run_id,
        sample_id=sample_id,
        pipeline_version=pipeline_version,
    )
    warnings: list[str] = []

    try:
        descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        existing = read_holder(lock_path)
        reason = _stale_reason(existing, lock_path) if reclaim_stale else None
        if reason is None:
            raise RunAlreadyRunning(
                existing,
                lock_path,
                _held_message(existing, lock_path),
            ) from None
        # Reclaiming is a plain unlink followed by a fresh O_EXCL create rather than an
        # overwrite, so two processes that both judge the lock stale still cannot both win.
        lock_path.unlink(missing_ok=True)
        try:
            descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            raise RunAlreadyRunning(
                read_holder(lock_path),
                lock_path,
                "another process reclaimed this run envelope first",
            ) from None
        warnings.append(f"Reclaimed a stale run lock: {reason}.")

    _write_lock(descriptor, holder)
    try:
        yield warnings
    finally:
        # Only remove a lock this process still owns. A lock reclaimed from us — because we
        # were wrongly judged dead — belongs to whoever holds it now.
        current = read_holder(lock_path)
        if current is not None and (current.pid, current.hostname) == (
            holder.pid,
            holder.hostname,
        ):
            lock_path.unlink(missing_ok=True)


def _held_message(holder: LockHolder | None, lock_path: Path) -> str:
    if holder is None:
        return (
            f"this run envelope is locked by {lock_path}, and the lock could not be read. "
            "Remove it only after confirming no run is in progress."
        )
    if holder.hostname != socket.gethostname():
        return (
            f"this run envelope is locked from another host: {holder.describe()}. "
            "A remote lock is never reclaimed automatically, because a crashed remote run "
            "and a running one look identical from here. Confirm the run has stopped, then "
            f"remove {lock_path}."
        )
    return (
        f"this run envelope is already in use: {holder.describe()}. "
        f"If that process has stopped, remove {lock_path}."
    )

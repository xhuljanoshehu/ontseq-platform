"""Deciding *what* to process and *whether it is ready* — without processing anything.

A watch folder is mostly two judgements, and both are easy to get subtly wrong:

**Is this directory finished being written?** A sequencer writes for hours. Picking a run up
while files are still arriving produces a truncated analysis that looks complete, which is
worse than not picking it up at all. Two mechanisms are offered and they are not equivalent:
an explicit marker file is authoritative and should be preferred wherever the producer can
write one; a quiescence window is a heuristic for producers that cannot, and it is a
heuristic — a slow network copy that stalls for longer than the window looks finished.

**Has this already been handled?** Without a memory, a failing sample is retried on every
poll forever, filling logs and burning a GPU on the same doomed run. The ledger records the
outcome of every attempt, and a failure is *not* retried by default: a deterministic failure
does not become a success by being repeated, and an operator noticing a stuck sample is the
intended outcome.

Everything here is pure judgement over the filesystem: no manifests, no pipeline, no
pydantic. That keeps the part with the awkward edge cases — clock skew, partial writes,
empty directories — testable on its own.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
import time
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

#: Kept beside the *output*, never inside the watch folder: a drop directory may be owned by
#: the instrument, mounted read-only, or wiped by whoever fills it.
LEDGER_FILENAME = ".ontseq-watch.json"

#: Directory names that are never work: dotfiles, and the places a run writes to.
IGNORED_PREFIXES = (".", "_")

#: Mirrors the sample_id contract in the manifest, so a directory that cannot yield a legal
#: sample identifier is rejected at discovery rather than deep inside validation.
SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")


class Readiness(StrEnum):
    """Why a candidate may or may not be processed yet."""

    READY = "ready"
    EMPTY = "empty"
    MARKER_MISSING = "marker_missing"
    STILL_CHANGING = "still_changing"
    UNREADABLE = "unreadable"


class Outcome(StrEnum):
    """What happened on an attempt. Recorded so it is never attempted blindly again."""

    COMPLETED = "completed"
    FAILED = "failed"
    #: Another process held the run envelope. Not a failure of this sample — try later.
    LOCKED = "locked"
    #: Rejected before any pipeline work, e.g. an unusable directory name.
    REJECTED = "rejected"


@dataclass(frozen=True)
class Candidate:
    """One drop directory, with the readiness judgement already made."""

    path: Path
    name: str
    readiness: Readiness
    detail: str
    #: Newest modification time seen inside, or ``None`` when nothing was readable.
    newest_mtime: float | None = None
    file_count: int = 0

    @property
    def is_ready(self) -> bool:
        return self.readiness is Readiness.READY


def _walk_files(root: Path) -> Iterator[Path]:
    """Yield every regular file below ``root``, ignoring what disappears mid-walk.

    A directory being written into changes under the walk. Treating a vanished entry as an
    error would make readiness depend on timing rather than on content.
    """
    for directory, _, filenames in os.walk(root):
        for filename in filenames:
            candidate = Path(directory) / filename
            try:
                if candidate.is_file():
                    yield candidate
            except OSError:
                continue


def inspect_directory(
    path: Path,
    *,
    ready_marker: str | None,
    quiet_seconds: float,
    now: float,
) -> tuple[Readiness, str, float | None, int]:
    """Judge one directory. Returns the readiness, a reason, the newest mtime and a count.

    When ``ready_marker`` is set it is the *only* thing consulted, because an explicit
    signal from the producer beats inferring completion from timestamps. Without one,
    quiescence is used and its heuristic nature is stated in the reason.
    """
    newest: float | None = None
    count = 0
    try:
        for item in _walk_files(path):
            count += 1
            try:
                mtime = item.stat().st_mtime
            except OSError:
                continue
            newest = mtime if newest is None else max(newest, mtime)
    except OSError as error:
        return Readiness.UNREADABLE, f"could not be read: {error}", None, 0

    if count == 0:
        return Readiness.EMPTY, "contains no files", None, 0

    if ready_marker is not None:
        if (path / ready_marker).exists():
            return Readiness.READY, f"the producer wrote {ready_marker}", newest, count
        return (
            Readiness.MARKER_MISSING,
            f"waiting for the producer to write {ready_marker}",
            newest,
            count,
        )

    if newest is None:
        return Readiness.UNREADABLE, "no modification time could be read", None, count

    age = now - newest
    if age < quiet_seconds:
        return (
            Readiness.STILL_CHANGING,
            f"last modified {age:.0f}s ago, below the {quiet_seconds:.0f}s quiet window",
            newest,
            count,
        )
    return (
        Readiness.READY,
        f"unchanged for {age:.0f}s (quiescence is a heuristic, not a completion signal)",
        newest,
        count,
    )


def discover(
    watch_dir: Path,
    *,
    ready_marker: str | None = None,
    quiet_seconds: float = 300.0,
    now: float | None = None,
) -> list[Candidate]:
    """Return every drop directory under ``watch_dir``, ready or not, sorted by name.

    Not-ready candidates are returned rather than filtered out, so a caller can report *why*
    a sample is sitting untouched. A watch folder that silently ignores things is a watch
    folder nobody trusts.
    """
    if not watch_dir.is_dir():
        raise NotADirectoryError(f"watch folder does not exist: {watch_dir}")
    moment = time.time() if now is None else now
    candidates: list[Candidate] = []
    for entry in sorted(watch_dir.iterdir(), key=lambda item: item.name):
        if not entry.is_dir() or entry.name.startswith(IGNORED_PREFIXES):
            continue
        readiness, detail, newest, count = inspect_directory(
            entry, ready_marker=ready_marker, quiet_seconds=quiet_seconds, now=moment
        )
        candidates.append(
            Candidate(
                path=entry,
                name=entry.name,
                readiness=readiness,
                detail=detail,
                newest_mtime=newest,
                file_count=count,
            )
        )
    return candidates


@dataclass(frozen=True)
class Attempt:
    """What a previous pass did with one drop directory."""

    name: str
    outcome: Outcome
    detail: str
    attempted_at: str
    run_id: str | None = None
    sample_id: str | None = None


class Ledger:
    """A durable memory of what has already been attempted.

    Written atomically, because a watcher killed mid-write must not come back to a corrupt
    ledger and treat every finished sample as new work.
    """

    def __init__(self, path: Path, attempts: dict[str, Attempt] | None = None) -> None:
        self.path = path
        self._attempts: dict[str, Attempt] = dict(attempts or {})

    @classmethod
    def load(cls, path: Path) -> Ledger:
        """Read the ledger, starting empty when it is absent or unreadable.

        A corrupt ledger is treated as no ledger. That re-attempts finished work, which
        resume makes cheap; the alternative — refusing to start — turns one bad file into a
        stopped watcher.
        """
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            entries = payload["attempts"]
        except (OSError, ValueError, KeyError, TypeError):
            return cls(path)
        attempts: dict[str, Attempt] = {}
        for name, item in entries.items() if isinstance(entries, dict) else []:
            try:
                attempts[str(name)] = Attempt(
                    name=str(name),
                    outcome=Outcome(item["outcome"]),
                    detail=str(item["detail"]),
                    attempted_at=str(item["attempted_at"]),
                    run_id=item.get("run_id"),
                    sample_id=item.get("sample_id"),
                )
            except (KeyError, TypeError, ValueError):
                continue
        return cls(path, attempts)

    def get(self, name: str) -> Attempt | None:
        return self._attempts.get(name)

    def record(self, attempt: Attempt) -> None:
        self._attempts[attempt.name] = attempt
        self.save()

    def attempts(self) -> Iterable[Attempt]:
        return tuple(self._attempts[name] for name in sorted(self._attempts))

    def save(self) -> None:
        payload = {
            "schema_version": "0.1.0",
            "attempts": {
                name: {
                    "outcome": attempt.outcome.value,
                    "detail": attempt.detail,
                    "attempted_at": attempt.attempted_at,
                    "run_id": attempt.run_id,
                    "sample_id": attempt.sample_id,
                }
                for name, attempt in sorted(self._attempts.items())
            },
        }
        self.path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, staged = tempfile.mkstemp(
            dir=self.path.parent, prefix=f".{self.path.name}.", suffix=".tmp"
        )
        try:
            with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2, sort_keys=True)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(staged, self.path)
        except BaseException:
            Path(staged).unlink(missing_ok=True)
            raise


def should_attempt(
    candidate: Candidate, previous: Attempt | None, *, retry_failed: bool = False
) -> tuple[bool, str]:
    """Decide whether to run this candidate now, and say why not when the answer is no.

    A previous ``LOCKED`` always retries: it means another process was working on the
    envelope, which says nothing about this sample. A previous ``FAILED`` does not, unless
    asked: a deterministic failure repeated every poll is noise that hides the one sample
    somebody needs to look at.
    """
    if not candidate.is_ready:
        return False, candidate.detail
    if previous is None:
        return True, "not attempted before"
    if previous.outcome is Outcome.COMPLETED:
        return False, f"already completed at {previous.attempted_at}"
    if previous.outcome is Outcome.LOCKED:
        return True, "previously blocked by another run; retrying"
    if previous.outcome is Outcome.REJECTED:
        return False, f"rejected at {previous.attempted_at}: {previous.detail}"
    if retry_failed:
        return True, f"retrying an earlier failure from {previous.attempted_at}"
    return (
        False,
        f"failed at {previous.attempted_at} and is not retried automatically "
        f"({previous.detail}); pass --retry-failed once the cause is understood",
    )


def sample_id_from_directory(name: str) -> str | None:
    """Derive a sample identifier from a directory name, or ``None`` when it cannot be.

    Deliberately conservative: no cleaning, no truncation, no substitution. A name that does
    not already satisfy the manifest contract is rejected so that somebody writes the
    identifier down, rather than the watcher inventing one that then appears on a report.
    """
    return name if SAMPLE_ID.fullmatch(name) else None

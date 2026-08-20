"""Answering "what happened to these runs?" without anyone parsing JSON by hand.

A run leaves a complete record in ``provenance/run.json``, and a watcher leaves one in its
ledger. Both are machine-readable and neither is readable at a glance, which is a problem
the moment the pipeline runs unattended: the person who needs to know whether last night
went well should not have to write a script to find out.

Two states this reports are worth naming, because they are the ones a bare directory
listing cannot show:

**RUNNING versus INTERRUPTED.** A lock file means a run claimed this envelope. Whether that
run still exists is a different question, and the answer is not always knowable — for a lock
taken on another host it is genuinely unknown, and this says so rather than guessing. An
INTERRUPTED envelope is where a run died; the next attempt will reclaim the lock and resume,
which is worth knowing before somebody deletes the directory to "start clean".

**UNFINISHED.** An envelope with neither a lock nor a run report is one where a run started
and never wrote a verdict. Rare, and exactly the case that a summary counting only PASS and
FAIL would hide by omission.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from .models import ModuleRunStatus
from .pipeline.envelope import sha256_file
from .pipeline.lock import LOCK_FILENAME, LockHolder, holder_is_running, read_holder
from .pipeline.review import RELEASE_RELATIVE, REVIEW_LOG, ReviewError, ReviewState
from .pipeline.review import current_state as review_state
from .pipeline.review import read_log as read_review_log
from .pipeline.state import RunReport
from .pipeline.watch import LEDGER_FILENAME, Ledger

RUN_REPORT_RELATIVE = "provenance/run.json"


class RunState(StrEnum):
    """What an envelope is, as opposed to what its last report said."""

    PASSED = "passed"
    FAILED = "failed"
    #: A lock is held by a process that still exists, or by an unreachable host.
    RUNNING = "running"
    #: A lock is held by a process that is gone. The next run reclaims it and resumes.
    INTERRUPTED = "interrupted"
    #: Neither lock nor report: a run started here and never recorded a verdict.
    UNFINISHED = "unfinished"
    #: A report exists but cannot be parsed, e.g. written by an incompatible version.
    UNREADABLE = "unreadable"


@dataclass(frozen=True)
class EnvelopeStatus:
    """One run envelope, as it stands right now."""

    run_id: str
    sample_id: str
    root: Path
    state: RunState
    holder: LockHolder | None = None
    #: ``None`` when the holder is on another host and liveness cannot be determined.
    holder_alive: bool | None = None
    report: RunReport | None = None
    detail: str = ""
    #: Whether anybody has signed this envelope off, and whether that still applies to the
    #: content on disk. Reported, never used to decide the exit code: a run that nobody has
    #: reviewed yet is not a fault, and a monitoring check that fires on every fresh run
    #: teaches people to ignore it.
    review: ReviewState | None = None
    review_detail: str = ""

    @property
    def unverified_stages(self) -> tuple[str, ...]:
        if self.report is None:
            return ()
        return tuple(item.value for item in self.report.unverified_stages)


def _review(root: Path) -> tuple[ReviewState | None, str]:
    """Resolve the envelope's sign-off trail, or ``(None, "")`` when it has none.

    Read here so an operator scanning many envelopes sees run health and sign-off in one
    pass. It never influences the exit code: `ontseq status` answers "did the runs work",
    and `ontseq review status` answers "may this leave the system". Folding the second into
    the first would make the check fire on every fresh run, which is how a monitoring signal
    becomes noise.
    """
    log = root / REVIEW_LOG
    if not log.is_file():
        return None, ""
    release = root / RELEASE_RELATIVE
    digest = sha256_file(release) if release.is_file() else None
    try:
        entries = read_review_log(log)
    except ReviewError as error:
        return ReviewState.UNREADABLE, str(error)
    return review_state(entries, digest)


def _classify(root: Path) -> EnvelopeStatus:
    run_id, sample_id = root.parent.name, root.name
    lock_path = root / LOCK_FILENAME
    report: RunReport | None = None
    detail = ""
    review, review_detail = _review(root)

    report_path = root / RUN_REPORT_RELATIVE
    if report_path.is_file():
        try:
            report = RunReport.model_validate_json(report_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            return EnvelopeStatus(
                run_id=run_id,
                sample_id=sample_id,
                root=root,
                state=RunState.UNREADABLE,
                detail=f"{RUN_REPORT_RELATIVE} could not be read: {error}",
                review=review,
                review_detail=review_detail,
            )

    if lock_path.exists():
        # A held lock outranks whatever the last report said: a run is in progress here, and
        # the report on disk is a snapshot from partway through it.
        holder = read_holder(lock_path)
        alive = holder_is_running(holder) if holder is not None else None
        if alive is False:
            detail = "a run died here; the next attempt reclaims the lock and resumes"
            state = RunState.INTERRUPTED
        else:
            detail = (
                "held from another host, so liveness cannot be checked from here"
                if alive is None and holder is not None
                else "a run is in progress"
            )
            state = RunState.RUNNING
        return EnvelopeStatus(
            run_id=run_id,
            sample_id=sample_id,
            root=root,
            state=state,
            holder=holder,
            holder_alive=alive,
            report=report,
            detail=detail,
            review=review,
            review_detail=review_detail,
        )

    if report is None:
        return EnvelopeStatus(
            run_id=run_id,
            sample_id=sample_id,
            root=root,
            state=RunState.UNFINISHED,
            detail="no run report was written; the run never reached a verdict",
            review=review,
            review_detail=review_detail,
        )
    return EnvelopeStatus(
        run_id=run_id,
        sample_id=sample_id,
        root=root,
        state=RunState.PASSED if report.passed else RunState.FAILED,
        report=report,
        detail=report.verdict_reason,
        review=review,
        review_detail=review_detail,
    )


def scan(output_dir: Path, *, run_id: str | None = None) -> list[EnvelopeStatus]:
    """Find every ``<run-id>/<sample-id>/`` envelope beneath an output directory.

    Envelopes are recognised by their shape rather than by a registry, so a directory
    copied in from elsewhere is reported like any other. Anything that does not look like an
    envelope is ignored silently — this is a reporting tool, not a validator of the output
    directory's tidiness.
    """
    if not output_dir.is_dir():
        raise NotADirectoryError(f"output directory does not exist: {output_dir}")
    found: list[EnvelopeStatus] = []
    for run_directory in sorted(output_dir.iterdir(), key=lambda item: item.name):
        if not run_directory.is_dir() or run_directory.name.startswith("."):
            continue
        if run_id is not None and run_directory.name != run_id:
            continue
        for sample_directory in sorted(run_directory.iterdir(), key=lambda item: item.name):
            if not sample_directory.is_dir() or sample_directory.name.startswith("."):
                continue
            if not _looks_like_envelope(sample_directory):
                continue
            found.append(_classify(sample_directory))
    return found


def _looks_like_envelope(path: Path) -> bool:
    """A run envelope always has a manifest directory, whatever else it got to."""
    return (path / "manifest").is_dir() or (path / RUN_REPORT_RELATIVE).is_file()


def render_text(statuses: list[EnvelopeStatus], *, verbose: bool = False) -> str:
    """Render a human summary. One line per envelope, plus stage detail when asked."""
    if not statuses:
        return "no run envelopes found"
    lines: list[str] = []
    for status in statuses:
        lines.append(
            f"{status.run_id}/{status.sample_id:<24} {status.state.value.upper():<12} "
            f"{status.detail}"
        )
        if status.holder is not None:
            lines.append(f"    lock: {status.holder.describe()}")
        if status.review is not None:
            lines.append(f"    review: {status.review.value.upper()} — {status.review_detail}")
        if status.unverified_stages:
            lines.append(
                "    UNVERIFIED ADAPTERS COMPLETED: " + ", ".join(status.unverified_stages)
            )
        if verbose and status.report is not None:
            for record in status.report.stages:
                marker = "resumed" if record.resumed else record.status.value
                lines.append(f"    {record.stage.value:<16} {marker:<10} {record.reason}")

    counts: dict[str, int] = {}
    for status in statuses:
        counts[status.state.value] = counts.get(status.state.value, 0) + 1
    summary = ", ".join(f"{count} {state}" for state, count in sorted(counts.items()))
    lines.append(f"total: {len(statuses)} envelope(s) — {summary}")
    return "\n".join(lines)


def render_json(statuses: list[EnvelopeStatus]) -> str:
    """Render the same information for a monitoring script."""
    payload = [
        {
            "run_id": status.run_id,
            "sample_id": status.sample_id,
            "state": status.state.value,
            "detail": status.detail,
            "holder": None
            if status.holder is None
            else {
                "pid": status.holder.pid,
                "hostname": status.holder.hostname,
                "acquired_at": status.holder.acquired_at,
                "alive": status.holder_alive,
            },
            "unverified_stages": list(status.unverified_stages),
            "review": None if status.review is None else status.review.value,
            "review_detail": status.review_detail,
            "stages": []
            if status.report is None
            else [
                {
                    "stage": record.stage.value,
                    "status": record.status.value,
                    "resumed": record.resumed,
                    "reason": record.reason,
                }
                for record in status.report.stages
            ],
        }
        for status in statuses
    ]
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def render_ledger(output_dir: Path) -> str:
    """Render the watch folder's memory, when one exists beside this output directory."""
    path = output_dir / LEDGER_FILENAME
    if not path.is_file():
        return ""
    attempts = list(Ledger.load(path).attempts())
    if not attempts:
        return f"watch ledger at {path.name}: no attempts recorded"
    lines = [f"watch ledger ({len(attempts)} attempt(s)):"]
    for attempt in attempts:
        lines.append(f"    {attempt.name:<28} {attempt.outcome.value.upper():<10} {attempt.detail}")
    return "\n".join(lines)


def exit_code(statuses: list[EnvelopeStatus]) -> int:
    """Map the worst state present onto an exit code a monitoring check can act on.

    Deliberately does not treat RUNNING as a problem: a monitoring check that fires while
    the pipeline is doing its job trains people to ignore it.
    """
    states = {status.state for status in statuses}
    if states & {RunState.FAILED, RunState.UNREADABLE}:
        return 2
    if states & {RunState.INTERRUPTED, RunState.UNFINISHED}:
        return 6
    return 0


def stage_counts(status: EnvelopeStatus) -> dict[str, int]:
    """Count stage outcomes, for a caller that wants a one-line health figure."""
    counts: dict[str, int] = {item.value: 0 for item in ModuleRunStatus}
    if status.report is not None:
        for record in status.report.stages:
            counts[record.status.value] += 1
    return counts

"""Wiring for ``ontseq review``: reading an envelope's trail and rendering it.

The judgement logic lives in :mod:`ontseq_platform.pipeline.review`, which carries no
dependencies and is unit tested on its own. This module knows only where things sit inside
a run envelope and how to present them.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from .pipeline.envelope import sha256_file
from .pipeline.review import (
    RELEASE_RELATIVE,
    REVIEW_LOG,
    Decision,
    ReviewEntry,
    ReviewError,
    ReviewState,
    accepted_reviewers,
    append_entry,
    current_state,
    read_log,
    verify_chain,
)


@dataclass(frozen=True)
class ReviewReport:
    """An envelope's review trail, resolved against the content currently on disk."""

    run_id: str
    sample_id: str
    root: Path
    state: ReviewState
    detail: str
    entries: tuple[ReviewEntry, ...]
    release_sha256: str | None
    reviewers: tuple[str, ...]


def release_digest(envelope_root: Path) -> str | None:
    """Digest of the release bundle, or ``None`` when the run produced none.

    A run that failed before the release stage has nothing to bind a review to, and that is
    reported rather than papered over with an empty-string digest that would compare equal
    to the next equally empty one.
    """
    path = envelope_root / RELEASE_RELATIVE
    return sha256_file(path) if path.is_file() else None


def inspect(envelope_root: Path) -> ReviewReport:
    """Resolve one envelope's review trail. Never writes."""
    if not envelope_root.is_dir():
        raise NotADirectoryError(f"run envelope does not exist: {envelope_root}")
    digest = release_digest(envelope_root)
    try:
        entries = read_log(envelope_root / REVIEW_LOG)
    except ReviewError as error:
        return ReviewReport(
            run_id=envelope_root.parent.name,
            sample_id=envelope_root.name,
            root=envelope_root,
            state=ReviewState.UNREADABLE,
            detail=str(error),
            entries=(),
            release_sha256=digest,
            reviewers=(),
        )
    state, detail = current_state(entries, digest)
    return ReviewReport(
        run_id=envelope_root.parent.name,
        sample_id=envelope_root.name,
        root=envelope_root,
        state=state,
        detail=detail,
        entries=tuple(entries),
        release_sha256=digest,
        reviewers=accepted_reviewers(entries, digest),
    )


def record(
    envelope_root: Path,
    *,
    decision: Decision,
    reviewer: str,
    note: str = "",
) -> ReviewEntry:
    """Append one judgement to an envelope's trail.

    Refuses when the run produced no release bundle: there would be nothing to bind the
    judgement to, and a review of unspecified content is not a review.
    """
    if not envelope_root.is_dir():
        raise NotADirectoryError(f"run envelope does not exist: {envelope_root}")
    digest = release_digest(envelope_root)
    if digest is None:
        raise ReviewError(
            f"{envelope_root} has no {RELEASE_RELATIVE}, so there is nothing to review. A "
            "run that did not reach the release stage cannot be signed off."
        )
    return append_entry(
        envelope_root / REVIEW_LOG,
        decision=decision,
        reviewer=reviewer,
        run_id=envelope_root.parent.name,
        sample_id=envelope_root.name,
        release_sha256=digest,
        note=note,
    )


def render_text(report: ReviewReport, *, verbose: bool = False) -> str:
    """Render for a person deciding whether a result may leave the system."""
    lines = [
        f"{report.run_id}/{report.sample_id}: {report.state.value.upper()}",
        f"  {report.detail}",
    ]
    if report.release_sha256 is None:
        lines.append("  no release bundle: this run produced nothing reviewable")
    else:
        lines.append(f"  release bundle {report.release_sha256[:16]}…")
    if report.reviewers:
        lines.append(f"  accepted by: {', '.join(report.reviewers)}")
    if report.entries:
        intact, detail = verify_chain(report.entries)
        lines.append(f"  chain: {'intact' if intact else 'BROKEN'} — {detail}")
    if verbose:
        for entry in report.entries:
            marker = "" if entry.release_sha256 == report.release_sha256 else "  (other content)"
            lines.append(f"    {entry.describe()}{marker}")
            if entry.note:
                lines.append(f"      note: {entry.note}")
    lines.append(
        "  identities are asserted, not authenticated; the chain is tamper-evident, not "
        "tamper-proof"
    )
    return "\n".join(lines)


def render_json(report: ReviewReport) -> str:
    """Render the same information for a release script."""
    payload = {
        "run_id": report.run_id,
        "sample_id": report.sample_id,
        "state": report.state.value,
        "detail": report.detail,
        "release_sha256": report.release_sha256,
        "accepted_by": list(report.reviewers),
        "entries": [
            {
                "sequence": entry.sequence,
                "decision": entry.decision.value,
                "reviewer": entry.reviewer,
                "identity_source": entry.identity_source,
                "recorded_at": entry.recorded_at,
                "release_sha256": entry.release_sha256,
                "note": entry.note,
                "entry_sha256": entry.entry_sha256,
            }
            for entry in report.entries
        ],
        "identity_is_authenticated": False,
        "chain_is_tamper_proof": False,
    }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"

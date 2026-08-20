"""Recording that a person examined a run and took responsibility for it.

A run produces evidence. Someone has to look at it and say so, and that statement has to
survive as a record — otherwise "this was reviewed" is a claim nobody can check later.

Four properties make the record worth having, and each one is a deliberate choice about
what *cannot* happen:

**A review is bound to content, not to a directory.** Every entry carries the SHA-256 of
the release bundle it was made against. Change anything the bundle covers and the review no
longer applies to what is on disk — it becomes ``STALE`` rather than silently continuing to
vouch for something else. This is the same rule the resume logic uses: identity by content,
never by name or timestamp.

**The log is append-only and its shape says so.** Each entry names the digest of the entry
before it, so removing, reordering or editing one breaks the chain at that point and every
point after it. A reader can detect that without any external record.

**Nothing is overwritten.** A later review does not replace an earlier one; it is appended,
and the history shows both. A reviewer who accepted and then rejected leaves both facts
behind, which is the entire purpose of an audit trail.

**The record never claims more than it is.** The reviewer identity is *asserted* — it comes
from the command line and nothing authenticated it. The chain is tamper-*evident*, not
tamper-proof: it has no key, so anyone who can rewrite the file can also recompute the whole
chain. It detects accidental corruption and casual editing, and that is all it detects. A
qualified electronic signature needs an authorised key, an identity provider and a records
policy, none of which exist here — so this module says ``asserted`` rather than pretending.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import Any

#: Where the review trail lives inside a run envelope.
REVIEW_DIRECTORY = "review"
REVIEW_LOG = "review/review.log.jsonl"
#: The artifact a review is bound to. It already covers the run report and every exportable
#: artifact by checksum, so binding to it binds to all of them.
RELEASE_RELATIVE = "release/release.json"

#: How the reviewer identity was established. Only one value is truthful today.
ASSERTED = "asserted"


class Decision(StrEnum):
    """What a reviewer concluded. Deliberately only two outcomes.

    Anything else — "needs repeating", "insufficient material" — is a rejection with a
    reason, and the reason is free text because a closed vocabulary invented here would be
    a clinical taxonomy nobody has agreed to.
    """

    ACCEPTED = "accepted"
    REJECTED = "rejected"


class ReviewState(StrEnum):
    """The state of an envelope's review trail, as it stands right now."""

    #: No review has been recorded.
    PENDING = "pending"
    #: The most recent review accepts exactly the content currently on disk.
    ACCEPTED = "accepted"
    #: The most recent review rejects exactly the content currently on disk.
    REJECTED = "rejected"
    #: Reviews exist, but the release bundle has changed since the last of them. The
    #: earlier judgement stands for what it saw and says nothing about what is here now.
    STALE = "stale"
    #: The log does not verify: an entry was removed, reordered or edited.
    BROKEN = "broken"
    #: The log exists but cannot be parsed at all.
    UNREADABLE = "unreadable"


class ReviewError(RuntimeError):
    """Raised when a review cannot be recorded."""


@dataclass(frozen=True)
class ReviewEntry:
    """One recorded judgement, and what it was a judgement about."""

    sequence: int
    decision: Decision
    reviewer: str
    identity_source: str
    recorded_at: str
    run_id: str
    sample_id: str
    #: Digest of ``release/release.json`` as it stood when this was recorded.
    release_sha256: str
    note: str
    previous_entry_sha256: str | None
    entry_sha256: str

    def describe(self) -> str:
        return (
            f"#{self.sequence} {self.decision.value.upper()} by {self.reviewer} "
            f"({self.identity_source}) at {self.recorded_at}"
        )


def _canonical(payload: Mapping[str, Any]) -> str:
    """Serialise an entry deterministically, so its digest is reproducible."""
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def entry_digest(payload: Mapping[str, Any]) -> str:
    """Digest of an entry's content, excluding the digest field itself."""
    body = {key: value for key, value in payload.items() if key != "entry_sha256"}
    return hashlib.sha256(_canonical(body).encode("utf-8")).hexdigest()


def _to_entry(payload: Mapping[str, Any]) -> ReviewEntry:
    """Build an entry from parsed JSON.

    Typed as ``Any`` rather than ``object`` because that is what it is: untrusted input
    whose shape is checked here by conversion, and whose failures the caller turns into a
    ``ReviewError`` naming the line.
    """
    return ReviewEntry(
        sequence=int(payload["sequence"]),
        decision=Decision(payload["decision"]),
        reviewer=str(payload["reviewer"]),
        identity_source=str(payload["identity_source"]),
        recorded_at=str(payload["recorded_at"]),
        run_id=str(payload["run_id"]),
        sample_id=str(payload["sample_id"]),
        release_sha256=str(payload["release_sha256"]),
        note=str(payload.get("note", "")),
        previous_entry_sha256=(
            None
            if payload.get("previous_entry_sha256") is None
            else str(payload["previous_entry_sha256"])
        ),
        entry_sha256=str(payload["entry_sha256"]),
    )


def _to_payload(entry: ReviewEntry) -> dict[str, Any]:
    return {
        "sequence": entry.sequence,
        "decision": entry.decision.value,
        "reviewer": entry.reviewer,
        "identity_source": entry.identity_source,
        "recorded_at": entry.recorded_at,
        "run_id": entry.run_id,
        "sample_id": entry.sample_id,
        "release_sha256": entry.release_sha256,
        "note": entry.note,
        "previous_entry_sha256": entry.previous_entry_sha256,
        "entry_sha256": entry.entry_sha256,
    }


def read_log(path: Path) -> list[ReviewEntry]:
    """Read every entry, in order. A missing log is an empty history, not an error."""
    if not path.is_file():
        return []
    entries: list[ReviewEntry] = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            entries.append(_to_entry(json.loads(line)))
        except (ValueError, KeyError, TypeError) as error:
            raise ReviewError(
                f"{path.name} line {number} is not a review entry: {error}"
            ) from error
    return entries


def verify_chain(entries: Sequence[ReviewEntry]) -> tuple[bool, str]:
    """Check that the log has not been edited, reordered or truncated in the middle.

    Detects accidental corruption and casual editing. It cannot detect a rewrite by someone
    who reruns the chain, because there is no key involved — see the module docstring.
    """
    previous: str | None = None
    for index, entry in enumerate(entries):
        if entry.sequence != index + 1:
            return False, f"entry {index + 1} is numbered {entry.sequence}; the log is out of order"
        if entry.previous_entry_sha256 != previous:
            return False, (
                f"entry {entry.sequence} does not follow the one before it; an entry was "
                "removed, reordered or edited"
            )
        recomputed = entry_digest(_to_payload(entry))
        if recomputed != entry.entry_sha256:
            return False, f"entry {entry.sequence} does not match its own digest; it was edited"
        previous = entry.entry_sha256
    return True, f"{len(entries)} entr(ies) verified"


def append_entry(
    log_path: Path,
    *,
    decision: Decision,
    reviewer: str,
    run_id: str,
    sample_id: str,
    release_sha256: str,
    note: str = "",
    identity_source: str = ASSERTED,
    now: datetime | None = None,
) -> ReviewEntry:
    """Append one judgement, refusing to write onto a log that no longer verifies.

    Appending to a broken chain would produce a record that looks continuous and is not,
    which is worse than refusing: the point of the chain is that a reader can trust its
    shape.
    """
    reviewer = reviewer.strip()
    if not reviewer:
        raise ReviewError("a review must name a reviewer; an anonymous record attests nothing")

    existing = read_log(log_path)
    if existing:
        intact, detail = verify_chain(existing)
        if not intact:
            raise ReviewError(f"refusing to append to a log that does not verify: {detail}")

    previous = existing[-1].entry_sha256 if existing else None
    payload: dict[str, Any] = {
        "sequence": len(existing) + 1,
        "decision": Decision(decision).value,
        "reviewer": reviewer,
        "identity_source": identity_source,
        "recorded_at": (now or datetime.now(UTC)).isoformat(),
        "run_id": run_id,
        "sample_id": sample_id,
        "release_sha256": release_sha256,
        "note": note,
        "previous_entry_sha256": previous,
    }
    payload["entry_sha256"] = entry_digest(payload)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    _append_atomically(log_path, _canonical(payload) + "\n")
    return _to_entry(payload)


def _append_atomically(path: Path, line: str) -> None:
    """Append by rewriting through a temporary file in the same directory.

    A plain append can leave a torn final line if the process dies mid-write, and a torn
    line breaks the chain permanently. Rewriting costs nothing at this size and means the
    log is either the old content or the old content plus one complete entry.
    """
    existing = path.read_text(encoding="utf-8") if path.is_file() else ""
    handle, temporary = tempfile.mkstemp(dir=str(path.parent), prefix=".review-", suffix=".tmp")
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            stream.write(existing + line)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def current_state(
    entries: Sequence[ReviewEntry], release_sha256: str | None
) -> tuple[ReviewState, str]:
    """Resolve the trail and the content on disk into one state and a reason."""
    if not entries:
        return ReviewState.PENDING, "no review has been recorded"
    intact, detail = verify_chain(entries)
    if not intact:
        return ReviewState.BROKEN, detail
    latest = entries[-1]
    if release_sha256 is None:
        return ReviewState.STALE, (
            "the release bundle this review was made against is no longer present"
        )
    if latest.release_sha256 != release_sha256:
        return ReviewState.STALE, (
            f"the release bundle changed after {latest.describe()}; that judgement stands "
            "for what it saw and says nothing about the content here now"
        )
    if latest.decision is Decision.REJECTED:
        reason = latest.note or "no reason given"
        return ReviewState.REJECTED, f"rejected by {latest.reviewer}: {reason}"
    return ReviewState.ACCEPTED, f"accepted by {latest.reviewer} at {latest.recorded_at}"


def accepted_reviewers(
    entries: Iterable[ReviewEntry], release_sha256: str | None
) -> tuple[str, ...]:
    """Distinct reviewers who accepted exactly the content currently on disk.

    Used for a four-eyes check. Counting acceptances rather than reviewers would let one
    person satisfy it twice, and counting acceptances of *other* content would let a review
    of a superseded bundle carry over to this one.
    """
    if release_sha256 is None:
        return ()
    seen: dict[str, None] = {}
    for entry in entries:
        if entry.decision is Decision.ACCEPTED and entry.release_sha256 == release_sha256:
            seen.setdefault(entry.reviewer, None)
    return tuple(seen)


def exit_code(state: ReviewState, *, reviewers: int = 0, required_reviewers: int = 0) -> int:
    """Map a review state onto an exit code a release script can act on.

    0 nothing stands in the way, 2 something is definitely wrong, 6 something is unfinished.
    The split matches ``ontseq status``, so one convention covers both commands.
    """
    if state in {ReviewState.REJECTED, ReviewState.BROKEN, ReviewState.UNREADABLE}:
        return 2
    if state in {ReviewState.PENDING, ReviewState.STALE}:
        return 6
    if reviewers < required_reviewers:
        return 6
    return 0

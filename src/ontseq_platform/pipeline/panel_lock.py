"""The panel lock as a runtime contract, not a document nobody reads.

A target BED is a list of coordinates. Everything that decides whether those coordinates
mean anything lives in the lock beside it: which assembly they belong to, whether the
design carries flanks, and whether anybody has confirmed it against the panel the
sequencer actually selected on.

The runner fingerprints the BED, so editing it re-runs the stage that used it. Nothing
reads the lock, which means a run can measure per-target depth against an unconfirmed
design and emit a report that looks exactly as finished as one against a confirmed panel.
That is the gap this module closes:

* :func:`load_panel_lock` refuses a lock missing the fields a reader needs to judge it;
* :func:`verify_panel_bed` proves the BED on disk is the one the lock describes;
* :func:`check_declared_role` fails closed when a manifest claims a role the lock denies;
* :func:`panel_usage_warnings` states, in words, what an unconfirmed panel does not
  establish.

Nothing here interprets biology, sets a threshold or decides whether a target is
adequately covered. It compares what the lock says with what the run assumes, and refuses
when the two disagree. Only :mod:`ontseq_platform.pipeline.envelope` is imported, so the
logic is testable without pydantic or a reference genome.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

from .envelope import sha256_file

#: The only status that means a person confirmed the design against its source.
CONFIRMED_STATUS = "confirmed"

#: Suffix the generator appends to a target whose label and coordinates disagree.
REVIEW_SUFFIX = "_REVIEW_REQUIRED"

_HEX = frozenset("0123456789abcdef")


class PanelLockError(ValueError):
    """The lock is unusable, or it contradicts what the run assumes.

    Raised rather than warned deliberately. A panel whose provenance cannot be established
    is not a degraded input that yields a weaker answer; it is an input whose coordinates
    mean something unknown, and measuring depth over it produces a number nobody can read.
    """


@dataclass(frozen=True)
class PanelLock:
    """What a panel lock states about itself. No field is inferred."""

    path: Path
    panel_version: str
    status: str
    genome_build: str
    role: str
    bed_path: str
    bed_sha256: str
    target_count: int
    #: Empty when the lock does not say what kind of thing it counts.
    target_type: str
    #: ``None`` when the lock does not state it, which is not the same as zero.
    unique_target_labels: int | None
    #: ``None`` means no validated gene count is claimed. This is the honest value for a
    #: panel reproduced from a laboratory record without independent curation.
    validated_gene_count: int | None
    promotion_blockers: tuple[str, ...]
    open_question_targets: tuple[str, ...]

    @property
    def confirmed(self) -> bool:
        """True only when the lock says so *and* lists nothing standing in the way."""
        return self.status == CONFIRMED_STATUS and not self.promotion_blockers


def _mapping(value: object, *, source: Path, key: str) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise PanelLockError(f"{source.name}: {key!r} must be a mapping")
    return value


def _text(document: Mapping[str, object], key: str, *, source: Path) -> str:
    value = document.get(key)
    if not isinstance(value, str) or not value.strip():
        raise PanelLockError(f"{source.name}: {key!r} must be a non-empty string")
    return value


def _count(document: Mapping[str, object], key: str, *, source: Path) -> int:
    value = document.get(key)
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise PanelLockError(f"{source.name}: {key!r} must be a non-negative integer")
    return value


def _optional_count(document: Mapping[str, object], key: str, *, source: Path) -> int | None:
    if document.get(key) is None:
        return None
    return _count(document, key, source=source)


def _strings(document: Mapping[str, object], key: str, *, source: Path) -> tuple[str, ...]:
    value = document.get(key)
    if value is None:
        return ()
    if not isinstance(value, list):
        raise PanelLockError(f"{source.name}: {key!r} must be a list")
    collected: list[str] = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            raise PanelLockError(f"{source.name}: every {key!r} entry must be a non-empty string")
        collected.append(item)
    return tuple(collected)


def _open_question_targets(document: Mapping[str, object], *, source: Path) -> tuple[str, ...]:
    value = document.get("open_questions")
    if value is None:
        return ()
    if not isinstance(value, list):
        raise PanelLockError(f"{source.name}: 'open_questions' must be a list")
    collected: list[str] = []
    for item in value:
        if isinstance(item, str) and item.strip():
            collected.append(item)
            continue
        if isinstance(item, dict):
            target = item.get("target")
            if isinstance(target, str) and target.strip():
                collected.append(target)
                continue
        raise PanelLockError(f"{source.name}: an open question must name a target")
    return tuple(collected)


def load_panel_lock(path: Path) -> PanelLock:
    """Read a panel lock, refusing one that cannot be judged.

    Identity fields are mandatory: without a version, a status, a build, a role and a BED
    digest there is nothing to check a run against. The descriptive fields added later are
    optional and their absence is preserved as ``None`` or an empty tuple, so a caller can
    tell "the lock does not say" apart from "the lock says zero".
    """
    try:
        raw: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except OSError as error:
        raise PanelLockError(f"panel lock cannot be read: {path}") from error
    except yaml.YAMLError as error:
        raise PanelLockError(f"{path.name}: not valid YAML") from error

    document = _mapping(raw, source=path, key="panel lock")
    bed = _mapping(document.get("bed"), source=path, key="bed")
    digest = _text(bed, "sha256", source=path)
    if len(digest) != 64 or not set(digest) <= _HEX:
        raise PanelLockError(f"{path.name}: 'bed.sha256' is not a lowercase SHA-256 digest")

    return PanelLock(
        path=path,
        panel_version=_text(document, "panel_version", source=path),
        status=_text(document, "status", source=path),
        genome_build=_text(document, "genome_build", source=path),
        role=_text(document, "role", source=path),
        bed_path=_text(bed, "path", source=path),
        bed_sha256=digest,
        target_count=_count(bed, "target_count", source=path),
        target_type=str(bed.get("target_type") or ""),
        unique_target_labels=_optional_count(bed, "unique_target_labels", source=path),
        validated_gene_count=_optional_count(bed, "validated_gene_count", source=path),
        promotion_blockers=_strings(document, "promotion_blockers", source=path),
        open_question_targets=_open_question_targets(document, source=path),
    )


def target_labels(bed_path: Path) -> tuple[str, ...]:
    """Return the fourth column of every interval, in file order.

    Comment and blank lines are skipped. A row without a label is an error rather than an
    empty string: an unnamed interval cannot be reported against.
    """
    labels: list[str] = []
    for number, line in enumerate(bed_path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 4 or not fields[3].strip():
            raise PanelLockError(f"{bed_path.name}: line {number} has no target label")
        labels.append(fields[3])
    return tuple(labels)


def verify_panel_bed(lock: PanelLock, bed_path: Path) -> None:
    """Prove the BED on disk is the one the lock describes, or refuse the run.

    Content, not the path: a lock beside a BED that was regenerated, re-sorted or hand
    edited describes a design that no longer exists.
    """
    if not bed_path.is_file():
        raise PanelLockError(f"target BED is missing: {bed_path}")
    digest = sha256_file(bed_path)
    if digest != lock.bed_sha256:
        raise PanelLockError(
            f"{bed_path.name} is not the design {lock.path.name} describes: the lock records "
            f"{lock.bed_sha256[:12]}... and the file hashes to {digest[:12]}.... Refusing to "
            "measure coverage against a panel whose provenance cannot be established"
        )
    labels = target_labels(bed_path)
    if len(labels) != lock.target_count:
        raise PanelLockError(
            f"{bed_path.name} holds {len(labels)} interval(s) but {lock.path.name} records "
            f"{lock.target_count}"
        )
    if lock.unique_target_labels is not None and len(set(labels)) != lock.unique_target_labels:
        raise PanelLockError(
            f"{bed_path.name} holds {len(set(labels))} distinct label(s) but "
            f"{lock.path.name} records {lock.unique_target_labels}"
        )


def check_declared_role(lock: PanelLock, declared_role: str) -> None:
    """Refuse a manifest whose declared BED role contradicts the lock.

    This is the failure that produces a plausible wrong number rather than an error. Depth
    over a buffered selection panel includes flanking sequence; read as an unbuffered
    analysis ROI it overstates what was covered, and nothing downstream can detect it.
    """
    if declared_role == lock.role:
        return
    raise PanelLockError(
        f"{lock.path.name} records role {lock.role!r} but the run declares "
        f"{declared_role!r}. These are not interchangeable: per-target depth over a "
        "buffered selection panel includes flanking sequence and must not be read as "
        "coverage of an unbuffered analysis ROI"
    )


def panel_usage_warnings(lock: PanelLock, *, labels: Sequence[str] = ()) -> tuple[str, ...]:
    """State what this panel does not establish, in words a reviewer can act on.

    Returned as warnings rather than raised: a research run against an unconfirmed design
    is legitimate and is how a design becomes confirmed. What is not legitimate is a report
    that does not say so.
    """
    warnings: list[str] = []
    if not lock.confirmed:
        warnings.append(
            f"Panel {lock.panel_version} is recorded as {lock.status!r} in "
            f"{lock.path.name}. Per-target coverage from this run is technical evidence "
            "about an unconfirmed design, not a statement about a validated panel."
        )
    for blocker in lock.promotion_blockers:
        warnings.append(f"Panel promotion blocker: {blocker}")
    if lock.validated_gene_count is None:
        counted = lock.target_type or "targets"
        warnings.append(
            f"The lock records {lock.target_count} {counted} and claims no validated gene "
            "count. The interval labels are reproduced from the laboratory source and were "
            "not independently curated."
        )
    if "buffered" in lock.role:
        warnings.append(
            f"The design is buffered ({lock.role}). Per-target depth therefore includes "
            "flanking sequence and is not coverage of the gene body alone."
        )
    for target in lock.open_question_targets:
        warnings.append(
            f"Target {target} carries an unresolved label/coordinate contradiction and "
            f"appears as {target}{REVIEW_SUFFIX}. Do not read it as evidence about "
            f"{target} until the source is resolved."
        )
    flagged = sorted({label for label in labels if label.endswith(REVIEW_SUFFIX)})
    for label in flagged:
        if label[: -len(REVIEW_SUFFIX)] not in lock.open_question_targets:
            warnings.append(
                f"The BED contains {label}, which the lock does not list as an open "
                "question. The derivative and the lock disagree about this target."
            )
    return tuple(warnings)

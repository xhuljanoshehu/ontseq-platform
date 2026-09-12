"""Guideline criteria (ELN 2022, WHO 2022, ICC 2022) and the gate that keeps drafts out.

The criteria table in ``configs/knowledge_bundles/GUIDELINE_CRITERIA_DRAFT_v0`` was drafted
by a language model from memory rather than transcribed from the guideline documents. That
is a legitimate way to *start* — a haematologist corrects a structured draft far faster than
they write one from nothing — and an illegitimate way to *finish*, because a criteria table
is indistinguishable from a verified one once it is loaded into a report.

So the distinction is carried in the data and enforced here. Every record declares
``verification``. :func:`load_reportable_criteria` returns only records a named reviewer has
marked ``verified``; anything still ``unverified_model_draft`` is refused, loudly, with the
record identifiers named. Review tooling that deliberately wants the draft calls
:func:`load_for_review` instead, which is the only way to obtain unverified content and is
named so that no call site reaches it by accident.

The second gate is about the assay rather than the text. ``assay_status`` records whether
ONTSeq can evaluate a criterion at all. Small-variant criteria — NPM1, FLT3-ITD, CEBPA,
TP53, the myelodysplasia-related genes — cannot be evaluated today because no variant caller
is wired in, and complex/monosomal karyotype can only ever be a lower bound because balanced
rearrangements outside the panel are invisible. :func:`risk_group_determinable` exists so
that a caller has to ask, and so that "not determinable" is a value the report can carry
rather than an absence a reader fills in with "normal".
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

VERIFIED = "verified"
DRAFT = "unverified_model_draft"
REJECTED = "rejected"

COMPUTABLE = "computable_with_current_assay"
NEEDS_SMALL_VARIANTS = "requires_small_variant_calling_not_implemented"
LOWER_BOUND_ONLY = "requires_complete_karyotype_assay_cannot_guarantee"
NOT_OBSERVABLE = "not_a_sequencing_observable"


class GuidelineCriteriaError(ValueError):
    """Raised when a criteria bundle is unusable, or unusable *for the requested purpose*."""


@dataclass(frozen=True)
class Criterion:
    record_id: str
    category: str
    display_name: str
    pattern_type: str
    detectable_by: tuple[str, ...]
    assay_status: str
    verification: str
    guideline_reference: str | None
    reviewer_note: str
    caveat: str
    genes: tuple[str, ...] = ()
    iscn: str | None = None

    @property
    def evaluable_today(self) -> bool:
        return self.assay_status == COMPUTABLE


@dataclass(frozen=True)
class CriteriaBundle:
    bundle_id: str
    schema_version: str
    reviewer: str | None
    review_date: str | None
    criteria: tuple[Criterion, ...]

    def by_verification(self, status: str) -> tuple[Criterion, ...]:
        return tuple(item for item in self.criteria if item.verification == status)


def _strings(raw: dict[str, object], key: str) -> tuple[str, ...]:
    value = raw.get(key, ())
    if not isinstance(value, list | tuple):
        raise GuidelineCriteriaError(f"criteria record field {key!r} must be a list")
    return tuple(str(item) for item in value)


def _criterion(raw: dict[str, object]) -> Criterion:
    try:
        return Criterion(
            record_id=str(raw["record_id"]),
            category=str(raw["category"]),
            display_name=str(raw["display_name"]),
            pattern_type=str(raw["pattern_type"]),
            detectable_by=_strings(raw, "detectable_by"),
            assay_status=str(raw["assay_status"]),
            verification=str(raw["verification"]),
            guideline_reference=(
                str(raw["guideline_reference"]) if raw.get("guideline_reference") else None
            ),
            reviewer_note=str(raw.get("reviewer_note", "")),
            caveat=str(raw.get("caveat", "")),
            genes=_strings(raw, "genes"),
            iscn=str(raw["iscn"]) if raw.get("iscn") else None,
        )
    except KeyError as error:
        raise GuidelineCriteriaError(f"criteria record is missing {error}") from error


def load_for_review(path: Path) -> CriteriaBundle:
    """Load every record, verified or not. For review tooling only.

    Named explicitly rather than offered as a flag, so that a call site cannot obtain
    unverified criteria without saying in its own text that it meant to.
    """
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GuidelineCriteriaError(f"criteria bundle is unreadable: {path}") from error
    records = document.get("records")
    if not isinstance(records, list) or not records:
        raise GuidelineCriteriaError(f"criteria bundle holds no records: {path}")
    provenance = document.get("provenance", {})
    return CriteriaBundle(
        bundle_id=str(document.get("bundle_id", "UNKNOWN")),
        schema_version=str(document.get("schema_version", "0")),
        reviewer=(str(provenance.get("reviewer")) if provenance.get("reviewer") else None),
        review_date=(str(provenance.get("review_date")) if provenance.get("review_date") else None),
        criteria=tuple(_criterion(item) for item in records),
    )


def load_reportable_criteria(path: Path) -> CriteriaBundle:
    """Load only criteria a named reviewer has verified. Refuse otherwise.

    Refusing names the offending records instead of returning a quietly shortened list: a
    silently smaller criteria table produces a silently wrong classification, which is the
    failure this whole module exists to prevent.
    """
    bundle = load_for_review(path)
    drafts = bundle.by_verification(DRAFT)
    if drafts:
        names = ", ".join(item.record_id for item in drafts[:5])
        more = f" and {len(drafts) - 5} more" if len(drafts) > 5 else ""
        raise GuidelineCriteriaError(
            f"{len(drafts)} criteria in {bundle.bundle_id} are still an unverified model "
            f"draft and must not reach a report: {names}{more}. A haematologist has to "
            "review them and set verification to 'verified' or 'rejected' first."
        )
    if bundle.reviewer is None:
        raise GuidelineCriteriaError(
            f"{bundle.bundle_id} names no reviewer; a verified bundle must record who verified it"
        )
    return CriteriaBundle(
        bundle_id=bundle.bundle_id,
        schema_version=bundle.schema_version,
        reviewer=bundle.reviewer,
        review_date=bundle.review_date,
        criteria=bundle.by_verification(VERIFIED),
    )


@dataclass(frozen=True)
class Determinability:
    """Whether a risk group may be stated at all, and what is missing if not."""

    determinable: bool
    blocking: tuple[Criterion, ...]

    def reason(self) -> str:
        if self.determinable:
            return "Every criterion this guideline needs can be evaluated by this assay."
        missing = ", ".join(sorted({item.display_name for item in self.blocking}))
        return (
            f"Not determinable: {len(self.blocking)} criteria cannot be evaluated by this "
            f"assay ({missing}). Absence of these findings is not evidence of their absence."
        )


def risk_group_determinable(criteria: Sequence[Criterion]) -> Determinability:
    """Report whether the assay can evaluate every criterion, never guessing from the rest.

    A risk group derived from the subset the assay happens to see is not a partial answer,
    it is a wrong one: the unevaluated criteria are overwhelmingly the adverse ones, so the
    result would be biased towards favourable.
    """
    blocking = tuple(item for item in criteria if not item.evaluable_today)
    return Determinability(determinable=not blocking, blocking=blocking)

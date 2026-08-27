from __future__ import annotations

from collections.abc import Iterable

from .models import (
    EventType,
    GenomicEvent,
    SvEvidencePolicy,
    SvObservability,
    SvValidationStatus,
)

_CONFIDENCE_ORDER = {"high": 0, "moderate": 1, "low": 2, "unclassified": 3}


def _max_support(event: GenomicEvent) -> int:
    return max((item.support_reads or 0 for item in event.evidence), default=0)


def _best_vaf(event: GenomicEvent) -> float | None:
    values = [
        item.variant_allele_fraction
        for item in event.evidence
        if item.variant_allele_fraction is not None
    ]
    return max(values) if values else None


def _has_precise_evidence(event: GenomicEvent) -> bool:
    return any(item.precise is True for item in event.evidence)


def _caller_count(event: GenomicEvent) -> int:
    return len({item.caller.strip().lower() for item in event.evidence if item.caller.strip()})


def _technically_clean(event: GenomicEvent) -> bool:
    """Return whether no attached caller evidence carries a non-PASS filter."""
    return all(not (item.filters and item.filters != ["PASS"]) for item in event.evidence)


def _default_policy() -> SvEvidencePolicy:
    """Backward-compatible default for library callers; production loads the YAML policy."""
    return SvEvidencePolicy(
        profile_id="sv-evidence-priority-technical-v1",
        status="technical_defaults_only",
        note="Unvalidated built-in fallback; production runs must load the versioned policy.",
    )


def classify_sv_event(event: GenomicEvent, policy: SvEvidencePolicy | None = None) -> GenomicEvent:
    """Assign an unvalidated technical confidence tier to an SV candidate.

    ``reportable`` deliberately remains false until assay-specific acceptance criteria have
    been established against independent truth data.
    """
    resolved = policy or _default_policy()
    support = _max_support(event)
    callers = _caller_count(event)
    vaf = _best_vaf(event)
    precise = _has_precise_evidence(event)
    clean = _technically_clean(event)

    score = 0
    reasons: list[str] = []

    if callers >= 2:
        score += resolved.caller_consensus_weight
        reasons.append(f"caller_concordance={callers}")
    elif callers == 1:
        score += resolved.single_caller_weight
        reasons.append("single_caller")

    if support >= resolved.support_high:
        score += resolved.support_high_weight
        reasons.append(f"support={support}>={resolved.support_high}")
    elif support >= resolved.support_moderate:
        score += resolved.support_moderate_weight
        reasons.append(f"support={support}>={resolved.support_moderate}")
    elif support >= resolved.support_minimum:
        score += resolved.support_minimum_weight
        reasons.append(f"support={support}>={resolved.support_minimum}")

    if precise:
        score += resolved.precise_breakpoint_weight
        reasons.append("precise_breakpoint")

    if vaf is not None:
        if vaf >= resolved.vaf_high:
            score += resolved.vaf_high_weight
            reasons.append(f"vaf={vaf:.3f}>={resolved.vaf_high:.3f}")
        elif vaf >= resolved.vaf_minimum:
            score += resolved.vaf_minimum_weight
            reasons.append(f"vaf={vaf:.3f}>={resolved.vaf_minimum:.3f}")
        else:
            reasons.append(f"vaf={vaf:.3f}<{resolved.vaf_minimum:.3f}")

    if not clean:
        score += resolved.non_pass_filter_weight
        reasons.append("non_pass_caller_filter")

    if event.observability == SvObservability.OBSERVED_ADEQUATELY:
        score += resolved.adequate_observability_weight
        reasons.append("both_breakpoints_observed_adequately")
    elif event.observability in {
        SvObservability.INSUFFICIENT_COVERAGE,
        SvObservability.OUTSIDE_TARGET,
    }:
        score += resolved.inadequate_observability_weight
        reasons.append(f"observability={event.observability.value}")

    if event.technical_flags:
        flag_count = len(set(event.technical_flags))
        penalty = min(
            resolved.maximum_context_penalty,
            abs(resolved.context_flag_weight) * flag_count,
        )
        score -= penalty
        reasons.append(f"artifact_context_penalty={penalty}")

    if event.known_rearrangement:
        score += resolved.known_aml_pattern_weight
        reasons.append(f"known_aml_pattern={event.known_rearrangement}")

    if (
        event.event_type in {EventType.TRANSLOCATION, EventType.INVERSION}
        and support >= resolved.support_moderate
    ):
        score += resolved.balanced_sv_weight
        reasons.append("balanced_sv_review_priority")

    confidence = (
        "high"
        if score >= resolved.high_score
        else "moderate"
        if score >= resolved.moderate_score
        else "low"
    )
    notes = [
        note
        for note in event.notes
        if not note.startswith("Technical SV prioritization (")
        and not note.startswith("Technical confidence is not clinical validation;")
    ]
    notes.append(
        f"Technical SV prioritization ({resolved.profile_id}; {resolved.status}): "
        + ", ".join(reasons)
        + f"; score={score}."
    )
    notes.append(
        "Technical confidence is not clinical validation; reportable remains false until "
        "assay-specific benchmark criteria pass."
    )
    validation_status = event.validation_status
    if confidence in {"high", "moderate"} and validation_status == SvValidationStatus.DETECTED:
        validation_status = SvValidationStatus.TECHNICALLY_SUPPORTED
    return event.model_copy(
        update={
            "confidence": confidence,
            "validation_status": validation_status,
            "reportable": False,
            "notes": notes,
        }
    )


def prioritize_sv_events(
    events: Iterable[GenomicEvent], policy: SvEvidencePolicy | None = None
) -> list[GenomicEvent]:
    """Classify and sort SV candidates into a deterministic review order."""
    classified = [classify_sv_event(event, policy) for event in events]
    return sorted(
        classified,
        key=lambda event: (
            _CONFIDENCE_ORDER[event.confidence],
            0 if event.event_type == EventType.TRANSLOCATION else 1,
            -_max_support(event),
            event.primary.chromosome,
            event.primary.start,
            event.event_id,
        ),
    )


def sv_review_queue(
    events: Iterable[GenomicEvent],
    *,
    include_low: bool = False,
    limit: int = 50,
    policy: SvEvidencePolicy | None = None,
) -> list[GenomicEvent]:
    """Return the compact human-review queue while preserving all raw events elsewhere."""
    if limit < 1:
        raise ValueError("review queue limit must be at least 1")
    allowed = {"high", "moderate", "low"} if include_low else {"high", "moderate"}
    return [event for event in prioritize_sv_events(events, policy) if event.confidence in allowed][
        :limit
    ]

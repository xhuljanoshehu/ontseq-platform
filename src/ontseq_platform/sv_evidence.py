from __future__ import annotations

from collections.abc import Iterable

from .models import EventType, GenomicEvent


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
    """Return whether no attached caller evidence carries a non-PASS filter.

    Sniffles normalization currently accepts only PASS records, but this remains useful once
    evidence from additional callers is merged into one event.
    """
    for item in event.evidence:
        if item.filters and item.filters != ["PASS"]:
            return False
    return True


def classify_sv_event(event: GenomicEvent) -> GenomicEvent:
    """Assign an *unvalidated technical* confidence tier to an SV candidate.

    This is deliberately not a clinical validation rule. It reduces the raw-candidate review
    burden by combining caller concordance and read-level evidence already present in the
    normalized event. ``reportable`` remains false until assay-specific benchmark acceptance
    criteria have been established against independent truth data.
    """
    support = _max_support(event)
    callers = _caller_count(event)
    vaf = _best_vaf(event)
    precise = _has_precise_evidence(event)
    clean = _technically_clean(event)

    score = 0
    reasons: list[str] = []

    if callers >= 2:
        score += 4
        reasons.append(f"caller_concordance={callers}")
    elif callers == 1:
        score += 1
        reasons.append("single_caller")

    if support >= 20:
        score += 4
        reasons.append(f"support={support}>=20")
    elif support >= 10:
        score += 2
        reasons.append(f"support={support}>=10")
    elif support >= 5:
        score += 1
        reasons.append(f"support={support}>=5")

    if precise:
        score += 1
        reasons.append("precise_breakpoint")

    if vaf is not None:
        if vaf >= 0.10:
            score += 2
            reasons.append(f"vaf={vaf:.3f}>=0.10")
        elif vaf >= 0.05:
            score += 1
            reasons.append(f"vaf={vaf:.3f}>=0.05")
        else:
            reasons.append(f"vaf={vaf:.3f}<0.05")

    if not clean:
        score -= 3
        reasons.append("non_pass_caller_filter")

    # Balanced rearrangements are the main use case of the adaptive-sampling arm. Give them
    # review priority, but do not make them reportable merely because of event type.
    if event.event_type in {EventType.TRANSLOCATION, EventType.INVERSION} and support >= 10:
        score += 1
        reasons.append("balanced_sv_review_priority")

    if score >= 8:
        confidence = "high"
    elif score >= 5:
        confidence = "moderate"
    else:
        confidence = "low"

    notes = list(event.notes)
    notes.append(
        "Technical SV prioritization (unvalidated): " + ", ".join(reasons) + f"; score={score}."
    )
    notes.append(
        "Technical confidence is not clinical validation; reportable remains false until "
        "assay-specific benchmark criteria pass."
    )
    return event.model_copy(
        update={
            "confidence": confidence,
            "reportable": False,
            "notes": notes,
        }
    )


def prioritize_sv_events(events: Iterable[GenomicEvent]) -> list[GenomicEvent]:
    """Classify and sort SV candidates into a deterministic review order."""
    classified = [classify_sv_event(event) for event in events]
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
) -> list[GenomicEvent]:
    """Return the compact human-review queue while preserving all raw events elsewhere."""
    if limit < 1:
        raise ValueError("review queue limit must be at least 1")
    allowed = {"high", "moderate", "low"} if include_low else {"high", "moderate"}
    return [event for event in prioritize_sv_events(events) if event.confidence in allowed][:limit]

from __future__ import annotations

from collections.abc import Iterable

from .models import EventType, FusionSupportStatus, GenomicEvent

_CONFIDENCE_ORDER = {"high": 0, "moderate": 1, "low": 2, "unclassified": 3}


def is_structural_variant(event: GenomicEvent) -> bool:
    return event.event_type in {
        EventType.DELETION,
        EventType.DUPLICATION,
        EventType.INVERSION,
        EventType.TRANSLOCATION,
        EventType.INSERTION,
        EventType.FUSION,
    } and not (
        event.event_type in {EventType.DELETION, EventType.DUPLICATION}
        and event.copy_number is not None
    )


def release_state(event: GenomicEvent) -> str:
    """Human-readable release state without turning ``False`` into a biological verdict."""
    return "REPORTABLE" if event.reportable else "BENCHMARK_REQUIRED"


def fusion_assessment(event: GenomicEvent) -> str:
    """Describe the evidence actually produced by the current research pipeline."""
    if event.fusion_status == FusionSupportStatus.VALIDATED:
        return "VALIDATED"
    if event.fusion_status == FusionSupportStatus.SUPPORTED:
        return "SUPPORTED"
    if event.fusion_status == FusionSupportStatus.CANDIDATE:
        return "KNOWLEDGE_MATCH_CANDIDATE"
    if event.fusion_evidence is not None:
        return "BREAKPOINT_EVIDENCE"
    return "NOT_APPLICABLE"


def review_priority(event: GenomicEvent) -> str:
    if event.known_rearrangement:
        return "HEMATOLOGY_REVIEW"
    if event.fusion_evidence is not None and event.confidence in {"high", "moderate"}:
        return "FUSION_REVIEW"
    if event.confidence in {"high", "moderate"}:
        return "TECHNICAL_REVIEW"
    return "BACKGROUND"


def gene_pair_label(event: GenomicEvent) -> str:
    if event.known_rearrangement and "::" in event.known_rearrangement:
        return event.known_rearrangement
    if event.fusion_evidence is not None:
        genes = [event.fusion_evidence.gene_a.gene, event.fusion_evidence.gene_b.gene]
        resolved = [gene for gene in genes if gene]
        if len(resolved) == 2:
            return "::".join(resolved)
    genes = [event.primary.gene, event.secondary.gene if event.secondary else None]
    resolved = [gene for gene in genes if gene]
    return "::".join(resolved) if resolved else "unannotated"


def caller_count(event: GenomicEvent) -> int:
    return len({item.caller.strip().lower() for item in event.evidence if item.caller.strip()})


def maximum_support(event: GenomicEvent) -> int:
    return max((item.support_reads or 0 for item in event.evidence), default=0)


def pathology_label(event: GenomicEvent) -> str:
    """Render source-attributed disease associations without implying a diagnosis."""
    if not event.known_pathologies:
        return "no curated pathology association"
    return "; ".join(
        f"{pathology.name} ({pathology.disease_id})" for pathology in event.known_pathologies
    )


def key_findings(events: Iterable[GenomicEvent], *, limit: int = 25) -> list[GenomicEvent]:
    candidates = [
        event
        for event in events
        if event.known_rearrangement
        or (event.fusion_evidence is not None and event.confidence in {"high", "moderate"})
    ]
    candidates.sort(
        key=lambda event: (
            0 if event.known_rearrangement else 1,
            _CONFIDENCE_ORDER[event.confidence],
            -caller_count(event),
            -maximum_support(event),
            event.event_id,
        )
    )
    return candidates[:limit]


def fusion_review_events(events: Iterable[GenomicEvent]) -> list[GenomicEvent]:
    candidates = [
        event
        for event in events
        if event.event_type == EventType.FUSION
        or event.fusion_evidence is not None
        or event.known_rearrangement is not None
    ]
    return sorted(
        candidates,
        key=lambda event: (
            0 if event.known_rearrangement else 1,
            _CONFIDENCE_ORDER[event.confidence],
            -caller_count(event),
            -maximum_support(event),
            event.event_id,
        ),
    )

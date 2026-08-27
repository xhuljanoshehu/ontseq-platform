from __future__ import annotations

import pytest

from ontseq_platform.models import Evidence, EventType, GenomicEvent, Locus
from ontseq_platform.sv_evidence import classify_sv_event, prioritize_sv_events, sv_review_queue


def _event(
    event_id: str,
    *,
    event_type: EventType = EventType.DELETION,
    support: int = 5,
    vaf: float | None = None,
    precise: bool | None = None,
    caller: str = "Sniffles2",
    secondary: bool = False,
) -> GenomicEvent:
    return GenomicEvent(
        event_id=event_id,
        event_type=event_type,
        primary=Locus(chromosome="chr2", start=133_011_912, end=133_011_913),
        secondary=(
            Locus(chromosome="chr21", start=9_827_580, end=9_827_581)
            if secondary
            else None
        ),
        length_bp=None if event_type == EventType.TRANSLOCATION else 121,
        evidence=[
            Evidence(
                caller=caller,
                caller_version="2.8.0",
                support_reads=support,
                variant_allele_fraction=vaf,
                filters=["PASS"],
                precise=precise,
            )
        ],
    )


def test_b418_like_translocation_is_shortlisted_but_never_reportable() -> None:
    event = _event(
        "SNIFFLES2-B418",
        event_type=EventType.TRANSLOCATION,
        support=55,
        secondary=True,
    )

    classified = classify_sv_event(event)

    assert classified.confidence == "moderate"
    assert classified.reportable is False
    assert any("support=55>=20" in note for note in classified.notes)
    assert any("not clinical validation" in note.lower() for note in classified.notes)


def test_weak_single_caller_event_stays_low_priority() -> None:
    classified = classify_sv_event(_event("weak", support=5))

    assert classified.confidence == "low"
    assert classified.reportable is False


def test_multicaller_concordance_can_raise_technical_confidence() -> None:
    event = _event("consensus", support=12, vaf=0.15, precise=True)
    event = event.model_copy(
        update={
            "evidence": event.evidence
            + [
                Evidence(
                    caller="cuteSV",
                    caller_version="2.1.3",
                    support_reads=11,
                    variant_allele_fraction=0.14,
                    filters=["PASS"],
                    precise=True,
                )
            ]
        }
    )

    classified = classify_sv_event(event)

    assert classified.confidence == "high"
    assert classified.reportable is False


def test_review_queue_excludes_low_priority_by_default() -> None:
    weak = _event("weak", support=5)
    moderate = _event(
        "b418-like",
        event_type=EventType.TRANSLOCATION,
        support=30,
        secondary=True,
    )

    queue = sv_review_queue([weak, moderate])

    assert [event.event_id for event in queue] == ["b418-like"]


def test_prioritization_puts_high_before_moderate_and_low() -> None:
    low = _event("low", support=5)
    moderate = _event(
        "moderate",
        event_type=EventType.TRANSLOCATION,
        support=30,
        secondary=True,
    )
    high = _event("high", support=20, vaf=0.2, precise=True)

    ordered = prioritize_sv_events([low, moderate, high])

    assert [event.event_id for event in ordered] == ["high", "moderate", "low"]


def test_review_queue_limit_must_be_positive() -> None:
    with pytest.raises(ValueError, match="at least 1"):
        sv_review_queue([_event("one")], limit=0)

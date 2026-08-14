from __future__ import annotations

import re

from .models import EventType, GenomicEvent, ISCNProposal


def _chromosome_number(chromosome: str) -> int:
    value = chromosome.removeprefix("chr")
    return {"X": 23, "Y": 24}.get(value, int(value) if value.isdigit() else 99)


def _band_range(event: GenomicEvent) -> str | None:
    start = event.primary.cytoband_start
    end = event.primary.cytoband_end or start
    if not start:
        return None
    return start if start == end else f"{start}{end}"


def _fragment(event: GenomicEvent) -> str | None:
    chrom = event.primary.chromosome.removeprefix("chr")
    if event.event_type == EventType.CHROMOSOME_GAIN:
        return f"+{chrom}"
    if event.event_type == EventType.CHROMOSOME_LOSS:
        return f"-{chrom}"
    band = _band_range(event)
    if event.event_type == EventType.DELETION and band:
        return f"del({chrom})({band})"
    if event.event_type == EventType.DUPLICATION and band:
        return f"dup({chrom})({band})"
    if event.event_type == EventType.INVERSION and band:
        return f"inv({chrom})({band})"
    if event.event_type == EventType.TRANSLOCATION and event.secondary:
        other = event.secondary.chromosome.removeprefix("chr")
        band_a = event.primary.cytoband_start
        band_b = event.secondary.cytoband_start
        if band_a and band_b:
            return f"t({chrom};{other})({band_a};{band_b})"
    return None


def build_iscn_proposal(
    events: list[GenomicEvent], chromosome_count: int = 46, sex_chromosomes: str = "XX"
) -> ISCNProposal:
    """Render a deliberately limited, unvalidated ISCN-like proposal.

    The function supports a small auditable subset only. It does not claim full ISCN 2024
    semantic conformance and must not be used for automatic clinical release.
    """

    if not re.fullmatch(r"(?:XX|XY|X|XXY|XYY)", sex_chromosomes):
        raise ValueError("Unsupported sex chromosome complement for subset renderer")

    ordered = sorted(
        events,
        key=lambda item: (
            _chromosome_number(item.primary.chromosome),
            item.primary.start,
            item.event_type.value,
        ),
    )
    fragments: list[str] = []
    event_ids: list[str] = []
    warnings: list[str] = [
        "Automatically generated proposal; expert cytogenetic review is mandatory."
    ]

    for event in ordered:
        if not event.reportable:
            continue
        fragment = _fragment(event)
        if fragment:
            fragments.append(fragment)
            event_ids.append(event.event_id)
        else:
            warnings.append(
                f"Event {event.event_id} is outside the implemented ISCN subset and was omitted."
            )

    base = f"{chromosome_count},{sex_chromosomes}"
    notation = base if not fragments else f"{base},{','.join(fragments)}"
    return ISCNProposal(notation=notation, source_event_ids=event_ids, warnings=warnings)

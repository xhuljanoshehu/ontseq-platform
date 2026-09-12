from __future__ import annotations

from collections.abc import Iterable

from .models import (
    EventType,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    SvConsensusPolicy,
    SvConsensusReport,
    SvValidationStatus,
)


def _chromosome_key(chromosome: str) -> tuple[int, str]:
    label = chromosome.removeprefix("chr")
    if label.isdigit():
        return int(label), ""
    return {"X": 23, "Y": 24}.get(label, 99), label


def _same_chromosome(left: str, right: str) -> bool:
    return left.removeprefix("chr") == right.removeprefix("chr")


def _locus_key(locus: Locus) -> tuple[int, str, int, int]:
    rank, suffix = _chromosome_key(locus.chromosome)
    return rank, suffix, locus.start, locus.end


def _ordered_breakends(event: GenomicEvent) -> tuple[Locus, Locus]:
    assert event.secondary is not None
    pair = (event.primary, event.secondary)
    return pair if _locus_key(pair[0]) <= _locus_key(pair[1]) else (pair[1], pair[0])


def _distance(left: Locus, right: Locus) -> int:
    return max(abs(left.start - right.start), abs(left.end - right.end))


def _reciprocal_overlap(left: Locus, right: Locus) -> float:
    intersection = max(0, min(left.end, right.end) - max(left.start, right.start))
    if intersection == 0:
        return 0.0
    return min(intersection / (left.end - left.start), intersection / (right.end - right.start))


def _length_ratio_difference(left: GenomicEvent, right: GenomicEvent) -> float:
    if left.length_bp is None or right.length_bp is None:
        return 0.0
    return abs(left.length_bp - right.length_bp) / max(left.length_bp, right.length_bp)


def _orientations(event: GenomicEvent) -> set[str]:
    swapped = (
        event.event_type == EventType.TRANSLOCATION
        and event.secondary is not None
        and _ordered_breakends(event)[0] != event.primary
    )
    result: set[str] = set()
    for evidence in event.evidence:
        orientation = evidence.supporting_read_strands
        if orientation is None:
            continue
        result.add(orientation[::-1] if swapped and len(orientation) == 2 else orientation)
    return result


def events_match(
    left: GenomicEvent, right: GenomicEvent, policy: SvConsensusPolicy
) -> tuple[bool, int | None]:
    """Return whether two normalized calls represent the same candidate event."""
    if left.event_type != right.event_type:
        return False, None
    if policy.require_orientation_when_available:
        left_orientation = _orientations(left)
        right_orientation = _orientations(right)
        if (
            left_orientation
            and right_orientation
            and left_orientation.isdisjoint(right_orientation)
        ):
            return False, None

    if left.event_type == EventType.TRANSLOCATION:
        if left.secondary is None or right.secondary is None:
            return False, None
        left_a, left_b = _ordered_breakends(left)
        right_a, right_b = _ordered_breakends(right)
        if not _same_chromosome(left_a.chromosome, right_a.chromosome) or not _same_chromosome(
            left_b.chromosome, right_b.chromosome
        ):
            return False, None
        distance = max(_distance(left_a, right_a), _distance(left_b, right_b))
        return distance <= policy.maximum_breakpoint_distance_bp, distance

    if not _same_chromosome(left.primary.chromosome, right.primary.chromosome):
        return False, None
    distance = _distance(left.primary, right.primary)
    if distance > policy.maximum_breakpoint_distance_bp:
        return False, distance
    if _length_ratio_difference(left, right) > policy.maximum_length_ratio_difference:
        return False, distance
    if left.event_type == EventType.INSERTION:
        return True, distance
    overlap = _reciprocal_overlap(left.primary, right.primary)
    return overlap >= policy.minimum_reciprocal_overlap, distance


def _callers(event: GenomicEvent) -> set[str]:
    return {
        evidence.caller.strip().lower() for evidence in event.evidence if evidence.caller.strip()
    }


def _canonical_loci(event: GenomicEvent) -> tuple[Locus, Locus | None]:
    if event.event_type == EventType.TRANSLOCATION and event.secondary is not None:
        return _ordered_breakends(event)
    return event.primary, event.secondary


def _event_sort_key(event: GenomicEvent) -> tuple[object, ...]:
    primary, secondary = _canonical_loci(event)
    return (
        event.event_type.value,
        _locus_key(primary),
        _locus_key(secondary) if secondary is not None else (),
        event.event_id,
    )


def _merge_cluster(
    cluster: list[GenomicEvent], *, cluster_number: int, maximum_distance: int
) -> GenomicEvent:
    representative = sorted(
        cluster,
        key=lambda event: (
            -len(_callers(event)),
            -max((item.support_reads or 0 for item in event.evidence), default=0),
            event.event_id,
        ),
    )[0]
    primary, secondary = _canonical_loci(representative)
    evidence = [item for event in cluster for item in event.evidence]
    source_ids = sorted(
        {source for event in cluster for source in (event.source_event_ids or [event.event_id])}
    )
    callers = _callers(representative.model_copy(update={"evidence": evidence}))
    notes = [note for event in cluster for note in event.notes]
    notes.append(
        f"Consolidated {len(source_ids)} normalized representation(s) using an unvalidated "
        f"maximum breakpoint distance of {maximum_distance} bp."
    )
    if len(callers) >= 2:
        notes.append("Independent caller concordance is technical evidence, not ground truth.")
    return representative.model_copy(
        update={
            "event_id": f"SVCLUSTER-{cluster_number:06d}",
            "primary": primary,
            "secondary": secondary,
            "evidence": evidence,
            "source_event_ids": source_ids,
            "breakpoint_distance_bp": maximum_distance if len(cluster) > 1 else 0,
            "validation_status": (
                SvValidationStatus.TECHNICALLY_SUPPORTED
                if len(callers) >= 2
                else SvValidationStatus.DETECTED
            ),
            "reportable": False,
            "notes": notes,
        }
    )


def consolidate_sv_events(
    events: Iterable[GenomicEvent], policy: SvConsensusPolicy
) -> list[GenomicEvent]:
    """Greedily cluster equivalent normalized calls in deterministic genomic order."""
    ordered = sorted(events, key=_event_sort_key)
    clusters: list[list[GenomicEvent]] = []
    distances: list[int] = []
    for event in ordered:
        placed = False
        for index, cluster in enumerate(clusters):
            shared_caller = _callers(event).intersection(_callers(cluster[0]))
            if not policy.merge_within_caller and shared_caller:
                continue
            matches, distance = events_match(cluster[0], event, policy)
            if matches:
                cluster.append(event)
                distances[index] = max(distances[index], distance or 0)
                placed = True
                break
        if not placed:
            clusters.append([event])
            distances.append(0)
    return [
        _merge_cluster(cluster, cluster_number=index, maximum_distance=distances[index - 1])
        for index, cluster in enumerate(clusters, start=1)
    ]


def build_consensus_report(
    *,
    sample_id: str,
    genome_build: GenomeBuild,
    events: Iterable[GenomicEvent],
    policy: SvConsensusPolicy,
) -> SvConsensusReport:
    materialized = list(events)
    consolidated = consolidate_sv_events(materialized, policy)
    callers = sorted(
        {item.caller for event in materialized for item in event.evidence}, key=str.lower
    )
    return SvConsensusReport(
        sample_id=sample_id,
        genome_build=genome_build,
        status=ModuleRunStatus.COMPLETED if consolidated else ModuleRunStatus.NO_CALL,
        policy=policy,
        events=consolidated,
        input_event_count=len(materialized),
        consolidated_event_count=len(consolidated),
        caller_names=callers,
        warnings=[policy.note] if policy.status != "validated" else [],
        limitations=[
            "Caller consensus is evidence, not an independently validated truth set.",
            "Breakpoint matching thresholds are technical defaults until benchmark calibration.",
        ],
    )

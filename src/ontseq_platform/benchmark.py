from __future__ import annotations

from dataclasses import dataclass

from .models import (
    BenchmarkCase,
    BenchmarkKind,
    BenchmarkMatch,
    BenchmarkMetrics,
    BenchmarkReport,
    EventType,
    GenomicEvent,
    Locus,
)

CNV_TYPES = {
    EventType.CHROMOSOME_GAIN,
    EventType.CHROMOSOME_LOSS,
    EventType.DELETION,
    EventType.DUPLICATION,
}
SV_TYPES = {
    EventType.DELETION,
    EventType.DUPLICATION,
    EventType.INVERSION,
    EventType.INSERTION,
    EventType.TRANSLOCATION,
    EventType.FUSION,
}


@dataclass(frozen=True)
class CandidateMatch:
    truth_index: int
    query_index: int
    match: BenchmarkMatch


def _maximum_cardinality_matches(
    candidates: list[CandidateMatch], case: BenchmarkCase
) -> list[CandidateMatch]:
    """Select a deterministic maximum-cardinality bipartite matching.

    Candidate scores order preferences, while augmenting paths ensure a locally attractive
    edge cannot reduce the total number of true-positive pairs.
    """
    adjacency: dict[int, list[CandidateMatch]] = {}
    for candidate in candidates:
        adjacency.setdefault(candidate.truth_index, []).append(candidate)
    for options in adjacency.values():
        options.sort(
            key=lambda item: (
                -item.match.score,
                item.match.query_event_id,
            )
        )

    matched_by_query: dict[int, CandidateMatch] = {}

    def augment(truth_index: int, visited_queries: set[int]) -> bool:
        for candidate in adjacency.get(truth_index, []):
            if candidate.query_index in visited_queries:
                continue
            visited_queries.add(candidate.query_index)
            incumbent = matched_by_query.get(candidate.query_index)
            if incumbent is None or augment(incumbent.truth_index, visited_queries):
                matched_by_query[candidate.query_index] = candidate
                return True
        return False

    truth_order = sorted(
        adjacency,
        key=lambda index: (
            len(adjacency[index]),
            case.truth_events[index].event_id,
        ),
    )
    for truth_index in truth_order:
        augment(truth_index, set())
    return sorted(
        matched_by_query.values(),
        key=lambda item: (item.match.truth_event_id, item.match.query_event_id),
    )


def _chromosome(value: str) -> str:
    return value.removeprefix("chr")


def _reciprocal_overlap(left: Locus, right: Locus) -> float:
    if _chromosome(left.chromosome) != _chromosome(right.chromosome):
        return 0.0
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start))
    if overlap == 0:
        return 0.0
    return min(overlap / (left.end - left.start), overlap / (right.end - right.start))


def _maximum_distance(left: Locus, right: Locus) -> int | None:
    if _chromosome(left.chromosome) != _chromosome(right.chromosome):
        return None
    return max(abs(left.start - right.start), abs(left.end - right.end))


def _cnv_candidate(
    truth: GenomicEvent, query: GenomicEvent, case: BenchmarkCase
) -> BenchmarkMatch | None:
    if truth.event_type != query.event_type:
        return None
    overlap = _reciprocal_overlap(truth.primary, query.primary)
    if overlap < case.thresholds.minimum_reciprocal_overlap:
        return None
    tolerance = case.thresholds.copy_number_tolerance
    if tolerance is not None:
        if truth.copy_number is None or query.copy_number is None:
            return None
        if abs(truth.copy_number - query.copy_number) > tolerance:
            return None
    return BenchmarkMatch(
        truth_event_id=truth.event_id,
        query_event_id=query.event_id,
        score=overlap,
        reciprocal_overlap=overlap,
    )


def _paired_distance(truth: GenomicEvent, query: GenomicEvent) -> int | None:
    if truth.secondary is None or query.secondary is None:
        return None
    direct_left = _maximum_distance(truth.primary, query.primary)
    direct_right = _maximum_distance(truth.secondary, query.secondary)
    direct = None if direct_left is None or direct_right is None else max(direct_left, direct_right)
    swapped_left = _maximum_distance(truth.primary, query.secondary)
    swapped_right = _maximum_distance(truth.secondary, query.primary)
    swapped = (
        None if swapped_left is None or swapped_right is None else max(swapped_left, swapped_right)
    )
    available = [value for value in (direct, swapped) if value is not None]
    return min(available) if available else None


def _sv_candidate(
    truth: GenomicEvent, query: GenomicEvent, case: BenchmarkCase
) -> BenchmarkMatch | None:
    if truth.event_type != query.event_type:
        return None
    if truth.secondary is not None or query.secondary is not None:
        distance = _paired_distance(truth, query)
    else:
        distance = _maximum_distance(truth.primary, query.primary)
    maximum = case.thresholds.maximum_breakpoint_distance_bp
    if distance is None or distance > maximum:
        return None
    score = 1.0 if maximum == 0 else max(0.0, 1 - (distance / (maximum + 1)))
    return BenchmarkMatch(
        truth_event_id=truth.event_id,
        query_event_id=query.event_id,
        score=score,
        maximum_breakpoint_distance_bp=distance,
    )


def _ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _validate_case(case: BenchmarkCase) -> None:
    allowed = CNV_TYPES if case.kind == BenchmarkKind.CNV else SV_TYPES
    for label, events in (("truth", case.truth_events), ("query", case.query_events)):
        identifiers = [event.event_id for event in events]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError(f"{label} event IDs must be unique")
        unsupported = sorted(
            {event.event_type.value for event in events if event.event_type not in allowed}
        )
        if unsupported:
            raise ValueError(
                f"{label} contains event types unsupported for {case.kind.value}: "
                + ", ".join(unsupported)
            )


def benchmark_case(case: BenchmarkCase) -> BenchmarkReport:
    _validate_case(case)
    candidate_matches: list[CandidateMatch] = []
    matcher = _cnv_candidate if case.kind == BenchmarkKind.CNV else _sv_candidate
    for truth_index, truth in enumerate(case.truth_events):
        for query_index, query in enumerate(case.query_events):
            match = matcher(truth, query, case)
            if match is not None:
                candidate_matches.append(
                    CandidateMatch(
                        truth_index=truth_index,
                        query_index=query_index,
                        match=match,
                    )
                )
    selected = _maximum_cardinality_matches(candidate_matches, case)
    used_truth = {candidate.truth_index for candidate in selected}
    used_query = {candidate.query_index for candidate in selected}
    matches = [candidate.match for candidate in selected]

    true_positive = len(matches)
    false_positive = len(case.query_events) - true_positive
    false_negative = len(case.truth_events) - true_positive
    precision = _ratio(true_positive, true_positive + false_positive)
    recall = _ratio(true_positive, true_positive + false_negative)
    f1 = (
        None
        if precision is None or recall is None or precision + recall == 0
        else 2 * precision * recall / (precision + recall)
    )
    return BenchmarkReport(
        case_id=case.case_id,
        kind=case.kind,
        genome_build=case.genome_build,
        thresholds=case.thresholds,
        strata=case.strata,
        metrics=BenchmarkMetrics(
            true_positive=true_positive,
            false_positive=false_positive,
            false_negative=false_negative,
            precision=precision,
            recall=recall,
            f1=f1,
        ),
        matches=matches,
        unmatched_truth_event_ids=[
            event.event_id
            for index, event in enumerate(case.truth_events)
            if index not in used_truth
        ],
        unmatched_query_event_ids=[
            event.event_id
            for index, event in enumerate(case.query_events)
            if index not in used_query
        ],
        warnings=[
            "This normalized-event benchmark is an engineering control, not clinical validation.",
            "Results depend on the locked truth set, representation and matching thresholds.",
            "Use Truvari or another validated format-aware comparator for production "
            "VCF benchmarks.",
        ],
    )

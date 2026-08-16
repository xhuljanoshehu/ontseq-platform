"""Dependency-free half-open genomic interval algebra.

Every coordinate in this module is zero-based and half-open (``[start, end)``), matching
BED semantics. The module deliberately avoids pydantic, numpy and any other third-party
dependency so that the numerical core of the CNV subsystem stays independently testable
and free of contract-layer concerns.

An *interval set* is a mapping from contig name to a list of ``(start, end)`` tuples. A
set is *normalized* when, per contig, the tuples are sorted, non-empty, non-overlapping
and non-adjacent. Public functions return normalized sets and accept arbitrary input.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

Interval = tuple[int, int]
IntervalSet = dict[str, list[Interval]]


def canonical_contig(name: str) -> str:
    """Return a contig name without the ``chr`` prefix.

    Truth files, caller outputs and reference locks disagree about the prefix. Comparing
    ``chr7`` against ``7`` as different contigs would silently destroy recall, so every
    comparison in this package routes contig names through this function.
    """
    return name.removeprefix("chr")


def normalize(intervals: Iterable[Interval]) -> list[Interval]:
    """Sort, drop empty spans and merge overlapping or touching intervals."""
    ordered = sorted((start, end) for start, end in intervals if end > start)
    if not ordered:
        return []
    merged: list[Interval] = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            if end > last_end:
                merged[-1] = (last_start, end)
        else:
            merged.append((start, end))
    return merged


def normalize_set(source: Mapping[str, Sequence[Interval]]) -> IntervalSet:
    """Normalize every contig of an interval set, dropping contigs that become empty."""
    result: IntervalSet = {}
    for contig, spans in source.items():
        merged = normalize(spans)
        if merged:
            result[canonical_contig(contig)] = merged
    return result


def total_length(source: Mapping[str, Sequence[Interval]]) -> int:
    """Return the summed length of every interval, without normalizing first.

    Callers that need overlap-free totals must normalize beforehand. Keeping this
    function naive is intentional: double counting overlapping input must stay visible
    rather than being silently repaired.
    """
    return sum(end - start for spans in source.values() for start, end in spans)


def intersect(left: Sequence[Interval], right: Sequence[Interval]) -> list[Interval]:
    """Intersect two interval lists using a linear sweep over normalized input."""
    first = normalize(left)
    second = normalize(right)
    result: list[Interval] = []
    index = 0
    for start, end in first:
        while index < len(second) and second[index][1] <= start:
            index += 1
        cursor = index
        while cursor < len(second) and second[cursor][0] < end:
            overlap_start = max(start, second[cursor][0])
            overlap_end = min(end, second[cursor][1])
            if overlap_end > overlap_start:
                result.append((overlap_start, overlap_end))
            cursor += 1
    return result


def subtract(left: Sequence[Interval], right: Sequence[Interval]) -> list[Interval]:
    """Return the part of ``left`` that is not covered by ``right``."""
    remaining = normalize(left)
    removals = normalize(right)
    if not removals:
        return remaining
    result: list[Interval] = []
    for start, end in remaining:
        cursor = start
        for removal_start, removal_end in removals:
            if removal_end <= cursor:
                continue
            if removal_start >= end:
                break
            if removal_start > cursor:
                result.append((cursor, min(removal_start, end)))
            cursor = max(cursor, removal_end)
            if cursor >= end:
                break
        if cursor < end:
            result.append((cursor, end))
    return normalize(result)


def intersect_set(
    left: Mapping[str, Sequence[Interval]], right: Mapping[str, Sequence[Interval]]
) -> IntervalSet:
    """Intersect two interval sets contig by contig."""
    normalized_right = normalize_set(right)
    result: IntervalSet = {}
    for contig, spans in normalize_set(left).items():
        other = normalized_right.get(contig)
        if not other:
            continue
        overlap = intersect(spans, other)
        if overlap:
            result[contig] = overlap
    return result


def subtract_set(
    left: Mapping[str, Sequence[Interval]], right: Mapping[str, Sequence[Interval]]
) -> IntervalSet:
    """Remove every interval of ``right`` from ``left``, contig by contig."""
    normalized_right = normalize_set(right)
    result: IntervalSet = {}
    for contig, spans in normalize_set(left).items():
        remainder = subtract(spans, normalized_right.get(contig, []))
        if remainder:
            result[contig] = remainder
    return result


def union_set(*sources: Mapping[str, Sequence[Interval]]) -> IntervalSet:
    """Merge any number of interval sets into one normalized set."""
    combined: dict[str, list[Interval]] = {}
    for source in sources:
        for contig, spans in source.items():
            combined.setdefault(canonical_contig(contig), []).extend(spans)
    return normalize_set(combined)


def overlap_length(span: Interval, contig: str, other: Mapping[str, Sequence[Interval]]) -> int:
    """Return how many bases of a single span are covered by an interval set."""
    spans = other.get(canonical_contig(contig))
    if not spans:
        return 0
    return sum(end - start for start, end in intersect([span], spans))


def contig_lengths_to_set(lengths: Mapping[str, int]) -> IntervalSet:
    """Build a whole-genome interval set from a contig-length mapping."""
    return {
        canonical_contig(contig): [(0, length)] for contig, length in lengths.items() if length > 0
    }

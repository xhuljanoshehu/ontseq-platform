"""Maximum matching must be total and correct on dense/ambiguous synthetic call sets."""

from __future__ import annotations

import itertools
import random

from test_benchmark import _event

from ontseq_platform.benchmark import CandidateMatch, _maximum_cardinality_matches
from ontseq_platform.models import (
    BenchmarkCase,
    BenchmarkKind,
    BenchmarkMatch,
    EventType,
    GenomeBuild,
)


def _case(n: int) -> BenchmarkCase:
    return BenchmarkCase(
        case_id="SYNTHETIC_MATCHING",
        kind=BenchmarkKind.SV,
        genome_build=GenomeBuild.GRCH38,
        truth_events=[
            _event(f"T{i:05}", EventType.DELETION, "chr1", i * 100, i * 100 + 50) for i in range(n)
        ],
        query_events=[
            _event(f"Q{i:05}", EventType.DELETION, "chr1", i * 100, i * 100 + 50) for i in range(n)
        ],
    )


def _edges(case: BenchmarkCase, pairs: list[tuple[int, int]]) -> list[CandidateMatch]:
    return [
        CandidateMatch(
            t,
            q,
            BenchmarkMatch(
                truth_event_id=case.truth_events[t].event_id,
                query_event_id=case.query_events[q].event_id,
                score=1.0,
            ),
        )
        for t, q in pairs
    ]


def test_long_augmenting_path_does_not_depend_on_python_recursion_limit() -> None:
    n = 1200
    case = _case(n)
    pairs = [(t, q) for t in range(n - 1) for q in (t, t + 1)] + [(n - 1, 0), (n - 1, n - 1)]
    found = _maximum_cardinality_matches(_edges(case, pairs), case)
    assert len(found) == n
    assert len({x.truth_index for x in found}) == n
    assert len({x.query_index for x in found}) == n


def test_all_three_by_three_graphs_match_independent_exhaustive_oracle() -> None:
    case = _case(3)
    universe = list(itertools.product(range(3), repeat=2))
    for mask in range(1 << len(universe)):
        pairs = [pair for i, pair in enumerate(universe) if mask & (1 << i)]
        legal = set(pairs)
        expected = max(
            sum((t, q) in legal for t, q in enumerate(order))
            for order in itertools.permutations(range(3))
        )
        found = _maximum_cardinality_matches(_edges(case, pairs), case)
        assert len(found) == expected
        assert len({x.query_index for x in found}) == len(found)
        assert len({x.truth_index for x in found}) == len(found)
        assert all((x.truth_index, x.query_index) in legal for x in found)


def test_candidate_permutation_does_not_change_deterministic_match() -> None:
    case = _case(5)
    candidates = _edges(case, list(itertools.product(range(5), repeat=2)))
    expected = _maximum_cardinality_matches(candidates, case)
    random.Random(45).shuffle(candidates)
    assert _maximum_cardinality_matches(candidates, case) == expected

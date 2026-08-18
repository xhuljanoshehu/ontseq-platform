from __future__ import annotations

from collections import defaultdict

from pydantic import Field, model_validator

from .fusion import FusionInterpretationReport
from .models import StrictModel


class FusionRedundancyGroup(StrictModel):
    source_event_ids: list[str] = Field(min_length=2)
    unordered_breakpoint_pair: tuple[str, str]
    auto_deduplicated: bool = False

    @model_validator(mode="after")
    def event_ids_are_unique(self) -> FusionRedundancyGroup:
        if len(self.source_event_ids) != len(set(self.source_event_ids)):
            raise ValueError("redundancy group source event ids must be unique")
        return self


class FusionRedundancyReport(StrictModel):
    schema_version: str = "0.1.0"
    groups: list[FusionRedundancyGroup] = Field(default_factory=list)
    candidate_count: int = Field(ge=0)
    potentially_redundant_candidate_count: int = Field(ge=0)
    auto_deduplicated: bool = False
    warnings: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def counts_are_consistent(self) -> FusionRedundancyReport:
        grouped_ids = {
            event_id for group in self.groups for event_id in group.source_event_ids
        }
        if self.potentially_redundant_candidate_count != len(grouped_ids):
            raise ValueError("fusion redundancy candidate accounting is inconsistent")
        if self.potentially_redundant_candidate_count > self.candidate_count:
            raise ValueError("redundant candidate count cannot exceed candidate count")
        return self


def _canonical_chromosome(chromosome: str) -> str:
    return chromosome[3:] if chromosome.startswith("chr") else chromosome


def _locus_key(chromosome: str, position_0based: int) -> str:
    return f"{_canonical_chromosome(chromosome)}:{position_0based}"


def analyze_fusion_redundancy(
    report: FusionInterpretationReport,
) -> FusionRedundancyReport:
    """Flag duplicate/reciprocal breakpoint pairs without collapsing source evidence.

    Two candidates are grouped when they share the same unordered pair of exact 0-based
    breakpoint loci. This catches exact duplicates and reciprocal BND representations while
    deliberately avoiding heuristic breakpoint windows. No candidate is removed or merged.
    """

    grouped: dict[tuple[str, str], list[str]] = defaultdict(list)
    for candidate in report.candidates:
        first = _locus_key(candidate.primary.chromosome, candidate.primary.position_0based)
        second = _locus_key(candidate.secondary.chromosome, candidate.secondary.position_0based)
        key = tuple(sorted((first, second)))
        grouped[key].append(candidate.source_event_id)

    groups = [
        FusionRedundancyGroup(
            source_event_ids=event_ids,
            unordered_breakpoint_pair=key,
        )
        for key, event_ids in sorted(grouped.items())
        if len(event_ids) > 1
    ]
    redundant_ids = {
        event_id for group in groups for event_id in group.source_event_ids
    }
    warnings: list[str] = []
    if groups:
        warnings.append(
            "Multiple candidate records share the same exact unordered breakpoint pair. "
            "They are preserved as separate source evidence and must not be counted as "
            "independent fusion events without review or validated deduplication."
        )
    return FusionRedundancyReport(
        groups=groups,
        candidate_count=len(report.candidates),
        potentially_redundant_candidate_count=len(redundant_ids),
        warnings=warnings,
    )

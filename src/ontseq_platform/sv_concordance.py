from __future__ import annotations

from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from .models import EventType, Evidence, Locus, StrictModel


class SVConcordanceStatus(StrEnum):
    EXACT_MATCH = "exact_match"
    NEAR_MATCH = "near_match"
    TOPOLOGY_CONFLICT = "topology_conflict"


class SVCallerObservation(StrictModel):
    observation_id: str = Field(min_length=1)
    caller: str = Field(min_length=1)
    caller_version: str = Field(min_length=1)
    source_event_id: str = Field(min_length=1)
    event_type: EventType
    primary: Locus
    secondary: Locus | None = None
    evidence: Evidence
    research_only: Literal[True] = True
    reportable: Literal[False] = False

    @model_validator(mode="after")
    def paired_event_requirements_and_provenance(self) -> SVCallerObservation:
        paired_types = {EventType.TRANSLOCATION, EventType.FUSION}
        if self.event_type in paired_types and self.secondary is None:
            raise ValueError("paired SV observation requires a secondary locus")
        if self.evidence.caller != self.caller:
            raise ValueError("SV observation caller and evidence caller must match")
        if self.evidence.caller_version != self.caller_version:
            raise ValueError("SV observation caller version and evidence version must match")
        return self


class SVConcordancePolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    maximum_breakpoint_distance_bp: int = Field(ge=0)
    semantics: Literal["software_comparison_only"] = "software_comparison_only"
    clinically_validated: Literal[False] = False
    note: str = Field(min_length=1)


class SVConcordancePair(StrictModel):
    left_observation_id: str
    right_observation_id: str
    left_caller: str
    right_caller: str
    left_event_type: EventType
    right_event_type: EventType
    status: SVConcordanceStatus
    maximum_breakpoint_distance_bp: int = Field(ge=0)
    breakpoint_order_swapped: bool = False
    evidence_semantics: Literal["support_only_not_truth"] = "support_only_not_truth"
    research_only: Literal[True] = True
    reportable: Literal[False] = False

    @model_validator(mode="after")
    def status_is_consistent(self) -> SVConcordancePair:
        if self.status == SVConcordanceStatus.EXACT_MATCH:
            if self.left_event_type != self.right_event_type:
                raise ValueError("exact concordance requires the same event type")
            if self.maximum_breakpoint_distance_bp != 0:
                raise ValueError("exact concordance requires zero breakpoint distance")
        if self.status == SVConcordanceStatus.NEAR_MATCH:
            if self.left_event_type != self.right_event_type:
                raise ValueError("near concordance requires the same event type")
            if self.maximum_breakpoint_distance_bp == 0:
                raise ValueError("zero-distance concordance must be classified as exact")
        if self.status == SVConcordanceStatus.TOPOLOGY_CONFLICT:
            if self.left_event_type == self.right_event_type:
                raise ValueError("topology conflict requires different event types")
        return self


class SVConcordanceReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    policy: SVConcordancePolicy
    pairs: list[SVConcordancePair] = Field(default_factory=list)
    unmatched_left_observation_ids: list[str] = Field(default_factory=list)
    unmatched_right_observation_ids: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    conclusion_semantics: Literal["caller_concordance_is_not_truth"] = (
        "caller_concordance_is_not_truth"
    )
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def observations_are_not_reused(self) -> SVConcordanceReport:
        left_ids = [pair.left_observation_id for pair in self.pairs]
        right_ids = [pair.right_observation_id for pair in self.pairs]
        if len(left_ids) != len(set(left_ids)):
            raise ValueError("left caller observation is reused across concordance pairs")
        if len(right_ids) != len(set(right_ids)):
            raise ValueError("right caller observation is reused across concordance pairs")
        if set(left_ids) & set(self.unmatched_left_observation_ids):
            raise ValueError("matched left observation cannot also be unmatched")
        if set(right_ids) & set(self.unmatched_right_observation_ids):
            raise ValueError("matched right observation cannot also be unmatched")
        return self


class _GeometryMatch(StrictModel):
    maximum_distance_bp: int = Field(ge=0)
    swapped: bool = False


def _canonical_chromosome(chromosome: str) -> str:
    return chromosome[3:] if chromosome.startswith("chr") else chromosome


def _locus_distance(left: Locus, right: Locus) -> int | None:
    if _canonical_chromosome(left.chromosome) != _canonical_chromosome(right.chromosome):
        return None
    return max(abs(left.start - right.start), abs(left.end - right.end))


def _geometry_match(
    left: SVCallerObservation,
    right: SVCallerObservation,
) -> _GeometryMatch | None:
    left_paired = left.secondary is not None
    right_paired = right.secondary is not None
    if left_paired != right_paired:
        return None

    primary_distance = _locus_distance(left.primary, right.primary)
    if not left_paired:
        if primary_distance is None:
            return None
        return _GeometryMatch(maximum_distance_bp=primary_distance)

    if left.secondary is None or right.secondary is None:
        return None

    direct_secondary_distance = _locus_distance(left.secondary, right.secondary)
    direct: _GeometryMatch | None = None
    if primary_distance is not None and direct_secondary_distance is not None:
        direct = _GeometryMatch(
            maximum_distance_bp=max(primary_distance, direct_secondary_distance),
            swapped=False,
        )

    swapped_primary_distance = _locus_distance(left.primary, right.secondary)
    swapped_secondary_distance = _locus_distance(left.secondary, right.primary)
    swapped: _GeometryMatch | None = None
    if swapped_primary_distance is not None and swapped_secondary_distance is not None:
        swapped = _GeometryMatch(
            maximum_distance_bp=max(swapped_primary_distance, swapped_secondary_distance),
            swapped=True,
        )

    if direct is None:
        return swapped
    if swapped is None:
        return direct
    if swapped.maximum_distance_bp < direct.maximum_distance_bp:
        return swapped
    return direct


def _pair_status(
    left: SVCallerObservation,
    right: SVCallerObservation,
    geometry: _GeometryMatch,
) -> SVConcordanceStatus:
    if left.event_type != right.event_type:
        return SVConcordanceStatus.TOPOLOGY_CONFLICT
    if geometry.maximum_distance_bp == 0:
        return SVConcordanceStatus.EXACT_MATCH
    return SVConcordanceStatus.NEAR_MATCH


def _priority(status: SVConcordanceStatus) -> int:
    return {
        SVConcordanceStatus.EXACT_MATCH: 0,
        SVConcordanceStatus.NEAR_MATCH: 1,
        SVConcordanceStatus.TOPOLOGY_CONFLICT: 2,
    }[status]


def compare_sv_caller_observations(
    left: list[SVCallerObservation],
    right: list[SVCallerObservation],
    policy: SVConcordancePolicy,
) -> SVConcordanceReport:
    """Compare normalized caller observations without producing biological truth.

    Matching is one-to-one and uses only an explicit software-comparison breakpoint tolerance.
    Exact and near concordance are evidence labels, not confirmation. Different event types at
    concordant geometry are preserved as topology conflicts rather than coerced into agreement.
    """

    if len({item.observation_id for item in left}) != len(left):
        raise ValueError("left caller observation ids must be unique")
    if len({item.observation_id for item in right}) != len(right):
        raise ValueError("right caller observation ids must be unique")

    candidate_pairs: list[
        tuple[int, int, int, int, SVConcordanceStatus, _GeometryMatch]
    ] = []
    for left_index, left_item in enumerate(left):
        for right_index, right_item in enumerate(right):
            geometry = _geometry_match(left_item, right_item)
            if geometry is None:
                continue
            if geometry.maximum_distance_bp > policy.maximum_breakpoint_distance_bp:
                continue
            status = _pair_status(left_item, right_item, geometry)
            candidate_pairs.append(
                (
                    _priority(status),
                    geometry.maximum_distance_bp,
                    left_index,
                    right_index,
                    status,
                    geometry,
                )
            )

    candidate_pairs.sort(key=lambda item: item[:4])
    used_left: set[int] = set()
    used_right: set[int] = set()
    pairs: list[SVConcordancePair] = []
    for _, _, left_index, right_index, status, geometry in candidate_pairs:
        if left_index in used_left or right_index in used_right:
            continue
        left_item = left[left_index]
        right_item = right[right_index]
        used_left.add(left_index)
        used_right.add(right_index)
        pairs.append(
            SVConcordancePair(
                left_observation_id=left_item.observation_id,
                right_observation_id=right_item.observation_id,
                left_caller=left_item.caller,
                right_caller=right_item.caller,
                left_event_type=left_item.event_type,
                right_event_type=right_item.event_type,
                status=status,
                maximum_breakpoint_distance_bp=geometry.maximum_distance_bp,
                breakpoint_order_swapped=geometry.swapped,
            )
        )

    unmatched_left = [
        item.observation_id for index, item in enumerate(left) if index not in used_left
    ]
    unmatched_right = [
        item.observation_id for index, item in enumerate(right) if index not in used_right
    ]

    warnings = [
        "Caller concordance is software evidence only and does not establish biological truth, "
        "functional fusion status, clinical validity or reportability."
    ]
    if any(pair.status == SVConcordanceStatus.TOPOLOGY_CONFLICT for pair in pairs):
        warnings.append(
            "At least one geometrically concordant caller pair has conflicting event types; "
            "the conflict is preserved for review."
        )

    return SVConcordanceReport(
        policy=policy,
        pairs=pairs,
        unmatched_left_observation_ids=unmatched_left,
        unmatched_right_observation_ids=unmatched_right,
        warnings=warnings,
    )

from __future__ import annotations

from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from .fusion import (
    FusionBreakpoint,
    FusionClassification,
    FusionGenePair,
    FusionInterpretationReport,
    ObservabilityStatus,
)
from .fusion_redundancy import FusionRedundancyReport, analyze_fusion_redundancy
from .models import GenomeBuild, ModuleRunStatus, ReviewStatus, StrictModel


class FusionReviewDisposition(StrEnum):
    CANDIDATE_EVIDENCE = "candidate_evidence"
    NO_CALL = "no_call"
    NOT_RUN = "not_run"
    FAILED = "failed"


class FusionReviewerCallerEvidence(StrictModel):
    caller: str = Field(min_length=1)
    caller_version: str = Field(min_length=1)
    support_reads: int | None = Field(default=None, ge=0)
    local_coverage: float | None = Field(default=None, ge=0)
    variant_allele_fraction: float | None = Field(default=None, ge=0, le=1)
    quality: float | None = Field(default=None, ge=0)
    filters: list[str] = Field(default_factory=list)
    precise: bool | None = None


class FusionReviewerCandidate(StrictModel):
    candidate_id: str = Field(min_length=1)
    source_event_id: str = Field(min_length=1)
    classification: FusionClassification
    primary: FusionBreakpoint
    secondary: FusionBreakpoint
    gene_pairs: list[FusionGenePair] = Field(default_factory=list)
    caller_evidence: list[FusionReviewerCallerEvidence] = Field(default_factory=list)
    known_pair_present: bool
    orientation_evidence_available: bool
    transcript_orientation_resolved: bool
    both_breakpoints_observable: bool
    potentially_redundant: bool = False
    redundancy_group_event_ids: list[str] = Field(default_factory=list)
    review_status: ReviewStatus = ReviewStatus.REVIEW_REQUIRED
    research_only: Literal[True] = True
    reportable: Literal[False] = False
    limitations: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def summary_fields_match_candidate(self) -> FusionReviewerCandidate:
        if self.known_pair_present != any(pair.known_pair for pair in self.gene_pairs):
            raise ValueError("reviewer known-pair summary is inconsistent")
        resolved = bool(self.gene_pairs) and all(
            pair.orientation_resolved for pair in self.gene_pairs
        )
        if self.transcript_orientation_resolved != resolved:
            raise ValueError("reviewer transcript-orientation summary is inconsistent")
        observable = (
            self.primary.observability == ObservabilityStatus.OBSERVABLE
            and self.secondary.observability == ObservabilityStatus.OBSERVABLE
        )
        if self.both_breakpoints_observable != observable:
            raise ValueError("reviewer breakpoint-observability summary is inconsistent")
        if self.potentially_redundant:
            if len(self.redundancy_group_event_ids) < 2:
                raise ValueError("redundant reviewer candidate requires a redundancy group")
            if self.source_event_id not in self.redundancy_group_event_ids:
                raise ValueError("reviewer candidate is absent from its redundancy group")
        elif self.redundancy_group_event_ids:
            raise ValueError("non-redundant reviewer candidate must not carry a redundancy group")
        return self


class FusionReviewerReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    source_status: ModuleRunStatus
    disposition: FusionReviewDisposition
    review_status: ReviewStatus = ReviewStatus.REVIEW_REQUIRED
    candidates: list[FusionReviewerCandidate] = Field(default_factory=list)
    source_translocation_count: int = Field(ge=0)
    unresolved_source_event_ids: list[str] = Field(default_factory=list)
    missing_breakend_descriptor_event_ids: list[str] = Field(default_factory=list)
    redundancy_group_count: int = Field(default=0, ge=0)
    absence_interpretation: Literal["not_established"] = "not_established"
    warnings: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def disposition_matches_source_status(self) -> FusionReviewerReport:
        expected = {
            ModuleRunStatus.COMPLETED: FusionReviewDisposition.CANDIDATE_EVIDENCE,
            ModuleRunStatus.NO_CALL: FusionReviewDisposition.NO_CALL,
            ModuleRunStatus.NOT_RUN: FusionReviewDisposition.NOT_RUN,
            ModuleRunStatus.FAILED: FusionReviewDisposition.FAILED,
        }[self.source_status]
        if self.disposition != expected:
            raise ValueError("fusion reviewer disposition does not match source status")
        if self.disposition == FusionReviewDisposition.CANDIDATE_EVIDENCE:
            if not self.candidates:
                raise ValueError("completed fusion review requires candidate evidence")
        elif self.candidates:
            raise ValueError("non-completed fusion review must not contain candidates")
        if len({item.candidate_id for item in self.candidates}) != len(self.candidates):
            raise ValueError("fusion reviewer candidate ids must be unique")
        if len({item.source_event_id for item in self.candidates}) != len(self.candidates):
            raise ValueError("fusion reviewer source event ids must be unique")
        return self


def _caller_evidence(candidate_evidence: list[object]) -> list[FusionReviewerCallerEvidence]:
    result: list[FusionReviewerCallerEvidence] = []
    for evidence in candidate_evidence:
        caller = getattr(evidence, "caller")
        caller_version = getattr(evidence, "caller_version")
        result.append(
            FusionReviewerCallerEvidence(
                caller=caller,
                caller_version=caller_version,
                support_reads=getattr(evidence, "support_reads"),
                local_coverage=getattr(evidence, "local_coverage"),
                variant_allele_fraction=getattr(evidence, "variant_allele_fraction"),
                quality=getattr(evidence, "quality"),
                filters=list(getattr(evidence, "filters")),
                precise=getattr(evidence, "precise"),
            )
        )
    return result


def _disposition(status: ModuleRunStatus) -> FusionReviewDisposition:
    return {
        ModuleRunStatus.COMPLETED: FusionReviewDisposition.CANDIDATE_EVIDENCE,
        ModuleRunStatus.NO_CALL: FusionReviewDisposition.NO_CALL,
        ModuleRunStatus.NOT_RUN: FusionReviewDisposition.NOT_RUN,
        ModuleRunStatus.FAILED: FusionReviewDisposition.FAILED,
    }[status]


def _unique_messages(messages: list[str]) -> list[str]:
    return list(dict.fromkeys(messages))


def build_fusion_reviewer_report(
    report: FusionInterpretationReport,
    *,
    redundancy: FusionRedundancyReport | None = None,
) -> FusionReviewerReport:
    """Build a privacy-safe reviewer view without creating a clinical conclusion.

    The reviewer contract separates source status, candidate evidence, observability and
    redundancy. It never converts caller agreement, a known gene pair or the absence of a
    candidate into clinical truth.
    """

    redundancy_report = redundancy or analyze_fusion_redundancy(report)
    if redundancy_report.candidate_count != len(report.candidates):
        raise ValueError("fusion redundancy report does not match candidate count")

    candidate_event_ids = {candidate.source_event_id for candidate in report.candidates}
    redundancy_by_event: dict[str, list[str]] = {}
    for group in redundancy_report.groups:
        if any(event_id not in candidate_event_ids for event_id in group.source_event_ids):
            raise ValueError("fusion redundancy report contains an unknown source event")
        for event_id in group.source_event_ids:
            redundancy_by_event[event_id] = list(group.source_event_ids)

    reviewer_candidates: list[FusionReviewerCandidate] = []
    for candidate in report.candidates:
        redundancy_group = redundancy_by_event.get(candidate.source_event_id, [])
        reviewer_candidates.append(
            FusionReviewerCandidate(
                candidate_id=candidate.candidate_id,
                source_event_id=candidate.source_event_id,
                classification=candidate.classification,
                primary=candidate.primary,
                secondary=candidate.secondary,
                gene_pairs=candidate.gene_pairs,
                caller_evidence=_caller_evidence(candidate.evidence),
                known_pair_present=any(pair.known_pair for pair in candidate.gene_pairs),
                orientation_evidence_available=candidate.breakend_descriptor is not None,
                transcript_orientation_resolved=bool(candidate.gene_pairs)
                and all(pair.orientation_resolved for pair in candidate.gene_pairs),
                both_breakpoints_observable=(
                    candidate.primary.observability == ObservabilityStatus.OBSERVABLE
                    and candidate.secondary.observability == ObservabilityStatus.OBSERVABLE
                ),
                potentially_redundant=bool(redundancy_group),
                redundancy_group_event_ids=redundancy_group,
                limitations=candidate.limitations,
            )
        )

    warnings = [*report.warnings, *redundancy_report.warnings]
    if report.status == ModuleRunStatus.NO_CALL:
        warnings.append(
            "NO_CALL means that the fusion module did not emit assessable candidate evidence; "
            "it is not a validated biological negative result."
        )
    elif report.status == ModuleRunStatus.FAILED:
        warnings.append(
            "The fusion module failed; absence of candidate evidence is uninterpretable."
        )
    elif report.status == ModuleRunStatus.NOT_RUN:
        warnings.append(
            "The fusion module was not run; no inference about fusion status is permitted."
        )

    return FusionReviewerReport(
        sample_id=report.sample_id,
        genome_build=report.genome_build,
        source_status=report.status,
        disposition=_disposition(report.status),
        candidates=reviewer_candidates,
        source_translocation_count=report.source_translocation_count,
        unresolved_source_event_ids=report.unresolved_source_event_ids,
        missing_breakend_descriptor_event_ids=report.missing_breakend_descriptor_event_ids,
        redundancy_group_count=len(redundancy_report.groups),
        warnings=_unique_messages(warnings),
    )

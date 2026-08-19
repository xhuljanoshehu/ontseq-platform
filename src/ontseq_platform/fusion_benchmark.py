from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Literal

import yaml
from pydantic import Field, model_validator

from .breakends import SnifflesJunctionOrientation
from .fusion import FusionClassification, FusionInterpretationReport, ObservabilityStatus
from .models import StrictModel


class PartnerExpectation(StrEnum):
    EXACT_PAIR = "exact_pair"
    ANY_PARTNER = "any_partner"
    NONE = "none"


class SyntheticFusionFixture(StrictModel):
    """Expected software behavior for a deliberately non-biological synthetic case.

    Family labels such as BCR::ABL1 are mnemonic benchmark semantics only. This contract
    deliberately contains no real genomic coordinates and must not be used as a clinical
    knowledge base or reportability whitelist.
    """

    fixture_id: str = Field(pattern=r"^[a-z0-9][a-z0-9_.-]+$")
    family_label: str = Field(min_length=1)
    intended_use: Literal["software_benchmark_only"] = "software_benchmark_only"
    coordinate_semantics: Literal["synthetic_nonbiological"] = "synthetic_nonbiological"
    clinical_truth: Literal[False] = False
    expected_classification: FusionClassification
    partner_expectation: PartnerExpectation = PartnerExpectation.NONE
    expected_gene_pair: tuple[str, str] | None = None
    anchor_gene: str | None = None
    expected_primary_observability: ObservabilityStatus | None = None
    expected_secondary_observability: ObservabilityStatus | None = None
    accepted_junction_orientations: list[SnifflesJunctionOrientation] = Field(default_factory=list)
    expected_candidate_count: int = Field(default=1, ge=0)
    expected_reportable: Literal[False] = False
    expected_research_only: Literal[True] = True
    expected_transcript_orientation_resolved: Literal[False] = False
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def partner_contract_is_consistent(self) -> SyntheticFusionFixture:
        if self.partner_expectation == PartnerExpectation.EXACT_PAIR:
            if self.expected_gene_pair is None or self.anchor_gene is not None:
                raise ValueError(
                    "exact-pair fixtures require expected_gene_pair and no anchor_gene"
                )
            first, second = self.expected_gene_pair
            if not first or not second or first.upper() == second.upper():
                raise ValueError("expected_gene_pair requires two distinct non-empty genes")
        elif self.partner_expectation == PartnerExpectation.ANY_PARTNER:
            if not self.anchor_gene or self.expected_gene_pair is not None:
                raise ValueError("any-partner fixtures require anchor_gene and no exact pair")
        elif self.expected_gene_pair is not None or self.anchor_gene is not None:
            raise ValueError("partner-free fixtures must not define gene expectations")
        return self


class SyntheticFusionBenchmarkSuite(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    suite_id: str = Field(pattern=r"^[a-z0-9][a-z0-9_.-]+$")
    status: Literal["synthetic_software_benchmark_only"] = "synthetic_software_benchmark_only"
    description: str = Field(min_length=1)
    fixtures: list[SyntheticFusionFixture] = Field(min_length=1)

    @model_validator(mode="after")
    def fixture_ids_are_unique(self) -> SyntheticFusionBenchmarkSuite:
        ids = [fixture.fixture_id for fixture in self.fixtures]
        if len(ids) != len(set(ids)):
            raise ValueError("synthetic fusion benchmark fixture ids must be unique")
        return self


class FusionBenchmarkCaseResult(StrictModel):
    fixture_id: str
    passed: bool
    failures: list[str] = Field(default_factory=list)


def load_synthetic_fusion_benchmark(path: Path) -> SyntheticFusionBenchmarkSuite:
    if not path.is_file():
        raise ValueError("synthetic fusion benchmark YAML is missing or unreadable")
    with path.open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError("synthetic fusion benchmark YAML must contain a mapping")
    return SyntheticFusionBenchmarkSuite.model_validate(payload)


def _canonical_pair(first: str, second: str) -> tuple[str, str]:
    return tuple(sorted((first.upper(), second.upper())))  # type: ignore[return-value]


def evaluate_fusion_benchmark_case(
    fixture: SyntheticFusionFixture,
    report: FusionInterpretationReport,
) -> FusionBenchmarkCaseResult:
    """Evaluate software behavior only; never infer analytical or clinical validity."""

    failures: list[str] = []
    candidates = report.candidates
    if len(candidates) != fixture.expected_candidate_count:
        failures.append(
            "candidate_count "
            f"expected={fixture.expected_candidate_count} observed={len(candidates)}"
        )

    for candidate in candidates:
        if candidate.reportable is not fixture.expected_reportable:
            failures.append(f"{candidate.candidate_id}: reportable state mismatch")
        if candidate.research_only is not fixture.expected_research_only:
            failures.append(f"{candidate.candidate_id}: research_only state mismatch")
        if candidate.classification != fixture.expected_classification:
            failures.append(
                f"{candidate.candidate_id}: classification expected="
                f"{fixture.expected_classification.value} observed={candidate.classification.value}"
            )
        if fixture.expected_primary_observability is not None and (
            candidate.primary.observability != fixture.expected_primary_observability
        ):
            failures.append(f"{candidate.candidate_id}: primary observability mismatch")
        if fixture.expected_secondary_observability is not None and (
            candidate.secondary.observability != fixture.expected_secondary_observability
        ):
            failures.append(f"{candidate.candidate_id}: secondary observability mismatch")
        if any(
            pair.orientation_resolved != fixture.expected_transcript_orientation_resolved
            for pair in candidate.gene_pairs
        ):
            failures.append(f"{candidate.candidate_id}: transcript orientation state mismatch")

        observed_pairs = {
            _canonical_pair(pair.gene_a, pair.gene_b) for pair in candidate.gene_pairs
        }
        if fixture.partner_expectation == PartnerExpectation.EXACT_PAIR:
            assert fixture.expected_gene_pair is not None
            expected_pair = _canonical_pair(*fixture.expected_gene_pair)
            if expected_pair not in observed_pairs:
                failures.append(f"{candidate.candidate_id}: expected gene pair not observed")
        elif fixture.partner_expectation == PartnerExpectation.ANY_PARTNER:
            assert fixture.anchor_gene is not None
            anchor = fixture.anchor_gene.upper()
            if not any(anchor in pair for pair in observed_pairs):
                failures.append(f"{candidate.candidate_id}: anchor gene not observed")

        if fixture.accepted_junction_orientations:
            descriptor = candidate.breakend_descriptor
            if descriptor is None:
                failures.append(f"{candidate.candidate_id}: BND descriptor missing")
            elif (
                descriptor.sniffles_junction_orientation
                not in fixture.accepted_junction_orientations
            ):
                failures.append(f"{candidate.candidate_id}: junction orientation not accepted")

    return FusionBenchmarkCaseResult(
        fixture_id=fixture.fixture_id,
        passed=not failures,
        failures=failures,
    )

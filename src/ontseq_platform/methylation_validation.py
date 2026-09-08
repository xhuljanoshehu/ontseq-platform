"""Prospective, donor-independent methylation read-group recovery contracts.

No result in this module estimates DNA mass, cellularity, tumour purity or an LoD.
Registration is an operator declaration and content lock, not a trusted timestamp.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
import unicodedata
from collections.abc import Mapping, Sequence
from datetime import datetime
from enum import StrEnum
from typing import Any, Literal, cast

from pydantic import Field, field_validator, model_validator

from .methylation_mixture import (
    MethylationMixtureLevelResult,
    MethylationMixturePolicy,
    MethylationReferenceMarker,
    _Counts,
    _draw_estimate,
    _point_estimate,
    _quantile,
    _reference_marker_digest,
)
from .models import GenomeBuild, ModuleRunStatus, StrictModel

SHA256 = r"^[0-9a-f]{64}$"
SampleRole = Literal["calibration_a", "calibration_b", "test_a", "test_b"]
SAMPLE_ROLES: tuple[SampleRole, ...] = ("calibration_a", "calibration_b", "test_a", "test_b")
_UNKNOWN = frozenset({"", "unknown", "pending", "unavailable", "not_reported", "na", "n/a"})


def _is_unknown_metadata(value: str) -> bool:
    return value.strip().casefold() in _UNKNOWN


def content_sha256(value: StrictModel | Mapping[str, Any]) -> str:
    """Hash canonical JSON content, independent of whitespace and key order."""
    payload = value.model_dump(mode="json") if isinstance(value, StrictModel) else dict(value)
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


class ValidationDecision(StrEnum):
    PASS = "PASS"
    FAIL = "FAIL"
    NO_CALL = "NO_CALL"


class RecoveryAcceptance(StrictModel):
    read_group_budget: int = Field(ge=2)
    maximum_absolute_bias: float = Field(ge=0, le=1)
    maximum_mae: float = Field(ge=0, le=1)
    maximum_rmse: float = Field(ge=0, le=1)
    minimum_conditional_interval_coverage: float = Field(ge=0, le=1)
    maximum_no_call_rate: float = Field(ge=0, lt=1)
    maximum_per_fraction_absolute_bias: float = Field(ge=0, le=1)
    minimum_per_fraction_conditional_interval_coverage: float = Field(ge=0, le=1)
    maximum_per_fraction_no_call_rate: float = Field(ge=0, lt=1)


class MethylationValidationMatrix(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    matrix_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$")
    status: Literal["research_acceptance_candidate_unvalidated"] = (
        "research_acceptance_candidate_unvalidated"
    )
    source_a_fractions: list[float] = Field(min_length=3)
    read_group_budgets: list[int] = Field(min_length=1)
    replicate_seeds: list[int] = Field(min_length=10)
    budget_unit: Literal["retained_whole_read_groups_not_genomic_coverage"] = (
        "retained_whole_read_groups_not_genomic_coverage"
    )
    acceptance: list[RecoveryAcceptance] = Field(min_length=1)
    estimator_policy: MethylationMixturePolicy
    rationale: str = Field(min_length=20)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_grid(self) -> MethylationValidationMatrix:
        if self.source_a_fractions != sorted(set(self.source_a_fractions)):
            raise ValueError("Fractions must be unique and sorted")
        if self.source_a_fractions[0] != 0 or self.source_a_fractions[-1] != 1:
            raise ValueError("Validation includes both pure-source endpoints")
        if any(not math.isfinite(x) or x < 0 or x > 1 for x in self.source_a_fractions):
            raise ValueError("Fractions must be finite probabilities")
        if self.read_group_budgets != sorted(set(self.read_group_budgets)):
            raise ValueError("Read-group budgets must be unique and sorted")
        if any(x < 2 for x in self.read_group_budgets):
            raise ValueError("Read-group budgets must be at least two")
        if len(set(self.replicate_seeds)) != len(self.replicate_seeds) or any(
            x < 0 or x > 2_147_483_647 for x in self.replicate_seeds
        ):
            raise ValueError("Replicate seeds must be distinct nonnegative 32-bit integers")
        if [x.read_group_budget for x in self.acceptance] != self.read_group_budgets:
            raise ValueError("Exactly one ordered acceptance policy is required per budget")
        if self.estimator_policy.confidence_level != 0.95:
            raise ValueError("This matrix evaluates nominal 95 percent conditional intervals")
        if self.estimator_policy.source_a_fractions != self.source_a_fractions:
            raise ValueError("Estimator and matrix fractions disagree")
        if self.estimator_policy.replicates != len(self.replicate_seeds):
            raise ValueError("Estimator and matrix replicate counts disagree")
        if self.estimator_policy.total_mixture_read_groups is not None:
            raise ValueError("Matrix budgets replace the single estimator read-group budget")
        for depth in self.read_group_budgets:
            if any(0 < f < 1 and round(depth * f) in {0, depth} for f in self.source_a_fractions):
                raise ValueError("Every budget must represent every interior mixture fraction")
        return self

    @property
    def expected_levels(self) -> int:
        return (
            len(self.source_a_fractions) * len(self.read_group_budgets) * len(self.replicate_seeds)
        )


class ValidationTechnicalMetadata(StrictModel):
    sequencing_platform: str
    flow_cell_product_code: str
    library_kit: str
    basecaller_name: str
    basecaller_version: str
    basecaller_model: str
    caller_name: str
    caller_version: str
    reference_id: str
    reference_genome_build: GenomeBuild
    reference_sha256: str = Field(pattern=SHA256)
    modification_context: Literal["CpG_5mC_vs_unmethylated", "CpG_5mC_vs_unmodified_C"] = (
        "CpG_5mC_vs_unmethylated"
    )
    modification_model: str | None = None
    cytosine_modifications: Literal["5mC", "5mC_5hmC"] | None = None
    modkit_version: str | None = None
    probability_transform: Literal["none"] = "none"
    adapter_name: str
    adapter_version: str
    # Hash the adapter's full normalisation/filter policy, including modkit options.
    adapter_policy_sha256: str = Field(pattern=SHA256)


class ValidationSample(StrictModel):
    accession: str = Field(min_length=3)
    dataset_version: str = Field(min_length=1)
    biological_sample_id: str = Field(min_length=3)
    donor_id: str = Field(min_length=3)
    biological_source: str = Field(min_length=3)
    cell_type: str = Field(min_length=3)
    material_type: str = Field(min_length=3)
    input_sha256: str = Field(pattern=SHA256)
    source_modbam_sha256: str | None = Field(default=None, pattern=SHA256)
    provenance_reference: str = Field(min_length=3)
    technical: ValidationTechnicalMetadata

    @field_validator(
        "biological_sample_id", "donor_id", "biological_source", "cell_type", mode="before"
    )
    @classmethod
    def identity_has_no_controls(cls, value: Any) -> Any:
        # Invisible controls must not create apparently identical but byte-distinct
        # samples, donors or biological types at the independence/source gates.
        # Ordinary whitespace normalization remains the inherited StrictModel behavior.
        if isinstance(value, str) and any(
            unicodedata.category(character) in {"Cc", "Cf"} for character in value
        ):
            raise ValueError("Biological identity fields must not contain control characters")
        return value


class ValidationCohort(StrictModel):
    samples: dict[SampleRole, ValidationSample]
    independence_evidence: str = Field(min_length=1)
    access_basis: Literal["public", "institutionally_authorized"]
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def four_roles(self) -> ValidationCohort:
        if set(self.samples) != set(SAMPLE_ROLES):
            raise ValueError("Exactly calibration_a, calibration_b, test_a and test_b are required")
        return self


def cohort_eligibility_reasons(cohort: ValidationCohort) -> list[str]:
    """Declarations are necessary evidence, not independently verified donor identity."""
    reasons: list[str] = []
    samples = cohort.samples
    if _is_unknown_metadata(cohort.independence_evidence):
        reasons.append("Biological independence evidence is missing.")
    for role, sample in samples.items():
        for field, value in sample.model_dump(exclude={"technical"}).items():
            if isinstance(value, str) and _is_unknown_metadata(value):
                reasons.append(f"{role}.{field} is unknown.")
        for field, value in sample.technical.model_dump().items():
            if isinstance(value, str) and _is_unknown_metadata(value):
                reasons.append(f"{role}.technical.{field} is unknown.")
        if sample.technical.caller_name == "dorado" and (
            sample.technical.modification_model is None
            or sample.technical.cytosine_modifications is None
            or sample.technical.modification_context != "CpG_5mC_vs_unmodified_C"
        ):
            reasons.append(f"{role}: Dorado modification semantics are incomplete.")
        if sample.technical.adapter_name == "modkit_extract_full_tsv" and (
            sample.source_modbam_sha256 is None or sample.technical.modkit_version is None
        ):
            reasons.append(f"{role}: modkit version or source modBAM fingerprint is missing.")
    if len({sample.biological_sample_id for sample in samples.values()}) != 4:
        reasons.append(
            "Four distinct biological samples are required; a read split is insufficient."
        )
    if len({sample.input_sha256 for sample in samples.values()}) != 4:
        reasons.append("Biological inputs must not share byte-identical fingerprints.")
    lineage_hashes: set[str] = set()
    for sample in samples.values():
        current = {sample.input_sha256}
        if sample.source_modbam_sha256 is not None:
            current.add(sample.source_modbam_sha256)
        if current & lineage_hashes:
            reasons.append("Input/source-modBAM provenance overlaps between biological roles.")
            break
        lineage_hashes.update(current)
    calibration_donors = {samples[x].donor_id for x in ("calibration_a", "calibration_b")}
    test_donors = {samples[x].donor_id for x in ("test_a", "test_b")}
    if calibration_donors & test_donors:
        reasons.append("Calibration and test donors overlap; donor-independent recovery is absent.")
    if len({content_sha256(s.technical) for s in samples.values()}) != 1:
        reasons.append("All four sources require identical declared technical protocols.")
    if len({s.material_type for s in samples.values()}) != 1:
        reasons.append("Biological material types are incompatible.")
    for suffix in ("a", "b"):
        calibration = samples[cast(SampleRole, f"calibration_{suffix}")]
        test = samples[cast(SampleRole, f"test_{suffix}")]
        if (calibration.biological_source, calibration.cell_type) != (
            test.biological_source,
            test.cell_type,
        ):
            reasons.append(f"Calibration/test source {suffix} biological type is incompatible.")
    if samples["calibration_a"].biological_source == samples["calibration_b"].biological_source:
        reasons.append("Source A and source B are not declared biologically distinct.")
    return reasons


class MethylationValidationRegistration(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    registration_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,69}$")
    matrix: MethylationValidationMatrix
    cohort: ValidationCohort
    code_sha256: str = Field(pattern=SHA256)
    software_version: str = Field(min_length=1)
    registered_at: datetime
    test_outcomes_unseen: bool
    matrix_sha256: str = Field(pattern=SHA256)
    cohort_sha256: str = Field(pattern=SHA256)
    lock_sha256: str = Field(pattern=SHA256)
    sensitive_output: Literal[True] = True
    timestamp_basis: Literal["operator_declaration_not_trusted_timestamp"] = (
        "operator_declaration_not_trusted_timestamp"
    )

    @model_validator(mode="after")
    def registration_lock(self) -> MethylationValidationRegistration:
        if self.registered_at.utcoffset() is None:
            raise ValueError("Registration timestamp requires a timezone")
        if content_sha256(self.matrix) != self.matrix_sha256:
            raise ValueError("Matrix content does not match registration checksum")
        if content_sha256(self.cohort) != self.cohort_sha256:
            raise ValueError("Cohort content does not match registration checksum")
        if (
            content_sha256(self.model_dump(mode="json", exclude={"lock_sha256"}))
            != self.lock_sha256
        ):
            raise ValueError("Registration content does not match its lock")
        return self


def preregister_validation(
    matrix: MethylationValidationMatrix,
    cohort: ValidationCohort,
    *,
    registration_id: str,
    registered_at: datetime,
    test_outcomes_unseen: bool,
    code_sha256: str,
    software_version: str,
) -> MethylationValidationRegistration:
    values: dict[str, Any] = {
        "schema_version": "0.1.0",
        "registration_id": registration_id,
        "matrix": matrix.model_dump(mode="json"),
        "cohort": cohort.model_dump(mode="json"),
        "code_sha256": code_sha256,
        "software_version": software_version,
        "registered_at": registered_at.isoformat().replace("+00:00", "Z"),
        "test_outcomes_unseen": test_outcomes_unseen,
        "matrix_sha256": content_sha256(matrix),
        "cohort_sha256": content_sha256(cohort),
        "timestamp_basis": "operator_declaration_not_trusted_timestamp",
        "sensitive_output": True,
    }
    values["lock_sha256"] = content_sha256(values)
    return MethylationValidationRegistration.model_validate(values)


def validation_level_seed(
    matrix: MethylationValidationMatrix, budget: int, fraction: float, replicate: int
) -> int:
    """Depth/fraction/replicate-specific seeded draws fixed by the registered matrix."""
    payload = f"{content_sha256(matrix)}:{budget}:{fraction:.12g}:{replicate}"
    return int(hashlib.sha256(payload.encode()).hexdigest()[:15], 16)


class ValidationLevelEvidence(StrictModel):
    read_group_budget: int = Field(ge=2)
    replicate_seed: int = Field(ge=0)
    result: MethylationMixtureLevelResult


class MethylationValidationEvidence(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    registration_sha256: str = Field(pattern=SHA256)
    code_sha256: str = Field(pattern=SHA256)
    software_version: str = Field(min_length=1)
    experiment_started_at: datetime
    experiment_completed_at: datetime
    observed_input_sha256: dict[SampleRole, str]
    read_group_counts: dict[SampleRole, int]
    reference_markers: list[MethylationReferenceMarker]
    reference_marker_sha256: str = Field(pattern=SHA256)
    calibration_selection_sha256: str = Field(pattern=SHA256)
    levels: list[ValidationLevelEvidence]
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def evidence_lock(self) -> MethylationValidationEvidence:
        if any(
            x.utcoffset() is None
            for x in (self.experiment_started_at, self.experiment_completed_at)
        ):
            raise ValueError("Experiment timestamps require timezones")
        if self.experiment_started_at > self.experiment_completed_at:
            raise ValueError("Experiment completion precedes its start")
        if set(self.observed_input_sha256) != set(SAMPLE_ROLES):
            raise ValueError("Observed input fingerprints require all four roles")
        if set(self.read_group_counts) != set(SAMPLE_ROLES) or any(
            x < 1 for x in self.read_group_counts.values()
        ):
            raise ValueError("Positive read-group counts require all four roles")
        if _reference_marker_digest(self.reference_markers) != self.reference_marker_sha256:
            raise ValueError("Reference markers do not match their checksum")
        return self


class RecoveryMetrics(StrictModel):
    expected_levels: int = Field(ge=1)
    observed_levels: int = Field(ge=0)
    completed_levels: int = Field(ge=0)
    explicit_no_call_levels: int = Field(ge=0)
    missing_levels: int = Field(ge=0)
    no_call_rate: float = Field(ge=0, le=1)
    mean_bias: float | None = None
    mean_absolute_error: float | None = None
    root_mean_squared_error: float | None = None
    conditional_interval_coverage: float | None = None
    mean_interval_width: float | None = None
    accuracy_scope: Literal["completed_levels_only"] = "completed_levels_only"
    no_call_denominator: Literal["all_preregistered_levels_including_missing"] = (
        "all_preregistered_levels_including_missing"
    )


class RecoveryStratum(StrictModel):
    read_group_budget: int
    source_a_fraction: float | None = None
    metrics: RecoveryMetrics
    decision: ValidationDecision
    reasons: list[str]


_LIMITATIONS = [
    "Research use only. Acceptance limits are prospective research candidates, not validated "
    "clinical thresholds.",
    "The estimand is the fraction of retained whole read groups from source A. No DNA mass, "
    "cell fraction, tumour purity, subclone fraction or detection limit is inferred.",
    "Intervals are conditional on fixed calibration and markers; calibration, donor, marker "
    "selection and within-read/cross-marker correlation uncertainty remain unmodelled.",
    "Seeded mixtures reuse finite test pools; repetitions are not independent donors or "
    "independent physical dilution experiments. Empirical coverage is descriptive.",
    "Registration timestamps and donor/protocol metadata are operator declarations with "
    "content checksums, not independently authenticated biological or temporal evidence.",
    "SNV, CNV, ploidy and subclones require separate evidence. A 10-20 percent detection "
    "limit remains an untested hypothesis.",
]


class MethylationValidationReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    registration: MethylationValidationRegistration
    evidence: MethylationValidationEvidence | None = None
    execution_no_call_reason: str | None = Field(default=None, min_length=1, max_length=500)
    decision: ValidationDecision
    reasons: list[str]
    overall_metrics: RecoveryMetrics
    strata: list[RecoveryStratum]
    limitations: list[str] = Field(default_factory=lambda: list(_LIMITATIONS))
    research_only: Literal[True] = True
    sensitive_output: Literal[True] = True
    validation_scope: Literal["donor_independent_in_silico_read_group_recovery_only"] = (
        "donor_independent_in_silico_read_group_recovery_only"
    )
    estimand: Literal["fraction_of_retained_whole_read_groups_from_source_a"] = (
        "fraction_of_retained_whole_read_groups_from_source_a"
    )

    @model_validator(mode="after")
    def report_is_recomputed(self) -> MethylationValidationReport:
        expected = _evaluation_values(
            self.registration, self.evidence, execution_no_call_reason=self.execution_no_call_reason
        )
        for name in ("decision", "reasons", "overall_metrics", "strata"):
            actual = getattr(self, name)
            if actual != expected[name]:
                raise ValueError(f"Validation report {name} does not match recomputed evidence")
        if self.limitations != _LIMITATIONS:
            raise ValueError("Scientific limitations must remain complete and unchanged")
        return self


def evaluate_methylation_validation(
    registration: MethylationValidationRegistration,
    evidence: MethylationValidationEvidence | None = None,
    *,
    execution_no_call_reason: str | None = None,
) -> MethylationValidationReport:
    # Revalidate nested objects: model_copy(update=...) is not a validation operation.
    registration = MethylationValidationRegistration.model_validate_json(
        registration.model_dump_json()
    )
    if evidence is not None:
        evidence = MethylationValidationEvidence.model_validate_json(evidence.model_dump_json())
    # All content is checked and all derived fields are computed immediately below.
    # Avoid performing every Monte Carlo replay twice during report construction.
    return MethylationValidationReport.model_construct(
        registration=registration,
        evidence=evidence,
        execution_no_call_reason=execution_no_call_reason,
        **_evaluation_values(
            registration, evidence, execution_no_call_reason=execution_no_call_reason
        ),
    )


def _evaluation_values(
    registration: MethylationValidationRegistration,
    evidence: MethylationValidationEvidence | None,
    *,
    execution_no_call_reason: str | None = None,
) -> dict[str, Any]:
    matrix = registration.matrix
    reasons = cohort_eligibility_reasons(registration.cohort)
    if execution_no_call_reason is not None:
        if evidence is not None:
            raise ValueError(
                "An execution NO_CALL reason cannot replace submitted experiment evidence"
            )
        if not isinstance(execution_no_call_reason, str) or not (
            1 <= len(execution_no_call_reason.strip()) <= 500
        ):
            raise ValueError("An execution NO_CALL reason must contain 1 to 500 characters")
        reasons.append(execution_no_call_reason)
    if not registration.test_outcomes_unseen:
        reasons.append("Test outcomes were examined before registration.")
    records: list[ValidationLevelEvidence] = []
    if evidence is None:
        reasons.append("No independent experiment evidence is available.")
    else:
        if evidence.registration_sha256 != registration.lock_sha256:
            reasons.append("Evidence is not bound to this preregistration lock.")
        if (evidence.code_sha256, evidence.software_version) != (
            registration.code_sha256,
            registration.software_version,
        ):
            reasons.append("Executed code/version differs from the preregistered implementation.")
        if registration.registered_at >= evidence.experiment_started_at:
            reasons.append("Registration did not precede the declared experiment start.")
        expected_hashes = {
            role: sample.input_sha256 for role, sample in registration.cohort.samples.items()
        }
        if evidence.observed_input_sha256 != expected_hashes:
            reasons.append("Observed source fingerprints differ from the preregistered cohort.")
        _validate_evidence_levels(matrix, evidence)
        records = evidence.levels
    if len(records) != matrix.expected_levels:
        reasons.append("The preregistered fraction/depth/replicate grid is incomplete.")
    overall = _metrics([x.result for x in records], matrix.expected_levels)
    strata: list[RecoveryStratum] = []
    for acceptance in matrix.acceptance:
        selected = [
            x.result for x in records if x.read_group_budget == acceptance.read_group_budget
        ]
        expected_per_fraction = len(matrix.replicate_seeds)
        expected_per_depth = expected_per_fraction * len(matrix.source_a_fractions)
        metrics = _metrics(selected, expected_per_depth)
        failed = _metric_failures(metrics, acceptance, per_fraction=False)
        strata.append(
            RecoveryStratum(
                read_group_budget=acceptance.read_group_budget,
                metrics=metrics,
                decision=_stratum_decision(metrics, failed),
                reasons=failed,
            )
        )
        for fraction in matrix.source_a_fractions:
            metrics = _metrics(
                [x for x in selected if x.target_source_a_fraction == fraction],
                expected_per_fraction,
            )
            failed = _metric_failures(metrics, acceptance, per_fraction=True)
            strata.append(
                RecoveryStratum(
                    read_group_budget=acceptance.read_group_budget,
                    source_a_fraction=fraction,
                    metrics=metrics,
                    decision=_stratum_decision(metrics, failed),
                    reasons=failed,
                )
            )
    if reasons or overall.completed_levels == 0:
        decision = ValidationDecision.NO_CALL
        if overall.completed_levels == 0:
            reasons.append("No preregistered level yielded an evaluable quantitative estimate.")
    elif any(x.decision != ValidationDecision.PASS for x in strata):
        decision = ValidationDecision.FAIL
        reasons.append("At least one preregistered depth or fraction stratum fails acceptance.")
    else:
        decision = ValidationDecision.PASS
    return {"decision": decision, "reasons": reasons, "overall_metrics": overall, "strata": strata}


def _metrics(levels: Sequence[MethylationMixtureLevelResult], expected: int) -> RecoveryMetrics:
    completed = [x for x in levels if x.status == ModuleRunStatus.COMPLETED]
    errors = [
        cast(float, x.estimated_source_a_fraction) - x.realized_source_a_read_group_fraction
        for x in completed
    ]
    n = len(completed)
    return RecoveryMetrics(
        expected_levels=expected,
        observed_levels=len(levels),
        completed_levels=n,
        explicit_no_call_levels=len(levels) - n,
        missing_levels=expected - len(levels),
        no_call_rate=(expected - n) / expected,
        mean_bias=sum(errors) / n if n else None,
        mean_absolute_error=sum(abs(x) for x in errors) / n if n else None,
        root_mean_squared_error=math.sqrt(sum(x * x for x in errors) / n) if n else None,
        conditional_interval_coverage=sum(
            cast(float, x.confidence_interval_lower)
            <= x.realized_source_a_read_group_fraction
            <= cast(float, x.confidence_interval_upper)
            for x in completed
        )
        / n
        if n
        else None,
        mean_interval_width=sum(
            cast(float, x.confidence_interval_upper) - cast(float, x.confidence_interval_lower)
            for x in completed
        )
        / n
        if n
        else None,
    )


def _metric_failures(
    metrics: RecoveryMetrics, limits: RecoveryAcceptance, *, per_fraction: bool
) -> list[str]:
    reasons: list[str] = []
    if metrics.missing_levels:
        reasons.append("Preregistered levels are missing.")
    max_no_call = (
        limits.maximum_per_fraction_no_call_rate if per_fraction else limits.maximum_no_call_rate
    )
    if metrics.no_call_rate > max_no_call:
        reasons.append(f"NO_CALL rate {metrics.no_call_rate:.4f} exceeds {max_no_call:.4f}.")
    if metrics.completed_levels == 0:
        reasons.append("Accuracy and interval coverage cannot be estimated.")
        return reasons
    bias_limit = (
        limits.maximum_per_fraction_absolute_bias if per_fraction else limits.maximum_absolute_bias
    )
    if abs(cast(float, metrics.mean_bias)) > bias_limit:
        reasons.append(f"Absolute bias exceeds {bias_limit:.4f}.")
    if not per_fraction:
        if cast(float, metrics.mean_absolute_error) > limits.maximum_mae:
            reasons.append(f"MAE exceeds {limits.maximum_mae:.4f}.")
        if cast(float, metrics.root_mean_squared_error) > limits.maximum_rmse:
            reasons.append(f"RMSE exceeds {limits.maximum_rmse:.4f}.")
    min_coverage = (
        limits.minimum_per_fraction_conditional_interval_coverage
        if per_fraction
        else limits.minimum_conditional_interval_coverage
    )
    if cast(float, metrics.conditional_interval_coverage) < min_coverage:
        reasons.append(f"Conditional 95 percent interval coverage is below {min_coverage:.4f}.")
    return reasons


def _stratum_decision(metrics: RecoveryMetrics, failures: Sequence[str]) -> ValidationDecision:
    if metrics.missing_levels or metrics.completed_levels == 0:
        return ValidationDecision.NO_CALL
    return ValidationDecision.FAIL if failures else ValidationDecision.PASS


def _validate_evidence_levels(
    matrix: MethylationValidationMatrix, evidence: MethylationValidationEvidence
) -> None:
    policy = matrix.estimator_policy
    markers = evidence.reference_markers
    marker_ids = [x.marker_id for x in markers]
    if marker_ids != sorted(set(marker_ids)):
        raise ValueError("Marker lock must be unique and canonically ordered")
    if policy.maximum_markers is not None and len(markers) > policy.maximum_markers:
        raise ValueError("Marker lock exceeds the registered policy")
    for marker in markers:
        if (marker.source_a_total_read_groups, marker.source_b_total_read_groups) != (
            evidence.read_group_counts["calibration_a"],
            evidence.read_group_counts["calibration_b"],
        ):
            raise ValueError("Marker counts disagree with independent calibration source counts")
        if (
            marker.source_a_valid_calls < policy.minimum_reference_valid_calls
            or marker.source_b_valid_calls < policy.minimum_reference_valid_calls
            or abs(marker.delta_beta) + 1e-12 < policy.minimum_absolute_delta_beta
        ):
            raise ValueError("Marker lock violates registered calibration filters")
    seen: set[tuple[int, float, int]] = set()
    level_ids: set[str] = set()
    for record in evidence.levels:
        level = record.result
        key = (record.read_group_budget, level.target_source_a_fraction, level.replicate)
        if key in seen or level.level_id in level_ids:
            raise ValueError("Duplicate grid cell or level identifier")
        seen.add(key)
        level_ids.add(level.level_id)
        if (
            record.read_group_budget not in matrix.read_group_budgets
            or level.target_source_a_fraction not in matrix.source_a_fractions
            or not 1 <= level.replicate <= len(matrix.replicate_seeds)
        ):
            raise ValueError("Evidence contains an unregistered grid cell")
        if record.replicate_seed != matrix.replicate_seeds[level.replicate - 1]:
            raise ValueError("Replicate seed differs from preregistered matrix")
        if level.seed != validation_level_seed(matrix, *key):
            raise ValueError("Mixture seed differs from preregistered matrix")
        if level.total_read_groups != record.read_group_budget:
            raise ValueError("Mixture read depth differs from preregistered whole-read budget")
        if level.source_a_read_groups != round(record.read_group_budget * key[1]):
            raise ValueError("Mixture source count disagrees with the requested fraction")
        if (
            level.source_a_read_groups > evidence.read_group_counts["test_a"]
            or level.source_b_read_groups > evidence.read_group_counts["test_b"]
        ):
            raise ValueError("Mixture exceeds the independent test pool")
        _validate_level_numerics(level, markers, policy)


def _validate_level_numerics(
    level: MethylationMixtureLevelResult,
    markers: Sequence[MethylationReferenceMarker],
    policy: MethylationMixturePolicy,
) -> None:
    observations = level.marker_observations
    if [x.marker_id for x in observations] != [x.marker_id for x in markers]:
        raise ValueError("Level observations disagree with the calibration marker lock")
    used = sum(x.valid_calls >= policy.minimum_mixture_valid_calls for x in observations)
    if (
        level.markers_available != len(markers)
        or level.markers_used != used
        or level.valid_mixture_calls != sum(x.valid_calls for x in observations)
        or level.confidence_level != policy.confidence_level
    ):
        raise ValueError("Level counts or confidence level disagree with registered evidence")
    counts = {
        x.marker_id: _Counts(
            modified=x.modified_calls,
            unmethylated=x.unmethylated_calls,
            ambiguous=x.ambiguous_calls,
        )
        for x in observations
    }
    raw, point, fit, _ = _point_estimate(
        markers,
        counts,
        minimum_mixture_valid_calls=policy.minimum_mixture_valid_calls,
        mixture_total_read_groups=level.total_read_groups,
    )
    expected_completed = (
        len(markers) >= policy.minimum_markers
        and used >= policy.minimum_markers
        and raw is not None
        and point is not None
        and fit is not None
    )
    expected_completed = bool(
        expected_completed
        and (
            -policy.maximum_extrapolation <= cast(float, raw) <= 1 + policy.maximum_extrapolation
            and cast(float, fit) <= policy.maximum_standardized_model_fit_rmse
        )
    )
    successful: list[float] = []
    lower: float | None = None
    upper: float | None = None
    if expected_completed:
        rng = random.Random(level.seed + 9_000_001)
        draws = [
            _draw_estimate(
                markers,
                counts,
                minimum_mixture_valid_calls=policy.minimum_mixture_valid_calls,
                mixture_total_read_groups=level.total_read_groups,
                fitted_source_a_fraction=cast(float, point),
                rng=rng,
            )
            for _ in range(policy.uncertainty_draws)
        ]
        successful = sorted(x for x in draws if x is not None)
        required = math.ceil(
            policy.uncertainty_draws * policy.minimum_successful_uncertainty_fraction
        )
        expected_completed = len(successful) >= required
        if expected_completed:
            tail = (1 - policy.confidence_level) / 2
            lower, upper = _quantile(successful, tail), _quantile(successful, 1 - tail)
            expected_completed = upper - lower <= policy.maximum_confidence_interval_width
    if (level.status == ModuleRunStatus.COMPLETED) != expected_completed:
        raise ValueError("Level status cannot be reproduced from counts, policy and seed")
    if level.successful_uncertainty_draws != len(successful):
        raise ValueError("Uncertainty draw count cannot be reproduced")
    if expected_completed:
        for actual, expected in zip(
            (
                level.unconstrained_source_a_fraction,
                level.estimated_source_a_fraction,
                level.model_fit_standardized_rmse,
                level.confidence_interval_lower,
                level.confidence_interval_upper,
            ),
            (raw, point, fit, lower, upper),
            strict=True,
        ):
            if not math.isclose(
                cast(float, actual), cast(float, expected), rel_tol=0, abs_tol=1e-12
            ):
                raise ValueError("Level estimate or interval cannot be reproduced from evidence")

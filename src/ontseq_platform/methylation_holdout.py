"""Prospective recovery in held-out read groups from two fixed biological sources.

Read-pool disjointness is not donor independence. This contract supports a scoped
technical recovery decision while donor/population transfer remains unevaluated.
"""

from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime
from typing import Any, Literal, cast

from pydantic import Field, model_validator

from .methylation_mixture import (
    MethylationReferenceMarker,
    _calibration_read_group_count,
    _reference_marker_digest,
)
from .methylation_validation import (
    SAMPLE_ROLES,
    SHA256,
    MethylationValidationEvidence,
    MethylationValidationMatrix,
    RecoveryMetrics,
    RecoveryStratum,
    SampleRole,
    ValidationDecision,
    ValidationLevelEvidence,
    ValidationSample,
    _is_unknown_metadata,
    _metric_failures,
    _metrics,
    _stratum_decision,
    _validate_evidence_levels,
    content_sha256,
)
from .models import StrictModel

HoldoutSourceRole = Literal["source_a", "source_b"]
HOLDOUT_SOURCE_ROLES: tuple[HoldoutSourceRole, ...] = ("source_a", "source_b")


class PairedHoldoutCohort(StrictModel):
    samples: dict[HoldoutSourceRole, ValidationSample]
    source_distinction_evidence: str = Field(min_length=1)
    access_basis: Literal["public", "institutionally_authorized"]
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def exactly_two_sources(self) -> PairedHoldoutCohort:
        if set(self.samples) != set(HOLDOUT_SOURCE_ROLES):
            raise ValueError("Exactly source_a and source_b are required")
        return self


class PairedHoldoutSplitPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    algorithm: Literal["sorted_read_identifiers_seeded_shuffle_v1"] = (
        "sorted_read_identifiers_seeded_shuffle_v1"
    )
    calibration_read_group_fraction: float = Field(gt=0, lt=1)
    source_a_split_seed: int = Field(ge=0, le=2_147_483_647)
    source_b_split_seed: int = Field(ge=0, le=2_147_483_647)
    split_unit: Literal["whole_read_group"] = "whole_read_group"

    @model_validator(mode="after")
    def separate_seeds(self) -> PairedHoldoutSplitPolicy:
        if self.source_a_split_seed == self.source_b_split_seed:
            raise ValueError("The two source split seeds must be explicitly distinct")
        return self


def expected_calibration_read_groups(total: int, fraction: float) -> int:
    """Use the existing deterministic round-and-clamp rule, with both pools nonempty."""
    if total < 2 or not 0 < fraction < 1:
        raise ValueError("Each source needs at least two read groups and a proper split fraction")
    return _calibration_read_group_count(total, fraction)


def holdout_split_manifest_sha256(
    *,
    registration_sha256: str,
    observed_input_sha256: Mapping[HoldoutSourceRole, str],
    source_read_group_counts: Mapping[HoldoutSourceRole, int],
    calibration_read_group_counts: Mapping[HoldoutSourceRole, int],
    test_read_group_counts: Mapping[HoldoutSourceRole, int],
    pool_selection_sha256: Mapping[SampleRole, str],
) -> str:
    """Bind pool-selection attestations and counts to their registered source inputs."""
    return content_sha256(
        {
            "registration_sha256": registration_sha256,
            "observed_input_sha256": dict(observed_input_sha256),
            "source_read_group_counts": dict(source_read_group_counts),
            "calibration_read_group_counts": dict(calibration_read_group_counts),
            "test_read_group_counts": dict(test_read_group_counts),
            "pool_selection_sha256": dict(pool_selection_sha256),
        }
    )


def holdout_cohort_eligibility_reasons(cohort: PairedHoldoutCohort) -> list[str]:
    """Require source/protocol provenance; donor independence is deliberately not inferred."""
    reasons: list[str] = []
    if _is_unknown_metadata(cohort.source_distinction_evidence):
        reasons.append("Evidence for the biological distinction of source A and B is missing.")
    for role, sample in cohort.samples.items():
        for field, value in sample.model_dump(exclude={"technical", "donor_id"}).items():
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
    source_a, source_b = (cohort.samples[role] for role in HOLDOUT_SOURCE_ROLES)
    if source_a.biological_source == source_b.biological_source:
        reasons.append("Source A and source B are not declared biologically distinct.")
    if source_a.biological_sample_id == source_b.biological_sample_id:
        reasons.append("The two input sources share a biological sample identifier.")
    lineage_a = {source_a.input_sha256, source_a.source_modbam_sha256} - {None}
    lineage_b = {source_b.input_sha256, source_b.source_modbam_sha256} - {None}
    if lineage_a & lineage_b:
        reasons.append("The two biological inputs share input or parent modBAM fingerprints.")
    if content_sha256(source_a.technical) != content_sha256(source_b.technical):
        reasons.append("Source A and source B have incompatible declared technical protocols.")
    if source_a.material_type != source_b.material_type:
        reasons.append("The two sources have incompatible biological material types.")
    return reasons


class PairedHoldoutRegistration(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    study_design: Literal["paired_source_read_holdout"] = "paired_source_read_holdout"
    registration_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,69}$")
    matrix: MethylationValidationMatrix
    cohort: PairedHoldoutCohort
    split_policy: PairedHoldoutSplitPolicy
    code_sha256: str = Field(pattern=SHA256)
    software_version: str = Field(min_length=1)
    registered_at: datetime
    test_outcomes_unseen: bool
    matrix_sha256: str = Field(pattern=SHA256)
    cohort_sha256: str = Field(pattern=SHA256)
    split_policy_sha256: str = Field(pattern=SHA256)
    lock_sha256: str = Field(pattern=SHA256)
    timestamp_basis: Literal["operator_declaration_not_trusted_timestamp"] = (
        "operator_declaration_not_trusted_timestamp"
    )
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def registration_lock(self) -> PairedHoldoutRegistration:
        if self.registered_at.utcoffset() is None:
            raise ValueError("Registration timestamp requires a timezone")
        if self.split_policy.calibration_read_group_fraction != (
            self.matrix.estimator_policy.calibration_read_group_fraction
        ):
            raise ValueError(
                "Explicit split fraction disagrees with the registered estimator policy"
            )
        for name in ("matrix", "cohort", "split_policy"):
            if content_sha256(getattr(self, name)) != getattr(self, f"{name}_sha256"):
                raise ValueError(f"Registered {name} does not match its checksum")
        if (
            content_sha256(self.model_dump(mode="json", exclude={"lock_sha256"}))
            != self.lock_sha256
        ):
            raise ValueError("Registration content does not match its lock")
        return self


def preregister_paired_holdout(
    matrix: MethylationValidationMatrix,
    cohort: PairedHoldoutCohort,
    split_policy: PairedHoldoutSplitPolicy,
    *,
    registration_id: str,
    registered_at: datetime,
    test_outcomes_unseen: bool,
    code_sha256: str,
    software_version: str,
) -> PairedHoldoutRegistration:
    values: dict[str, Any] = {
        "schema_version": "0.1.0",
        "study_design": "paired_source_read_holdout",
        "registration_id": registration_id,
        "matrix": matrix.model_dump(mode="json"),
        "cohort": cohort.model_dump(mode="json"),
        "split_policy": split_policy.model_dump(mode="json"),
        "code_sha256": code_sha256,
        "software_version": software_version,
        "registered_at": registered_at.isoformat().replace("+00:00", "Z"),
        "test_outcomes_unseen": test_outcomes_unseen,
        "matrix_sha256": content_sha256(matrix),
        "cohort_sha256": content_sha256(cohort),
        "split_policy_sha256": content_sha256(split_policy),
        "timestamp_basis": "operator_declaration_not_trusted_timestamp",
        "sensitive_output": True,
    }
    values["lock_sha256"] = content_sha256(values)
    return PairedHoldoutRegistration.model_validate(values)


class PairedHoldoutEvidence(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    study_design: Literal["paired_source_read_holdout"] = "paired_source_read_holdout"
    registration_sha256: str = Field(pattern=SHA256)
    code_sha256: str = Field(pattern=SHA256)
    software_version: str = Field(min_length=1)
    experiment_started_at: datetime
    experiment_completed_at: datetime
    observed_input_sha256: dict[HoldoutSourceRole, str]
    source_read_group_counts: dict[HoldoutSourceRole, int]
    calibration_read_group_counts: dict[HoldoutSourceRole, int]
    test_read_group_counts: dict[HoldoutSourceRole, int]
    # These keys denote four COMPUTATIONAL POOLS, never four biological samples.
    pool_selection_sha256: dict[SampleRole, str]
    split_manifest_sha256: str = Field(pattern=SHA256)
    reference_markers: list[MethylationReferenceMarker]
    reference_marker_sha256: str = Field(pattern=SHA256)
    levels: list[ValidationLevelEvidence]
    pool_role_definition: Literal["calibration_and_test_read_partitions_of_two_sources"] = (
        "calibration_and_test_read_partitions_of_two_sources"
    )
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def evidence_is_coherent(self) -> PairedHoldoutEvidence:
        for timestamp in (self.experiment_started_at, self.experiment_completed_at):
            if timestamp.utcoffset() is None:
                raise ValueError("Experiment timestamps require timezones")
        if self.experiment_started_at > self.experiment_completed_at:
            raise ValueError("Experiment completion precedes its start")
        for mapping in (
            self.observed_input_sha256,
            self.source_read_group_counts,
            self.calibration_read_group_counts,
            self.test_read_group_counts,
        ):
            if set(mapping) != set(HOLDOUT_SOURCE_ROLES):
                raise ValueError("Both biological source keys are required")
        for role in HOLDOUT_SOURCE_ROLES:
            total = self.source_read_group_counts[role]
            calibration = self.calibration_read_group_counts[role]
            test = self.test_read_group_counts[role]
            if total < 2 or min(calibration, test) < 1 or calibration + test != total:
                raise ValueError("Calibration and test pool counts must partition each source")
        if set(self.pool_selection_sha256) != set(SAMPLE_ROLES):
            raise ValueError("Four computational pool selection digests are required")
        for digest in (*self.observed_input_sha256.values(), *self.pool_selection_sha256.values()):
            if len(digest) != 64 or any(x not in "0123456789abcdef" for x in digest):
                raise ValueError("Source and pool fingerprints must be lowercase SHA-256")
        if _reference_marker_digest(self.reference_markers) != self.reference_marker_sha256:
            raise ValueError("Reference markers do not match their checksum")
        if self.split_manifest_sha256 != holdout_split_manifest_sha256(
            registration_sha256=self.registration_sha256,
            observed_input_sha256=self.observed_input_sha256,
            source_read_group_counts=self.source_read_group_counts,
            calibration_read_group_counts=self.calibration_read_group_counts,
            test_read_group_counts=self.test_read_group_counts,
            pool_selection_sha256=self.pool_selection_sha256,
        ):
            raise ValueError(
                "Source/pool contents do not match the registered split manifest digest"
            )
        return self


_LIMITATIONS = [
    "Research use only. PASS means recovery of retained whole read-group fractions in the "
    "registered two-source holdout experiment, not biological population or clinical validation.",
    "Calibration and test pools are disjoint read groups from the same two biological sources. "
    "They are not four independent samples or independent calibration/test donors.",
    "Donor-independent validation is NOT_EVALUATED. DNA mass, cell fraction, tumour purity, "
    "subclone fraction and a detection limit are not inferred.",
    "The nominal 95 percent intervals condition on fixed calibration and markers. Calibration, "
    "donor, marker-selection and within-read/cross-marker correlation uncertainty are unmodelled.",
    "Seeded mixtures reuse finite held-out pools. Repetitions are not independent biological "
    "samples or physical dilution experiments; empirical interval coverage is descriptive.",
    "Metadata, unseen-outcome declarations and timestamps are operator attestations with "
    "content locks, not independently authenticated biological or temporal evidence.",
    "SNV, CNV, ploidy and subclones require separate evidence. A 10-20 percent detection "
    "limit remains an untested hypothesis.",
]


class PairedHoldoutRecoveryReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    registration: PairedHoldoutRegistration
    evidence: PairedHoldoutEvidence | None = None
    execution_no_call_reason: str | None = Field(default=None, min_length=1, max_length=500)
    decision: ValidationDecision
    reasons: list[str]
    overall_metrics: RecoveryMetrics
    strata: list[RecoveryStratum]
    validation_scope: Literal["held_out_read_group_recovery_within_two_sources"] = (
        "held_out_read_group_recovery_within_two_sources"
    )
    donor_independent_validation: Literal["NOT_EVALUATED"] = "NOT_EVALUATED"
    estimand: Literal["fraction_of_retained_whole_read_groups_from_source_a"] = (
        "fraction_of_retained_whole_read_groups_from_source_a"
    )
    limitations: list[str] = Field(default_factory=lambda: list(_LIMITATIONS))
    research_only: Literal[True] = True
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def report_is_recomputed(self) -> PairedHoldoutRecoveryReport:
        expected = _evaluation_values(
            self.registration, self.evidence, execution_no_call_reason=self.execution_no_call_reason
        )
        for name in ("decision", "reasons", "overall_metrics", "strata"):
            if getattr(self, name) != expected[name]:
                raise ValueError(f"Holdout report {name} disagrees with recomputed evidence")
        if self.limitations != _LIMITATIONS:
            raise ValueError("The complete scoped research limitations must remain unchanged")
        return self


def evaluate_paired_holdout(
    registration: PairedHoldoutRegistration,
    evidence: PairedHoldoutEvidence | None = None,
    *,
    execution_no_call_reason: str | None = None,
) -> PairedHoldoutRecoveryReport:
    registration = PairedHoldoutRegistration.model_validate_json(registration.model_dump_json())
    if evidence is not None:
        evidence = PairedHoldoutEvidence.model_validate_json(evidence.model_dump_json())
    return PairedHoldoutRecoveryReport.model_construct(
        registration=registration,
        evidence=evidence,
        execution_no_call_reason=execution_no_call_reason,
        **_evaluation_values(
            registration, evidence, execution_no_call_reason=execution_no_call_reason
        ),
    )


def _numerical_pool_view(evidence: PairedHoldoutEvidence) -> MethylationValidationEvidence:
    """Internal numeric compatibility view; it introduces no biological sample identities.

    The existing replay helper only uses marker locks, computational pool counts and
    levels. The two input hashes therefore repeat across their two actual partitions.
    This internal object is never interpreted or exported as an independent study.
    """
    role_source = {role: cast(HoldoutSourceRole, f"source_{role[-1]}") for role in SAMPLE_ROLES}
    counts = {
        role: (
            evidence.calibration_read_group_counts
            if role.startswith("calibration")
            else evidence.test_read_group_counts
        )[source]
        for role, source in role_source.items()
    }
    return MethylationValidationEvidence(
        registration_sha256=evidence.registration_sha256,
        code_sha256=evidence.code_sha256,
        software_version=evidence.software_version,
        experiment_started_at=evidence.experiment_started_at,
        experiment_completed_at=evidence.experiment_completed_at,
        observed_input_sha256={
            role: evidence.observed_input_sha256[source] for role, source in role_source.items()
        },
        read_group_counts=counts,
        reference_markers=evidence.reference_markers,
        reference_marker_sha256=evidence.reference_marker_sha256,
        calibration_selection_sha256=content_sha256(
            {
                role: evidence.pool_selection_sha256[role]
                for role in ("calibration_a", "calibration_b")
            }
        ),
        levels=evidence.levels,
    )


def _evaluation_values(
    registration: PairedHoldoutRegistration,
    evidence: PairedHoldoutEvidence | None,
    *,
    execution_no_call_reason: str | None = None,
) -> dict[str, Any]:
    matrix = registration.matrix
    reasons = holdout_cohort_eligibility_reasons(registration.cohort)
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
        reasons.append("Held-out test outcomes were examined before registration.")
    records: list[ValidationLevelEvidence] = []
    if evidence is None:
        reasons.append("No registered two-source holdout experiment evidence is available.")
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
            reasons.append("Observed source fingerprints differ from the registered two sources.")
        for role in HOLDOUT_SOURCE_ROLES:
            expected_count = expected_calibration_read_groups(
                evidence.source_read_group_counts[role],
                registration.split_policy.calibration_read_group_fraction,
            )
            if evidence.calibration_read_group_counts[role] != expected_count:
                raise ValueError("Calibration pool size disagrees with the registered split policy")
        _validate_evidence_levels(matrix, _numerical_pool_view(evidence))
        records = evidence.levels
    if len(records) != matrix.expected_levels:
        reasons.append("The preregistered fraction/depth/replicate grid is incomplete.")
    overall = _metrics([x.result for x in records], matrix.expected_levels)
    strata: list[RecoveryStratum] = []
    for acceptance in matrix.acceptance:
        selected = [
            x.result for x in records if x.read_group_budget == acceptance.read_group_budget
        ]
        for fraction in (None, *matrix.source_a_fractions):
            per_fraction = fraction is not None
            expected = len(matrix.replicate_seeds) * (
                1 if per_fraction else len(matrix.source_a_fractions)
            )
            levels = (
                [x for x in selected if x.target_source_a_fraction == fraction]
                if per_fraction
                else selected
            )
            metrics = _metrics(levels, expected)
            failed = _metric_failures(metrics, acceptance, per_fraction=per_fraction)
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
            reasons.append("No registered level yielded an evaluable quantitative estimate.")
    elif any(x.decision != ValidationDecision.PASS for x in strata):
        decision = ValidationDecision.FAIL
        reasons.append(
            "At least one registered depth or fraction stratum fails recovery acceptance."
        )
    else:
        decision = ValidationDecision.PASS
    return {"decision": decision, "reasons": reasons, "overall_metrics": overall, "strata": strata}

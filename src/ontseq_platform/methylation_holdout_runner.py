"""Execute a registered whole-read holdout experiment from two biological sources."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from typing import Literal

from pydantic import model_validator

from . import __version__
from .methylation_holdout import (
    HOLDOUT_SOURCE_ROLES,
    HoldoutSourceRole,
    PairedHoldoutEvidence,
    PairedHoldoutRecoveryReport,
    PairedHoldoutRegistration,
    evaluate_paired_holdout,
    holdout_cohort_eligibility_reasons,
    holdout_split_manifest_sha256,
)
from .methylation_mixture import (
    _aggregate,
    _NanopolishCall,
    _reference_marker_digest,
    _reference_markers,
    _selection_digest,
    _split_reads,
)
from .methylation_validation import SampleRole
from .methylation_validation_runner import (
    ValidationLocalInput,
    generate_validation_levels,
    parse_validation_source,
    validation_software_sha256,
)
from .models import StrictModel
from .reference import sha256_file


class PairedHoldoutLocalInputs(StrictModel):
    samples: dict[HoldoutSourceRole, ValidationLocalInput]
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def exactly_two_sources(self) -> PairedHoldoutLocalInputs:
        if set(self.samples) != set(HOLDOUT_SOURCE_ROLES):
            raise ValueError("Local inputs must bind exactly source_a and source_b")
        if len({item.input_format for item in self.samples.values()}) != 1:
            raise ValueError("One registered holdout experiment cannot mix adapter formats")
        return self


@dataclass(frozen=True)
class _ReadPool:
    """A computational read partition, never an additional biological sample."""

    calls_by_read: Mapping[str, tuple[_NanopolishCall, ...]]


def execute_registered_paired_holdout(
    registration: PairedHoldoutRegistration,
    inputs: PairedHoldoutLocalInputs,
) -> PairedHoldoutRecoveryReport:
    registration = PairedHoldoutRegistration.model_validate_json(registration.model_dump_json())
    inputs = PairedHoldoutLocalInputs.model_validate_json(inputs.model_dump_json())
    if (
        holdout_cohort_eligibility_reasons(registration.cohort)
        or not registration.test_outcomes_unseen
    ):
        return evaluate_paired_holdout(registration)
    if validation_software_sha256() != registration.code_sha256:
        raise ValueError("Software/runtime changed after holdout registration")
    started = datetime.now(UTC)
    if started <= registration.registered_at:
        raise ValueError("Holdout experiment must begin after its declared registration")
    sources = {
        role: parse_validation_source(
            role,
            registration.cohort.samples[role],
            inputs.samples[role],
            registration.matrix.estimator_policy,
        )
        for role in HOLDOUT_SOURCE_ROLES
    }
    if set(sources["source_a"].calls_by_read) & set(sources["source_b"].calls_by_read):
        raise ValueError("Read identifiers overlap between the two declared biological sources")
    calibration: dict[HoldoutSourceRole, _ReadPool] = {}
    test: dict[HoldoutSourceRole, _ReadPool] = {}
    pool_digests: dict[SampleRole, str] = {}
    split = registration.split_policy
    roles: tuple[tuple[HoldoutSourceRole, SampleRole, SampleRole, int], ...] = (
        ("source_a", "calibration_a", "test_a", split.source_a_split_seed),
        ("source_b", "calibration_b", "test_b", split.source_b_split_seed),
    )
    for role, calibration_role, test_role, seed in roles:
        all_calls = sources[role].calls_by_read
        if len(all_calls) < 2:
            return evaluate_paired_holdout(
                registration,
                execution_no_call_reason="Each source needs at least two retained read groups.",
            )
        calibration_names, test_names = _split_reads(
            tuple(all_calls),
            fraction=split.calibration_read_group_fraction,
            seed=seed,
        )
        calibration[role] = _ReadPool({name: all_calls[name] for name in calibration_names})
        test[role] = _ReadPool({name: all_calls[name] for name in test_names})
        # Context names identify computational pools, not biological sample identities.
        pool_digests[calibration_role] = _selection_digest(
            calibration_names, (), context=f"paired-holdout:{calibration_role}"
        )
        pool_digests[test_role] = _selection_digest(
            test_names, (), context=f"paired-holdout:{test_role}"
        )
    if any(
        len(pool.calls_by_read) < max(registration.matrix.read_group_budgets)
        for pool in test.values()
    ):
        report = evaluate_paired_holdout(
            registration,
            execution_no_call_reason=(
                "The fixed holdout pools cannot supply the largest registered read-group budget."
            ),
        )
    else:
        markers = _reference_markers(
            _aggregate(calibration["source_a"], tuple(calibration["source_a"].calls_by_read)),
            _aggregate(calibration["source_b"], tuple(calibration["source_b"].calls_by_read)),
            source_a_total_read_groups=len(calibration["source_a"].calls_by_read),
            source_b_total_read_groups=len(calibration["source_b"].calls_by_read),
            policy=registration.matrix.estimator_policy,
        )
        levels = generate_validation_levels(
            registration.matrix,
            test["source_a"],
            test["source_b"],
            markers,
        )
        observed_hashes = {
            role: registration.cohort.samples[role].input_sha256 for role in HOLDOUT_SOURCE_ROLES
        }
        source_counts = {role: len(sources[role].calls_by_read) for role in HOLDOUT_SOURCE_ROLES}
        calibration_counts = {
            role: len(calibration[role].calls_by_read) for role in HOLDOUT_SOURCE_ROLES
        }
        test_counts = {role: len(test[role].calls_by_read) for role in HOLDOUT_SOURCE_ROLES}
        evidence = PairedHoldoutEvidence(
            registration_sha256=registration.lock_sha256,
            code_sha256=validation_software_sha256(),
            software_version=__version__,
            experiment_started_at=started,
            experiment_completed_at=datetime.now(UTC),
            observed_input_sha256=observed_hashes,
            source_read_group_counts=source_counts,
            calibration_read_group_counts=calibration_counts,
            test_read_group_counts=test_counts,
            pool_selection_sha256=pool_digests,
            split_manifest_sha256=holdout_split_manifest_sha256(
                registration_sha256=registration.lock_sha256,
                observed_input_sha256=observed_hashes,
                source_read_group_counts=source_counts,
                calibration_read_group_counts=calibration_counts,
                test_read_group_counts=test_counts,
                pool_selection_sha256=pool_digests,
            ),
            reference_markers=markers,
            reference_marker_sha256=_reference_marker_digest(markers),
            levels=levels,
        )
        report = evaluate_paired_holdout(registration, evidence)
    for role in HOLDOUT_SOURCE_ROLES:
        if (
            sha256_file(inputs.samples[role].calls_path)
            != registration.cohort.samples[role].input_sha256
        ):
            raise ValueError(f"{role}: input changed during the holdout experiment")
    if validation_software_sha256() != registration.code_sha256:
        raise ValueError("Software/runtime changed during the holdout experiment")
    return report

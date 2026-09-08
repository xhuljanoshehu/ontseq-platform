from __future__ import annotations

import hashlib
import unittest
from datetime import UTC, datetime, timedelta

from pydantic import ValidationError

from ontseq_platform.methylation_holdout import (
    HOLDOUT_SOURCE_ROLES,
    PairedHoldoutCohort,
    PairedHoldoutEvidence,
    PairedHoldoutRecoveryReport,
    PairedHoldoutRegistration,
    PairedHoldoutSplitPolicy,
    evaluate_paired_holdout,
    expected_calibration_read_groups,
    holdout_cohort_eligibility_reasons,
    holdout_split_manifest_sha256,
    preregister_paired_holdout,
)
from ontseq_platform.methylation_mixture import (
    MethylationMixturePolicy,
    MethylationSourceSummary,
    NanopolishCallState,
    _aggregate,
    _estimate_level,
    _MarkerKey,
    _NanopolishCall,
    _ParsedSource,
    _reference_marker_digest,
    _reference_markers,
    _selection_digest,
    _split_reads,
)
from ontseq_platform.methylation_validation import (
    MethylationValidationMatrix,
    RecoveryAcceptance,
    ValidationDecision,
    ValidationLevelEvidence,
    ValidationSample,
    ValidationTechnicalMetadata,
    content_sha256,
    validation_level_seed,
)
from ontseq_platform.models import FileFingerprint, GenomeBuild


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def _registration(*, minimum_markers: int = 2) -> PairedHoldoutRegistration:
    policy = MethylationMixturePolicy(
        profile_id="synthetic-holdout-policy",
        source_a_fractions=[0, 0.5, 1],
        replicates=10,
        minimum_reference_valid_calls=2,
        minimum_mixture_valid_calls=1,
        minimum_markers=minimum_markers,
        maximum_markers=10,
        uncertainty_draws=100,
        maximum_confidence_interval_width=1.0,
        maximum_standardized_model_fit_rmse=100.0,
        maximum_extrapolation=1.0,
        note="Synthetic technical regression only; no human data.",
    )
    matrix = MethylationValidationMatrix(
        matrix_id="synthetic-holdout-matrix",
        source_a_fractions=[0, 0.5, 1],
        read_group_budgets=[20],
        replicate_seeds=list(range(10)),
        estimator_policy=policy,
        acceptance=[
            RecoveryAcceptance(
                read_group_budget=20,
                maximum_absolute_bias=0.05,
                maximum_mae=0.1,
                maximum_rmse=0.15,
                minimum_conditional_interval_coverage=0.9,
                maximum_no_call_rate=0.1,
                maximum_per_fraction_absolute_bias=0.1,
                minimum_per_fraction_conditional_interval_coverage=0.8,
                maximum_per_fraction_no_call_rate=0.2,
            )
        ],
        rationale="Small deterministic software test; no biological evidence.",
    )
    technical = ValidationTechnicalMetadata(
        sequencing_platform="synthetic-ONT",
        flow_cell_product_code="synthetic-flowcell",
        library_kit="synthetic-kit",
        basecaller_name="synthetic-basecaller",
        basecaller_version="1.0.0",
        basecaller_model="synthetic-model",
        caller_name="nanopolish",
        caller_version="0.14.0",
        reference_id="synthetic-reference",
        reference_genome_build=GenomeBuild.GRCH38,
        reference_sha256=_sha("synthetic-reference"),
        adapter_name="nanopolish_call_methylation_tsv",
        adapter_version="nanopolish-call-table-v2",
        adapter_policy_sha256=content_sha256(policy),
    )
    cohort = PairedHoldoutCohort(
        samples={
            role: ValidationSample(
                accession=f"SYNTHETIC_{role}",
                dataset_version="synthetic-v1",
                biological_sample_id=f"SYNTHETIC_SAMPLE_{role}",
                donor_id="unknown",
                biological_source=f"synthetic-biology-{role[-1]}",
                cell_type=f"synthetic-cell-{role[-1]}",
                material_type="synthetic-DNA",
                input_sha256=_sha(role),
                provenance_reference="Synthetic software fixture only.",
                technical=technical,
            )
            for role in HOLDOUT_SOURCE_ROLES
        },
        source_distinction_evidence="Two synthetic labels; no actual biological sources.",
        access_basis="public",
    )
    return preregister_paired_holdout(
        matrix,
        cohort,
        PairedHoldoutSplitPolicy(
            calibration_read_group_fraction=0.25,
            source_a_split_seed=101,
            source_b_split_seed=202,
        ),
        registration_id="synthetic-holdout-registration",
        registered_at=datetime(2026, 9, 5, tzinfo=UTC),
        test_outcomes_unseen=True,
        code_sha256=_sha("synthetic-code"),
        software_version="synthetic-v1",
    )


def _source(
    role: str, seed: int, *, reverse_test: bool = False
) -> tuple[_ParsedSource, tuple[str, ...], tuple[str, ...]]:
    names = tuple(f"SYNTHETIC_{role}_{read:03d}" for read in range(160))
    calibration, test = _split_reads(names, fraction=0.25, seed=seed)
    calls = {}
    for is_test, pool in ((False, calibration), (True, test)):
        rate = 4 if role == "source_a" else 1
        if is_test and reverse_test:
            rate = 5 - rate
        for index, name in enumerate(pool):
            calls[name] = tuple(
                _NanopolishCall(
                    marker=_MarkerKey(
                        "chr1", "+", 100 + marker * 10, 100 + marker * 10, 1, _sha("synthetic-CG")
                    ),
                    state=(
                        NanopolishCallState.METHYLATED
                        if index % 5 < rate
                        else NanopolishCallState.UNMETHYLATED
                    ),
                )
                for marker in range(3)
            )
    source = _ParsedSource(
        summary=MethylationSourceSummary(
            source_id=f"synthetic-{role}",
            genome_build=GenomeBuild.GRCH38,
            fingerprint=FileFingerprint(size_bytes=480, sha256=_sha(role)),
            decoded_size_bytes=480,
            rows_total=480,
            canonical_call_rows=480,
            informative_call_rows=480,
            ambiguous_call_rows=0,
            skipped_noncanonical_rows=0,
            read_group_count=160,
        ),
        calls_by_read=calls,
    )
    return source, calibration, test


def _evidence(
    registration: PairedHoldoutRegistration, *, reverse_test: bool = False
) -> PairedHoldoutEvidence:
    source_a, calibration_a, test_a = _source("source_a", 101, reverse_test=reverse_test)
    source_b, calibration_b, test_b = _source("source_b", 202, reverse_test=reverse_test)
    markers = _reference_markers(
        _aggregate(source_a, calibration_a),
        _aggregate(source_b, calibration_b),
        source_a_total_read_groups=40,
        source_b_total_read_groups=40,
        policy=registration.matrix.estimator_policy,
    )
    levels = []
    for fraction in registration.matrix.source_a_fractions:
        for replicate, seed in enumerate(registration.matrix.replicate_seeds, 1):
            count_a = round(20 * fraction)
            level = _estimate_level(
                level_id=f"SYNTHETIC.F{fraction}.R{replicate}",
                replicate=replicate,
                seed=validation_level_seed(registration.matrix, 20, fraction, replicate),
                target_fraction=fraction,
                source_a_reads=test_a[:count_a],
                source_b_reads=test_b[: 20 - count_a],
                source_a=source_a,
                source_b=source_b,
                reference_markers=markers,
                policy=registration.matrix.estimator_policy,
            )
            levels.append(
                ValidationLevelEvidence(
                    read_group_budget=20,
                    replicate_seed=seed,
                    result=level,
                )
            )
    split_values = {
        "registration_sha256": registration.lock_sha256,
        "observed_input_sha256": {role: _sha(role) for role in HOLDOUT_SOURCE_ROLES},
        "source_read_group_counts": {role: 160 for role in HOLDOUT_SOURCE_ROLES},
        "calibration_read_group_counts": {role: 40 for role in HOLDOUT_SOURCE_ROLES},
        "test_read_group_counts": {role: 120 for role in HOLDOUT_SOURCE_ROLES},
        "pool_selection_sha256": {
            "calibration_a": _selection_digest(calibration_a, (), context="calibration_a"),
            "calibration_b": _selection_digest((), calibration_b, context="calibration_b"),
            "test_a": _selection_digest(test_a, (), context="test_a"),
            "test_b": _selection_digest((), test_b, context="test_b"),
        },
    }
    return PairedHoldoutEvidence(
        **split_values,
        split_manifest_sha256=holdout_split_manifest_sha256(**split_values),
        code_sha256=registration.code_sha256,
        software_version=registration.software_version,
        experiment_started_at=registration.registered_at + timedelta(seconds=1),
        experiment_completed_at=registration.registered_at + timedelta(seconds=2),
        reference_markers=markers,
        reference_marker_sha256=_reference_marker_digest(markers),
        levels=levels,
    )


def _refresh_split_manifest(evidence: PairedHoldoutEvidence) -> None:
    evidence.split_manifest_sha256 = holdout_split_manifest_sha256(
        registration_sha256=evidence.registration_sha256,
        observed_input_sha256=evidence.observed_input_sha256,
        source_read_group_counts=evidence.source_read_group_counts,
        calibration_read_group_counts=evidence.calibration_read_group_counts,
        test_read_group_counts=evidence.test_read_group_counts,
        pool_selection_sha256=evidence.pool_selection_sha256,
    )


class PairedHoldoutTests(unittest.TestCase):
    def test_missing_evidence_is_scoped_no_call(self) -> None:
        report = evaluate_paired_holdout(_registration())
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.no_call_rate, 1)
        self.assertEqual(report.donor_independent_validation, "NOT_EVALUATED")
        self.assertEqual(report.validation_scope, "held_out_read_group_recovery_within_two_sources")
        self.assertTrue(report.research_only)

    def test_two_sources_pass_without_donor_independence_claim(self) -> None:
        registration = _registration()
        self.assertEqual(holdout_cohort_eligibility_reasons(registration.cohort), [])
        report = evaluate_paired_holdout(registration, _evidence(registration))
        self.assertEqual(report.decision, ValidationDecision.PASS)
        self.assertEqual(report.donor_independent_validation, "NOT_EVALUATED")
        self.assertEqual(set(report.registration.cohort.samples), {"source_a", "source_b"})
        self.assertEqual(
            PairedHoldoutRecoveryReport.model_validate_json(report.model_dump_json()), report
        )

    def test_split_is_reproducible_disjoint_and_whole_read(self) -> None:
        source, calibration, test = _source("source_a", 101)
        self.assertEqual(set(calibration) | set(test), set(source.calls_by_read))
        self.assertFalse(set(calibration) & set(test))
        repeated = _split_reads(tuple(reversed(source.calls_by_read)), fraction=0.25, seed=101)
        self.assertEqual(repeated, (calibration, test))
        self.assertEqual(expected_calibration_read_groups(160, 0.25), 40)
        self.assertEqual(expected_calibration_read_groups(2, 0.01), 1)
        self.assertEqual(expected_calibration_read_groups(2, 0.99), 1)
        self.assertTrue(all(len(source.calls_by_read[x]) == 3 for x in calibration + test))

    def test_test_biology_does_not_select_markers_but_reversed_recovery_fails(self) -> None:
        registration = _registration()
        baseline, altered = _evidence(registration), _evidence(registration, reverse_test=True)
        self.assertEqual(baseline.reference_marker_sha256, altered.reference_marker_sha256)
        self.assertEqual(baseline.reference_markers, altered.reference_markers)
        report = evaluate_paired_holdout(registration, altered)
        self.assertEqual(report.decision, ValidationDecision.FAIL)
        self.assertGreater(report.overall_metrics.mean_absolute_error or 0, 0.6)

    def test_insufficient_markers_and_incomplete_grid_cannot_pass(self) -> None:
        registration = _registration(minimum_markers=5)
        report = evaluate_paired_holdout(registration, _evidence(registration))
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.explicit_no_call_levels, 30)
        registration = _registration()
        evidence = _evidence(registration)
        evidence.levels.pop()
        report = evaluate_paired_holdout(registration, evidence)
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.missing_levels, 1)
        self.assertAlmostEqual(report.overall_metrics.no_call_rate, 1 / 30)

    def test_split_seed_registration_and_pool_counts_cannot_be_changed(self) -> None:
        registration = _registration()
        registration.split_policy.source_a_split_seed = 999
        with self.assertRaises(ValidationError):
            evaluate_paired_holdout(registration)
        registration = _registration()
        evidence = _evidence(registration)
        evidence.calibration_read_group_counts["source_a"] = 41
        evidence.test_read_group_counts["source_a"] = 119
        _refresh_split_manifest(evidence)
        with self.assertRaisesRegex(ValueError, "split policy"):
            evaluate_paired_holdout(registration, evidence)
        evidence = _evidence(registration)
        evidence.test_read_group_counts["source_a"] = 119
        with self.assertRaisesRegex(ValueError, "partition"):
            evaluate_paired_holdout(registration, evidence)

    def test_numerical_tampering_reuses_registered_replay_checks(self) -> None:
        for field in ("seed", "interval", "duplicate"):
            with self.subTest(field=field):
                registration = _registration()
                evidence = _evidence(registration)
                if field == "seed":
                    evidence.levels[0].result.seed += 1
                elif field == "interval":
                    evidence.levels[0].result.confidence_interval_upper = 1.0
                else:
                    evidence.levels.append(evidence.levels[0])
                with self.assertRaises(ValueError):
                    evaluate_paired_holdout(registration, evidence)

    def test_input_provenance_or_late_registration_prevents_pass(self) -> None:
        for field in ("input", "code", "time"):
            with self.subTest(field=field):
                registration = _registration()
                evidence = _evidence(registration)
                if field == "input":
                    evidence.observed_input_sha256["source_a"] = _sha("wrong-input")
                    _refresh_split_manifest(evidence)
                elif field == "code":
                    evidence.code_sha256 = _sha("wrong-code")
                else:
                    evidence.experiment_started_at = registration.registered_at
                self.assertEqual(
                    evaluate_paired_holdout(registration, evidence).decision,
                    ValidationDecision.NO_CALL,
                )

    def test_saved_pool_selection_digest_tampering_is_rejected(self) -> None:
        registration = _registration()
        evidence = _evidence(registration)
        evidence.pool_selection_sha256["test_a"] = _sha("tampered-test-pool")
        with self.assertRaisesRegex(ValidationError, "split manifest digest"):
            evaluate_paired_holdout(registration, evidence)

    def test_scope_and_derived_report_decision_are_not_editable(self) -> None:
        registration = _registration()
        report = evaluate_paired_holdout(registration, _evidence(registration))
        payload = report.model_dump(mode="json")
        payload["donor_independent_validation"] = "PASS"
        with self.assertRaises(ValidationError):
            PairedHoldoutRecoveryReport.model_validate(payload)
        payload = evaluate_paired_holdout(registration).model_dump(mode="json")
        payload["decision"] = "PASS"
        with self.assertRaises(ValidationError):
            PairedHoldoutRecoveryReport.model_validate(payload)

    def test_identical_biological_inputs_and_parent_bam_are_ineligible(self) -> None:
        cohort = _registration().cohort
        cohort.samples["source_a"].source_modbam_sha256 = _sha("shared-parent")
        cohort.samples["source_b"].source_modbam_sha256 = _sha("shared-parent")
        self.assertIn("parent modBAM", " ".join(holdout_cohort_eligibility_reasons(cohort)))
        cohort = _registration().cohort
        cohort.samples["source_b"].biological_source = cohort.samples["source_a"].biological_source
        self.assertIn("biologically distinct", " ".join(holdout_cohort_eligibility_reasons(cohort)))

    def test_blank_distinction_and_unknown_metadata_are_ineligible(self) -> None:
        cohort = _registration().cohort
        cohort.source_distinction_evidence = "   "
        cohort.samples["source_a"].technical.basecaller_model = " UnKnOwN "
        reasons = " ".join(holdout_cohort_eligibility_reasons(cohort))
        self.assertIn("biological distinction of source A and B is missing", reasons)
        self.assertIn("source_a.technical.basecaller_model is unknown", reasons)
        with self.assertRaises(ValidationError):
            PairedHoldoutCohort.model_validate(cohort.model_dump(mode="json"))

    def test_invisible_source_identity_suffix_cannot_register(self) -> None:
        registration = _registration()
        cohort = registration.cohort
        cohort.samples["source_b"].biological_source = (
            cohort.samples["source_a"].biological_source + "\u200b"
        )
        with self.assertRaisesRegex(ValidationError, "control characters"):
            preregister_paired_holdout(
                registration.matrix,
                cohort,
                registration.split_policy,
                registration_id="synthetic-control-rejection",
                registered_at=registration.registered_at,
                test_outcomes_unseen=True,
                code_sha256=registration.code_sha256,
                software_version=registration.software_version,
            )

    def test_execution_insufficiency_is_reported_without_overriding_results(self) -> None:
        registration = _registration()
        report = evaluate_paired_holdout(
            registration, execution_no_call_reason="Held-out source cannot supply requested budget."
        )
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(
            PairedHoldoutRecoveryReport.model_validate_json(report.model_dump_json()), report
        )
        with self.assertRaises(ValueError):
            evaluate_paired_holdout(
                registration,
                _evidence(registration),
                execution_no_call_reason="Suppress completed evidence.",
            )


if __name__ == "__main__":
    unittest.main()

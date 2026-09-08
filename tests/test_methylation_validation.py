from __future__ import annotations

import hashlib
import unittest
from datetime import UTC, datetime, timedelta
from pathlib import Path

from pydantic import ValidationError

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
)
from ontseq_platform.methylation_validation import (
    SAMPLE_ROLES,
    MethylationValidationEvidence,
    MethylationValidationMatrix,
    MethylationValidationRegistration,
    MethylationValidationReport,
    RecoveryAcceptance,
    ValidationCohort,
    ValidationDecision,
    ValidationLevelEvidence,
    ValidationSample,
    ValidationTechnicalMetadata,
    cohort_eligibility_reasons,
    content_sha256,
    evaluate_methylation_validation,
    preregister_validation,
    validation_level_seed,
)
from ontseq_platform.models import FileFingerprint, GenomeBuild


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _matrix(*, minimum_markers: int = 2) -> MethylationValidationMatrix:
    fractions = [0.0, 0.5, 1.0]
    policy = MethylationMixturePolicy(
        profile_id="synthetic-independent-policy",
        source_a_fractions=fractions,
        replicates=10,
        minimum_reference_valid_calls=2,
        minimum_mixture_valid_calls=1,
        minimum_markers=minimum_markers,
        maximum_markers=10,
        uncertainty_draws=100,
        maximum_confidence_interval_width=1.0,
        maximum_standardized_model_fit_rmse=100.0,
        maximum_extrapolation=1.0,
        note="Synthetic software tests only; no real biological validation evidence.",
    )
    return MethylationValidationMatrix(
        matrix_id="synthetic-independent-matrix",
        source_a_fractions=fractions,
        read_group_budgets=[20],
        replicate_seeds=list(range(10)),
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
        estimator_policy=policy,
        rationale="Small deterministic fixture for software checks, not scientific evidence.",
    )


def _registration(
    matrix: MethylationValidationMatrix | None = None, *, outcomes_unseen: bool = True
) -> MethylationValidationRegistration:
    matrix = matrix or _matrix()
    metadata = ValidationTechnicalMetadata(
        sequencing_platform="Oxford Nanopore",
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
        adapter_policy_sha256=content_sha256(matrix.estimator_policy),
    )
    cohort = ValidationCohort(
        samples={
            role: ValidationSample(
                accession=f"SYNTHETIC_{role}",
                dataset_version="synthetic-v1",
                biological_sample_id=f"SYNTHETIC_SAMPLE_{role}",
                donor_id=f"SYNTHETIC_DONOR_{role}",
                biological_source=f"synthetic-source-{role[-1]}",
                cell_type=f"synthetic-cell-{role[-1]}",
                material_type="synthetic-native-DNA",
                input_sha256=_sha(role),
                provenance_reference="synthetic fixture, no human data",
                technical=metadata,
            )
            for role in SAMPLE_ROLES
        },
        independence_evidence="Four synthetic software fixtures; no real biological donors.",
        access_basis="public",
    )
    return preregister_validation(
        matrix,
        cohort,
        registration_id="synthetic-registration",
        registered_at=datetime(2026, 9, 5, tzinfo=UTC),
        test_outcomes_unseen=outcomes_unseen,
        code_sha256=_sha("synthetic-code"),
        software_version="synthetic-v1",
    )


def _source(role: str, *, rate: int) -> _ParsedSource:
    calls = {}
    for read in range(40):
        calls[f"SYNTHETIC_{role}_{read}"] = tuple(
            _NanopolishCall(
                marker=_MarkerKey(
                    "chr1", "+", 100 + marker * 10, 100 + marker * 10, 1, _sha("synthetic-CG")
                ),
                state=(
                    NanopolishCallState.METHYLATED
                    if read % 5 < rate
                    else NanopolishCallState.UNMETHYLATED
                ),
            )
            for marker in range(3)
        )
    return _ParsedSource(
        summary=MethylationSourceSummary(
            source_id=f"synthetic-{role}",
            genome_build=GenomeBuild.GRCH38,
            fingerprint=FileFingerprint(size_bytes=120, sha256=_sha(role)),
            decoded_size_bytes=120,
            rows_total=120,
            canonical_call_rows=120,
            informative_call_rows=120,
            ambiguous_call_rows=0,
            skipped_noncanonical_rows=0,
            read_group_count=40,
        ),
        calls_by_read=calls,
    )


def _evidence(
    registration: MethylationValidationRegistration, *, reversed_test: bool = False
) -> MethylationValidationEvidence:
    matrix = registration.matrix
    source_a, source_b = _source("calibration_a", rate=4), _source("calibration_b", rate=1)
    test_a = _source("test_a", rate=1 if reversed_test else 4)
    test_b = _source("test_b", rate=4 if reversed_test else 1)
    markers = _reference_markers(
        _aggregate(source_a, tuple(source_a.calls_by_read)),
        _aggregate(source_b, tuple(source_b.calls_by_read)),
        source_a_total_read_groups=40,
        source_b_total_read_groups=40,
        policy=matrix.estimator_policy,
    )
    levels = []
    for fraction in matrix.source_a_fractions:
        for replicate, seed in enumerate(matrix.replicate_seeds, 1):
            count_a = round(20 * fraction)
            level = _estimate_level(
                level_id=f"synthetic-F{fraction}-R{replicate}",
                replicate=replicate,
                seed=validation_level_seed(matrix, 20, fraction, replicate),
                target_fraction=fraction,
                source_a_reads=tuple(test_a.calls_by_read)[:count_a],
                source_b_reads=tuple(test_b.calls_by_read)[: 20 - count_a],
                source_a=test_a,
                source_b=test_b,
                reference_markers=markers,
                policy=matrix.estimator_policy,
            )
            levels.append(
                ValidationLevelEvidence(
                    read_group_budget=20,
                    replicate_seed=seed,
                    result=level,
                )
            )
    return MethylationValidationEvidence(
        registration_sha256=registration.lock_sha256,
        code_sha256=registration.code_sha256,
        software_version=registration.software_version,
        experiment_started_at=registration.registered_at + timedelta(seconds=1),
        experiment_completed_at=registration.registered_at + timedelta(seconds=2),
        observed_input_sha256={role: _sha(role) for role in SAMPLE_ROLES},
        read_group_counts={role: 40 for role in SAMPLE_ROLES},
        reference_markers=markers,
        reference_marker_sha256=_reference_marker_digest(markers),
        calibration_selection_sha256=_sha("synthetic-calibration-selection"),
        levels=levels,
    )


class MethylationValidationTests(unittest.TestCase):
    def test_repository_matrix_is_explicit_prospective_candidate(self) -> None:
        path = Path(__file__).parents[1] / "configs/methylation/independent_validation.v1.json"
        matrix = MethylationValidationMatrix.model_validate_json(path.read_text())
        self.assertEqual(matrix.expected_levels, 270)
        self.assertEqual(matrix.read_group_budgets, [1000, 5000, 10000])
        self.assertEqual(matrix.status, "research_acceptance_candidate_unvalidated")

    def test_missing_evidence_is_unconditional_no_call(self) -> None:
        report = evaluate_methylation_validation(_registration())
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.missing_levels, 30)
        self.assertEqual(report.overall_metrics.no_call_rate, 1.0)
        self.assertIsNone(report.overall_metrics.mean_absolute_error)
        self.assertTrue(report.registration.sensitive_output)
        self.assertTrue(report.research_only)

    def test_synthetic_recovery_and_json_round_trip(self) -> None:
        registration = _registration()
        report = evaluate_methylation_validation(registration, _evidence(registration))
        self.assertEqual(report.decision, ValidationDecision.PASS)
        self.assertLess(report.overall_metrics.mean_absolute_error or 0, 1e-10)
        self.assertEqual(report.overall_metrics.no_call_rate, 0)
        self.assertEqual(
            MethylationValidationReport.model_validate_json(report.model_dump_json()), report
        )

    def test_execution_insufficiency_is_an_explicit_replayable_no_call(self) -> None:
        registration = _registration()
        report = evaluate_methylation_validation(
            registration, execution_no_call_reason="Insufficient retained test read groups."
        )
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertIn("Insufficient retained test read groups.", report.reasons)
        self.assertEqual(
            MethylationValidationReport.model_validate_json(report.model_dump_json()), report
        )
        with self.assertRaises(ValueError):
            evaluate_methylation_validation(
                registration, _evidence(registration), execution_no_call_reason="Hide results."
            )

    def test_reversed_biology_fails_even_with_zero_aggregate_bias(self) -> None:
        registration = _registration()
        report = evaluate_methylation_validation(
            registration, _evidence(registration, reversed_test=True)
        )
        self.assertEqual(report.decision, ValidationDecision.FAIL)
        self.assertLess(abs(report.overall_metrics.mean_bias or 0), 1e-10)
        self.assertGreater(report.overall_metrics.mean_absolute_error or 0, 0.6)
        self.assertTrue(any(x.decision == ValidationDecision.FAIL for x in report.strata))

    def test_incomplete_grid_cannot_pass_and_missing_levels_count(self) -> None:
        registration = _registration()
        evidence = _evidence(registration)
        evidence.levels.pop()
        report = evaluate_methylation_validation(registration, evidence)
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.missing_levels, 1)
        self.assertAlmostEqual(report.overall_metrics.no_call_rate, 1 / 30)
        self.assertLess(report.overall_metrics.mean_absolute_error or 0, 1e-10)

    def test_insufficient_markers_are_explicit_no_call(self) -> None:
        registration = _registration(_matrix(minimum_markers=5))
        report = evaluate_methylation_validation(registration, _evidence(registration))
        self.assertEqual(report.decision, ValidationDecision.NO_CALL)
        self.assertEqual(report.overall_metrics.explicit_no_call_levels, 30)
        self.assertIsNone(report.overall_metrics.conditional_interval_coverage)

    def test_read_split_same_donor_and_unknown_protocol_are_ineligible(self) -> None:
        cohort = _registration().cohort
        cohort.samples["test_a"].donor_id = cohort.samples["calibration_a"].donor_id
        cohort.samples["test_b"].biological_sample_id = cohort.samples[
            "calibration_b"
        ].biological_sample_id
        cohort.samples["test_b"].technical.caller_version = "unknown"
        reasons = " ".join(cohort_eligibility_reasons(cohort))
        self.assertIn("donors overlap", reasons)
        self.assertIn("read split", reasons)
        self.assertIn("unknown", reasons)

    def test_biological_identities_reject_invisible_controls_before_normalization(self) -> None:
        sample = _registration().cohort.samples["test_a"].model_dump(mode="json")
        for field in ("biological_sample_id", "donor_id", "biological_source", "cell_type"):
            for control in ("\x00", "\n", "\x7f", "\u200b"):
                with self.subTest(field=field, control=repr(control)):
                    payload = {**sample, field: "SYNTHETIC" + control + "IDENTITY"}
                    with self.assertRaisesRegex(ValidationError, "control characters"):
                        ValidationSample.model_validate(payload)

    def test_biological_identity_whitespace_normalization_and_case_are_preserved(self) -> None:
        sample = _registration().cohort.samples["test_a"].model_dump(mode="json")
        sample["donor_id"] = "  Synthetic_Donor_A  "
        self.assertEqual(ValidationSample.model_validate(sample).donor_id, "Synthetic_Donor_A")

    def test_unknown_metadata_helpers_normalize_mutated_values(self) -> None:
        cohort = _registration().cohort
        cohort.independence_evidence = "   "
        cohort.samples["test_a"].technical.basecaller_model = " UnKnOwN "
        reasons = " ".join(cohort_eligibility_reasons(cohort))
        self.assertIn("independence evidence is missing", reasons)
        self.assertIn("test_a.technical.basecaller_model is unknown", reasons)
        with self.assertRaises(ValidationError):
            ValidationCohort.model_validate(cohort.model_dump(mode="json"))

    def test_posthoc_registration_and_wrong_code_or_inputs_cannot_pass(self) -> None:
        for corruption in ("seen", "time", "code", "input"):
            with self.subTest(corruption=corruption):
                registration = _registration(outcomes_unseen=corruption != "seen")
                evidence = _evidence(registration)
                if corruption == "time":
                    evidence.experiment_started_at = registration.registered_at
                elif corruption == "code":
                    evidence.code_sha256 = _sha("changed-code")
                elif corruption == "input":
                    evidence.observed_input_sha256["test_a"] = _sha("changed-input")
                report = evaluate_methylation_validation(registration, evidence)
                self.assertEqual(report.decision, ValidationDecision.NO_CALL)

    def test_reused_parent_modbam_is_ineligible_even_with_distinct_exports(self) -> None:
        cohort = _registration().cohort
        cohort.samples["calibration_a"].source_modbam_sha256 = _sha("synthetic-parent-modbam")
        cohort.samples["test_a"].source_modbam_sha256 = _sha("synthetic-parent-modbam")
        self.assertIn("provenance overlaps", " ".join(cohort_eligibility_reasons(cohort)))
        cohort.samples["test_a"].source_modbam_sha256 = cohort.samples["calibration_b"].input_sha256
        self.assertIn("provenance overlaps", " ".join(cohort_eligibility_reasons(cohort)))

    def test_registration_matrix_and_cohort_tampering_are_rejected(self) -> None:
        for field in ("matrix", "cohort", "code"):
            with self.subTest(field=field):
                registration = _registration()
                if field == "matrix":
                    registration.matrix.acceptance[0].maximum_mae = 1
                elif field == "cohort":
                    registration.cohort.samples["test_a"].donor_id = "changed-donor"
                else:
                    registration.code_sha256 = _sha("changed-code")
                with self.assertRaises(ValidationError):
                    evaluate_methylation_validation(registration)

    def test_seed_interval_point_and_duplicate_tampering_are_rejected(self) -> None:
        for corruption in ("seed", "interval", "point", "duplicate"):
            with self.subTest(corruption=corruption):
                registration = _registration()
                evidence = _evidence(registration)
                level = evidence.levels[0].result
                if corruption == "seed":
                    level.seed += 1
                elif corruption == "interval":
                    level.confidence_interval_upper = 1.0
                elif corruption == "point":
                    level.estimated_source_a_fraction = 0.25
                else:
                    evidence.levels.append(evidence.levels[0])
                with self.assertRaises(ValueError):
                    evaluate_methylation_validation(registration, evidence)

    def test_derived_metrics_or_decision_cannot_be_overwritten_on_reload(self) -> None:
        registration = _registration()
        report = evaluate_methylation_validation(registration, _evidence(registration))
        payload = report.model_dump(mode="json")
        payload["overall_metrics"]["mean_absolute_error"] = 0.2
        with self.assertRaises(ValidationError):
            MethylationValidationReport.model_validate(payload)
        payload = evaluate_methylation_validation(registration).model_dump(mode="json")
        payload["decision"] = "PASS"
        with self.assertRaises(ValidationError):
            MethylationValidationReport.model_validate(payload)


if __name__ == "__main__":
    unittest.main()

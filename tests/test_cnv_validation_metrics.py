from __future__ import annotations

import hashlib
import json
import unittest
from datetime import UTC, datetime

from pydantic import ValidationError

from ontseq_platform.cnv_validation_contracts import (
    CnvAcceptanceMetric,
    CnvAcceptanceQuestion,
    CnvAssessabilityMask,
    CnvCallerLock,
    CnvDataBasis,
    CnvNegativeUniverse,
    CnvQuantitativeTruth,
    CnvRepeatKind,
    CnvStratificationDimension,
    CnvStratificationPlan,
    CnvTruthSource,
    CnvValidationCohort,
    CnvValidationMatrix,
    CnvValidationSpecimen,
)
from ontseq_platform.cnv_validation_evidence import (
    CnvContributionStatus,
    CnvEvidenceRecordKind,
    CnvFullEvidenceRecord,
    CnvNativeArtifactReference,
    CnvNormalizedEventRecord,
    CnvRunOutcomeState,
    CnvValidationEvidenceManifest,
    canonical_evidence_sha256,
    seal_cnv_validation_evidence,
)
from ontseq_platform.cnv_validation_metrics import (
    CnvLaneEventAssignment,
    CnvMetricEvaluationState,
    CnvMetricResult,
    CnvNegativeUnitAssessment,
    CnvAcceptanceDecision,
    CnvAcceptanceResult,
    CnvReproducibilityPair,
    CnvValidationMetricReport,
    aggregate_cnv_event_metrics,
    aggregate_cnv_quantitative_metrics,
    aggregate_cnv_reproducibility_metrics,
    aggregate_cnv_run_state_metrics,
    aggregate_cnv_validation,
    assign_cnv_validation_lanes,
    cnv_numeric_band,
)
from ontseq_platform.cnv_validation_registration import (
    CnvValidationRegistration,
    preregister_cnv_validation,
)
from ontseq_platform.models import (
    BenchmarkThresholds,
    EventType,
    GenomeBuild,
    GenomicEvent,
    Locus,
)


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _matrix() -> CnvValidationMatrix:
    return CnvValidationMatrix(
        matrix_id="synthetic-cnv-metrics-v1",
        callers=[
            CnvCallerLock(
                caller_id="qdnaseq_ace",
                caller_version="synthetic-1",
                adapter_policy_sha256=_sha("qdnaseq-policy"),
                execution_identity_sha256=_sha("qdnaseq-runtime"),
            )
        ],
        genome_builds=[GenomeBuild.GRCH38],
        data_bases=[CnvDataBasis.LCWGS_GENOME_WIDE],
        stratification=CnvStratificationPlan(
            primary_dimensions=[
                CnvStratificationDimension.COVERAGE,
                CnvStratificationDimension.TUMOR_FRACTION,
                CnvStratificationDimension.EVENT_CLASS,
                CnvStratificationDimension.EVENT_SIZE,
            ],
            secondary_dimensions=[
                CnvStratificationDimension.CALLER,
                CnvStratificationDimension.GENOME_BUILD,
                CnvStratificationDimension.DATA_BASIS,
                CnvStratificationDimension.BIN_SIZE,
            ],
            coverage_cutpoints_x=[2.0, 5.0, 10.0],
            tumor_fraction_cutpoints=[0.1, 0.2, 0.5],
            event_size_cutpoints_bp=[1_000_000, 5_000_000, 10_000_000],
        ),
        event_classes=[EventType.DELETION, EventType.DUPLICATION],
        qdnaseq_bin_sizes_kbp=[100, 500, 1000],
        matching_thresholds=BenchmarkThresholds(
            minimum_reciprocal_overlap=0.5,
            copy_number_tolerance=0.25,
        ),
        acceptance=[
            CnvAcceptanceQuestion(
                question_id="primary-sensitivity",
                metric=CnvAcceptanceMetric.SENSITIVITY,
                minimum_evaluable_denominator=1,
                minimum_acceptable=0.9,
            )
        ],
        rationale="Synthetic metric-assignment contract without biological performance claim.",
    )


def _truth_event() -> GenomicEvent:
    return GenomicEvent(
        event_id="truth-del-1",
        event_type=EventType.DELETION,
        primary=Locus(chromosome="7", start=1_000_000, end=6_000_000),
        copy_number=1.0,
    )


def _specimen() -> CnvValidationSpecimen:
    mask_sha = _sha("mask")
    return CnvValidationSpecimen(
        specimen_id="SYNTHETIC_CNV_001",
        biological_specimen_id="SYNTHETIC_BIO_001",
        repeat_kind=CnvRepeatKind.INDEPENDENT,
        access_basis="public",
        material_type="synthetic-DNA",
        genome_build=GenomeBuild.GRCH38,
        data_basis=CnvDataBasis.LCWGS_GENOME_WIDE,
        reference_id="synthetic-grch38",
        reference_sha256=_sha("reference"),
        input_sha256=_sha("input"),
        coverage_x=5.25,
        coverage_definition="synthetic mean autosomal depth",
        tumor_fraction=0.21,
        tumor_fraction_method="synthetic orthogonal fraction",
        tumor_fraction_timepoint="same synthetic aliquot",
        truth_sources=[
            CnvTruthSource(
                method_name="synthetic-karyotype",
                method_version="v1",
                resource_id="synthetic-truth-v1",
                resource_sha256=_sha("truth"),
                provenance_reference="synthetic fixture only",
            )
        ],
        assessability_mask=CnvAssessabilityMask(
            resource_id="synthetic-mask-v1",
            resource_sha256=mask_sha,
            unit="regions",
            definition="Synthetic locked assessability mask for metrics contract tests.",
        ),
        truth_events=[_truth_event()],
        negative_universe=CnvNegativeUniverse(
            universe_id="synthetic-negative-bins",
            unit="genomic_bins",
            assessable_units=100,
            resource_sha256=_sha("negative-universe"),
            assessability_mask_sha256=mask_sha,
            definition="Synthetic negative-bin universe for metrics contract tests.",
        ),
    )


def _registration() -> CnvValidationRegistration:
    return preregister_cnv_validation(
        _matrix(),
        CnvValidationCohort(specimens=[_specimen()]),
        registration_id="synthetic-cnv-metrics-registration",
        registered_at=datetime(2026, 9, 18, tzinfo=UTC),
        code_sha256=_sha("code"),
        software_version="0.8.2",
    )


def _parameters() -> dict[str, object]:
    return {"bin_size_kbp": 500, "ace_penalty": 0.6}


def _artifact() -> CnvNativeArtifactReference:
    return CnvNativeArtifactReference(
        artifact_id="segments-500",
        role="caller_segments",
        relative_path="cnv/qdnaseq/500/segments.tsv",
        sha256=_sha("segments"),
        size_bytes=123,
        media_type="text/tab-separated-values",
    )


def _full_record(
    registration: CnvValidationRegistration,
    *,
    record_id: str,
    native_record_id: str,
    record_kind: CnvEvidenceRecordKind,
    run_outcome: CnvRunOutcomeState = CnvRunOutcomeState.OBSERVED,
    contribution_status: CnvContributionStatus = CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
    native_artifact_ids: list[str] | None = None,
    outcome_reason: str | None = None,
    event_type: EventType | None = None,
    primary: Locus | None = None,
    copy_number: float | None = None,
    caller_version: str = "synthetic-1",
    truth_support_resource_ids: list[str] | None = None,
    specimen_id: str | None = None,
    replicate_id: str = "replicate-1",
    fit_group_id: str | None = None,
    selected_fit: bool | None = None,
    cellularity: float | None = None,
    ploidy: float | None = None,
    fit_error: float | None = None,
) -> CnvFullEvidenceRecord:
    specimen = (
        registration.cohort.specimens[0]
        if specimen_id is None
        else next(item for item in registration.cohort.specimens if item.specimen_id == specimen_id)
    )
    parameters = _parameters()
    return CnvFullEvidenceRecord(
        record_id=record_id,
        registration_sha256=registration.lock_sha256,
        specimen_id=specimen.specimen_id,
        biological_specimen_id=specimen.biological_specimen_id,
        caller_id="qdnaseq_ace",
        caller_version=caller_version,
        adapter_policy_sha256=_sha("qdnaseq-policy"),
        execution_identity_sha256=_sha("qdnaseq-runtime"),
        genome_build=specimen.genome_build,
        data_basis=specimen.data_basis,
        reference_id=specimen.reference_id,
        reference_sha256=specimen.reference_sha256,
        input_sha256=specimen.input_sha256,
        coverage_x=specimen.coverage_x,
        coverage_definition=specimen.coverage_definition,
        tumor_fraction=specimen.tumor_fraction,
        tumor_fraction_method=specimen.tumor_fraction_method,
        tumor_fraction_timepoint=specimen.tumor_fraction_timepoint,
        bin_size_kbp=500,
        repeat_kind=specimen.repeat_kind,
        replicate_id=replicate_id,
        repeat_group_id=specimen.repeat_group_id,
        record_kind=record_kind,
        native_record_id=native_record_id,
        run_outcome=run_outcome,
        contribution_status=contribution_status,
        outcome_reason=outcome_reason,
        truth_support_resource_ids=truth_support_resource_ids or [],
        caller_parameters=parameters,
        caller_parameters_sha256=canonical_evidence_sha256(parameters),
        dependency_versions={"QDNAseq": "synthetic-1", "ACE": "synthetic-1"},
        native_artifact_ids=native_artifact_ids or [],
        fit_group_id=fit_group_id,
        selected_fit=selected_fit,
        cellularity=cellularity,
        ploidy=ploidy,
        fit_error=fit_error,
        event_type=event_type,
        primary=primary,
        copy_number=copy_number,
    )


def _run_summary(
    registration: CnvValidationRegistration,
    **updates: object,
) -> CnvFullEvidenceRecord:
    values: dict[str, object] = {
        "record_id": "run-summary-500",
        "native_record_id": "run-summary-500",
        "record_kind": CnvEvidenceRecordKind.RUN_SUMMARY,
    }
    values.update(updates)
    return _full_record(registration, **values)


def _event_source(
    registration: CnvValidationRegistration,
    **updates: object,
) -> CnvFullEvidenceRecord:
    values: dict[str, object] = {
        "record_id": "segment-500-1",
        "native_record_id": "segment-500-1",
        "record_kind": CnvEvidenceRecordKind.CALLER_SEGMENT,
        "native_artifact_ids": ["segments-500"],
        "event_type": EventType.DELETION,
        "primary": Locus(chromosome="7", start=1_100_000, end=6_100_000),
        "copy_number": 1.1,
    }
    values.update(updates)
    return _full_record(registration, **values)


def _normalized_event(
    registration: CnvValidationRegistration,
    *,
    specimen_id: str | None = None,
    replicate_id: str = "replicate-1",
    **updates: object,
) -> CnvNormalizedEventRecord:
    specimen = (
        registration.cohort.specimens[0]
        if specimen_id is None
        else next(item for item in registration.cohort.specimens if item.specimen_id == specimen_id)
    )
    values: dict[str, object] = {
        "normalized_event_id": "normalized-del-1",
        "registration_sha256": registration.lock_sha256,
        "specimen_id": specimen.specimen_id,
        "biological_specimen_id": specimen.biological_specimen_id,
        "caller_id": "qdnaseq_ace",
        "caller_version": "synthetic-1",
        "adapter_policy_sha256": _sha("qdnaseq-policy"),
        "execution_identity_sha256": _sha("qdnaseq-runtime"),
        "normalization_policy_sha256": _sha("normalization"),
        "genome_build": specimen.genome_build,
        "data_basis": specimen.data_basis,
        "reference_id": specimen.reference_id,
        "reference_sha256": specimen.reference_sha256,
        "input_sha256": specimen.input_sha256,
        "coverage_x": specimen.coverage_x,
        "coverage_definition": specimen.coverage_definition,
        "tumor_fraction": specimen.tumor_fraction,
        "tumor_fraction_method": specimen.tumor_fraction_method,
        "tumor_fraction_timepoint": specimen.tumor_fraction_timepoint,
        "bin_size_kbp": 500,
        "repeat_kind": specimen.repeat_kind,
        "replicate_id": replicate_id,
        "repeat_group_id": specimen.repeat_group_id,
        "source_full_evidence_ids": ["segment-500-1"],
        "event_type": EventType.DELETION,
        "primary": Locus(chromosome="7", start=1_100_000, end=6_100_000),
        "normalized_copy_number": 1.1,
        "contribution_status": CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
    }
    values.update(updates)
    return CnvNormalizedEventRecord.model_validate(values)


def _manifest(
    registration: CnvValidationRegistration,
    *,
    full_evidence: list[CnvFullEvidenceRecord] | None = None,
    normalized_events: list[CnvNormalizedEventRecord] | None = None,
) -> CnvValidationEvidenceManifest:
    records = full_evidence or [_run_summary(registration), _event_source(registration)]
    artifacts = [_artifact()] if any(item.native_artifact_ids for item in records) else []
    return seal_cnv_validation_evidence(
        manifest_id="SYNTHETIC_CNV_METRICS_EVIDENCE",
        registration_sha256=registration.lock_sha256,
        native_artifacts=artifacts,
        full_evidence=records,
        normalized_events=(
            [_normalized_event(registration)] if normalized_events is None else normalized_events
        ),
    )


class CnvLaneAssignmentTests(unittest.TestCase):
    def test_manifest_must_bind_exact_registration_lock(self) -> None:
        registration = _registration()
        evidence = _manifest(registration)
        payload = evidence.model_dump(mode="json")
        different_registration = _sha("different-registration")
        payload["registration_sha256"] = different_registration
        for record in payload["full_evidence"]:
            record["registration_sha256"] = different_registration
        for event in payload["normalized_events"]:
            event["registration_sha256"] = different_registration
        payload["manifest_sha256"] = canonical_evidence_sha256(
            {key: value for key, value in payload.items() if key != "manifest_sha256"}
        )
        mismatched = CnvValidationEvidenceManifest.model_validate(payload)
        with self.assertRaises(ValueError):
            assign_cnv_validation_lanes(registration, mismatched)

    def test_registered_caller_identity_is_fail_closed(self) -> None:
        registration = _registration()
        evidence = _manifest(
            registration,
            full_evidence=[
                _run_summary(registration, caller_version="synthetic-2"),
                _event_source(registration, caller_version="synthetic-2"),
            ],
            normalized_events=[_normalized_event(registration, caller_version="synthetic-2")],
        )
        with self.assertRaises(ValueError):
            assign_cnv_validation_lanes(registration, evidence)

    def test_lane_requires_exactly_one_run_summary(self) -> None:
        registration = _registration()
        without_summary = _manifest(
            registration,
            full_evidence=[_event_source(registration)],
        )
        with self.assertRaises(ValueError):
            assign_cnv_validation_lanes(registration, without_summary)

        duplicate_summary = _manifest(
            registration,
            full_evidence=[
                _run_summary(registration),
                _run_summary(
                    registration,
                    record_id="run-summary-500-b",
                    native_record_id="run-summary-500-b",
                ),
                _event_source(registration),
            ],
        )
        with self.assertRaises(ValueError):
            assign_cnv_validation_lanes(registration, duplicate_summary)

    def test_observed_lane_uses_registered_matching_thresholds_once(self) -> None:
        registration = _registration()
        assignments = assign_cnv_validation_lanes(
            registration,
            _manifest(registration),
        )
        self.assertEqual(len(assignments), 1)
        assignment = assignments[0]
        self.assertIsInstance(assignment, CnvLaneEventAssignment)
        self.assertEqual(assignment.run_outcome, CnvRunOutcomeState.OBSERVED)
        self.assertEqual(assignment.true_positive, 1)
        self.assertEqual(assignment.false_positive, 0)
        self.assertEqual(assignment.false_negative, 0)
        self.assertEqual(len(assignment.matches), 1)
        self.assertEqual(assignment.matches[0].truth_event_id, "truth-del-1")
        self.assertEqual(
            assignment.matches[0].normalized_event_id,
            "normalized-del-1",
        )
        self.assertEqual(
            assignment.minimum_reciprocal_overlap,
            registration.matrix.matching_thresholds.minimum_reciprocal_overlap,
        )

    def test_secondary_event_cannot_steal_primary_truth_match(self) -> None:
        registration = _registration()
        secondary_source = _event_source(
            registration,
            record_id="segment-500-secondary",
            native_record_id="segment-500-secondary",
            contribution_status=CnvContributionStatus.SECONDARY,
            primary=Locus(chromosome="7", start=1_000_000, end=6_000_000),
            copy_number=1.0,
        )
        secondary_event = _normalized_event(
            registration,
            normalized_event_id="normalized-del-secondary",
            source_full_evidence_ids=["segment-500-secondary"],
            contribution_status=CnvContributionStatus.SECONDARY,
            primary=Locus(chromosome="7", start=1_000_000, end=6_000_000),
            normalized_copy_number=1.0,
        )
        evidence = _manifest(
            registration,
            full_evidence=[
                _run_summary(registration),
                _event_source(registration),
                secondary_source,
            ],
            normalized_events=[
                _normalized_event(registration),
                secondary_event,
            ],
        )
        assignment = assign_cnv_validation_lanes(registration, evidence)[0]
        self.assertEqual(assignment.true_positive, 1)
        self.assertEqual(assignment.false_positive, 0)
        self.assertEqual(
            [item.normalized_event_id for item in assignment.matches],
            ["normalized-del-1"],
        )

    def test_terminal_lane_does_not_fabricate_event_assignments(self) -> None:
        registration = _registration()
        for outcome, contribution in (
            (CnvRunOutcomeState.FAILED, CnvContributionStatus.FAILED),
            (CnvRunOutcomeState.NO_CALL, CnvContributionStatus.NO_CALL),
            (
                CnvRunOutcomeState.NOT_ASSESSABLE,
                CnvContributionStatus.NOT_ASSESSABLE,
            ),
        ):
            with self.subTest(outcome=outcome):
                evidence = _manifest(
                    registration,
                    full_evidence=[
                        _run_summary(
                            registration,
                            run_outcome=outcome,
                            contribution_status=contribution,
                            outcome_reason="Synthetic terminal lane state.",
                        )
                    ],
                    normalized_events=[],
                )
                assignment = assign_cnv_validation_lanes(registration, evidence)[0]
                self.assertEqual(assignment.run_outcome, outcome)
                self.assertIsNone(assignment.true_positive)
                self.assertEqual(assignment.matches, [])

    def test_biological_negative_is_not_treated_as_successful_caller_lane(self) -> None:
        registration = _registration()
        run = _run_summary(
            registration,
            run_outcome=CnvRunOutcomeState.BIOLOGICAL_NEGATIVE_ORTHOGONAL_TRUTH,
            contribution_status=CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
            outcome_reason="Synthetic orthogonal biological negative truth.",
            truth_support_resource_ids=["synthetic-truth-v1"],
        )
        evidence = _manifest(
            registration,
            full_evidence=[run],
            normalized_events=[],
        )
        assignment = assign_cnv_validation_lanes(registration, evidence)[0]
        self.assertEqual(
            assignment.run_outcome,
            CnvRunOutcomeState.BIOLOGICAL_NEGATIVE_ORTHOGONAL_TRUTH,
        )
        self.assertIsNone(assignment.true_positive)
        self.assertEqual(assignment.matches, [])


def _registration_with_truth_event(event: GenomicEvent) -> CnvValidationRegistration:
    specimen = _specimen().model_copy(update={"truth_events": [event]})
    return preregister_cnv_validation(
        _matrix(),
        CnvValidationCohort(specimens=[specimen]),
        registration_id="synthetic-cnv-metrics-registration-custom-truth",
        registered_at=datetime(2026, 9, 18, tzinfo=UTC),
        code_sha256=_sha("code"),
        software_version="0.8.2",
    )


def _find_metric(
    metrics: list[CnvMetricResult],
    metric: CnvAcceptanceMetric,
    **stratum: str | int | float | bool,
) -> CnvMetricResult:
    matches = [
        item
        for item in metrics
        if item.metric == metric
        and all(item.stratum_key.get(key) == value for key, value in stratum.items())
    ]
    if len(matches) != 1:
        raise AssertionError(
            f"Expected exactly one {metric.value} metric for {stratum}, found {len(matches)}"
        )
    return matches[0]


class CnvProspectiveMetricAggregationTests(unittest.TestCase):
    def test_numeric_band_boundaries_are_left_closed_right_open(self) -> None:
        cutpoints = [2.0, 5.0, 10.0]
        self.assertEqual(cnv_numeric_band(1.999, cutpoints), "<2")
        self.assertEqual(cnv_numeric_band(2.0, cutpoints), "[2,5)")
        self.assertEqual(cnv_numeric_band(5.0, cutpoints), "[5,10)")
        self.assertEqual(cnv_numeric_band(10.0, cutpoints), ">=10")

    def test_overall_and_fully_stratified_event_metrics_are_emitted(self) -> None:
        registration = _registration()
        evidence = _manifest(registration)
        metrics = aggregate_cnv_event_metrics(registration, evidence)

        overall_sensitivity = _find_metric(
            metrics,
            CnvAcceptanceMetric.SENSITIVITY,
            scope="overall",
        )
        overall_precision = _find_metric(
            metrics,
            CnvAcceptanceMetric.PRECISION,
            scope="overall",
        )
        overall_f1 = _find_metric(
            metrics,
            CnvAcceptanceMetric.F1,
            scope="overall",
        )
        self.assertEqual(overall_sensitivity.value, 1.0)
        self.assertEqual(overall_precision.value, 1.0)
        self.assertEqual(overall_f1.value, 1.0)

        stratified = _find_metric(
            metrics,
            CnvAcceptanceMetric.SENSITIVITY,
            caller="qdnaseq_ace",
            genome_build="GRCh38",
            data_basis="lcwgs_genome_wide",
            coverage="[5,10)",
            tumor_fraction="[0.2,0.5)",
            event_class="deletion",
            event_size="[5000000,10000000)",
            bin_size=500,
        )
        self.assertEqual(stratified.value, 1.0)
        self.assertEqual(stratified.evaluable_denominator, 1)

    def test_registered_empty_event_stratum_is_not_evaluable_not_zero(self) -> None:
        registration = _registration()
        metrics = aggregate_cnv_event_metrics(registration, _manifest(registration))
        empty = _find_metric(
            metrics,
            CnvAcceptanceMetric.SENSITIVITY,
            caller="qdnaseq_ace",
            genome_build="GRCh38",
            data_basis="lcwgs_genome_wide",
            coverage="[5,10)",
            tumor_fraction="[0.2,0.5)",
            event_class="duplication",
            event_size="[5000000,10000000)",
            bin_size=500,
        )
        self.assertEqual(empty.state, CnvMetricEvaluationState.NOT_EVALUABLE)
        self.assertIsNone(empty.value)
        self.assertEqual(empty.evaluable_denominator, 0)

    def test_matching_is_not_repeated_inside_event_size_strata(self) -> None:
        truth = GenomicEvent(
            event_id="truth-del-cross-boundary",
            event_type=EventType.DELETION,
            primary=Locus(chromosome="7", start=1_000_000, end=5_900_000),
            copy_number=1.0,
        )
        registration = _registration_with_truth_event(truth)
        source = _event_source(
            registration,
            primary=Locus(chromosome="7", start=900_000, end=6_000_000),
            copy_number=1.0,
        )
        event = _normalized_event(
            registration,
            primary=Locus(chromosome="7", start=900_000, end=6_000_000),
            normalized_copy_number=1.0,
        )
        evidence = _manifest(
            registration,
            full_evidence=[_run_summary(registration), source],
            normalized_events=[event],
        )
        metrics = aggregate_cnv_event_metrics(registration, evidence)

        recall_side = _find_metric(
            metrics,
            CnvAcceptanceMetric.SENSITIVITY,
            caller="qdnaseq_ace",
            genome_build="GRCh38",
            data_basis="lcwgs_genome_wide",
            coverage="[5,10)",
            tumor_fraction="[0.2,0.5)",
            event_class="deletion",
            event_size="[1000000,5000000)",
            bin_size=500,
        )
        precision_side = _find_metric(
            metrics,
            CnvAcceptanceMetric.PRECISION,
            caller="qdnaseq_ace",
            genome_build="GRCh38",
            data_basis="lcwgs_genome_wide",
            coverage="[5,10)",
            tumor_fraction="[0.2,0.5)",
            event_class="deletion",
            event_size="[5000000,10000000)",
            bin_size=500,
        )
        self.assertEqual(recall_side.value, 1.0)
        self.assertEqual(precision_side.value, 1.0)
        self.assertEqual(recall_side.component_counts["tp_truth"], 1)
        self.assertEqual(recall_side.component_counts["fn"], 0)
        self.assertEqual(precision_side.component_counts["tp_query"], 1)
        self.assertEqual(precision_side.component_counts["fp"], 0)

    def test_metric_provenance_keeps_source_ids_without_mutating_evidence(self) -> None:
        registration = _registration()
        evidence = _manifest(registration)
        locked_sha = evidence.manifest_sha256
        metrics = aggregate_cnv_event_metrics(registration, evidence)
        result = _find_metric(
            metrics,
            CnvAcceptanceMetric.SENSITIVITY,
            scope="overall",
        )
        self.assertEqual(evidence.manifest_sha256, locked_sha)
        self.assertEqual(result.registration_sha256, registration.lock_sha256)
        self.assertEqual(result.evidence_manifest_sha256, evidence.manifest_sha256)
        self.assertEqual(result.truth_event_ids, ["truth-del-1"])
        self.assertEqual(result.normalized_event_ids, ["normalized-del-1"])
        self.assertIn("run-summary-500", result.full_evidence_ids)
        self.assertIn("segment-500-1", result.full_evidence_ids)
        self.assertEqual(len(result.lane_ids), 1)

    def test_false_positive_burden_uses_observed_lane_denominator(self) -> None:
        registration = _registration()
        metrics = aggregate_cnv_event_metrics(registration, _manifest(registration))
        result = _find_metric(
            metrics,
            CnvAcceptanceMetric.FALSE_POSITIVE_BURDEN,
            scope="overall",
        )
        self.assertEqual(result.value, 0.0)
        self.assertEqual(result.evaluable_denominator, 1)
        self.assertEqual(result.component_counts["fp"], 0)


class CnvRunStateAndSpecificityTests(unittest.TestCase):
    def test_run_state_rates_keep_terminal_outcomes_separate(self) -> None:
        registration = _registration()
        cases = (
            (
                CnvRunOutcomeState.NO_CALL,
                CnvContributionStatus.NO_CALL,
                CnvAcceptanceMetric.NO_CALL_RATE,
            ),
            (
                CnvRunOutcomeState.FAILED,
                CnvContributionStatus.FAILED,
                CnvAcceptanceMetric.TECHNICAL_FAILURE_RATE,
            ),
            (
                CnvRunOutcomeState.NOT_ASSESSABLE,
                CnvContributionStatus.NOT_ASSESSABLE,
                CnvAcceptanceMetric.NOT_ASSESSABLE_RATE,
            ),
        )
        for outcome, contribution, expected_metric in cases:
            with self.subTest(outcome=outcome):
                evidence = _manifest(
                    registration,
                    full_evidence=[
                        _run_summary(
                            registration,
                            run_outcome=outcome,
                            contribution_status=contribution,
                            outcome_reason="Synthetic terminal lane state.",
                        )
                    ],
                    normalized_events=[],
                )
                metrics = aggregate_cnv_run_state_metrics(registration, evidence)
                result = _find_metric(metrics, expected_metric, scope="overall")
                self.assertEqual(result.value, 1.0)
                self.assertEqual(result.evaluable_denominator, 1)

    def test_run_state_denominator_is_executed_lanes_not_matrix_cross_product(self) -> None:
        registration = _registration()
        metrics = aggregate_cnv_run_state_metrics(registration, _manifest(registration))
        result = _find_metric(
            metrics,
            CnvAcceptanceMetric.NO_CALL_RATE,
            scope="overall",
        )
        self.assertEqual(result.value, 0.0)
        self.assertEqual(result.evaluable_denominator, 1)
        self.assertEqual(result.component_counts["executed_lanes"], 1)

        missing_resolution = _find_metric(
            metrics,
            CnvAcceptanceMetric.NO_CALL_RATE,
            caller="qdnaseq_ace",
            genome_build="GRCh38",
            data_basis="lcwgs_genome_wide",
            coverage="[5,10)",
            tumor_fraction="[0.2,0.5)",
            bin_size=100,
        )
        self.assertEqual(
            missing_resolution.state,
            CnvMetricEvaluationState.NOT_EVALUABLE,
        )
        self.assertIsNone(missing_resolution.value)
        self.assertEqual(missing_resolution.evaluable_denominator, 0)

    def test_specificity_requires_explicit_negative_unit_assessment(self) -> None:
        registration = _registration()
        evidence = _manifest(registration)
        metrics = aggregate_cnv_run_state_metrics(registration, evidence)
        result = _find_metric(
            metrics,
            CnvAcceptanceMetric.SPECIFICITY,
            scope="overall",
        )
        self.assertEqual(result.state, CnvMetricEvaluationState.NOT_EVALUABLE)
        self.assertIsNone(result.value)

    def test_negative_unit_assessment_is_bound_to_registered_universe(self) -> None:
        registration = _registration()
        evidence = _manifest(registration)
        assignment = assign_cnv_validation_lanes(registration, evidence)[0]
        specimen = registration.cohort.specimens[0]
        universe = specimen.negative_universe
        assert universe is not None

        assessment = CnvNegativeUnitAssessment(
            assessment_id="negative-assessment-1",
            registration_sha256=registration.lock_sha256,
            evidence_manifest_sha256=evidence.manifest_sha256,
            lane_id=assignment.lane_id,
            specimen_id=specimen.specimen_id,
            universe_id=universe.universe_id,
            universe_resource_sha256=universe.resource_sha256,
            assessability_mask_sha256=specimen.assessability_mask.resource_sha256,
            assessed_units=100,
            false_positive_units=4,
            full_evidence_ids=[assignment.run_summary_full_evidence_id],
        )
        metrics = aggregate_cnv_run_state_metrics(
            registration,
            evidence,
            negative_assessments=[assessment],
        )
        specificity = _find_metric(
            metrics,
            CnvAcceptanceMetric.SPECIFICITY,
            scope="overall",
        )
        self.assertEqual(specificity.value, 0.96)
        self.assertEqual(specificity.evaluable_denominator, 100)
        self.assertEqual(specificity.component_counts["true_negative_units"], 96)
        self.assertEqual(specificity.component_counts["false_positive_units"], 4)

        wrong_universe = assessment.model_copy(
            update={"universe_resource_sha256": _sha("wrong-universe")}
        )
        with self.assertRaises(ValueError):
            aggregate_cnv_run_state_metrics(
                registration,
                evidence,
                negative_assessments=[wrong_universe],
            )

    def test_negative_unit_assessment_rejects_invalid_unit_counts(self) -> None:
        registration = _registration()
        evidence = _manifest(registration)
        assignment = assign_cnv_validation_lanes(registration, evidence)[0]
        specimen = registration.cohort.specimens[0]
        universe = specimen.negative_universe
        assert universe is not None

        with self.assertRaises(ValueError):
            CnvNegativeUnitAssessment(
                assessment_id="negative-assessment-invalid",
                registration_sha256=registration.lock_sha256,
                evidence_manifest_sha256=evidence.manifest_sha256,
                lane_id=assignment.lane_id,
                specimen_id=specimen.specimen_id,
                universe_id=universe.universe_id,
                universe_resource_sha256=universe.resource_sha256,
                assessability_mask_sha256=specimen.assessability_mask.resource_sha256,
                assessed_units=100,
                false_positive_units=101,
                full_evidence_ids=[assignment.run_summary_full_evidence_id],
            )

        incomplete = CnvNegativeUnitAssessment(
            assessment_id="negative-assessment-incomplete",
            registration_sha256=registration.lock_sha256,
            evidence_manifest_sha256=evidence.manifest_sha256,
            lane_id=assignment.lane_id,
            specimen_id=specimen.specimen_id,
            universe_id=universe.universe_id,
            universe_resource_sha256=universe.resource_sha256,
            assessability_mask_sha256=specimen.assessability_mask.resource_sha256,
            assessed_units=99,
            false_positive_units=0,
            full_evidence_ids=[assignment.run_summary_full_evidence_id],
        )
        with self.assertRaises(ValueError):
            aggregate_cnv_run_state_metrics(
                registration,
                evidence,
                negative_assessments=[incomplete],
            )

    def test_specificity_is_not_event_class_specific_without_registered_unit_classes(self) -> None:
        registration = _registration()
        evidence = _manifest(registration)
        assignment = assign_cnv_validation_lanes(registration, evidence)[0]
        specimen = registration.cohort.specimens[0]
        universe = specimen.negative_universe
        assert universe is not None
        assessment = CnvNegativeUnitAssessment(
            assessment_id="negative-assessment-1",
            registration_sha256=registration.lock_sha256,
            evidence_manifest_sha256=evidence.manifest_sha256,
            lane_id=assignment.lane_id,
            specimen_id=specimen.specimen_id,
            universe_id=universe.universe_id,
            universe_resource_sha256=universe.resource_sha256,
            assessability_mask_sha256=specimen.assessability_mask.resource_sha256,
            assessed_units=100,
            false_positive_units=0,
            full_evidence_ids=[assignment.run_summary_full_evidence_id],
        )
        metrics = aggregate_cnv_run_state_metrics(
            registration,
            evidence,
            negative_assessments=[assessment],
        )
        event_specific = [
            item
            for item in metrics
            if item.metric == CnvAcceptanceMetric.SPECIFICITY and "event_class" in item.stratum_key
        ]
        self.assertEqual(event_specific, [])


def _registration_with_quantitative_truth() -> CnvValidationRegistration:
    specimen = _specimen().model_copy(
        update={
            "quantitative_truth": CnvQuantitativeTruth(
                truth_source_resource_id="synthetic-truth-v1",
                cellularity=0.40,
                ploidy=2.20,
                note="Synthetic orthogonal quantitative truth.",
            )
        }
    )
    return preregister_cnv_validation(
        _matrix(),
        CnvValidationCohort(specimens=[specimen]),
        registration_id="synthetic-cnv-quantitative-registration",
        registered_at=datetime(2026, 9, 18, tzinfo=UTC),
        code_sha256=_sha("code"),
        software_version="0.8.2",
    )


def _selected_fit(
    registration: CnvValidationRegistration,
    *,
    record_id: str = "fit-500-selected",
    specimen_id: str | None = None,
    replicate_id: str = "replicate-1",
    selected_fit: bool = True,
    contribution_status: CnvContributionStatus = CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
    cellularity: float = 0.35,
    ploidy: float = 2.30,
) -> CnvFullEvidenceRecord:
    return _full_record(
        registration,
        record_id=record_id,
        native_record_id=record_id,
        record_kind=CnvEvidenceRecordKind.CALLER_FIT,
        contribution_status=contribution_status,
        specimen_id=specimen_id,
        replicate_id=replicate_id,
        fit_group_id="ace-fit-grid-500",
        selected_fit=selected_fit,
        cellularity=cellularity,
        ploidy=ploidy,
        fit_error=0.1,
    )


def _repeat_registration(
    *, second_biological_id: str = "SYNTHETIC_BIO_REPEAT"
) -> CnvValidationRegistration:
    first = _specimen().model_copy(
        update={
            "specimen_id": "SYNTHETIC_REPEAT_1",
            "biological_specimen_id": "SYNTHETIC_BIO_REPEAT",
            "repeat_kind": CnvRepeatKind.BETWEEN_RUN,
            "repeat_group_id": "repeat-group-1",
            "input_sha256": _sha("repeat-input-1"),
        }
    )
    second = _specimen().model_copy(
        update={
            "specimen_id": "SYNTHETIC_REPEAT_2",
            "biological_specimen_id": second_biological_id,
            "repeat_kind": CnvRepeatKind.BETWEEN_RUN,
            "repeat_group_id": "repeat-group-1",
            "input_sha256": _sha("repeat-input-2"),
        }
    )
    return preregister_cnv_validation(
        _matrix(),
        CnvValidationCohort(specimens=[first, second]),
        registration_id="synthetic-cnv-repeat-registration",
        registered_at=datetime(2026, 9, 18, tzinfo=UTC),
        code_sha256=_sha("code"),
        software_version="0.8.2",
    )


def _repeat_manifest(
    registration: CnvValidationRegistration,
    *,
    with_events: bool = True,
) -> CnvValidationEvidenceManifest:
    records: list[CnvFullEvidenceRecord] = []
    events: list[CnvNormalizedEventRecord] = []
    for index, specimen in enumerate(registration.cohort.specimens, start=1):
        replicate_id = f"repeat-{index}"
        records.append(
            _run_summary(
                registration,
                record_id=f"run-summary-repeat-{index}",
                native_record_id=f"run-summary-repeat-{index}",
                specimen_id=specimen.specimen_id,
                replicate_id=replicate_id,
            )
        )
        if with_events:
            source_id = f"segment-repeat-{index}"
            records.append(
                _event_source(
                    registration,
                    record_id=source_id,
                    native_record_id=source_id,
                    specimen_id=specimen.specimen_id,
                    replicate_id=replicate_id,
                    native_artifact_ids=["segments-500"],
                )
            )
            events.append(
                _normalized_event(
                    registration,
                    specimen_id=specimen.specimen_id,
                    replicate_id=replicate_id,
                    normalized_event_id=f"normalized-repeat-{index}",
                    source_full_evidence_ids=[source_id],
                )
            )
    return _manifest(
        registration,
        full_evidence=records,
        normalized_events=events,
    )


class CnvQuantitativeAndReproducibilityTests(unittest.TestCase):
    def test_copy_number_error_reports_mae_and_signed_bias(self) -> None:
        registration = _registration()
        evidence = _manifest(registration)
        metrics = aggregate_cnv_quantitative_metrics(registration, evidence)
        result = _find_metric(
            metrics,
            CnvAcceptanceMetric.COPY_NUMBER_ERROR,
            scope="overall",
        )
        self.assertAlmostEqual(result.value or 0.0, 0.1)
        self.assertAlmostEqual(result.summary_values["signed_bias"], 0.1)
        self.assertEqual(result.evaluable_denominator, 1)
        self.assertEqual(result.component_counts["eligible_pairs"], 1)

    def test_cellularity_and_ploidy_require_explicit_quantitative_truth(self) -> None:
        registration = _registration()
        evidence = _manifest(
            registration,
            full_evidence=[
                _run_summary(registration),
                _event_source(registration),
                _selected_fit(registration),
            ],
        )
        metrics = aggregate_cnv_quantitative_metrics(registration, evidence)
        self.assertEqual(
            _find_metric(
                metrics,
                CnvAcceptanceMetric.CELLULARITY_ERROR,
                scope="overall",
            ).state,
            CnvMetricEvaluationState.NOT_EVALUABLE,
        )
        self.assertEqual(
            _find_metric(
                metrics,
                CnvAcceptanceMetric.PLOIDY_ERROR,
                scope="overall",
            ).state,
            CnvMetricEvaluationState.NOT_EVALUABLE,
        )

    def test_selected_fit_quantitative_errors_do_not_promote_alternative_fits(self) -> None:
        registration = _registration_with_quantitative_truth()
        selected = _selected_fit(registration)
        alternative = _selected_fit(
            registration,
            record_id="fit-500-alternative",
            selected_fit=False,
            contribution_status=CnvContributionStatus.SECONDARY,
            cellularity=0.55,
            ploidy=3.10,
        )
        evidence = _manifest(
            registration,
            full_evidence=[
                _run_summary(registration),
                _event_source(registration),
                selected,
                alternative,
            ],
        )
        metrics = aggregate_cnv_quantitative_metrics(registration, evidence)
        cellularity = _find_metric(
            metrics,
            CnvAcceptanceMetric.CELLULARITY_ERROR,
            scope="overall",
        )
        ploidy = _find_metric(
            metrics,
            CnvAcceptanceMetric.PLOIDY_ERROR,
            scope="overall",
        )
        self.assertAlmostEqual(cellularity.value or 0.0, 0.05)
        self.assertAlmostEqual(ploidy.value or 0.0, 0.10)
        self.assertIn(selected.record_id, cellularity.full_evidence_ids)
        self.assertNotIn(alternative.record_id, cellularity.full_evidence_ids)
        self.assertIn(alternative.record_id, [item.record_id for item in evidence.full_evidence])

    def test_multiple_selected_fits_in_one_lane_fail_closed(self) -> None:
        registration = _registration_with_quantitative_truth()
        evidence = _manifest(
            registration,
            full_evidence=[
                _run_summary(registration),
                _event_source(registration),
                _selected_fit(registration),
                _selected_fit(registration, record_id="fit-500-selected-2"),
            ],
        )
        with self.assertRaises(ValueError):
            aggregate_cnv_quantitative_metrics(registration, evidence)

    def test_reproducibility_pairs_same_biological_specimen_and_matches_events(self) -> None:
        registration = _repeat_registration()
        evidence = _repeat_manifest(registration)
        metrics, pairs = aggregate_cnv_reproducibility_metrics(registration, evidence)
        result = _find_metric(
            metrics,
            CnvAcceptanceMetric.REPRODUCIBILITY,
            scope="overall",
        )
        self.assertEqual(result.value, 1.0)
        self.assertEqual(result.evaluable_denominator, 1)
        self.assertEqual(len(pairs), 1)
        self.assertIsInstance(pairs[0], CnvReproducibilityPair)
        self.assertEqual(pairs[0].concordance, 1.0)
        self.assertEqual(
            pairs[0].matched_normalized_event_pairs,
            [["normalized-repeat-1", "normalized-repeat-2"]],
        )

    def test_repeat_group_cannot_mix_biological_specimens(self) -> None:
        registration = _repeat_registration(second_biological_id="DIFFERENT_BIO")
        evidence = _repeat_manifest(registration)
        with self.assertRaises(ValueError):
            aggregate_cnv_reproducibility_metrics(registration, evidence)

    def test_two_empty_repeat_call_sets_are_not_artificial_perfect_concordance(self) -> None:
        registration = _repeat_registration()
        evidence = _repeat_manifest(registration, with_events=False)
        metrics, pairs = aggregate_cnv_reproducibility_metrics(registration, evidence)
        result = _find_metric(
            metrics,
            CnvAcceptanceMetric.REPRODUCIBILITY,
            scope="overall",
        )
        self.assertEqual(result.state, CnvMetricEvaluationState.NOT_EVALUABLE)
        self.assertIsNone(result.value)
        self.assertEqual(len(pairs), 1)
        self.assertIsNone(pairs[0].concordance)


def _registration_with_acceptance(
    question: CnvAcceptanceQuestion,
) -> CnvValidationRegistration:
    matrix = _matrix().model_copy(update={"acceptance": [question]})
    return preregister_cnv_validation(
        matrix,
        CnvValidationCohort(specimens=[_specimen()]),
        registration_id="synthetic-cnv-acceptance-registration",
        registered_at=datetime(2026, 9, 18, tzinfo=UTC),
        code_sha256=_sha("code"),
        software_version="0.8.2",
    )


def _event_sensitivity_stratum() -> dict[str, str | int | float | bool]:
    return {
        "caller": "qdnaseq_ace",
        "genome_build": "GRCh38",
        "data_basis": "lcwgs_genome_wide",
        "coverage": "[5,10)",
        "tumor_fraction": "[0.2,0.5)",
        "event_class": "deletion",
        "event_size": "[5000000,10000000)",
        "bin_size": 500,
    }


class CnvAcceptanceAndReportTests(unittest.TestCase):
    def test_acceptance_filters_match_metric_stratum_exactly(self) -> None:
        question = CnvAcceptanceQuestion(
            question_id="overall-sensitivity-pass",
            metric=CnvAcceptanceMetric.SENSITIVITY,
            minimum_evaluable_denominator=1,
            minimum_acceptable=0.9,
            filters={"scope": "overall"},
        )
        registration = _registration_with_acceptance(question)
        report = aggregate_cnv_validation(
            registration,
            _manifest(registration),
            report_id="synthetic-cnv-report-pass",
        )
        self.assertIsInstance(report, CnvValidationMetricReport)
        self.assertEqual(len(report.acceptance_results), 1)
        result = report.acceptance_results[0]
        self.assertIsInstance(result, CnvAcceptanceResult)
        self.assertEqual(result.decision, CnvAcceptanceDecision.PASS)
        self.assertEqual(result.metric_value, 1.0)
        self.assertIsNotNone(result.metric_id)

    def test_acceptance_below_registered_denominator_is_not_evaluable(self) -> None:
        question = CnvAcceptanceQuestion(
            question_id="overall-sensitivity-underpowered",
            metric=CnvAcceptanceMetric.SENSITIVITY,
            minimum_evaluable_denominator=2,
            minimum_acceptable=0.9,
            filters={"scope": "overall"},
        )
        registration = _registration_with_acceptance(question)
        report = aggregate_cnv_validation(
            registration,
            _manifest(registration),
            report_id="synthetic-cnv-report-underpowered",
        )
        result = report.acceptance_results[0]
        self.assertEqual(result.decision, CnvAcceptanceDecision.NOT_EVALUABLE)
        self.assertEqual(result.evaluable_denominator, 1)
        self.assertIsNone(result.accepted)

    def test_acceptance_bounds_produce_fail_without_mutating_metric(self) -> None:
        question = CnvAcceptanceQuestion(
            question_id="overall-sensitivity-fail",
            metric=CnvAcceptanceMetric.SENSITIVITY,
            minimum_evaluable_denominator=1,
            maximum_acceptable=0.9,
            filters={"scope": "overall"},
        )
        registration = _registration_with_acceptance(question)
        report = aggregate_cnv_validation(
            registration,
            _manifest(registration),
            report_id="synthetic-cnv-report-fail",
        )
        result = report.acceptance_results[0]
        self.assertEqual(result.decision, CnvAcceptanceDecision.FAIL)
        self.assertFalse(result.accepted)
        self.assertEqual(result.metric_value, 1.0)

    def test_missing_exact_stratum_is_explicit_not_evaluable(self) -> None:
        question = CnvAcceptanceQuestion(
            question_id="missing-stratum",
            metric=CnvAcceptanceMetric.SENSITIVITY,
            minimum_evaluable_denominator=1,
            minimum_acceptable=0.9,
            filters={"caller": "not-a-registered-caller"},
        )
        registration = _registration_with_acceptance(question)
        report = aggregate_cnv_validation(
            registration,
            _manifest(registration),
            report_id="synthetic-cnv-report-missing-stratum",
        )
        result = report.acceptance_results[0]
        self.assertEqual(result.decision, CnvAcceptanceDecision.NOT_EVALUABLE)
        self.assertIsNone(result.metric_id)
        self.assertIsNone(result.metric_value)

    def test_all_no_call_lane_maps_event_acceptance_to_no_call(self) -> None:
        question = CnvAcceptanceQuestion(
            question_id="event-sensitivity-no-call",
            metric=CnvAcceptanceMetric.SENSITIVITY,
            minimum_evaluable_denominator=1,
            minimum_acceptable=0.9,
            filters=_event_sensitivity_stratum(),
        )
        registration = _registration_with_acceptance(question)
        evidence = _manifest(
            registration,
            full_evidence=[
                _run_summary(
                    registration,
                    run_outcome=CnvRunOutcomeState.NO_CALL,
                    contribution_status=CnvContributionStatus.NO_CALL,
                    outcome_reason="Synthetic no-call lane.",
                )
            ],
            normalized_events=[],
        )
        report = aggregate_cnv_validation(
            registration,
            evidence,
            report_id="synthetic-cnv-report-no-call",
        )
        result = report.acceptance_results[0]
        self.assertEqual(result.decision, CnvAcceptanceDecision.NO_CALL)
        self.assertIsNone(result.accepted)

    def test_report_binds_registration_evidence_thresholds_and_metric_provenance(self) -> None:
        question = CnvAcceptanceQuestion(
            question_id="overall-sensitivity-provenance",
            metric=CnvAcceptanceMetric.SENSITIVITY,
            minimum_evaluable_denominator=1,
            minimum_acceptable=0.9,
            filters={"scope": "overall"},
        )
        registration = _registration_with_acceptance(question)
        evidence = _manifest(registration)
        report = aggregate_cnv_validation(
            registration,
            evidence,
            report_id="synthetic-cnv-report-provenance",
        )
        self.assertEqual(report.registration_sha256, registration.lock_sha256)
        self.assertEqual(report.evidence_manifest_sha256, evidence.manifest_sha256)
        self.assertEqual(
            report.matching_thresholds,
            registration.matrix.matching_thresholds,
        )
        self.assertTrue(report.assignments)
        self.assertTrue(report.metrics)
        self.assertTrue(
            all(
                item.registration_sha256 == registration.lock_sha256
                and item.evidence_manifest_sha256 == evidence.manifest_sha256
                for item in report.metrics
            )
        )
        metric_ids = [item.metric_id for item in report.metrics]
        self.assertEqual(len(metric_ids), len(set(metric_ids)))

    def test_report_lock_is_deterministic_and_rejects_metric_assignment_acceptance_tamper(
        self,
    ) -> None:
        question = CnvAcceptanceQuestion(
            question_id="overall-sensitivity-lock",
            metric=CnvAcceptanceMetric.SENSITIVITY,
            minimum_evaluable_denominator=1,
            minimum_acceptable=0.9,
            filters={"scope": "overall"},
        )
        registration = _registration_with_acceptance(question)
        evidence = _manifest(registration)
        first = aggregate_cnv_validation(
            registration,
            evidence,
            report_id="synthetic-cnv-report-lock",
        )
        second = aggregate_cnv_validation(
            registration,
            evidence,
            report_id="synthetic-cnv-report-lock",
        )
        self.assertEqual(first.report_sha256, second.report_sha256)

        metric_payload = first.model_dump(mode="json")
        metric_payload["metrics"][0]["reason"] = "tampered metric"
        with self.assertRaises(ValidationError):
            CnvValidationMetricReport.model_validate(metric_payload)

        assignment_payload = first.model_dump(mode="json")
        assignment_payload["assignments"][0]["minimum_reciprocal_overlap"] = 0.123
        with self.assertRaises(ValidationError):
            CnvValidationMetricReport.model_validate(assignment_payload)

        acceptance_payload = first.model_dump(mode="json")
        acceptance_payload["acceptance_results"][0]["decision"] = "FAIL"
        with self.assertRaises(ValidationError):
            CnvValidationMetricReport.model_validate(acceptance_payload)

    def test_report_contains_no_timestamp_or_production_winner_surface(self) -> None:
        question = CnvAcceptanceQuestion(
            question_id="overall-sensitivity-no-winner",
            metric=CnvAcceptanceMetric.SENSITIVITY,
            minimum_evaluable_denominator=1,
            minimum_acceptable=0.9,
            filters={"scope": "overall"},
        )
        registration = _registration_with_acceptance(question)
        report = aggregate_cnv_validation(
            registration,
            _manifest(registration),
            report_id="synthetic-cnv-report-no-winner",
        )
        serialized = json.dumps(report.model_dump(mode="json"), sort_keys=True)
        self.assertNotIn("created_at", serialized)
        self.assertNotIn("registered_at", serialized)
        self.assertNotIn("winner", serialized.lower())
        self.assertNotIn("selected_caller", serialized.lower())


if __name__ == "__main__":
    unittest.main()

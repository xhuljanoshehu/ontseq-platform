from __future__ import annotations

import hashlib
import unittest
from datetime import UTC, datetime

from ontseq_platform.cnv_validation_contracts import (
    CnvAcceptanceMetric,
    CnvAcceptanceQuestion,
    CnvAssessabilityMask,
    CnvCallerLock,
    CnvDataBasis,
    CnvNegativeUniverse,
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
    assign_cnv_validation_lanes,
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
) -> CnvFullEvidenceRecord:
    specimen = registration.cohort.specimens[0]
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
        replicate_id="replicate-1",
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
    **updates: object,
) -> CnvNormalizedEventRecord:
    specimen = registration.cohort.specimens[0]
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
        "replicate_id": "replicate-1",
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
            normalized_events=[
                _normalized_event(registration, caller_version="synthetic-2")
            ],
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


if __name__ == "__main__":
    unittest.main()

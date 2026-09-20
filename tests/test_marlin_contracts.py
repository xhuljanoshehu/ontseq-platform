from datetime import UTC, datetime

import pytest
from pydantic import ValidationError

from ontseq_platform.marlin_contracts import (
    MarlinArtifactLock,
    MarlinBridgeLock,
    MarlinClassificationDecision,
    MarlinFeatureSummary,
    MarlinGroupedScore,
    MarlinModelUnitScore,
    MarlinPredictionReport,
    MarlinProbeObservation,
    MarlinRuntimeCompatibilityProfile,
    MarlinSourceKind,
)
from ontseq_platform.models import FileFingerprint, GenomeBuild, ModuleRunStatus

SHA0 = "0" * 64
SHA1 = "1" * 64
SHA2 = "2" * 64
SHA3 = "3" * 64
SHA4 = "4" * 64
SHA5 = "5" * 64
SHA6 = "6" * 64


def _feature_summary(*, observed: int = 10, explicit_na: int = 2) -> MarlinFeatureSummary:
    return MarlinFeatureSummary(
        observed_model_feature_count=observed,
        explicit_na_feature_count=explicit_na,
        absent_feature_count=357340 - observed - explicit_na,
        non_model_probe_count=1,
        observed_fraction=observed / 357340,
        feature_vector_sha256=SHA0,
        feature_artifact_sha256=SHA1,
    )


def _artifact_lock() -> MarlinArtifactLock:
    return MarlinArtifactLock(
        lock_id="MARLIN_V1_GRCH37_TEST",
        model_version="1.0.0",
        model_source_uri="https://doi.org/10.5281/zenodo.15565404",
        model_sha256=SHA0,
        code_source_uri="https://github.com/hovestadt/MARLIN",
        code_version_or_commit="442aa603415a54f62e7367794f9a31c6bc20fc2d",
        code_manifest_sha256=SHA1,
        feature_source_uri="marlin_v1.features.RData",
        feature_sha256=SHA2,
        canonical_feature_list_sha256=SHA3,
        class_annotation_source_uri="marlin_v1.class_annotations.xlsx",
        class_annotation_sha256=SHA4,
        probe_resource_source_uri="marlin_v1.probes_hg19.bed.gz",
        probe_resource_sha256=SHA5,
        genome_build=GenomeBuild.GRCH37,
        R_version="4.1.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version_if_used="3.11.0",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 13, tzinfo=UTC),
    )


def _bridge_lock() -> MarlinBridgeLock:
    return MarlinBridgeLock(
        bridge_id="MARLIN_MODKIT_BRIDGE_TEST",
        status="validated_same_specimen_bridge",
        genome_build=GenomeBuild.GRCH37,
        adapter_version="marlin-modkit-bridge-v1",
        modkit_version="0.6.4",
        probe_resource_sha256=SHA0,
        feature_artifact_sha256=SHA1,
        validation_precomputed_input_sha256=SHA2,
        validation_modkit_bedmethyl_sha256=SHA3,
        validation_precomputed_feature_vector_sha256=SHA4,
        validation_modkit_feature_vector_sha256=SHA5,
        compared_feature_count=1000,
        concordant_feature_count=990,
        feature_agreement_fraction=0.99,
        evidence_reference="same-specimen-bridge-study-001",
        validated_at=datetime(2026, 9, 13, tzinfo=UTC),
    )


def test_feature_summary_requires_marlin_v1_feature_count() -> None:
    with pytest.raises(ValidationError):
        MarlinFeatureSummary(
            expected_feature_count=42,
            observed_model_feature_count=10,
            explicit_na_feature_count=0,
            absent_feature_count=32,
            non_model_probe_count=0,
            observed_fraction=10 / 42,
            feature_vector_sha256=SHA0,
            feature_artifact_sha256=SHA1,
        )


def test_feature_summary_counts_partition_model_features() -> None:
    with pytest.raises(ValidationError, match="partition"):
        MarlinFeatureSummary(
            observed_model_feature_count=10,
            explicit_na_feature_count=2,
            absent_feature_count=3,
            non_model_probe_count=0,
            observed_fraction=10 / 357340,
            feature_vector_sha256=SHA0,
            feature_artifact_sha256=SHA1,
        )


def test_feature_summary_observed_fraction_must_match_count() -> None:
    with pytest.raises(ValidationError, match="observed_fraction"):
        MarlinFeatureSummary(
            observed_model_feature_count=10,
            explicit_na_feature_count=0,
            absent_feature_count=357330,
            non_model_probe_count=0,
            observed_fraction=0.5,
            feature_vector_sha256=SHA0,
            feature_artifact_sha256=SHA1,
        )


def test_decision_values_are_unambiguous() -> None:
    assert MarlinClassificationDecision.HIGH_CONFIDENCE.value == "HIGH_CONFIDENCE"
    assert MarlinClassificationDecision.UNKNOWN.value == "UNKNOWN"
    assert MarlinSourceKind.PRECOMPUTED_METHYLATION.value == "PRECOMPUTED_METHYLATION"
    assert MarlinSourceKind.MODKIT_DERIVED.value == "MODKIT_DERIVED"


def test_probe_observation_rejects_noncanonical_interval() -> None:
    with pytest.raises(ValidationError):
        MarlinProbeObservation(
            chromosome="chr1", start=10, end=10, methylation_fraction=0.4, probe_id="cg1"
        )


def test_probe_observation_rejects_nonfinite_fraction() -> None:
    with pytest.raises(ValidationError):
        MarlinProbeObservation(
            chromosome="chr1",
            start=10,
            end=11,
            methylation_fraction=float("nan"),
            probe_id="cg1",
        )


def test_artifact_lock_is_grch37_and_exact_marlin_shape() -> None:
    lock = _artifact_lock()
    assert lock.expected_feature_count == 357340
    assert lock.expected_model_unit_count == 42
    with pytest.raises(ValidationError):
        MarlinArtifactLock(**{**lock.model_dump(), "genome_build": GenomeBuild.GRCH38})


def test_bridge_lock_records_same_specimen_evidence_and_modkit_identity() -> None:
    lock = _bridge_lock()
    assert lock.status == "validated_same_specimen_bridge"
    assert lock.modkit_version == "0.6.4"
    assert lock.feature_agreement_fraction == pytest.approx(0.99)
    assert lock.research_only is True


def test_bridge_lock_rejects_non_grch37_or_inconsistent_agreement() -> None:
    lock = _bridge_lock()
    with pytest.raises(ValidationError, match="GRCh37"):
        MarlinBridgeLock.model_validate({**lock.model_dump(), "genome_build": GenomeBuild.GRCH38})
    with pytest.raises(ValidationError, match="agreement"):
        MarlinBridgeLock.model_validate({**lock.model_dump(), "feature_agreement_fraction": 0.5})
    with pytest.raises(ValidationError, match="concordant"):
        MarlinBridgeLock.model_validate(
            {
                **lock.model_dump(),
                "compared_feature_count": 10,
                "concordant_feature_count": 11,
                "feature_agreement_fraction": 1.0,
            }
        )


def test_runtime_profile_requires_42_scores_and_positive_tolerance() -> None:
    with pytest.raises(ValidationError):
        MarlinRuntimeCompatibilityProfile(
            profile_id="p",
            reference_runtime_lock_id="lock",
            feature_vector_sha256=SHA0,
            reference_scores=[1 / 41] * 41,
            absolute_score_tolerance=1e-6,
            score_sum_tolerance=1e-5,
            top_model_unit_index=0,
            execution_backend="cpu",
            created_at=datetime(2026, 9, 13, tzinfo=UTC),
        )


def test_prediction_report_requires_runtime_profile_id() -> None:
    with pytest.raises(ValidationError, match="runtime_profile_id"):
        MarlinPredictionReport(
            sample_id="AL_001",
            status=ModuleRunStatus.NO_CALL,
            decision=MarlinClassificationDecision.UNKNOWN,
            source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
            genome_build=GenomeBuild.GRCH37,
            input_fingerprint=FileFingerprint(size_bytes=100, sha256=SHA6),
            feature_summary=_feature_summary(observed=0, explicit_na=0),
            raw_model_scores=[],
            class_scores=[],
            family_scores=[],
            lineage_scores=[],
            artifact_lock_id="MARLIN_V1_GRCH37_TEST",
        )


def _raw_uniform_scores() -> list[MarlinModelUnitScore]:
    return [
        MarlinModelUnitScore(model_id=i + 1, label=f"unit-{i + 1}", score=1 / 42) for i in range(42)
    ]


def test_prediction_report_completed_unknown_for_valid_low_confidence() -> None:
    grouped = [MarlinGroupedScore(label="Class A", score=0.69)]
    report = MarlinPredictionReport(
        sample_id="AL_001",
        status=ModuleRunStatus.COMPLETED,
        decision=MarlinClassificationDecision.UNKNOWN,
        source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
        genome_build=GenomeBuild.GRCH37,
        input_fingerprint=FileFingerprint(size_bytes=100, sha256=SHA6),
        feature_summary=_feature_summary(),
        raw_model_scores=_raw_uniform_scores(),
        class_scores=grouped,
        family_scores=[MarlinGroupedScore(label="Family A", score=0.8)],
        lineage_scores=[MarlinGroupedScore(label="Lineage A", score=0.9)],
        top_class="Class A",
        top_class_score=0.69,
        confidence_threshold=0.8,
        artifact_lock_id="MARLIN_V1_GRCH37_TEST",
        runtime_profile_id="PROFILE",
    )
    assert report.status is ModuleRunStatus.COMPLETED
    assert report.decision is MarlinClassificationDecision.UNKNOWN


def test_prediction_report_no_call_cannot_carry_high_confidence() -> None:
    with pytest.raises(ValidationError, match="NO_CALL"):
        MarlinPredictionReport(
            sample_id="AL_001",
            status=ModuleRunStatus.NO_CALL,
            decision=MarlinClassificationDecision.HIGH_CONFIDENCE,
            source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
            genome_build=GenomeBuild.GRCH37,
            input_fingerprint=FileFingerprint(size_bytes=100, sha256=SHA6),
            feature_summary=_feature_summary(observed=0, explicit_na=0),
            raw_model_scores=[],
            class_scores=[],
            family_scores=[],
            lineage_scores=[],
            top_class=None,
            top_class_score=None,
            confidence_threshold=0.8,
            artifact_lock_id="MARLIN_V1_GRCH37_TEST",
            runtime_profile_id="PROFILE",
        )


def test_prediction_report_no_call_requires_zero_evidence_payload() -> None:
    with pytest.raises(ValidationError, match="NO_CALL"):
        MarlinPredictionReport(
            sample_id="AL_001",
            status=ModuleRunStatus.NO_CALL,
            decision=MarlinClassificationDecision.UNKNOWN,
            source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
            genome_build=GenomeBuild.GRCH37,
            input_fingerprint=FileFingerprint(size_bytes=100, sha256=SHA6),
            feature_summary=_feature_summary(observed=0, explicit_na=0),
            raw_model_scores=_raw_uniform_scores(),
            class_scores=[],
            family_scores=[],
            lineage_scores=[],
            top_class=None,
            top_class_score=None,
            confidence_threshold=0.8,
            artifact_lock_id="MARLIN_V1_GRCH37_TEST",
            runtime_profile_id="PROFILE",
        )


def test_prediction_report_requires_raw_model_unit_order() -> None:
    raw = _raw_uniform_scores()
    raw[0], raw[1] = raw[1], raw[0]
    with pytest.raises(ValidationError, match="order"):
        MarlinPredictionReport(
            sample_id="AL_001",
            status=ModuleRunStatus.COMPLETED,
            decision=MarlinClassificationDecision.UNKNOWN,
            source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
            genome_build=GenomeBuild.GRCH37,
            input_fingerprint=FileFingerprint(size_bytes=100, sha256=SHA6),
            feature_summary=_feature_summary(),
            raw_model_scores=raw,
            class_scores=[MarlinGroupedScore(label="Class A", score=0.69)],
            family_scores=[MarlinGroupedScore(label="Family A", score=0.8)],
            lineage_scores=[MarlinGroupedScore(label="Lineage A", score=0.9)],
            top_class="Class A",
            top_class_score=0.69,
            confidence_threshold=0.8,
            artifact_lock_id="MARLIN_V1_GRCH37_TEST",
            runtime_profile_id="PROFILE",
        )


def test_prediction_report_ties_require_canonical_label_tiebreak() -> None:
    with pytest.raises(ValidationError, match="top classification"):
        MarlinPredictionReport(
            sample_id="AL_001",
            status=ModuleRunStatus.COMPLETED,
            decision=MarlinClassificationDecision.UNKNOWN,
            source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
            genome_build=GenomeBuild.GRCH37,
            input_fingerprint=FileFingerprint(size_bytes=100, sha256=SHA6),
            feature_summary=_feature_summary(),
            raw_model_scores=_raw_uniform_scores(),
            class_scores=[
                MarlinGroupedScore(label="Class Z", score=0.5),
                MarlinGroupedScore(label="Class A", score=0.5),
            ],
            family_scores=[MarlinGroupedScore(label="Family A", score=1.0)],
            lineage_scores=[MarlinGroupedScore(label="Lineage A", score=1.0)],
            top_class="Class Z",
            top_class_score=0.5,
            confidence_threshold=0.8,
            artifact_lock_id="MARLIN_V1_GRCH37_TEST",
            runtime_profile_id="PROFILE",
        )


def _valid_dual_runtime_report_payload() -> dict[str, object]:
    from ontseq_platform.marlin_contracts import MarlinArtifactSetIdentity, MarlinRuntimeIdentity

    artifact_set = MarlinArtifactSetIdentity(
        artifact_set_sha256=SHA6,
        model_version="1.0.0",
        model_sha256=SHA0,
        code_version_or_commit="442aa603415a54f62e7367794f9a31c6bc20fc2d",
        code_manifest_sha256=SHA1,
        feature_sha256=SHA2,
        canonical_feature_list_sha256=SHA3,
        class_annotation_sha256=SHA4,
        probe_resource_sha256=SHA5,
        genome_build=GenomeBuild.GRCH37,
    )
    reference_runtime = MarlinRuntimeIdentity(
        runtime_id="REF_RUNTIME",
        R_version="4.1.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version="3.10.21",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 15, tzinfo=UTC),
    )
    candidate_runtime = reference_runtime.model_copy(
        update={"runtime_id": "CANDIDATE_RUNTIME", "R_version": "4.2.3"}
    )
    scores = [0.5, 0.5] + [0.0] * 40
    return {
        "comparison_id": "DUAL_RUNTIME_TEST",
        "artifact_set_identity": artifact_set,
        "reference_artifact_lock_id": "REFERENCE_LOCK",
        "candidate_artifact_lock_id": "CANDIDATE_LOCK",
        "feature_vector_sha256": SHA6,
        "reference_profile_id": "REFERENCE_PROFILE",
        "reference_runtime_identity": reference_runtime,
        "candidate_runtime_identity": candidate_runtime,
        "absolute_score_tolerance": 1e-7,
        "score_sum_tolerance": 1e-5,
        "reference_scores": scores,
        "candidate_scores": scores,
        "absolute_differences": [0.0] * 42,
        "max_absolute_difference": 0.0,
        "reference_top_model_unit_index": 0,
        "candidate_top_model_unit_index": 0,
        "all_scores_within_tolerance": True,
        "top_model_unit_matches": True,
        "softmax_invariants_pass": True,
        "verdict": "PASS",
        "created_at": datetime(2026, 9, 15, tzinfo=UTC),
    }


def test_dual_runtime_report_rejects_inconsistent_score_gate_evidence() -> None:
    from ontseq_platform.marlin_contracts import MarlinDualRuntimeCompatibilityReport

    payload = _valid_dual_runtime_report_payload()
    payload["candidate_scores"] = [0.4999998, 0.5000002] + [0.0] * 40
    payload["absolute_differences"] = [2e-7, 2e-7] + [0.0] * 40
    payload["max_absolute_difference"] = 2e-7
    payload["candidate_top_model_unit_index"] = 1
    payload["all_scores_within_tolerance"] = False
    payload["top_model_unit_matches"] = False
    with pytest.raises(ValidationError, match="inconsistent"):
        MarlinDualRuntimeCompatibilityReport.model_validate(payload)


def test_dual_runtime_report_rejects_fail_verdict_when_all_gates_pass() -> None:
    from ontseq_platform.marlin_contracts import MarlinDualRuntimeCompatibilityReport

    payload = _valid_dual_runtime_report_payload()
    payload["verdict"] = "FAIL"
    with pytest.raises(ValidationError, match="PASS"):
        MarlinDualRuntimeCompatibilityReport.model_validate(payload)


def test_runtime_identity_requires_timezone_aware_timestamp() -> None:
    from ontseq_platform.marlin_contracts import MarlinRuntimeIdentity

    with pytest.raises(ValidationError, match="timezone"):
        MarlinRuntimeIdentity(
            runtime_id="REF_RUNTIME",
            R_version="4.1.3",
            keras_version="2.13.0",
            tensorflow_version="2.13.0",
            python_version="3.10.21",
            execution_backend="cpu",
            created_at=datetime(2026, 9, 15),
        )

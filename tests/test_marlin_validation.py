from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import pytest

from ontseq_platform.marlin_contracts import (
    MarlinClassificationDecision,
    MarlinFeatureSummary,
    MarlinGroupedScore,
    MarlinModelUnitScore,
    MarlinPredictionReport,
    MarlinSourceKind,
)
from ontseq_platform.marlin_validation import (
    MarlinValidationManifest,
    MarlinValidationSample,
    execute_marlin_validation,
)
from ontseq_platform.models import FileFingerprint, GenomeBuild, ModuleRunStatus
from ontseq_platform.reference import sha256_file


def _sample(
    tmp_path: Path,
    *,
    sha: str | None = None,
    expected: str | None = "Class A",
) -> MarlinValidationSample:
    path = tmp_path / "sample.txt.gz"
    path.write_bytes(b"processed-cpg-fixture")
    return MarlinValidationSample(
        accession="GSM_TEST_001",
        sample_id="AL_TEST_001",
        public_source_uri="https://example.invalid/GSM_TEST_001",
        local_input_path=path,
        input_sha256=sha or sha256_file(path),
        genome_build=GenomeBuild.GRCH37,
        expected_class=expected,
    )


def _manifest(sample: MarlinValidationSample) -> MarlinValidationManifest:
    return MarlinValidationManifest(
        validation_id="GSE280090_TEST_V1",
        artifact_lock_id="LOCK",
        runtime_profile_id="PROFILE",
        samples=[sample],
        locked_at=datetime(2026, 9, 13, tzinfo=UTC),
    )


def _report(
    *,
    input_sha: str,
    feature_sha: str = "a" * 64,
    decision: MarlinClassificationDecision = MarlinClassificationDecision.HIGH_CONFIDENCE,
    top_class: str = "Class A",
) -> MarlinPredictionReport:
    observed = 1
    top_score = (
        1.0
        if decision is MarlinClassificationDecision.HIGH_CONFIDENCE
        else 0.7
    )
    return MarlinPredictionReport(
        sample_id="AL_TEST_001",
        status=ModuleRunStatus.COMPLETED,
        decision=decision,
        source_kind=MarlinSourceKind.PRECOMPUTED_METHYLATION,
        genome_build=GenomeBuild.GRCH37,
        input_fingerprint=FileFingerprint(size_bytes=21, sha256=input_sha),
        feature_summary=MarlinFeatureSummary(
            observed_model_feature_count=observed,
            explicit_na_feature_count=0,
            absent_feature_count=357340 - observed,
            non_model_probe_count=0,
            observed_fraction=observed / 357340,
            feature_vector_sha256=feature_sha,
            feature_artifact_sha256="b" * 64,
        ),
        raw_model_scores=[
            MarlinModelUnitScore(
                model_id=i,
                label=f"model-unit-{i}",
                score=top_score if i == 1 else (1.0 - top_score if i == 2 else 0.0),
            )
            for i in range(1, 43)
        ],
        class_scores=[
            MarlinGroupedScore(label=top_class, score=top_score),
            MarlinGroupedScore(label="Other", score=1.0 - top_score),
        ],
        family_scores=[
            MarlinGroupedScore(label="Family A", score=top_score),
            MarlinGroupedScore(label="Family Other", score=1.0 - top_score),
        ],
        lineage_scores=[
            MarlinGroupedScore(label="Lineage A", score=top_score),
            MarlinGroupedScore(label="Lineage Other", score=1.0 - top_score),
        ],
        top_class=top_class,
        top_class_score=top_score,
        artifact_lock_id="LOCK",
    )


def test_validation_refuses_input_checksum_mismatch(tmp_path: Path) -> None:
    sample = _sample(tmp_path, sha="0" * 64)

    def should_not_run(
        _: MarlinValidationSample,
    ) -> tuple[MarlinPredictionReport, str | None, float]:
        raise AssertionError("classification must not run before checksum preflight passes")

    with pytest.raises(ValueError, match="SHA-256"):
        execute_marlin_validation(_manifest(sample), classify_sample=should_not_run)


def test_repeated_identical_run_requires_same_feature_digest_and_decision(tmp_path: Path) -> None:
    sample = _sample(tmp_path)
    report = _report(input_sha=sample.input_sha256)
    calls = 0

    def classify(_: MarlinValidationSample) -> tuple[MarlinPredictionReport, str | None, float]:
        nonlocal calls
        calls += 1
        return report, "c" * 64, 0.25

    result = execute_marlin_validation(_manifest(sample), classify_sample=classify)
    assert calls == 2
    assert result.sample_count == 1
    assert result.deterministic_rerun_failures == 0
    assert result.high_confidence_count == 1
    assert result.concordant_count == 1
    assert result.results[0].feature_vector_sha256 == "a" * 64
    assert result.results[0].raw_score_sha256 == "c" * 64


def test_nondeterministic_second_run_is_failed_not_silently_pooled(tmp_path: Path) -> None:
    sample = _sample(tmp_path)
    first = _report(input_sha=sample.input_sha256, feature_sha="a" * 64)
    second = _report(input_sha=sample.input_sha256, feature_sha="f" * 64)
    reports = iter([(first, "c" * 64, 0.2), (second, "d" * 64, 0.2)])

    result = execute_marlin_validation(
        _manifest(sample), classify_sample=lambda _: next(reports)
    )
    assert result.deterministic_rerun_failures == 1
    assert result.failed_count == 1
    assert result.results[0].status == "FAILED"
    assert "determin" in result.results[0].reason.lower()


def test_unknown_is_counted_separately_from_no_call(tmp_path: Path) -> None:
    sample = _sample(tmp_path)
    report = _report(
        input_sha=sample.input_sha256,
        decision=MarlinClassificationDecision.UNKNOWN,
    )
    result = execute_marlin_validation(
        _manifest(sample), classify_sample=lambda _: (report, "c" * 64, 0.1)
    )
    assert result.unknown_count == 1
    assert result.no_call_count == 0


def test_manifest_rejects_non_grch37_or_threshold_drift(tmp_path: Path) -> None:
    sample = _sample(tmp_path)
    with pytest.raises(Exception):
        MarlinValidationManifest(
            validation_id="BAD",
            artifact_lock_id="LOCK",
            runtime_profile_id="PROFILE",
            confidence_threshold=0.7,
            samples=[sample],
            locked_at=datetime(2026, 9, 13, tzinfo=UTC),
        )
    with pytest.raises(Exception):
        MarlinValidationSample.model_validate(
            {**sample.model_dump(), "genome_build": GenomeBuild.GRCH38}
        )


def test_failed_execution_reason_does_not_export_local_exception_text(tmp_path: Path) -> None:
    sample = _sample(tmp_path)

    def fail(_: MarlinValidationSample):
        raise ValueError("secret local path /mnt/c/Users/example/patient-file.txt")

    result = execute_marlin_validation(_manifest(sample), classify_sample=fail)
    assert result.failed_count == 1
    assert "patient-file" not in result.results[0].reason
    assert result.results[0].reason == "Classification execution failed: ValueError"

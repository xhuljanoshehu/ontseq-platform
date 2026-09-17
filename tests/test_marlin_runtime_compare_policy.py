from __future__ import annotations

from datetime import UTC, datetime

import pytest

from ontseq_platform.marlin_contracts import (
    MarlinArtifactLock,
    MarlinModelUnitScore,
    MarlinRuntimeCompatibilityProfile,
)
from ontseq_platform.marlin_runtime import MarlinRuntimeProbeReport, MarlinRuntimeResult
from ontseq_platform.marlin_runtime_compare import (
    compare_marlin_runtime_results,
    require_same_marlin_artifact_set,
    runtime_identity_from_probe,
)
from ontseq_platform.models import GenomeBuild


def _lock(*, lock_id: str, r_version: str) -> MarlinArtifactLock:
    return MarlinArtifactLock(
        lock_id=lock_id,
        model_version="1.0.0",
        model_source_uri="model",
        model_sha256="0" * 64,
        code_source_uri="code",
        code_version_or_commit="442aa603415a54f62e7367794f9a31c6bc20fc2d",
        code_manifest_sha256="1" * 64,
        feature_source_uri="features",
        feature_sha256="2" * 64,
        canonical_feature_list_sha256="3" * 64,
        class_annotation_source_uri="classes",
        class_annotation_sha256="4" * 64,
        probe_resource_source_uri="probes",
        probe_resource_sha256="5" * 64,
        genome_build=GenomeBuild.GRCH37,
        R_version=r_version,
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version_if_used="3.10.21",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 17, tzinfo=UTC),
    )


def _identity(lock: MarlinArtifactLock, runtime_id: str):
    return runtime_identity_from_probe(
        lock,
        MarlinRuntimeProbeReport(
            R_version=lock.R_version,
            keras_version=lock.keras_version,
            tensorflow_version=lock.tensorflow_version,
            python_version=lock.python_version_if_used or "3.10.21",
            execution_backend=lock.execution_backend,
        ),
        runtime_id=runtime_id,
        created_at=datetime(2026, 9, 17, 6, 0, tzinfo=UTC),
    )


def test_direct_comparison_rejects_widened_runtime_policy() -> None:
    reference_lock = _lock(lock_id="REFERENCE_LOCK", r_version="4.1.3")
    candidate_lock = _lock(lock_id="CANDIDATE_LOCK", r_version="4.2.3")
    artifact_identity = require_same_marlin_artifact_set(reference_lock, candidate_lock)
    reference_scores = [0.6, 0.4] + [0.0] * 40
    candidate_scores = [0.60000015, 0.39999985] + [0.0] * 40
    profile = MarlinRuntimeCompatibilityProfile(
        profile_id="REFERENCE_PROFILE",
        reference_runtime_lock_id=reference_lock.lock_id,
        feature_vector_sha256="a" * 64,
        reference_scores=reference_scores,
        absolute_score_tolerance=2e-7,
        score_sum_tolerance=1e-5,
        top_model_unit_index=0,
        execution_backend="cpu",
        created_at=datetime(2026, 9, 17, 6, 0, tzinfo=UTC),
    )
    candidate_result = MarlinRuntimeResult(
        artifact_lock_id=candidate_lock.lock_id,
        feature_vector_sha256=profile.feature_vector_sha256,
        model_scores=[
            MarlinModelUnitScore(model_id=index + 1, label=f"unit-{index + 1}", score=score)
            for index, score in enumerate(candidate_scores)
        ],
        raw_score_sha256="b" * 64,
        execution_backend="cpu",
    )

    with pytest.raises(ValueError, match="absolute score tolerance"):
        compare_marlin_runtime_results(
            comparison_id="POLICY_TEST",
            artifact_set_identity=artifact_identity,
            reference_lock=reference_lock,
            candidate_lock=candidate_lock,
            reference_runtime_identity=_identity(reference_lock, "REFERENCE_RUNTIME"),
            candidate_runtime_identity=_identity(candidate_lock, "CANDIDATE_RUNTIME"),
            reference_profile=profile,
            candidate_result=candidate_result,
            created_at=datetime(2026, 9, 17, 6, 5, tzinfo=UTC),
        )

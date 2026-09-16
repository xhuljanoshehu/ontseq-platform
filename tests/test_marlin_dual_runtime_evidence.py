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
    derive_marlin_artifact_set_identity,
    verify_marlin_artifact_set_identity,
)
from ontseq_platform.models import GenomeBuild

SHA0 = "0" * 64
SHA1 = "1" * 64
SHA2 = "2" * 64
SHA3 = "3" * 64
SHA4 = "4" * 64
SHA5 = "5" * 64
SHA6 = "6" * 64
FEATURE_SHA = "a" * 64


def _lock(*, candidate: bool = False) -> MarlinArtifactLock:
    return MarlinArtifactLock(
        lock_id="CANDIDATE_LOCK" if candidate else "REFERENCE_LOCK",
        model_version="1.0.0",
        model_source_uri="model",
        model_sha256=SHA0,
        code_source_uri="code",
        code_version_or_commit="442aa603415a54f62e7367794f9a31c6bc20fc2d",
        code_manifest_sha256=SHA1,
        feature_source_uri="features",
        feature_sha256=SHA2,
        canonical_feature_list_sha256=SHA3,
        class_annotation_source_uri="classes",
        class_annotation_sha256=SHA4,
        probe_resource_source_uri="probes",
        probe_resource_sha256=SHA5,
        genome_build=GenomeBuild.GRCH37,
        R_version="4.2.3" if candidate else "4.1.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version_if_used="3.10.21",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 16, tzinfo=UTC),
    )


def _probe(lock: MarlinArtifactLock) -> MarlinRuntimeProbeReport:
    return MarlinRuntimeProbeReport(
        R_version=lock.R_version,
        keras_version=lock.keras_version,
        tensorflow_version=lock.tensorflow_version,
        python_version=lock.python_version_if_used or "3.10.21",
        execution_backend="cpu",
    )


def _profile(created_at: datetime) -> MarlinRuntimeCompatibilityProfile:
    return MarlinRuntimeCompatibilityProfile(
        profile_id="REFERENCE_PROFILE",
        reference_runtime_lock_id="REFERENCE_LOCK",
        feature_vector_sha256=FEATURE_SHA,
        reference_scores=[0.6, 0.4] + [0.0] * 40,
        absolute_score_tolerance=1e-7,
        score_sum_tolerance=1e-5,
        top_model_unit_index=0,
        execution_backend="cpu",
        created_at=created_at,
    )


def _candidate_result() -> MarlinRuntimeResult:
    scores = [0.6, 0.4] + [0.0] * 40
    return MarlinRuntimeResult(
        artifact_lock_id="CANDIDATE_LOCK",
        feature_vector_sha256=FEATURE_SHA,
        model_scores=[
            MarlinModelUnitScore(model_id=index + 1, label=f"unit-{index + 1}", score=score)
            for index, score in enumerate(scores)
        ],
        raw_score_sha256=SHA6,
        execution_backend="cpu",
    )


def test_artifact_set_identity_rejects_tampered_digest() -> None:
    identity = derive_marlin_artifact_set_identity(_lock())
    tampered = identity.model_copy(update={"artifact_set_sha256": "f" * 64})

    with pytest.raises(ValueError, match="digest"):
        verify_marlin_artifact_set_identity(tampered)


def test_marlin_v1_runtime_profile_policy_rejects_widened_tolerance() -> None:
    from ontseq_platform.marlin_runtime_freeze import verify_marlin_v1_runtime_profile_policy

    profile = _profile(datetime(2026, 9, 16, 4, 0, tzinfo=UTC)).model_copy(
        update={"absolute_score_tolerance": 2e-7}
    )

    with pytest.raises(ValueError, match="absolute score tolerance"):
        verify_marlin_v1_runtime_profile_policy(profile)


def test_marlin_v1_runtime_profile_policy_accepts_frozen_policy() -> None:
    from ontseq_platform.marlin_runtime_freeze import verify_marlin_v1_runtime_profile_policy

    verify_marlin_v1_runtime_profile_policy(_profile(datetime(2026, 9, 16, 4, 0, tzinfo=UTC)))

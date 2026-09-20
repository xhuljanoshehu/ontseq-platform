from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import pytest

from ontseq_platform.marlin_contracts import MarlinModelUnitScore
from ontseq_platform.marlin_runtime import (
    MarlinRuntimeResult,
    create_runtime_compatibility_profile,
    verify_runtime_compatibility,
)


def _runtime_result(
    *,
    lock_id: str = "MARLIN_V1_RUNTIME_TEST",
    backend: str = "cpu",
    scores: tuple[float, ...] | None = None,
) -> MarlinRuntimeResult:
    values = scores if scores is not None else tuple([0.6, 0.4] + [0.0] * 40)
    return MarlinRuntimeResult(
        artifact_lock_id=lock_id,
        feature_vector_sha256="e" * 64,
        model_scores=[
            MarlinModelUnitScore(model_id=index + 1, label=f"unit-{index + 1}", score=score)
            for index, score in enumerate(values)
        ],
        raw_score_sha256="d" * 64,
        execution_backend=backend,
    )


def _profile():
    return create_runtime_compatibility_profile(
        _runtime_result(),
        profile_id="MARLIN_CPU_REF_V1",
        absolute_score_tolerance=1e-7,
        score_sum_tolerance=1e-5,
        created_at=datetime(2026, 9, 13, tzinfo=UTC),
    )


def test_create_profile_requires_explicit_positive_tolerance() -> None:
    with pytest.raises(ValueError, match="positive"):
        create_runtime_compatibility_profile(
            _runtime_result(),
            profile_id="MARLIN_CPU_REF_V1",
            absolute_score_tolerance=0.0,
            score_sum_tolerance=1e-5,
            created_at=datetime(2026, 9, 13, tzinfo=UTC),
        )


def test_profile_freezes_reference_result() -> None:
    result = _runtime_result()
    profile = _profile()
    assert profile.reference_runtime_lock_id == result.artifact_lock_id
    assert profile.feature_vector_sha256 == result.feature_vector_sha256
    assert profile.reference_scores == [item.score for item in result.model_scores]
    assert profile.top_model_unit_index == 0
    assert profile.execution_backend == "cpu"


def test_exact_candidate_passes_frozen_profile() -> None:
    verify_runtime_compatibility(_profile(), _runtime_result())


def test_compatibility_requires_same_top_model_unit() -> None:
    reference = _runtime_result(scores=tuple([0.50000001, 0.49999999] + [0.0] * 40))
    profile = create_runtime_compatibility_profile(
        reference,
        profile_id="MARLIN_CPU_REF_V1",
        absolute_score_tolerance=1.0,
        score_sum_tolerance=1e-5,
        created_at=datetime(2026, 9, 13, tzinfo=UTC),
    )
    candidate = _runtime_result(scores=tuple([0.49999999, 0.50000001] + [0.0] * 40))
    with pytest.raises(ValueError, match="top model unit"):
        verify_runtime_compatibility(profile, candidate)


def test_compatibility_rejects_score_delta_above_frozen_tolerance() -> None:
    candidate = _runtime_result(scores=tuple([0.599, 0.401] + [0.0] * 40))
    with pytest.raises(ValueError, match="tolerance"):
        verify_runtime_compatibility(_profile(), candidate)


def test_compatibility_rejects_different_runtime_lock_or_backend() -> None:
    profile = _profile()
    with pytest.raises(ValueError, match="runtime lock"):
        verify_runtime_compatibility(profile, _runtime_result(lock_id="OTHER_LOCK"))
    with pytest.raises(ValueError, match="backend"):
        verify_runtime_compatibility(profile, _runtime_result(backend="gpu"))


def test_compatibility_rejects_different_feature_vector() -> None:
    candidate = _runtime_result().model_copy(update={"feature_vector_sha256": "f" * 64})
    with pytest.raises(ValueError, match="feature vector"):
        verify_runtime_compatibility(_profile(), candidate)


def test_marlin_technical_config_pins_published_semantics() -> None:
    import yaml

    config = Path("configs/methylation/marlin_v1.technical.yaml")
    data = yaml.safe_load(config.read_text(encoding="utf-8"))
    assert data["status"] == "technical_defaults_and_published_semantics_only"
    assert data["genome_build"] == "GRCh37"
    assert data["expected_feature_count"] == 357340
    assert data["expected_model_unit_count"] == 42
    assert data["confidence_threshold"] == 0.8
    assert "runtime_compatibility_tolerance" not in data

    dual = data["dual_runtime_compatibility"]
    assert dual == {
        "status": "engineering_gate",
        "absolute_score_tolerance": 1e-7,
        "score_sum_tolerance": 1e-5,
        "fixture_generator_version": "sha256-index-mod3-v1",
        "research_only": True,
    }

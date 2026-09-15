from __future__ import annotations

import hashlib
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path

from .execution import CommandRunner
from .marlin_contracts import (
    MarlinArtifactLock,
    MarlinFeatureSummary,
    MarlinRuntimeCompatibilityProfile,
)
from .marlin_features import MarlinFeatureVector, marlin_feature_vector_sha256
from .marlin_runtime import (
    MarlinRuntimeProbeReport,
    MarlinRuntimeResult,
    create_runtime_compatibility_profile,
    run_marlin_inference,
)

_EXPECTED_FEATURE_COUNT = 357_340
_FIXTURE_SALT = b"ONTSeq-MARLIN-runtime-fixture-v1\x00"
_FIXTURE_VALUE_MAP = (-1.0, 0.0, 1.0)

MARLIN_RUNTIME_FIXTURE_GENERATOR_VERSION = "sha256-index-mod3-v1"
MARLIN_RUNTIME_FIXTURE_ABSOLUTE_SCORE_TOLERANCE = 1e-7
MARLIN_RUNTIME_FIXTURE_SCORE_SUM_TOLERANCE = 1e-5


def generate_frozen_runtime_fixture(lock: MarlinArtifactLock) -> MarlinFeatureVector:
    """Generate the canonical non-biological MARLIN runtime fixture.

    The fixture is independent of biological data and of model outcomes. Each zero-based
    feature index is hashed under a versioned fixed salt; the first digest byte modulo three
    selects -1, 0 or +1. The exact vector digest is regression-locked in tests.
    """
    if lock.expected_feature_count != _EXPECTED_FEATURE_COUNT:
        raise ValueError("MARLIN runtime freeze requires exactly 357340 locked features")

    values = tuple(
        _FIXTURE_VALUE_MAP[hashlib.sha256(_FIXTURE_SALT + index.to_bytes(8, "big")).digest()[0] % 3]
        for index in range(_EXPECTED_FEATURE_COUNT)
    )
    observed_count = sum(value != 0.0 for value in values)
    absent_count = _EXPECTED_FEATURE_COUNT - observed_count
    vector_sha256 = marlin_feature_vector_sha256(values)
    summary = MarlinFeatureSummary(
        observed_model_feature_count=observed_count,
        explicit_na_feature_count=0,
        absent_feature_count=absent_count,
        non_model_probe_count=0,
        observed_fraction=observed_count / _EXPECTED_FEATURE_COUNT,
        feature_vector_sha256=vector_sha256,
        feature_artifact_sha256=lock.canonical_feature_list_sha256,
    )
    return MarlinFeatureVector(values=values, summary=summary)


def render_frozen_runtime_fixture(vector: MarlinFeatureVector) -> str:
    """Render the canonical fixture as one integer per LF-terminated line."""
    if len(vector.values) != _EXPECTED_FEATURE_COUNT:
        raise ValueError("MARLIN runtime fixture requires exactly 357340 feature values")
    observed_sha256 = marlin_feature_vector_sha256(vector.values)
    if observed_sha256 != vector.summary.feature_vector_sha256:
        raise ValueError("MARLIN runtime fixture digest does not match its feature summary")
    return "".join(f"{int(value)}\n" for value in vector.values)


def verify_runtime_probe_matches_lock(
    lock: MarlinArtifactLock,
    probe: MarlinRuntimeProbeReport,
) -> None:
    """Require the live model runtime to match the identity recorded by the artifact lock."""
    expected = {
        "R_version": lock.R_version,
        "keras_version": lock.keras_version,
        "tensorflow_version": lock.tensorflow_version,
        "python_version": lock.python_version_if_used,
        "execution_backend": lock.execution_backend,
    }
    observed = {
        "R_version": probe.R_version,
        "keras_version": probe.keras_version,
        "tensorflow_version": probe.tensorflow_version,
        "python_version": probe.python_version,
        "execution_backend": probe.execution_backend,
    }
    for field, expected_value in expected.items():
        observed_value = observed[field]
        if expected_value != observed_value:
            raise ValueError(
                f"MARLIN live runtime {field} differs from artifact lock: "
                f"expected {expected_value!r}, observed {observed_value!r}"
            )


def freeze_runtime_compatibility(
    lock: MarlinArtifactLock,
    *,
    model_path: Path,
    inference_script: Path,
    runner: CommandRunner,
    work_dir: Path,
    profile_id: str,
    created_at: datetime,
    rscript_path: str = "Rscript",
    timeout_seconds: int = 300,
    class_labels: Sequence[str] | None = None,
) -> tuple[MarlinFeatureVector, MarlinRuntimeResult, MarlinRuntimeCompatibilityProfile]:
    """Run the fixed non-biological vector once and freeze its 42-score reference profile."""
    vector = generate_frozen_runtime_fixture(lock)
    labels = (
        tuple(class_labels)
        if class_labels is not None
        else tuple(f"model-unit-{index}" for index in range(1, 43))
    )
    result = run_marlin_inference(
        vector,
        lock,
        model_path=model_path,
        class_labels=labels,
        inference_script=inference_script,
        runner=runner,
        work_dir=work_dir,
        rscript_path=rscript_path,
        timeout_seconds=timeout_seconds,
    )
    profile = create_runtime_compatibility_profile(
        result,
        profile_id=profile_id,
        absolute_score_tolerance=MARLIN_RUNTIME_FIXTURE_ABSOLUTE_SCORE_TOLERANCE,
        score_sum_tolerance=MARLIN_RUNTIME_FIXTURE_SCORE_SUM_TOLERANCE,
        created_at=created_at,
    )
    return vector, result, profile


__all__ = [
    "MARLIN_RUNTIME_FIXTURE_ABSOLUTE_SCORE_TOLERANCE",
    "MARLIN_RUNTIME_FIXTURE_GENERATOR_VERSION",
    "MARLIN_RUNTIME_FIXTURE_SCORE_SUM_TOLERANCE",
    "freeze_runtime_compatibility",
    "generate_frozen_runtime_fixture",
    "render_frozen_runtime_fixture",
    "verify_runtime_probe_matches_lock",
]

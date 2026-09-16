from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

from .execution import CommandRunner
from .marlin_contracts import (
    MarlinArtifactLock,
    MarlinArtifactSetIdentity,
    MarlinDualRuntimeCompatibilityReport,
    MarlinDualRuntimeVerdict,
    MarlinRuntimeCompatibilityProfile,
    MarlinRuntimeIdentity,
)
from .marlin_runtime import (
    MarlinRuntimeProbeReport,
    MarlinRuntimeResult,
    load_frozen_runtime_fixture,
    probe_marlin_runtime,
    run_marlin_inference,
)
from .reference import sha256_file


def derive_marlin_artifact_set_identity(lock: MarlinArtifactLock) -> MarlinArtifactSetIdentity:
    """Derive the immutable MARLIN classifier identity independent of execution runtime."""
    payload: dict[str, object] = {
        "model_version": lock.model_version,
        "model_sha256": lock.model_sha256,
        "code_version_or_commit": lock.code_version_or_commit,
        "code_manifest_sha256": lock.code_manifest_sha256,
        "feature_sha256": lock.feature_sha256,
        "canonical_feature_list_sha256": lock.canonical_feature_list_sha256,
        "class_annotation_sha256": lock.class_annotation_sha256,
        "probe_resource_sha256": lock.probe_resource_sha256,
        "genome_build": lock.genome_build.value,
        "expected_feature_count": lock.expected_feature_count,
        "expected_model_unit_count": lock.expected_model_unit_count,
        "preprocessing_contract_version": lock.preprocessing_contract_version,
        "runtime_contract_version": lock.runtime_contract_version,
    }
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    digest = hashlib.sha256(encoded).hexdigest()
    return MarlinArtifactSetIdentity(
        artifact_set_sha256=digest,
        model_version=lock.model_version,
        model_sha256=lock.model_sha256,
        code_version_or_commit=lock.code_version_or_commit,
        code_manifest_sha256=lock.code_manifest_sha256,
        feature_sha256=lock.feature_sha256,
        canonical_feature_list_sha256=lock.canonical_feature_list_sha256,
        class_annotation_sha256=lock.class_annotation_sha256,
        probe_resource_sha256=lock.probe_resource_sha256,
        genome_build=lock.genome_build,
        expected_feature_count=lock.expected_feature_count,
        expected_model_unit_count=lock.expected_model_unit_count,
        preprocessing_contract_version=lock.preprocessing_contract_version,
        runtime_contract_version=lock.runtime_contract_version,
    )


def require_same_marlin_artifact_set(
    reference_lock: MarlinArtifactLock,
    candidate_lock: MarlinArtifactLock,
) -> MarlinArtifactSetIdentity:
    """Require two runtime-specific execution locks to bind the same immutable artifacts."""
    reference = derive_marlin_artifact_set_identity(reference_lock)
    candidate = derive_marlin_artifact_set_identity(candidate_lock)
    if reference.artifact_set_sha256 != candidate.artifact_set_sha256:
        raise ValueError(
            "MARLIN reference and candidate execution locks describe different artifact sets"
        )
    return reference


def _require_runtime_identity_matches_lock(
    identity: MarlinRuntimeIdentity,
    lock: MarlinArtifactLock,
    *,
    role: str,
) -> None:
    expected = {
        "R_version": lock.R_version,
        "keras_version": lock.keras_version,
        "tensorflow_version": lock.tensorflow_version,
        "python_version": lock.python_version_if_used,
        "execution_backend": lock.execution_backend,
        "runtime_contract_version": lock.runtime_contract_version,
    }
    observed = {
        "R_version": identity.R_version,
        "keras_version": identity.keras_version,
        "tensorflow_version": identity.tensorflow_version,
        "python_version": identity.python_version,
        "execution_backend": identity.execution_backend,
        "runtime_contract_version": identity.runtime_contract_version,
    }
    for field, expected_value in expected.items():
        observed_value = observed[field]
        if observed_value != expected_value:
            raise ValueError(
                f"MARLIN {role} runtime identity {field} differs from execution lock: "
                f"expected {expected_value!r}, observed {observed_value!r}"
            )


def runtime_identity_from_probe(
    lock: MarlinArtifactLock,
    probe: MarlinRuntimeProbeReport,
    *,
    runtime_id: str,
    created_at: datetime,
) -> MarlinRuntimeIdentity:
    """Bind live-probed runtime identity to the exact runtime fields of an execution lock."""
    identity = MarlinRuntimeIdentity(
        runtime_id=runtime_id,
        R_version=probe.R_version,
        keras_version=probe.keras_version,
        tensorflow_version=probe.tensorflow_version,
        python_version=probe.python_version,
        execution_backend=probe.execution_backend,
        runtime_contract_version=lock.runtime_contract_version,
        created_at=created_at,
    )
    _require_runtime_identity_matches_lock(identity, lock, role="live")
    return identity


def compare_marlin_runtime_results(
    *,
    comparison_id: str,
    artifact_set_identity: MarlinArtifactSetIdentity,
    reference_lock: MarlinArtifactLock,
    candidate_lock: MarlinArtifactLock,
    reference_runtime_identity: MarlinRuntimeIdentity,
    candidate_runtime_identity: MarlinRuntimeIdentity,
    reference_profile: MarlinRuntimeCompatibilityProfile,
    candidate_result: MarlinRuntimeResult,
    created_at: datetime,
) -> MarlinDualRuntimeCompatibilityReport:
    """Compare one candidate runtime result with a frozen reference numerical oracle."""
    expected_artifact_set = require_same_marlin_artifact_set(reference_lock, candidate_lock)
    if artifact_set_identity != expected_artifact_set:
        raise ValueError("MARLIN supplied artifact-set identity does not match execution locks")

    _require_runtime_identity_matches_lock(
        reference_runtime_identity,
        reference_lock,
        role="reference",
    )
    _require_runtime_identity_matches_lock(
        candidate_runtime_identity,
        candidate_lock,
        role="candidate",
    )

    if reference_profile.reference_runtime_lock_id != reference_lock.lock_id:
        raise ValueError("MARLIN reference profile belongs to a different reference execution lock")
    if reference_profile.execution_backend != reference_lock.execution_backend:
        raise ValueError("MARLIN reference profile backend differs from reference execution lock")
    if candidate_result.artifact_lock_id != candidate_lock.lock_id:
        raise ValueError("MARLIN candidate result belongs to a different candidate execution lock")
    if candidate_result.execution_backend != candidate_lock.execution_backend:
        raise ValueError("MARLIN candidate result backend differs from candidate execution lock")
    if candidate_result.feature_vector_sha256 != reference_profile.feature_vector_sha256:
        raise ValueError("MARLIN candidate feature vector differs from frozen reference profile")

    reference_scores = tuple(reference_profile.reference_scores)
    candidate_scores = tuple(item.score for item in candidate_result.model_scores)
    differences = tuple(
        abs(reference - candidate)
        for reference, candidate in zip(reference_scores, candidate_scores, strict=True)
    )
    all_within = all(
        difference <= reference_profile.absolute_score_tolerance for difference in differences
    )
    reference_top = reference_profile.top_model_unit_index
    candidate_top = max(range(42), key=candidate_scores.__getitem__)
    top_matches = candidate_top == reference_top
    softmax_ok = (
        abs(sum(reference_scores) - 1.0) <= reference_profile.score_sum_tolerance
        and abs(sum(candidate_scores) - 1.0) <= reference_profile.score_sum_tolerance
    )
    verdict = (
        MarlinDualRuntimeVerdict.PASS
        if all_within and top_matches and softmax_ok
        else MarlinDualRuntimeVerdict.FAIL
    )

    return MarlinDualRuntimeCompatibilityReport(
        comparison_id=comparison_id,
        artifact_set_identity=artifact_set_identity,
        reference_artifact_lock_id=reference_lock.lock_id,
        candidate_artifact_lock_id=candidate_lock.lock_id,
        feature_vector_sha256=reference_profile.feature_vector_sha256,
        reference_profile_id=reference_profile.profile_id,
        reference_runtime_identity=reference_runtime_identity,
        candidate_runtime_identity=candidate_runtime_identity,
        absolute_score_tolerance=reference_profile.absolute_score_tolerance,
        score_sum_tolerance=reference_profile.score_sum_tolerance,
        reference_scores=list(reference_scores),
        candidate_scores=list(candidate_scores),
        absolute_differences=list(differences),
        max_absolute_difference=max(differences),
        reference_top_model_unit_index=reference_top,
        candidate_top_model_unit_index=candidate_top,
        all_scores_within_tolerance=all_within,
        top_model_unit_matches=top_matches,
        softmax_invariants_pass=softmax_ok,
        verdict=verdict,
        created_at=created_at,
    )


def execute_marlin_dual_runtime_comparison(
    *,
    comparison_id: str,
    reference_lock: MarlinArtifactLock,
    candidate_lock: MarlinArtifactLock,
    reference_runtime_identity: MarlinRuntimeIdentity,
    reference_profile: MarlinRuntimeCompatibilityProfile,
    runtime_fixture_path: Path,
    model_path: Path,
    inference_script: Path,
    candidate_probe_script: Path,
    runner: CommandRunner,
    work_dir: Path,
    candidate_rscript_path: str,
    created_at: datetime,
) -> MarlinDualRuntimeCompatibilityReport:
    """Qualify one live candidate runtime against a frozen MARLIN reference profile."""
    artifact_set_identity = require_same_marlin_artifact_set(reference_lock, candidate_lock)
    _require_runtime_identity_matches_lock(
        reference_runtime_identity,
        reference_lock,
        role="reference",
    )
    if reference_profile.reference_runtime_lock_id != reference_lock.lock_id:
        raise ValueError("MARLIN reference profile belongs to a different reference execution lock")
    if reference_profile.execution_backend != reference_lock.execution_backend:
        raise ValueError("MARLIN reference profile backend differs from reference execution lock")

    fixture = load_frozen_runtime_fixture(runtime_fixture_path, reference_lock)
    if fixture.summary.feature_vector_sha256 != reference_profile.feature_vector_sha256:
        raise ValueError("MARLIN runtime fixture feature vector differs from reference profile")

    model = Path(model_path)
    if not model.is_file():
        raise ValueError(f"MARLIN model file is missing: {model}")
    if sha256_file(model) != artifact_set_identity.model_sha256:
        raise ValueError("MARLIN model SHA-256 differs from shared artifact set")

    live_probe = probe_marlin_runtime(
        probe_script=candidate_probe_script,
        runner=runner,
        rscript_path=candidate_rscript_path,
        timeout_seconds=120,
    )
    candidate_runtime_identity = runtime_identity_from_probe(
        candidate_lock,
        live_probe,
        runtime_id=f"{candidate_lock.lock_id}:live",
        created_at=created_at,
    )
    candidate_result = run_marlin_inference(
        fixture,
        candidate_lock,
        model_path=model,
        class_labels=tuple(f"model-unit-{index}" for index in range(1, 43)),
        inference_script=inference_script,
        runner=runner,
        work_dir=Path(work_dir) / "candidate-runtime",
        rscript_path=candidate_rscript_path,
        timeout_seconds=300,
    )
    return compare_marlin_runtime_results(
        comparison_id=comparison_id,
        artifact_set_identity=artifact_set_identity,
        reference_lock=reference_lock,
        candidate_lock=candidate_lock,
        reference_runtime_identity=reference_runtime_identity,
        candidate_runtime_identity=candidate_runtime_identity,
        reference_profile=reference_profile,
        candidate_result=candidate_result,
        created_at=created_at,
    )


__all__ = [
    "compare_marlin_runtime_results",
    "derive_marlin_artifact_set_identity",
    "execute_marlin_dual_runtime_comparison",
    "require_same_marlin_artifact_set",
    "runtime_identity_from_probe",
]

from __future__ import annotations

import hashlib
import json
from datetime import datetime

from .marlin_contracts import (
    MarlinArtifactLock,
    MarlinArtifactSetIdentity,
    MarlinRuntimeIdentity,
)
from .marlin_runtime import MarlinRuntimeProbeReport


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


def runtime_identity_from_probe(
    lock: MarlinArtifactLock,
    probe: MarlinRuntimeProbeReport,
    *,
    runtime_id: str,
    created_at: datetime,
) -> MarlinRuntimeIdentity:
    """Bind live-probed runtime identity to the exact runtime fields of an execution lock."""
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
        if observed_value != expected_value:
            raise ValueError(
                f"MARLIN live runtime {field} differs from execution lock: "
                f"expected {expected_value!r}, observed {observed_value!r}"
            )

    return MarlinRuntimeIdentity(
        runtime_id=runtime_id,
        R_version=probe.R_version,
        keras_version=probe.keras_version,
        tensorflow_version=probe.tensorflow_version,
        python_version=probe.python_version,
        execution_backend=probe.execution_backend,
        runtime_contract_version=lock.runtime_contract_version,
        created_at=created_at,
    )


__all__ = [
    "derive_marlin_artifact_set_identity",
    "require_same_marlin_artifact_set",
    "runtime_identity_from_probe",
]

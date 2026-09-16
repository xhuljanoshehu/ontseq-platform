from __future__ import annotations

from datetime import UTC, datetime

import pytest
from pydantic import ValidationError

from ontseq_platform.marlin_contracts import MarlinArtifactLock
from ontseq_platform.marlin_runtime import MarlinRuntimeProbeReport
from ontseq_platform.marlin_runtime_compare import (
    derive_marlin_artifact_set_identity,
    require_same_marlin_artifact_set,
    runtime_identity_from_probe,
)
from ontseq_platform.models import GenomeBuild

SHA0 = "0" * 64
SHA1 = "1" * 64
SHA2 = "2" * 64
SHA3 = "3" * 64
SHA4 = "4" * 64
SHA5 = "5" * 64
SHA6 = "6" * 64


def _lock(**updates: object) -> MarlinArtifactLock:
    payload: dict[str, object] = {
        "lock_id": "MARLIN_REFERENCE_LOCK",
        "model_version": "1.0.0",
        "model_source_uri": "https://doi.org/10.5281/zenodo.15565404",
        "model_sha256": SHA0,
        "code_source_uri": "https://github.com/hovestadt/MARLIN",
        "code_version_or_commit": "442aa603415a54f62e7367794f9a31c6bc20fc2d",
        "code_manifest_sha256": SHA1,
        "feature_source_uri": "marlin_v1.features.RData",
        "feature_sha256": SHA2,
        "canonical_feature_list_sha256": SHA3,
        "class_annotation_source_uri": "marlin_v1.class_annotations.xlsx",
        "class_annotation_sha256": SHA4,
        "probe_resource_source_uri": "marlin_v1.probes_hg19.bed.gz",
        "probe_resource_sha256": SHA5,
        "genome_build": GenomeBuild.GRCH37,
        "R_version": "4.1.3",
        "keras_version": "2.13.0",
        "tensorflow_version": "2.13.0",
        "python_version_if_used": "3.10.21",
        "execution_backend": "cpu",
        "created_at": datetime(2026, 9, 15, 18, 0, tzinfo=UTC),
    }
    payload.update(updates)
    return MarlinArtifactLock.model_validate(payload)


def _probe(**updates: object) -> MarlinRuntimeProbeReport:
    payload: dict[str, object] = {
        "R_version": "4.1.3",
        "keras_version": "2.13.0",
        "tensorflow_version": "2.13.0",
        "python_version": "3.10.21",
        "execution_backend": "cpu",
    }
    payload.update(updates)
    return MarlinRuntimeProbeReport.model_validate(payload)


def test_artifact_set_digest_ignores_execution_runtime_and_source_uri_formatting() -> None:
    reference = _lock()
    candidate = _lock(
        lock_id="MARLIN_CANDIDATE_LOCK",
        model_source_uri="zenodo:15565404",
        code_source_uri="github:hovestadt/MARLIN@442aa603",
        feature_source_uri="local:marlin_v1.features.RData",
        class_annotation_source_uri="local:marlin_v1.class_annotations.xlsx",
        probe_resource_source_uri="local:marlin_v1.probes_hg19.bed.gz",
        R_version="4.2.3",
        python_version_if_used="3.10.22",
        execution_backend="gpu",
        created_at=datetime(2026, 9, 16, 5, 0, tzinfo=UTC),
    )

    reference_identity = derive_marlin_artifact_set_identity(reference)
    candidate_identity = derive_marlin_artifact_set_identity(candidate)

    assert reference_identity.artifact_set_sha256 == candidate_identity.artifact_set_sha256
    assert reference_identity == candidate_identity
    assert require_same_marlin_artifact_set(reference, candidate) == reference_identity


@pytest.mark.parametrize(
    ("field", "replacement"),
    [
        ("model_version", "1.0.1"),
        ("model_sha256", SHA6),
        ("code_version_or_commit", "different-commit"),
        ("code_manifest_sha256", SHA6),
        ("feature_sha256", SHA6),
        ("canonical_feature_list_sha256", SHA6),
        ("class_annotation_sha256", SHA6),
        ("probe_resource_sha256", SHA6),
    ],
)
def test_artifact_set_digest_changes_for_immutable_identity_fields(
    field: str, replacement: object
) -> None:
    reference = _lock()
    changed = _lock(**{field: replacement})

    assert (
        derive_marlin_artifact_set_identity(reference).artifact_set_sha256
        != derive_marlin_artifact_set_identity(changed).artifact_set_sha256
    )
    with pytest.raises(ValueError, match="artifact set"):
        require_same_marlin_artifact_set(reference, changed)


def test_artifact_lock_still_rejects_non_grch37_build() -> None:
    with pytest.raises(ValidationError, match="GRCh37"):
        _lock(genome_build=GenomeBuild.GRCH38)


def test_runtime_identity_is_constructed_only_from_matching_live_probe() -> None:
    lock = _lock()
    created_at = datetime(2026, 9, 16, 5, 10, tzinfo=UTC)

    identity = runtime_identity_from_probe(
        lock,
        _probe(),
        runtime_id="MARLIN_REFERENCE_RUNTIME",
        created_at=created_at,
    )

    assert identity.runtime_id == "MARLIN_REFERENCE_RUNTIME"
    assert identity.R_version == "4.1.3"
    assert identity.python_version == "3.10.21"
    assert identity.execution_backend == "cpu"
    assert identity.created_at == created_at
    assert identity.runtime_contract_version == lock.runtime_contract_version


@pytest.mark.parametrize(
    ("probe_field", "replacement"),
    [
        ("R_version", "4.2.3"),
        ("keras_version", "2.14.0"),
        ("tensorflow_version", "2.14.0"),
        ("python_version", "3.10.22"),
        ("execution_backend", "gpu"),
    ],
)
def test_runtime_identity_rejects_probe_that_disagrees_with_execution_lock(
    probe_field: str, replacement: object
) -> None:
    with pytest.raises(ValueError, match=probe_field):
        runtime_identity_from_probe(
            _lock(),
            _probe(**{probe_field: replacement}),
            runtime_id="MARLIN_REFERENCE_RUNTIME",
            created_at=datetime(2026, 9, 16, 5, 10, tzinfo=UTC),
        )

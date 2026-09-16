from __future__ import annotations

from datetime import UTC, datetime

import pytest
from pydantic import ValidationError

from ontseq_platform.marlin_contracts import (
    MarlinArtifactLock,
    MarlinDualRuntimeVerdict,
    MarlinModelUnitScore,
    MarlinRuntimeCompatibilityProfile,
)
from ontseq_platform.marlin_runtime import MarlinRuntimeProbeReport, MarlinRuntimeResult
from ontseq_platform.marlin_runtime_compare import (
    compare_marlin_runtime_results,
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
FEATURE_VECTOR_SHA = "a" * 64


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


def _candidate_lock(**updates: object) -> MarlinArtifactLock:
    payload: dict[str, object] = {
        "lock_id": "MARLIN_CANDIDATE_LOCK",
        "R_version": "4.2.3",
        "python_version_if_used": "3.10.21",
        "created_at": datetime(2026, 9, 16, 5, 0, tzinfo=UTC),
    }
    payload.update(updates)
    return _lock(**payload)


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


def _identity_for_lock(lock: MarlinArtifactLock, runtime_id: str):
    return runtime_identity_from_probe(
        lock,
        _probe(
            R_version=lock.R_version,
            keras_version=lock.keras_version,
            tensorflow_version=lock.tensorflow_version,
            python_version=lock.python_version_if_used,
            execution_backend=lock.execution_backend,
        ),
        runtime_id=runtime_id,
        created_at=datetime(2026, 9, 16, 5, 10, tzinfo=UTC),
    )


def _reference_scores() -> list[float]:
    return [0.6, 0.4] + [0.0] * 40


def _profile(
    scores: list[float] | None = None,
    *,
    absolute_score_tolerance: float = 1e-7,
    score_sum_tolerance: float = 1e-5,
) -> MarlinRuntimeCompatibilityProfile:
    values = scores or _reference_scores()
    return MarlinRuntimeCompatibilityProfile(
        profile_id="MARLIN_REFERENCE_PROFILE",
        reference_runtime_lock_id="MARLIN_REFERENCE_LOCK",
        feature_vector_sha256=FEATURE_VECTOR_SHA,
        reference_scores=values,
        absolute_score_tolerance=absolute_score_tolerance,
        score_sum_tolerance=score_sum_tolerance,
        top_model_unit_index=max(range(42), key=values.__getitem__),
        execution_backend="cpu",
        created_at=datetime(2026, 9, 16, 5, 15, tzinfo=UTC),
    )


def _runtime_result(
    candidate_lock: MarlinArtifactLock,
    scores: list[float],
    *,
    feature_vector_sha256: str = FEATURE_VECTOR_SHA,
    artifact_lock_id: str | None = None,
) -> MarlinRuntimeResult:
    return MarlinRuntimeResult(
        artifact_lock_id=artifact_lock_id or candidate_lock.lock_id,
        feature_vector_sha256=feature_vector_sha256,
        model_scores=[
            MarlinModelUnitScore(model_id=index + 1, label=f"unit-{index + 1}", score=score)
            for index, score in enumerate(scores)
        ],
        raw_score_sha256=SHA6,
        execution_backend=candidate_lock.execution_backend,
    )


def _compare(
    candidate_scores: list[float],
    *,
    reference_scores: list[float] | None = None,
    absolute_score_tolerance: float = 1e-7,
    score_sum_tolerance: float = 1e-5,
):
    reference_lock = _lock()
    candidate_lock = _candidate_lock()
    artifact_identity = require_same_marlin_artifact_set(reference_lock, candidate_lock)
    reference_identity = _identity_for_lock(reference_lock, "MARLIN_REFERENCE_RUNTIME")
    candidate_identity = _identity_for_lock(candidate_lock, "MARLIN_CANDIDATE_RUNTIME")
    profile = _profile(
        reference_scores,
        absolute_score_tolerance=absolute_score_tolerance,
        score_sum_tolerance=score_sum_tolerance,
    )
    result = _runtime_result(candidate_lock, candidate_scores)
    return compare_marlin_runtime_results(
        comparison_id="MARLIN_DUAL_RUNTIME_TEST",
        artifact_set_identity=artifact_identity,
        reference_lock=reference_lock,
        candidate_lock=candidate_lock,
        reference_runtime_identity=reference_identity,
        candidate_runtime_identity=candidate_identity,
        reference_profile=profile,
        candidate_result=result,
        created_at=datetime(2026, 9, 16, 5, 20, tzinfo=UTC),
    )


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


def test_dual_runtime_exact_tolerance_boundary_passes() -> None:
    candidate_scores = _reference_scores()
    candidate_scores[2] = 1e-7

    report = _compare(candidate_scores)

    assert report.verdict is MarlinDualRuntimeVerdict.PASS
    assert report.max_absolute_difference == pytest.approx(1e-7, rel=0, abs=0)
    assert report.all_scores_within_tolerance is True
    assert report.top_model_unit_matches is True
    assert report.softmax_invariants_pass is True


def test_dual_runtime_one_score_above_tolerance_fails() -> None:
    candidate_scores = _reference_scores()
    candidate_scores[2] = 1.0000001e-7

    report = _compare(candidate_scores)

    assert report.verdict is MarlinDualRuntimeVerdict.FAIL
    assert report.all_scores_within_tolerance is False


def test_dual_runtime_top_unit_mismatch_fails_with_scores_inside_tolerance() -> None:
    reference_scores = [0.5, 0.49999995] + [0.0] * 40
    candidate_scores = [0.49999995, 0.5] + [0.0] * 40

    report = _compare(candidate_scores, reference_scores=reference_scores)

    assert report.max_absolute_difference < 1e-7
    assert report.all_scores_within_tolerance is True
    assert report.top_model_unit_matches is False
    assert report.verdict is MarlinDualRuntimeVerdict.FAIL


def test_dual_runtime_softmax_policy_can_fail_below_runtime_parser_limit() -> None:
    candidate_scores = _reference_scores()
    candidate_scores[2] = 1e-7
    candidate_scores[3] = 1e-7

    report = _compare(candidate_scores, score_sum_tolerance=1e-7)

    assert report.all_scores_within_tolerance is True
    assert report.top_model_unit_matches is True
    assert report.softmax_invariants_pass is False
    assert report.verdict is MarlinDualRuntimeVerdict.FAIL


def test_dual_runtime_rejects_candidate_result_from_wrong_execution_lock() -> None:
    reference_lock = _lock()
    candidate_lock = _candidate_lock()
    artifact_identity = require_same_marlin_artifact_set(reference_lock, candidate_lock)
    with pytest.raises(ValueError, match="candidate.*lock"):
        compare_marlin_runtime_results(
            comparison_id="MARLIN_DUAL_RUNTIME_TEST",
            artifact_set_identity=artifact_identity,
            reference_lock=reference_lock,
            candidate_lock=candidate_lock,
            reference_runtime_identity=_identity_for_lock(
                reference_lock, "MARLIN_REFERENCE_RUNTIME"
            ),
            candidate_runtime_identity=_identity_for_lock(
                candidate_lock, "MARLIN_CANDIDATE_RUNTIME"
            ),
            reference_profile=_profile(),
            candidate_result=_runtime_result(
                candidate_lock,
                _reference_scores(),
                artifact_lock_id="WRONG_CANDIDATE_LOCK",
            ),
            created_at=datetime(2026, 9, 16, 5, 20, tzinfo=UTC),
        )


def test_dual_runtime_rejects_candidate_fixture_digest_mismatch() -> None:
    reference_lock = _lock()
    candidate_lock = _candidate_lock()
    artifact_identity = require_same_marlin_artifact_set(reference_lock, candidate_lock)
    with pytest.raises(ValueError, match="feature vector"):
        compare_marlin_runtime_results(
            comparison_id="MARLIN_DUAL_RUNTIME_TEST",
            artifact_set_identity=artifact_identity,
            reference_lock=reference_lock,
            candidate_lock=candidate_lock,
            reference_runtime_identity=_identity_for_lock(
                reference_lock, "MARLIN_REFERENCE_RUNTIME"
            ),
            candidate_runtime_identity=_identity_for_lock(
                candidate_lock, "MARLIN_CANDIDATE_RUNTIME"
            ),
            reference_profile=_profile(),
            candidate_result=_runtime_result(
                candidate_lock,
                _reference_scores(),
                feature_vector_sha256=SHA6,
            ),
            created_at=datetime(2026, 9, 16, 5, 20, tzinfo=UTC),
        )

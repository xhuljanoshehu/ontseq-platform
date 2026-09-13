from __future__ import annotations

import hashlib

import pytest

from ontseq_platform.marlin_contracts import MarlinProbeObservation
from ontseq_platform.marlin_features import (
    _build_feature_vector_for_ids,
    build_marlin_feature_vector,
)
from ontseq_platform.marlin_input import MarlinPrecomputedInput
from ontseq_platform.models import FileFingerprint, GenomeBuild

SHA = "a" * 64


def _obs(probe_id: str, value: float | None, start: int) -> MarlinProbeObservation:
    return MarlinProbeObservation(
        chromosome="chr1",
        start=start,
        end=start + 1,
        methylation_fraction=value,
        probe_id=probe_id,
    )


def test_vector_is_ordered_by_locked_features_not_input_rows() -> None:
    observations = {
        "cgB": _obs("cgB", 0.9, 2),
        "cgA": _obs("cgA", 0.49, 1),
        "cgC": _obs("cgC", None, 3),
    }
    vector = _build_feature_vector_for_ids(
        ("cgA", "cgB", "cgC", "cgD"),
        observations,
        feature_artifact_sha256=SHA,
    )
    assert vector.values == (-1.0, 1.0, 0.0, 0.0)
    assert vector.summary.observed_model_feature_count == 2
    assert vector.summary.explicit_na_feature_count == 1
    assert vector.summary.absent_feature_count == 1


def test_exact_half_is_positive() -> None:
    vector = _build_feature_vector_for_ids(
        ("cgA",), {"cgA": _obs("cgA", 0.5, 1)}, feature_artifact_sha256=SHA
    )
    assert vector.values == (1.0,)


def test_vector_digest_uses_canonical_single_byte_encoding() -> None:
    vector = _build_feature_vector_for_ids(
        ("cgA", "cgB", "cgC"),
        {"cgA": _obs("cgA", 0.1, 1), "cgB": _obs("cgB", None, 2)},
        feature_artifact_sha256=SHA,
    )
    expected = hashlib.sha256(b"\xff\x00\x00").hexdigest()
    assert vector.summary.feature_vector_sha256 == expected


def test_non_model_probe_is_counted_but_not_in_vector() -> None:
    vector = _build_feature_vector_for_ids(
        ("cgA",),
        {"cgA": _obs("cgA", 0.7, 1), "off-model": _obs("off-model", 0.2, 2)},
        feature_artifact_sha256=SHA,
    )
    assert vector.values == (1.0,)
    assert vector.summary.non_model_probe_count == 1


def test_feature_vector_is_invariant_to_observation_order() -> None:
    a = _obs("cgA", 0.2, 1)
    b = _obs("cgB", 0.8, 2)
    first = _build_feature_vector_for_ids(
        ("cgA", "cgB"), {"cgA": a, "cgB": b}, feature_artifact_sha256=SHA
    )
    second = _build_feature_vector_for_ids(
        ("cgA", "cgB"), {"cgB": b, "cgA": a}, feature_artifact_sha256=SHA
    )
    assert first.values == second.values
    assert first.summary.feature_vector_sha256 == second.summary.feature_vector_sha256


def test_duplicate_locked_feature_ids_are_rejected() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        _build_feature_vector_for_ids(
            ("cgA", "cgA"), {"cgA": _obs("cgA", 0.7, 1)}, feature_artifact_sha256=SHA
        )


def test_production_builder_requires_exact_marlin_v1_count() -> None:
    source = MarlinPrecomputedInput(
        genome_build=GenomeBuild.GRCH37,
        input_fingerprint=FileFingerprint(size_bytes=10, sha256="b" * 64),
        observations=[_obs("cgA", 0.7, 1)],
    )
    with pytest.raises(ValueError, match="357340"):
        build_marlin_feature_vector(source, ("cgA",), feature_artifact_sha256=SHA)


def test_production_builder_refuses_non_grch37_even_if_constructed_unsafely() -> None:
    source = MarlinPrecomputedInput.model_construct(
        schema_version="0.1.0",
        adapter_id="marlin_probe_bed_v1",
        source_kind="PRECOMPUTED_METHYLATION",
        genome_build=GenomeBuild.GRCH38,
        input_fingerprint=FileFingerprint(size_bytes=10, sha256="b" * 64),
        observations=[_obs("cgA", 0.7, 1)],
    )
    full_features = tuple(f"cg{i}" for i in range(357340))
    with pytest.raises(ValueError, match="GRCh37"):
        build_marlin_feature_vector(source, full_features, feature_artifact_sha256=SHA)

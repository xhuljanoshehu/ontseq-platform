from __future__ import annotations

import gzip
import hashlib
from datetime import UTC, datetime
from pathlib import Path

import pytest

from ontseq_platform.marlin_contracts import (
    MarlinBridgeLock,
    MarlinProbeObservation,
    MarlinSourceKind,
)
from ontseq_platform.marlin_features import (
    _build_feature_vector_for_ids,
    build_marlin_feature_vector,
    build_marlin_from_modkit,
)
from ontseq_platform.marlin_input import MarlinPrecomputedInput, parse_marlin_modkit_probe_input
from ontseq_platform.models import FileFingerprint, GenomeBuild
from ontseq_platform.reference import sha256_file

SHA = "a" * 64


def _obs(probe_id: str, value: float | None, start: int) -> MarlinProbeObservation:
    return MarlinProbeObservation(
        chromosome="chr1",
        start=start,
        end=start + 1,
        methylation_fraction=value,
        probe_id=probe_id,
    )


def _bedmethyl_row(
    chromosome: str,
    start: int,
    *,
    valid: int,
    modified: int,
    code: str = "C",
    other_mod: int = 0,
) -> str:
    canonical = valid - modified - other_mod
    fields = [
        chromosome,
        str(start),
        str(start + 1),
        code,
        "0",
        ".",
        str(start),
        str(start + 1),
        "0",
        str(valid),
        "0",
        str(modified),
        str(canonical),
        str(other_mod),
        "0",
        "0",
        "0",
        "0",
    ]
    return "\t".join(fields) + "\n"


def _write_probe_resource(path: Path) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write("chr1\t10\t11\tcgA\n")
        handle.write("chr1\t11\t12\tcgA\n")
        handle.write("chr1\t20\t21\tcgB\n")


def _bridge_lock(probe_resource: Path) -> MarlinBridgeLock:
    return MarlinBridgeLock(
        bridge_id="MARLIN_MODKIT_BRIDGE_TEST",
        status="validated_same_specimen_bridge",
        genome_build=GenomeBuild.GRCH37,
        adapter_version="marlin-modkit-bridge-v1",
        modkit_version="0.6.4",
        probe_resource_sha256=sha256_file(probe_resource),
        feature_artifact_sha256=SHA,
        validation_precomputed_input_sha256="1" * 64,
        validation_modkit_bedmethyl_sha256="2" * 64,
        validation_precomputed_feature_vector_sha256="3" * 64,
        validation_modkit_feature_vector_sha256="4" * 64,
        compared_feature_count=1000,
        concordant_feature_count=990,
        feature_agreement_fraction=0.99,
        evidence_reference="same-specimen-bridge-study-001",
        validated_at=datetime(2026, 9, 13, tzinfo=UTC),
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


def test_modkit_bridge_is_disabled_without_validated_bridge_lock(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="native MARLIN bridge is not validated"):
        build_marlin_from_modkit(
            modkit_probe_calls=tmp_path / "probe_calls.bedmethyl",
            probe_resource=tmp_path / "marlin_v1.probes_hg19.bed.gz",
            bridge_lock=None,
        )


def test_modkit_probe_adapter_reproduces_published_sum_mod_over_sum_valid(
    tmp_path: Path,
) -> None:
    probe_resource = tmp_path / "probes.bed.gz"
    _write_probe_resource(probe_resource)
    calls = tmp_path / "calls.bedmethyl"
    calls.write_text(
        _bedmethyl_row("chr1", 10, valid=10, modified=7)
        + _bedmethyl_row("chr1", 11, valid=2, modified=1)
        + _bedmethyl_row("chr1", 30, valid=4, modified=4),
        encoding="utf-8",
    )
    lock = _bridge_lock(probe_resource)

    parsed = parse_marlin_modkit_probe_input(
        calls,
        probe_resource=probe_resource,
        bridge_lock=lock,
        modkit_version="0.6.4",
    )

    assert parsed.source_kind is MarlinSourceKind.MODKIT_DERIVED
    assert parsed.pileup_semantics == "5mC+5hmC-combine-mods-v1"
    by_probe = {item.probe_id: item for item in parsed.observations}
    assert by_probe["cgA"].start == 10
    assert by_probe["cgA"].end == 12
    assert by_probe["cgA"].methylation_fraction == pytest.approx(8 / 12)
    assert by_probe["cgB"].methylation_fraction is None


def test_modkit_bridge_requires_combined_cytosine_bedmethyl(tmp_path: Path) -> None:
    probe_resource = tmp_path / "probes.bed.gz"
    _write_probe_resource(probe_resource)
    calls = tmp_path / "separate-5mc.bedmethyl"
    calls.write_text(
        _bedmethyl_row("chr1", 10, valid=10, modified=7, code="m"),
        encoding="utf-8",
    )
    lock = _bridge_lock(probe_resource)

    with pytest.raises(ValueError, match="combine-mods"):
        parse_marlin_modkit_probe_input(
            calls,
            probe_resource=probe_resource,
            bridge_lock=lock,
            modkit_version="0.6.4",
        )


def test_modkit_bridge_rejects_other_mod_counts_under_combined_semantics(tmp_path: Path) -> None:
    probe_resource = tmp_path / "probes.bed.gz"
    _write_probe_resource(probe_resource)
    calls = tmp_path / "invalid-combined.bedmethyl"
    calls.write_text(
        _bedmethyl_row("chr1", 10, valid=10, modified=7, code="C", other_mod=1),
        encoding="utf-8",
    )
    lock = _bridge_lock(probe_resource)

    with pytest.raises(ValueError, match="N_other_mod"):
        parse_marlin_modkit_probe_input(
            calls,
            probe_resource=probe_resource,
            bridge_lock=lock,
            modkit_version="0.6.4",
        )


def test_modkit_bridge_builds_model_vector_only_after_evidence_gate(tmp_path: Path) -> None:
    probe_resource = tmp_path / "probes.bed.gz"
    _write_probe_resource(probe_resource)
    calls = tmp_path / "calls.bedmethyl"
    calls.write_text(
        _bedmethyl_row("chr1", 10, valid=10, modified=7)
        + _bedmethyl_row("chr1", 11, valid=2, modified=1),
        encoding="utf-8",
    )
    lock = _bridge_lock(probe_resource)
    feature_ids = ("cgA", "cgB") + tuple(f"dummy-{i}" for i in range(357338))

    vector = build_marlin_from_modkit(
        modkit_probe_calls=calls,
        probe_resource=probe_resource,
        bridge_lock=lock,
        feature_ids=feature_ids,
        feature_artifact_sha256=SHA,
        modkit_version="0.6.4",
    )

    assert vector.values[:2] == (1.0, 0.0)
    assert vector.summary.observed_model_feature_count == 1
    assert vector.summary.explicit_na_feature_count == 1


def test_modkit_bridge_rejects_resource_or_modkit_identity_drift(tmp_path: Path) -> None:
    probe_resource = tmp_path / "probes.bed.gz"
    _write_probe_resource(probe_resource)
    calls = tmp_path / "calls.bedmethyl"
    calls.write_text(_bedmethyl_row("chr1", 10, valid=10, modified=7), encoding="utf-8")
    lock = _bridge_lock(probe_resource)

    with pytest.raises(ValueError, match="modkit version"):
        parse_marlin_modkit_probe_input(
            calls,
            probe_resource=probe_resource,
            bridge_lock=lock,
            modkit_version="0.7.0",
        )

    probe_resource.write_bytes(b"changed")
    with pytest.raises(ValueError, match="probe-resource SHA-256"):
        parse_marlin_modkit_probe_input(
            calls,
            probe_resource=probe_resource,
            bridge_lock=lock,
            modkit_version="0.6.4",
        )

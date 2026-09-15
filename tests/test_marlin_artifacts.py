from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import pytest

from ontseq_platform.marlin_artifacts import (
    MarlinArtifactPaths,
    create_marlin_artifact_lock,
    load_exported_feature_ids,
    verify_marlin_artifact_lock,
)
from ontseq_platform.models import GenomeBuild


def _write_full_feature_list(path: Path) -> None:
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        for index in range(357340):
            handle.write(f"cg{index:09d}\n")


def _artifact_paths(tmp_path: Path) -> MarlinArtifactPaths:
    model = tmp_path / "marlin_v1.model.hdf5"
    feature_rdata = tmp_path / "marlin_v1.features.RData"
    canonical_features = tmp_path / "marlin_v1.features.txt"
    class_annotations = tmp_path / "marlin_v1.class_annotations.xlsx"
    probe_bed = tmp_path / "marlin_v1.probes_hg19.bed.gz"
    model.write_bytes(b"model")
    feature_rdata.write_bytes(b"rdata")
    _write_full_feature_list(canonical_features)
    class_annotations.write_bytes(b"xlsx")
    probe_bed.write_bytes(b"bedgz")
    return MarlinArtifactPaths(
        model=model,
        feature_rdata=feature_rdata,
        canonical_feature_list=canonical_features,
        class_annotations=class_annotations,
        probe_bed=probe_bed,
    )


def _create_lock(paths: MarlinArtifactPaths, *, build: GenomeBuild = GenomeBuild.GRCH37):
    return create_marlin_artifact_lock(
        paths,
        lock_id="MARLIN_V1_TEST",
        code_version_or_commit="442aa603415a54f62e7367794f9a31c6bc20fc2d",
        code_manifest_sha256="a" * 64,
        genome_build=build,
        R_version="4.1.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version_if_used="3.11.0",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 13, tzinfo=UTC),
    )


def test_load_exported_features_requires_exact_count(tmp_path: Path) -> None:
    path = tmp_path / "features.txt"
    path.write_text("cg1\ncg2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="357340"):
        load_exported_feature_ids(path)


def test_load_exported_features_rejects_duplicate_before_count_check(tmp_path: Path) -> None:
    path = tmp_path / "features.txt"
    path.write_text("cg1\ncg1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate"):
        load_exported_feature_ids(path)


def test_load_exported_features_rejects_blank_or_padded_ids(tmp_path: Path) -> None:
    blank = tmp_path / "blank.txt"
    blank.write_text("cg1\n\ncg2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="blank"):
        load_exported_feature_ids(blank)
    padded = tmp_path / "padded.txt"
    padded.write_text(" cg1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="whitespace"):
        load_exported_feature_ids(padded)


def test_create_lock_hashes_every_artifact_and_fixed_shape(tmp_path: Path) -> None:
    paths = _artifact_paths(tmp_path)
    lock = _create_lock(paths)
    assert lock.genome_build is GenomeBuild.GRCH37
    assert lock.model_version == "1.0.0"
    assert lock.expected_feature_count == 357340
    assert lock.expected_model_unit_count == 42
    assert lock.model_sha256 != lock.feature_sha256
    assert lock.canonical_feature_list_sha256 != lock.feature_sha256
    verify_marlin_artifact_lock(lock, paths)


def test_create_lock_rejects_non_grch37_build(tmp_path: Path) -> None:
    paths = _artifact_paths(tmp_path)
    with pytest.raises(ValueError, match="GRCh37"):
        _create_lock(paths, build=GenomeBuild.GRCH38)


def test_verify_lock_detects_changed_model(tmp_path: Path) -> None:
    paths = _artifact_paths(tmp_path)
    lock = _create_lock(paths)
    paths.model.write_bytes(b"changed")
    with pytest.raises(ValueError, match="model SHA-256"):
        verify_marlin_artifact_lock(lock, paths)


def test_verify_lock_detects_changed_feature_list(tmp_path: Path) -> None:
    paths = _artifact_paths(tmp_path)
    lock = _create_lock(paths)
    with paths.canonical_feature_list.open("at", encoding="utf-8") as handle:
        handle.write("cg999999999\n")
    with pytest.raises(ValueError, match="feature"):
        verify_marlin_artifact_lock(lock, paths)


def test_missing_artifact_is_rejected(tmp_path: Path) -> None:
    paths = _artifact_paths(tmp_path)
    paths.model.unlink()
    with pytest.raises(ValueError, match="model"):
        _create_lock(paths)

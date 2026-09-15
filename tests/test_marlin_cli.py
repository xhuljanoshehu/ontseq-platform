from __future__ import annotations

import sys
from pathlib import Path

import pytest

from ontseq_platform import entrypoint


def test_global_help_lists_all_marlin_commands(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(sys, "argv", ["ontseq", "--help"])
    entrypoint.main()
    output = capsys.readouterr().out
    assert "marlin-lock" in output
    assert "marlin-features" in output
    assert "marlin-classify" in output
    assert "marlin-validate" in output


def test_marlin_classify_requires_explicit_grch37(monkeypatch: pytest.MonkeyPatch) -> None:
    from ontseq_platform.marlin_cli import _parser

    with pytest.raises(SystemExit):
        _parser().parse_args(["marlin-classify", "--input", "x", "--output", "y"])
    with pytest.raises(SystemExit):
        _parser().parse_args(
            [
                "marlin-classify",
                "--input",
                "x",
                "--output",
                "y",
                "--genome-build",
                "GRCh38",
            ]
        )


def test_marlin_validate_slot_is_discoverable() -> None:
    from ontseq_platform.marlin_cli import _parser

    args = _parser().parse_args(
        [
            "marlin-validate",
            "--manifest",
            "manifest.json",
            "--artifact-lock",
            "lock.json",
            "--runtime-profile",
            "profile.json",
            "--runtime-fixture",
            "runtime-fixture.txt",
            "--model",
            "model.hdf5",
            "--feature-rdata",
            "features.RData",
            "--feature-list",
            "features.txt",
            "--class-annotations",
            "classes.xlsx",
            "--probe-bed",
            "probes.bed.gz",
            "--inference-script",
            "infer.R",
            "--output",
            "report.json",
        ]
    )
    assert args.command == "marlin-validate"
    assert args.runtime_fixture == Path("runtime-fixture.txt")
    assert args.output == Path("report.json")


def _write_features(path: Path) -> tuple[str, ...]:
    ids = tuple(f"cg{i:09d}" for i in range(357340))
    path.write_text("\n".join(ids) + "\n", encoding="utf-8")
    return ids


def _lock_for_feature_list(feature_list: Path):
    from datetime import UTC, datetime

    from ontseq_platform.marlin_contracts import MarlinArtifactLock
    from ontseq_platform.models import GenomeBuild
    from ontseq_platform.reference import sha256_file

    return MarlinArtifactLock(
        lock_id="LOCK",
        model_version="1.0.0",
        model_source_uri="model",
        model_sha256="0" * 64,
        code_source_uri="code",
        code_version_or_commit="commit",
        code_manifest_sha256="1" * 64,
        feature_source_uri="features.RData",
        feature_sha256="2" * 64,
        canonical_feature_list_sha256=sha256_file(feature_list),
        class_annotation_source_uri="classes.xlsx",
        class_annotation_sha256="3" * 64,
        probe_resource_source_uri="probes.bed.gz",
        probe_resource_sha256="4" * 64,
        genome_build=GenomeBuild.GRCH37,
        R_version="4.1.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version_if_used="3.11",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 13, tzinfo=UTC),
    )


def test_marlin_features_writes_locked_vector_and_summary(tmp_path: Path) -> None:
    from ontseq_platform.marlin_cli import main

    feature_list = tmp_path / "features.txt"
    features = _write_features(feature_list)
    lock = _lock_for_feature_list(feature_list)
    lock_path = tmp_path / "lock.json"
    lock_path.write_text(lock.model_dump_json(indent=2), encoding="utf-8")
    source = tmp_path / "sample.bed"
    source.write_text(f"chr1\t1\t2\t0.75\t{features[0]}\n", encoding="utf-8")
    vector = tmp_path / "vector.txt"
    summary = tmp_path / "summary.json"

    main(
        [
            "marlin-features",
            "--input",
            str(source),
            "--genome-build",
            "GRCh37",
            "--artifact-lock",
            str(lock_path),
            "--feature-list",
            str(feature_list),
            "--output-vector",
            str(vector),
            "--output-summary",
            str(summary),
        ]
    )

    assert vector.read_text(encoding="utf-8").splitlines()[0] == "1"
    payload = __import__("json").loads(summary.read_text(encoding="utf-8"))
    assert payload["observed_model_feature_count"] == 1
    assert payload["feature_artifact_sha256"] == lock.canonical_feature_list_sha256


def test_marlin_classify_zero_evidence_returns_no_call_without_r_runtime(tmp_path: Path) -> None:
    from datetime import UTC, datetime

    from openpyxl import Workbook

    from ontseq_platform.marlin_artifacts import (
        MarlinArtifactPaths,
        create_marlin_artifact_lock,
    )
    from ontseq_platform.marlin_cli import main
    from ontseq_platform.marlin_contracts import MarlinRuntimeCompatibilityProfile
    from ontseq_platform.models import GenomeBuild

    feature_list = tmp_path / "features.txt"
    _write_features(feature_list)
    model = tmp_path / "model.hdf5"
    feature_rdata = tmp_path / "features.RData"
    classes = tmp_path / "classes.xlsx"
    probe_bed = tmp_path / "probes.bed.gz"
    model.write_bytes(b"model")
    feature_rdata.write_bytes(b"features")
    probe_bed.write_bytes(b"probes")
    wb = Workbook()
    ws = wb.active
    ws.append(["model_id", "class_name_current", "mcf", "lineage"])
    for i in range(1, 43):
        ws.append([i, f"Class {i}", f"Family {i}", "Lineage"])
    wb.save(classes)
    lock = create_marlin_artifact_lock(
        MarlinArtifactPaths(model, feature_rdata, feature_list, classes, probe_bed),
        lock_id="LOCK",
        code_version_or_commit="commit",
        code_manifest_sha256="1" * 64,
        genome_build=GenomeBuild.GRCH37,
        R_version="4.1.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version_if_used="3.11",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 13, tzinfo=UTC),
    )
    lock_path = tmp_path / "lock.json"
    lock_path.write_text(lock.model_dump_json(indent=2), encoding="utf-8")
    profile = MarlinRuntimeCompatibilityProfile(
        profile_id="PROFILE",
        reference_runtime_lock_id="LOCK",
        feature_vector_sha256="a" * 64,
        reference_scores=[1 / 42] * 42,
        absolute_score_tolerance=1e-7,
        score_sum_tolerance=1e-5,
        top_model_unit_index=0,
        execution_backend="cpu",
        created_at=datetime(2026, 9, 13, tzinfo=UTC),
    )
    profile_path = tmp_path / "profile.json"
    profile_path.write_text(profile.model_dump_json(indent=2), encoding="utf-8")
    source = tmp_path / "sample.bed"
    source.write_text("chr1\t1\t2\t0.75\toff-model\n", encoding="utf-8")
    output = tmp_path / "result.json"

    main(
        [
            "marlin-classify",
            "--sample-id",
            "SAMPLE_001",
            "--input",
            str(source),
            "--genome-build",
            "GRCh37",
            "--artifact-lock",
            str(lock_path),
            "--runtime-profile",
            str(profile_path),
            "--model",
            str(model),
            "--feature-rdata",
            str(feature_rdata),
            "--feature-list",
            str(feature_list),
            "--class-annotations",
            str(classes),
            "--probe-bed",
            str(probe_bed),
            "--inference-script",
            str(tmp_path / "missing.R"),
            "--output",
            str(output),
        ]
    )
    payload = __import__("json").loads(output.read_text(encoding="utf-8"))
    assert payload["status"] == "NO_CALL"
    assert payload["decision"] == "UNKNOWN"
    assert payload["raw_model_scores"] == []

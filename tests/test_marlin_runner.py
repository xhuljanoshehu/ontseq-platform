from __future__ import annotations

import hashlib
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

import pytest
from openpyxl import Workbook

from ontseq_platform.execution import CommandResult
from ontseq_platform.marlin_artifacts import (
    MarlinArtifactPaths,
    create_marlin_artifact_lock,
)
from ontseq_platform.marlin_contracts import MarlinRuntimeCompatibilityProfile
from ontseq_platform.marlin_runner import MarlinRunResources, run_precomputed_marlin_classification
from ontseq_platform.models import GenomeBuild, ModuleRunStatus


def _resources(tmp_path: Path):
    feature_list = tmp_path / "features.txt"
    with feature_list.open("wt", encoding="utf-8") as handle:
        for i in range(357340):
            handle.write(f"cg{i:09d}\n")
    model = tmp_path / "model.hdf5"
    feature_rdata = tmp_path / "features.RData"
    classes = tmp_path / "classes.xlsx"
    probe_bed = tmp_path / "probes.bed.gz"
    inference_script = tmp_path / "infer.R"
    model.write_bytes(b"model")
    feature_rdata.write_bytes(b"features")
    probe_bed.write_bytes(b"probes")
    inference_script.write_text("# not executed for zero-evidence fixture\n", encoding="utf-8")
    wb = Workbook()
    ws = wb.active
    ws.append(["model_id", "class_name_current", "mcf", "lineage"])
    for i in range(1, 43):
        ws.append([i, f"Class {i}", f"Family {i}", "Lineage"])
    wb.save(classes)
    paths = MarlinArtifactPaths(model, feature_rdata, feature_list, classes, probe_bed)
    lock = create_marlin_artifact_lock(
        paths,
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
    return MarlinRunResources(paths, inference_script), lock, profile


def test_zero_evidence_runner_skips_external_runtime(tmp_path: Path) -> None:
    resources, lock, profile = _resources(tmp_path)
    input_path = tmp_path / "sample.bed"
    input_path.write_text("chr1\t1\t2\t0.9\toff-model\n", encoding="utf-8")

    class ExplodingRunner:
        def run(self, *args, **kwargs):
            raise AssertionError("runtime must not execute for zero observed model features")

    report, runtime_result = run_precomputed_marlin_classification(
        sample_id="SAMPLE",
        input_path=input_path,
        genome_build=GenomeBuild.GRCH37,
        lock=lock,
        runtime_profile=profile,
        resources=resources,
        runner=ExplodingRunner(),
        work_dir=tmp_path / "work",
    )
    assert report.status is ModuleRunStatus.NO_CALL
    assert runtime_result is None


def test_observed_runner_rejects_scores_outside_frozen_runtime_profile(tmp_path: Path) -> None:
    resources, lock, _ = _resources(tmp_path)
    input_path = tmp_path / "sample.bed"
    input_path.write_text("chr1\t1\t2\t0.9\tcg000000000\n", encoding="utf-8")
    feature_sha256 = hashlib.sha256(b"\x01" + b"\x00" * 357339).hexdigest()
    profile = MarlinRuntimeCompatibilityProfile(
        profile_id="PROFILE",
        reference_runtime_lock_id="LOCK",
        feature_vector_sha256=feature_sha256,
        reference_scores=[1 / 42] * 42,
        absolute_score_tolerance=1e-7,
        score_sum_tolerance=1e-5,
        top_model_unit_index=0,
        execution_backend="cpu",
        created_at=datetime(2026, 9, 13, tzinfo=UTC),
    )

    class DriftedScoreRunner:
        def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
            scores = [0.6, 0.4] + [0.0] * 40
            payload = ["model_id\tscore"]
            payload.extend(f"{index}\t{score:.17g}" for index, score in enumerate(scores, 1))
            Path(argv[-1]).write_text("\n".join(payload) + "\n", encoding="utf-8")
            return CommandResult(
                argv=tuple(argv), returncode=0, stdout="runtime-log", stderr=""
            )

    with pytest.raises(ValueError, match="frozen|tolerance"):
        run_precomputed_marlin_classification(
            sample_id="SAMPLE",
            input_path=input_path,
            genome_build=GenomeBuild.GRCH37,
            lock=lock,
            runtime_profile=profile,
            resources=resources,
            runner=DriftedScoreRunner(),
            work_dir=tmp_path / "work",
        )

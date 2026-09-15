from __future__ import annotations

import hashlib
from collections import Counter
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

import pytest

from ontseq_platform.execution import CommandResult
from ontseq_platform.marlin_contracts import MarlinArtifactLock
from ontseq_platform.marlin_runtime import MarlinRuntimeProbeReport
from ontseq_platform.marlin_runtime_freeze import (
    MARLIN_RUNTIME_FIXTURE_ABSOLUTE_SCORE_TOLERANCE,
    MARLIN_RUNTIME_FIXTURE_GENERATOR_VERSION,
    MARLIN_RUNTIME_FIXTURE_SCORE_SUM_TOLERANCE,
    freeze_runtime_compatibility,
    generate_frozen_runtime_fixture,
    render_frozen_runtime_fixture,
    verify_runtime_probe_matches_lock,
)
from ontseq_platform.models import GenomeBuild


class FakeRunner:
    def __init__(self, payload: str) -> None:
        self.payload = payload
        self.argv: tuple[str, ...] | None = None

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        self.argv = tuple(argv)
        Path(self.argv[-1]).write_text(self.payload, encoding="utf-8")
        return CommandResult(
            argv=self.argv,
            returncode=0,
            stdout="runtime-log",
            stderr="",
        )


def _payload(scores: Sequence[float]) -> str:
    rows = ["model_id\tscore"]
    rows.extend(f"{index}\t{score:.17g}" for index, score in enumerate(scores, 1))
    return "\n".join(rows) + "\n"


def _lock(model_sha256: str) -> MarlinArtifactLock:
    return MarlinArtifactLock(
        lock_id="MARLIN_V1_FREEZE_TEST",
        model_version="1.0.0",
        model_source_uri="https://doi.org/10.5281/zenodo.15565404",
        model_sha256=model_sha256,
        code_source_uri="https://github.com/hovestadt/MARLIN",
        code_version_or_commit="442aa603415a54f62e7367794f9a31c6bc20fc2d",
        code_manifest_sha256="a" * 64,
        feature_source_uri="marlin_v1.features.RData",
        feature_sha256="b" * 64,
        canonical_feature_list_sha256="f" * 64,
        class_annotation_source_uri="marlin_v1.class_annotations.xlsx",
        class_annotation_sha256="c" * 64,
        probe_resource_source_uri="marlin_v1.probes_hg19.bed.gz",
        probe_resource_sha256="d" * 64,
        genome_build=GenomeBuild.GRCH37,
        R_version="4.2.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version_if_used="3.10.21",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 15, tzinfo=UTC),
    )


def _probe() -> MarlinRuntimeProbeReport:
    return MarlinRuntimeProbeReport(
        R_version="4.2.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version="3.10.21",
        execution_backend="cpu",
    )


def test_runtime_fixture_generator_is_deterministic_and_contract_locked() -> None:
    lock = _lock("0" * 64)

    first = generate_frozen_runtime_fixture(lock)
    second = generate_frozen_runtime_fixture(lock)

    assert MARLIN_RUNTIME_FIXTURE_GENERATOR_VERSION == "sha256-index-mod3-v1"
    assert len(first.values) == 357340
    assert first.values == second.values
    assert first.values[:12] == (
        -1.0,
        1.0,
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
    )
    assert first.summary.feature_vector_sha256 == (
        "b7f7c11d52f777d5e9c45f5d61912d42c21fb2a6dc403a7303d8eb0bf5140cde"
    )
    assert Counter(first.values) == Counter({-1.0: 119915, 0.0: 118324, 1.0: 119101})
    assert first.summary.observed_model_feature_count == 239016
    assert first.summary.absent_feature_count == 118324
    assert first.summary.explicit_na_feature_count == 0
    assert first.summary.feature_artifact_sha256 == lock.canonical_feature_list_sha256


def test_runtime_fixture_text_is_byte_reproducible() -> None:
    vector = generate_frozen_runtime_fixture(_lock("0" * 64))

    first = render_frozen_runtime_fixture(vector)
    second = render_frozen_runtime_fixture(vector)

    assert first == second
    assert first.endswith("\n")
    assert len(first.encode("utf-8")) == 834595
    assert hashlib.sha256(first.encode("utf-8")).hexdigest() == (
        "c41b18446193c5ca93e64e054762356203a020b9a53d9330b31e979d61b27d30"
    )


def test_live_runtime_probe_must_match_artifact_lock() -> None:
    lock = _lock("0" * 64)
    verify_runtime_probe_matches_lock(lock, _probe())

    for field, value in (
        ("R_version", "4.3.0"),
        ("keras_version", "2.14.0"),
        ("tensorflow_version", "2.14.0"),
        ("python_version", "3.10.20"),
        ("execution_backend", "gpu"),
    ):
        with pytest.raises(ValueError, match=field):
            verify_runtime_probe_matches_lock(lock, _probe().model_copy(update={field: value}))


def test_freeze_runtime_compatibility_runs_model_once_and_binds_profile(tmp_path: Path) -> None:
    model = tmp_path / "model.hdf5"
    model.write_bytes(b"model")
    script = tmp_path / "infer.R"
    script.write_text("# synthetic inference fixture\n", encoding="utf-8")
    lock = _lock(hashlib.sha256(b"model").hexdigest())
    scores = [0.6, 0.4] + [0.0] * 40
    runner = FakeRunner(_payload(scores))
    created_at = datetime(2026, 9, 15, 8, 0, tzinfo=UTC)

    vector, result, profile = freeze_runtime_compatibility(
        lock,
        model_path=model,
        inference_script=script,
        runner=runner,
        work_dir=tmp_path / "work",
        profile_id="MARLIN_V1_CPU_FROZEN",
        created_at=created_at,
        rscript_path="Rscript",
    )

    assert runner.argv is not None
    assert runner.argv[0] == "Rscript"
    assert profile.reference_runtime_lock_id == lock.lock_id
    assert profile.feature_vector_sha256 == vector.summary.feature_vector_sha256
    assert [item.score for item in result.model_scores] == scores
    assert profile.reference_scores == scores
    assert profile.top_model_unit_index == 0
    assert profile.execution_backend == "cpu"
    assert profile.created_at == created_at
    assert profile.absolute_score_tolerance == MARLIN_RUNTIME_FIXTURE_ABSOLUTE_SCORE_TOLERANCE
    assert profile.score_sum_tolerance == MARLIN_RUNTIME_FIXTURE_SCORE_SUM_TOLERANCE
    assert MARLIN_RUNTIME_FIXTURE_ABSOLUTE_SCORE_TOLERANCE == 1e-7
    assert MARLIN_RUNTIME_FIXTURE_SCORE_SUM_TOLERANCE == 1e-5

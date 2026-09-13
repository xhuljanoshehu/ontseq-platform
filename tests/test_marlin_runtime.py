from __future__ import annotations

import hashlib
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

import pytest

from ontseq_platform.execution import CommandResult, ToolExecutionError
from ontseq_platform.marlin_contracts import MarlinArtifactLock, MarlinFeatureSummary
from ontseq_platform.marlin_features import MarlinFeatureVector
from ontseq_platform.marlin_runtime import run_marlin_inference
from ontseq_platform.models import GenomeBuild


class FakeRunner:
    def __init__(
        self,
        payload: str | None,
        *,
        returncode: int = 0,
        stderr: str = "",
    ) -> None:
        self.payload = payload
        self.returncode = returncode
        self.stderr = stderr
        self.argv: tuple[str, ...] | None = None

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        self.argv = tuple(argv)
        if self.payload is not None:
            Path(self.argv[-1]).write_text(self.payload, encoding="utf-8")
        return CommandResult(
            argv=self.argv,
            returncode=self.returncode,
            stdout="runtime-log",
            stderr=self.stderr,
        )


def _vector() -> MarlinFeatureVector:
    values = (1.0,) + (0.0,) * 357339
    digest = hashlib.sha256(b"\x01" + b"\x00" * 357339).hexdigest()
    return MarlinFeatureVector(
        values=values,
        summary=MarlinFeatureSummary(
            observed_model_feature_count=1,
            explicit_na_feature_count=0,
            absent_feature_count=357339,
            non_model_probe_count=0,
            observed_fraction=1 / 357340,
            feature_vector_sha256=digest,
            feature_artifact_sha256="f" * 64,
        ),
    )


def _lock(model_sha256: str) -> MarlinArtifactLock:
    return MarlinArtifactLock(
        lock_id="MARLIN_V1_RUNTIME_TEST",
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
        R_version="4.1.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version_if_used="3.11.0",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 13, tzinfo=UTC),
    )


def _labels() -> tuple[str, ...]:
    return tuple(f"unit-{index}" for index in range(1, 43))


def _payload(scores: Sequence[float], ids: Sequence[int] | None = None) -> str:
    ids = tuple(ids) if ids is not None else tuple(range(1, len(scores) + 1))
    rows = ["model_id\tscore"]
    rows.extend(f"{model_id}\t{score:.17g}" for model_id, score in zip(ids, scores, strict=True))
    return "\n".join(rows) + "\n"


def _runtime_paths(tmp_path: Path) -> tuple[Path, Path, MarlinArtifactLock]:
    model = tmp_path / "model.hdf5"
    model.write_bytes(b"model")
    script = tmp_path / "infer.R"
    script.write_text("# fixture\n", encoding="utf-8")
    model_sha = hashlib.sha256(b"model").hexdigest()
    return model, script, _lock(model_sha)


def _run(tmp_path: Path, runner: FakeRunner, scores: Sequence[float] | None = None):
    model, script, lock = _runtime_paths(tmp_path)
    if scores is not None:
        runner.payload = _payload(scores)
    return run_marlin_inference(
        _vector(),
        lock,
        model_path=model,
        class_labels=_labels(),
        inference_script=script,
        runner=runner,
        work_dir=tmp_path / "work",
        rscript_path="Rscript",
        timeout_seconds=30,
    )


def test_runtime_accepts_exact_42_score_softmax(tmp_path: Path) -> None:
    runner = FakeRunner(None)
    result = _run(tmp_path, runner, [1 / 42] * 42)
    assert len(result.model_scores) == 42
    assert result.model_scores[0].label == "unit-1"
    assert result.artifact_lock_id == "MARLIN_V1_RUNTIME_TEST"
    assert runner.argv is not None
    assert runner.argv[0] == "Rscript"


def test_runtime_rejects_41_scores(tmp_path: Path) -> None:
    runner = FakeRunner(_payload([1 / 41] * 41))
    with pytest.raises(ValueError, match="exactly 42"):
        _run(tmp_path, runner)


def test_runtime_rejects_nonfinite_score(tmp_path: Path) -> None:
    scores = [1 / 42] * 42
    payload = _payload(scores).replace(f"1\t{scores[0]:.17g}", "1\tnan", 1)
    with pytest.raises(ValueError, match="finite"):
        _run(tmp_path, FakeRunner(payload))


def test_runtime_rejects_softmax_sum_outside_tolerance(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="softmax"):
        _run(tmp_path, FakeRunner(_payload([0.01] * 42)))


def test_runtime_rejects_duplicate_or_out_of_order_model_ids(tmp_path: Path) -> None:
    scores = [1 / 42] * 42
    ids = list(range(1, 43))
    ids[1] = 1
    with pytest.raises(ValueError, match="model IDs"):
        _run(tmp_path, FakeRunner(_payload(scores, ids)))
    ids = list(range(1, 43))
    ids[0], ids[1] = ids[1], ids[0]
    with pytest.raises(ValueError, match="model IDs"):
        _run(tmp_path, FakeRunner(_payload(scores, ids)))


def test_runtime_rejects_invalid_class_labels(tmp_path: Path) -> None:
    model, script, lock = _runtime_paths(tmp_path)
    runner = FakeRunner(_payload([1 / 42] * 42))
    with pytest.raises(ValueError, match="42 unique"):
        run_marlin_inference(
            _vector(),
            lock,
            model_path=model,
            class_labels=("same",) * 42,
            inference_script=script,
            runner=runner,
            work_dir=tmp_path / "work",
        )


def test_runtime_rejects_nonzero_process_exit(tmp_path: Path) -> None:
    with pytest.raises(ToolExecutionError, match="failed"):
        _run(tmp_path, FakeRunner(None, returncode=3, stderr="keras failed"))


def test_runtime_rejects_missing_output_file(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="no score output"):
        _run(tmp_path, FakeRunner(None))


def test_runtime_rejects_model_hash_mismatch(tmp_path: Path) -> None:
    model, script, _ = _runtime_paths(tmp_path)
    with pytest.raises(ValueError, match="model SHA-256"):
        run_marlin_inference(
            _vector(),
            _lock("0" * 64),
            model_path=model,
            class_labels=_labels(),
            inference_script=script,
            runner=FakeRunner(_payload([1 / 42] * 42)),
            work_dir=tmp_path / "work",
        )


def test_runtime_rejects_feature_digest_mismatch(tmp_path: Path) -> None:
    vector = _vector()
    vector = MarlinFeatureVector(
        values=vector.values,
        summary=vector.summary.model_copy(update={"feature_vector_sha256": "0" * 64}),
    )
    model, script, lock = _runtime_paths(tmp_path)
    with pytest.raises(ValueError, match="feature-vector SHA-256"):
        run_marlin_inference(
            vector,
            lock,
            model_path=model,
            class_labels=_labels(),
            inference_script=script,
            runner=FakeRunner(_payload([1 / 42] * 42)),
            work_dir=tmp_path / "work",
        )

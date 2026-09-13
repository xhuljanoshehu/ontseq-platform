from __future__ import annotations

import hashlib
import math
import shutil
import struct
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .execution import CommandRunner, ToolExecutionError
from .marlin_contracts import MarlinArtifactLock, MarlinModelUnitScore
from .marlin_features import MarlinFeatureVector, marlin_feature_vector_sha256
from .models import StrictModel
from .reference import sha256_file

_EXPECTED_FEATURE_COUNT = 357340
_EXPECTED_MODEL_UNITS = 42
_SOFTMAX_SUM_TOLERANCE = 1e-5


class MarlinRuntimeResult(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    artifact_lock_id: str = Field(min_length=1)
    feature_vector_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    model_scores: list[MarlinModelUnitScore]
    raw_score_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    execution_backend: str = Field(min_length=1)

    @model_validator(mode="after")
    def result_has_exact_softmax_shape(self) -> MarlinRuntimeResult:
        if len(self.model_scores) != _EXPECTED_MODEL_UNITS:
            raise ValueError("MARLIN runtime result requires exactly 42 model scores")
        ids = [item.model_id for item in self.model_scores]
        labels = [item.label for item in self.model_scores]
        if ids != list(range(1, 43)) or len(set(labels)) != 42:
            raise ValueError("MARLIN runtime model IDs must be ordered 1..42 with unique labels")
        if abs(sum(item.score for item in self.model_scores) - 1.0) > _SOFTMAX_SUM_TOLERANCE:
            raise ValueError("MARLIN runtime scores violate softmax-sum tolerance")
        return self


def _require_file(path: Path, *, label: str) -> Path:
    candidate = Path(path)
    try:
        resolved = candidate.resolve(strict=True)
    except FileNotFoundError as exc:
        raise ValueError(f"MARLIN {label} file is missing: {candidate}") from exc
    if not resolved.is_file():
        raise ValueError(f"MARLIN {label} path is not a file: {candidate}")
    return resolved


def _validate_class_labels(class_labels: Sequence[str]) -> tuple[str, ...]:
    labels = tuple(class_labels)
    if len(labels) != 42 or len(set(labels)) != 42:
        raise ValueError("MARLIN runtime requires exactly 42 unique class labels")
    if any(not label or label != label.strip() for label in labels):
        raise ValueError("MARLIN class labels must be non-empty and whitespace-normalized")
    return labels


def _score_digest(scores: Sequence[MarlinModelUnitScore]) -> str:
    payload = b"".join(struct.pack(">d", item.score) for item in scores)
    return hashlib.sha256(payload).hexdigest()


def _parse_score_output(path: Path, labels: Sequence[str]) -> list[MarlinModelUnitScore]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "model_id\tscore":
        raise ValueError("MARLIN runtime score output has an invalid header")
    rows = lines[1:]
    if len(rows) != _EXPECTED_MODEL_UNITS:
        raise ValueError(f"MARLIN runtime must emit exactly 42 scores; found {len(rows)}")

    scores: list[MarlinModelUnitScore] = []
    for expected_id, row in enumerate(rows, 1):
        fields = row.split("\t")
        if len(fields) != 2:
            raise ValueError(
                f"MARLIN runtime row {expected_id} must contain 2 tab-separated fields"
            )
        try:
            model_id = int(fields[0])
        except ValueError as exc:
            raise ValueError(f"MARLIN runtime row {expected_id} has invalid model ID") from exc
        if model_id != expected_id:
            raise ValueError("MARLIN runtime model IDs must be in locked order 1..42")
        try:
            score = float(fields[1])
        except ValueError as exc:
            raise ValueError(f"MARLIN runtime row {expected_id} has invalid score") from exc
        if not math.isfinite(score):
            raise ValueError(f"MARLIN runtime row {expected_id} score must be finite")
        if not 0 <= score <= 1:
            raise ValueError(f"MARLIN runtime row {expected_id} score must be between 0 and 1")
        scores.append(
            MarlinModelUnitScore(
                model_id=model_id, label=labels[model_id - 1], score=score
            )
        )

    if abs(sum(item.score for item in scores) - 1.0) > _SOFTMAX_SUM_TOLERANCE:
        raise ValueError("MARLIN runtime scores violate softmax-sum tolerance")
    return scores


def run_marlin_inference(
    feature_vector: MarlinFeatureVector,
    lock: MarlinArtifactLock,
    *,
    model_path: Path,
    class_labels: Sequence[str],
    inference_script: Path,
    runner: CommandRunner,
    work_dir: Path,
    rscript_path: str = "Rscript",
    timeout_seconds: int = 300,
) -> MarlinRuntimeResult:
    """Run only the locked neural-network inference over an ONTSeq-built feature vector."""
    labels = _validate_class_labels(class_labels)
    if len(feature_vector.values) != _EXPECTED_FEATURE_COUNT:
        raise ValueError("MARLIN runtime requires exactly 357340 feature values")
    observed_vector_sha = marlin_feature_vector_sha256(feature_vector.values)
    if observed_vector_sha != feature_vector.summary.feature_vector_sha256:
        raise ValueError("MARLIN feature-vector SHA-256 does not match the feature summary")
    if feature_vector.summary.feature_artifact_sha256 != lock.canonical_feature_list_sha256:
        raise ValueError("MARLIN feature vector was not built from the locked feature artifact")

    model = _require_file(model_path, label="model")
    script = _require_file(inference_script, label="inference script")
    model_sha_before = sha256_file(model)
    if model_sha_before != lock.model_sha256:
        raise ValueError("MARLIN model SHA-256 does not match the artifact lock")

    root = Path(work_dir)
    root.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".marlin-runtime-", dir=root))
    vector_path = stage / "feature-vector.txt"
    score_path = stage / "scores.tsv"
    try:
        vector_path.write_text(
            "".join(f"{int(value)}\n" for value in feature_vector.values),
            encoding="utf-8",
            newline="\n",
        )
        argv = (
            rscript_path,
            str(script),
            str(vector_path),
            str(model),
            str(score_path),
        )
        result = runner.run(argv, timeout_seconds=timeout_seconds)
        if result.returncode != 0:
            detail = result.stderr.strip() or "no stderr"
            raise ToolExecutionError(
                f"MARLIN inference failed with exit code {result.returncode}: {detail}"
            )
        if not score_path.is_file():
            raise ValueError("MARLIN inference succeeded but produced no score output file")
        scores = _parse_score_output(score_path, labels)
        model_sha_after = sha256_file(model)
        if model_sha_after != model_sha_before:
            raise ValueError("MARLIN model changed while inference was running")
        return MarlinRuntimeResult(
            artifact_lock_id=lock.lock_id,
            feature_vector_sha256=observed_vector_sha,
            model_scores=scores,
            raw_score_sha256=_score_digest(scores),
            execution_backend=lock.execution_backend,
        )
    finally:
        shutil.rmtree(stage, ignore_errors=True)

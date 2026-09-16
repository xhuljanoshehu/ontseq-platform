from __future__ import annotations

import hashlib
import math
import shutil
import struct
import tempfile
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .execution import CommandRunner, ToolExecutionError
from .marlin_contracts import (
    MarlinArtifactLock,
    MarlinFeatureSummary,
    MarlinModelUnitScore,
    MarlinRuntimeCompatibilityProfile,
)
from .marlin_features import MarlinFeatureVector, marlin_feature_vector_sha256
from .models import StrictModel
from .reference import sha256_file

_EXPECTED_FEATURE_COUNT = 357340
_EXPECTED_MODEL_UNITS = 42
_SOFTMAX_SUM_TOLERANCE = 1e-5
_RUNTIME_PROBE_FIELDS = frozenset(
    {
        "R_version",
        "keras_version",
        "tensorflow_version",
        "python_version",
        "execution_backend",
    }
)


class MarlinRuntimeProbeReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    R_version: str = Field(min_length=1)
    keras_version: str = Field(min_length=1)
    tensorflow_version: str = Field(min_length=1)
    python_version: str = Field(pattern=r"^\d+\.\d+\.\d+$")
    execution_backend: Literal["cpu", "gpu"]
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def runtime_identity_is_normalized(self) -> MarlinRuntimeProbeReport:
        version_fields = (
            self.R_version,
            self.keras_version,
            self.tensorflow_version,
            self.python_version,
        )
        if any(value != value.strip() for value in version_fields):
            raise ValueError("MARLIN runtime probe version fields must be whitespace-normalized")
        return self


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


def _parse_runtime_probe_output(path: Path) -> MarlinRuntimeProbeReport:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "key\tvalue":
        raise ValueError("MARLIN runtime probe output has an invalid header")

    observed: dict[str, str] = {}
    for line_number, row in enumerate(lines[1:], 2):
        fields = row.split("\t")
        if len(fields) != 2 or not fields[0] or not fields[1]:
            raise ValueError(
                f"MARLIN runtime probe row {line_number} must contain non-empty key/value fields"
            )
        key, value = fields
        if key in observed:
            raise ValueError(f"MARLIN runtime probe contains duplicate field {key!r}")
        observed[key] = value

    if set(observed) != _RUNTIME_PROBE_FIELDS:
        missing = sorted(_RUNTIME_PROBE_FIELDS - set(observed))
        extra = sorted(set(observed) - _RUNTIME_PROBE_FIELDS)
        raise ValueError(
            "MARLIN runtime probe fields must match the locked contract; "
            f"missing={missing}, extra={extra}"
        )
    return MarlinRuntimeProbeReport.model_validate(observed)


def probe_marlin_runtime(
    *,
    probe_script: Path,
    runner: CommandRunner,
    rscript_path: str = "Rscript",
    timeout_seconds: int = 120,
) -> MarlinRuntimeProbeReport:
    """Execute the live R/Keras/TensorFlow probe and normalize its runtime identity."""
    script = _require_file(probe_script, label="runtime probe script")
    with tempfile.TemporaryDirectory(prefix="ontseq-marlin-runtime-probe-") as directory:
        output = Path(directory) / "runtime-probe.tsv"
        result = runner.run(
            (rscript_path, str(script), str(output)),
            timeout_seconds=timeout_seconds,
        )
        if result.returncode != 0:
            detail = result.stderr.strip() or "no stderr"
            raise ToolExecutionError(
                f"MARLIN runtime probe failed with exit code {result.returncode}: {detail}"
            )
        if not output.is_file():
            raise ValueError("MARLIN runtime probe succeeded but produced no output file")
        return _parse_runtime_probe_output(output)


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
                model_id=model_id,
                label=labels[model_id - 1],
                score=score,
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


def create_runtime_compatibility_profile(
    reference_result: MarlinRuntimeResult,
    *,
    profile_id: str,
    absolute_score_tolerance: float,
    score_sum_tolerance: float,
    created_at: datetime,
) -> MarlinRuntimeCompatibilityProfile:
    """Freeze a runtime's fixed-vector output before biological validation.

    Tolerances are caller-supplied experimental decisions. This function deliberately does
    not infer them from a biological cohort or from the observed score differences.
    """
    if not math.isfinite(absolute_score_tolerance) or absolute_score_tolerance <= 0:
        raise ValueError("absolute_score_tolerance must be a positive finite value")
    if not math.isfinite(score_sum_tolerance) or score_sum_tolerance <= 0:
        raise ValueError("score_sum_tolerance must be a positive finite value")
    scores = [item.score for item in reference_result.model_scores]
    top_index = max(range(_EXPECTED_MODEL_UNITS), key=scores.__getitem__)
    return MarlinRuntimeCompatibilityProfile(
        profile_id=profile_id,
        reference_runtime_lock_id=reference_result.artifact_lock_id,
        feature_vector_sha256=reference_result.feature_vector_sha256,
        reference_scores=scores,
        absolute_score_tolerance=absolute_score_tolerance,
        score_sum_tolerance=score_sum_tolerance,
        top_model_unit_index=top_index,
        execution_backend=reference_result.execution_backend,
        created_at=created_at,
    )


def verify_runtime_compatibility(
    profile: MarlinRuntimeCompatibilityProfile,
    candidate: MarlinRuntimeResult,
) -> None:
    """Refuse a runtime result that drifts outside a previously frozen profile."""
    if candidate.artifact_lock_id != profile.reference_runtime_lock_id:
        raise ValueError("MARLIN candidate runtime lock differs from compatibility profile")
    if candidate.execution_backend != profile.execution_backend:
        raise ValueError("MARLIN candidate execution backend differs from compatibility profile")
    if candidate.feature_vector_sha256 != profile.feature_vector_sha256:
        raise ValueError("MARLIN candidate feature vector differs from compatibility profile")

    candidate_scores = [item.score for item in candidate.model_scores]
    candidate_top = max(range(_EXPECTED_MODEL_UNITS), key=candidate_scores.__getitem__)
    if candidate_top != profile.top_model_unit_index:
        raise ValueError("MARLIN candidate changed the top model unit")
    if abs(sum(candidate_scores) - 1.0) > profile.score_sum_tolerance:
        raise ValueError("MARLIN candidate violates frozen score-sum tolerance")

    for model_id, (reference_score, candidate_score) in enumerate(
        zip(profile.reference_scores, candidate_scores, strict=True), 1
    ):
        if abs(candidate_score - reference_score) > profile.absolute_score_tolerance:
            raise ValueError(
                f"MARLIN candidate score exceeds frozen absolute tolerance at model unit {model_id}"
            )


def load_frozen_runtime_fixture(
    fixture_path: Path,
    lock: MarlinArtifactLock,
) -> MarlinFeatureVector:
    """Load the canonical {-1,0,+1} runtime fixture without biological preprocessing."""
    fixture = _require_file(fixture_path, label="runtime fixture")
    before_sha256 = sha256_file(fixture)
    lines = fixture.read_text(encoding="utf-8").splitlines()
    if len(lines) != _EXPECTED_FEATURE_COUNT:
        raise ValueError("MARLIN runtime fixture requires exactly 357340 feature values")

    values: list[float] = []
    for line_number, text in enumerate(lines, 1):
        if text not in {"-1", "0", "1"}:
            raise ValueError(
                f"MARLIN runtime fixture line {line_number} must be exactly -1, 0 or 1"
            )
        values.append(float(text))

    after_sha256 = sha256_file(fixture)
    if after_sha256 != before_sha256:
        raise ValueError("MARLIN runtime fixture changed while it was being loaded")

    value_tuple = tuple(values)
    observed_count = sum(value != 0 for value in value_tuple)
    absent_count = _EXPECTED_FEATURE_COUNT - observed_count
    vector_sha256 = marlin_feature_vector_sha256(value_tuple)
    summary = MarlinFeatureSummary(
        observed_model_feature_count=observed_count,
        explicit_na_feature_count=0,
        absent_feature_count=absent_count,
        non_model_probe_count=0,
        observed_fraction=observed_count / _EXPECTED_FEATURE_COUNT,
        feature_vector_sha256=vector_sha256,
        feature_artifact_sha256=lock.canonical_feature_list_sha256,
    )
    return MarlinFeatureVector(values=value_tuple, summary=summary)


def verify_frozen_runtime_fixture(
    profile: MarlinRuntimeCompatibilityProfile,
    lock: MarlinArtifactLock,
    *,
    fixture_path: Path,
    model_path: Path,
    inference_script: Path,
    runner: CommandRunner,
    work_dir: Path,
    rscript_path: str = "Rscript",
    timeout_seconds: int = 300,
) -> MarlinRuntimeResult:
    """Re-execute the frozen non-biological fixture before external validation."""
    if profile.reference_runtime_lock_id != lock.lock_id:
        raise ValueError("MARLIN runtime profile belongs to a different artifact/runtime lock")
    if profile.execution_backend != lock.execution_backend:
        raise ValueError("MARLIN runtime profile execution backend differs from artifact lock")

    vector = load_frozen_runtime_fixture(fixture_path, lock)
    if vector.summary.feature_vector_sha256 != profile.feature_vector_sha256:
        raise ValueError(
            "MARLIN runtime fixture feature vector differs from frozen compatibility profile"
        )

    result = run_marlin_inference(
        vector,
        lock,
        model_path=model_path,
        class_labels=tuple(f"model-unit-{index}" for index in range(1, 43)),
        inference_script=inference_script,
        runner=runner,
        work_dir=work_dir,
        rscript_path=rscript_path,
        timeout_seconds=timeout_seconds,
    )
    verify_runtime_compatibility(profile, result)
    return result

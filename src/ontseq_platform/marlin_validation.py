from __future__ import annotations

import math
import os
import tempfile
import time
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .marlin_artifacts import MarlinArtifactPaths
from .marlin_contracts import (
    MarlinArtifactLock,
    MarlinClassificationDecision,
    MarlinPredictionReport,
    MarlinRuntimeCompatibilityProfile,
    MarlinSourceKind,
)
from .execution import SubprocessRunner
from .io import load_model
from .marlin_runner import MarlinRunResources, run_precomputed_marlin_classification
from .models import GenomeBuild, ModuleRunStatus, StrictModel
from .reference import sha256_file

_SHA256 = r"^[0-9a-f]{64}$"


class MarlinValidationSample(StrictModel):
    accession: str = Field(min_length=3)
    sample_id: str = Field(min_length=3)
    public_source_uri: str = Field(min_length=1)
    local_input_path: Path
    input_sha256: str = Field(pattern=_SHA256)
    genome_build: GenomeBuild
    expected_class: str | None = Field(default=None, min_length=1)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def v1_is_grch37(self) -> MarlinValidationSample:
        if self.genome_build is not GenomeBuild.GRCH37:
            raise ValueError("MARLIN v1 external validation currently requires GRCh37/hg19")
        return self


class MarlinValidationManifest(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    validation_id: str = Field(min_length=3)
    artifact_lock_id: str = Field(min_length=1)
    runtime_profile_id: str = Field(min_length=1)
    threshold_policy_id: Literal["MARLIN_V1_PUBLISHED_0.8"] = "MARLIN_V1_PUBLISHED_0.8"
    confidence_threshold: Literal[0.8] = 0.8
    samples: list[MarlinValidationSample] = Field(min_length=1)
    locked_at: datetime
    posthoc_tuning_allowed: Literal[False] = False
    research_only: Literal[True] = True
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def manifest_is_locked_and_unique(self) -> MarlinValidationManifest:
        if self.locked_at.utcoffset() is None:
            raise ValueError("MARLIN validation manifest timestamp requires a timezone")
        accessions = [item.accession for item in self.samples]
        sample_ids = [item.sample_id for item in self.samples]
        if len(accessions) != len(set(accessions)):
            raise ValueError("MARLIN validation accessions must be unique")
        if len(sample_ids) != len(set(sample_ids)):
            raise ValueError("MARLIN validation sample IDs must be unique")
        return self


class MarlinValidationSampleResult(StrictModel):
    accession: str
    sample_id: str
    status: Literal["COMPLETED", "NO_CALL", "FAILED"]
    decision: MarlinClassificationDecision | None = None
    top_class: str | None = None
    top_class_score: float | None = Field(default=None, ge=0, le=1.00001)
    expected_class: str | None = None
    concordant: bool | None = None
    feature_vector_sha256: str | None = Field(default=None, pattern=_SHA256)
    raw_score_sha256: str | None = Field(default=None, pattern=_SHA256)
    observed_model_feature_count: int | None = Field(default=None, ge=0, le=357340)
    observed_fraction: float | None = Field(default=None, ge=0, le=1)
    runtime_seconds: float = Field(ge=0)
    reason: str = ""
    research_only: Literal[True] = True


class MarlinValidationCohortReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    validation_id: str
    artifact_lock_id: str
    runtime_profile_id: str
    confidence_threshold: Literal[0.8] = 0.8
    sample_count: int = Field(ge=0)
    high_confidence_count: int = Field(ge=0)
    unknown_count: int = Field(ge=0)
    no_call_count: int = Field(ge=0)
    failed_count: int = Field(ge=0)
    concordant_count: int = Field(ge=0)
    discordant_count: int = Field(ge=0)
    deterministic_rerun_failures: int = Field(ge=0)
    results: list[MarlinValidationSampleResult]
    generated_at: datetime
    research_only: Literal[True] = True
    sensitive_output: Literal[True] = True

    @model_validator(mode="after")
    def counts_match_results(self) -> MarlinValidationCohortReport:
        if self.generated_at.utcoffset() is None:
            raise ValueError("MARLIN validation report timestamp requires a timezone")
        if self.sample_count != len(self.results):
            raise ValueError("MARLIN validation sample_count differs from results")
        outcome_count = (
            self.high_confidence_count
            + self.unknown_count
            + self.no_call_count
            + self.failed_count
        )
        if outcome_count != self.sample_count:
            raise ValueError("MARLIN validation status counts do not partition the cohort")
        comparable = sum(
            item.expected_class is not None and item.top_class is not None
            for item in self.results
        )
        if self.concordant_count + self.discordant_count != comparable:
            raise ValueError("MARLIN validation concordance counts differ from comparable results")
        return self


ClassificationCallback = Callable[
    [MarlinValidationSample], tuple[MarlinPredictionReport, str | None, float]
]


def _validate_registered_inputs(manifest: MarlinValidationManifest) -> None:
    for sample in manifest.samples:
        path = sample.local_input_path
        if not path.is_file():
            raise ValueError(f"MARLIN validation input is missing for {sample.accession}")
        observed = sha256_file(path)
        if observed != sample.input_sha256:
            raise ValueError(
                f"MARLIN validation input SHA-256 mismatch for {sample.accession}"
            )


def _validate_report_identity(
    sample: MarlinValidationSample,
    manifest: MarlinValidationManifest,
    report: MarlinPredictionReport,
) -> None:
    if report.sample_id != sample.sample_id:
        raise ValueError("MARLIN validation report sample ID differs from manifest")
    if report.genome_build is not sample.genome_build:
        raise ValueError("MARLIN validation report genome build differs from manifest")
    if report.source_kind is not MarlinSourceKind.PRECOMPUTED_METHYLATION:
        raise ValueError("MARLIN external validation requires PRECOMPUTED_METHYLATION")
    if report.artifact_lock_id != manifest.artifact_lock_id:
        raise ValueError("MARLIN validation report artifact lock differs from manifest")
    if report.input_fingerprint.sha256 != sample.input_sha256:
        raise ValueError("MARLIN validation report input SHA-256 differs from manifest")
    if report.confidence_threshold != manifest.confidence_threshold:
        raise ValueError("MARLIN validation report threshold differs from manifest")


def _reports_are_deterministic(
    first: MarlinPredictionReport,
    second: MarlinPredictionReport,
    first_raw_digest: str | None,
    second_raw_digest: str | None,
    *,
    absolute_score_tolerance: float,
) -> bool:
    if (
        first.feature_summary.feature_vector_sha256
        != second.feature_summary.feature_vector_sha256
        or first.status is not second.status
        or first.decision is not second.decision
        or first.top_class != second.top_class
    ):
        return False
    if first.top_class_score is None or second.top_class_score is None:
        if first.top_class_score != second.top_class_score:
            return False
    elif abs(first.top_class_score - second.top_class_score) > absolute_score_tolerance:
        return False
    if len(first.raw_model_scores) != len(second.raw_model_scores):
        return False
    if any(
        left.model_id != right.model_id
        or left.label != right.label
        or abs(left.score - right.score) > absolute_score_tolerance
        for left, right in zip(first.raw_model_scores, second.raw_model_scores, strict=True)
    ):
        return False
    if absolute_score_tolerance == 0 and first_raw_digest != second_raw_digest:
        return False
    return True


def _sample_result(
    sample: MarlinValidationSample,
    report: MarlinPredictionReport,
    raw_score_sha256: str | None,
    runtime_seconds: float,
) -> MarlinValidationSampleResult:
    concordant = (
        None
        if sample.expected_class is None or report.top_class is None
        else report.top_class == sample.expected_class
     )
    status = "NO_CALL" if report.status is ModuleRunStatus.NO_CALL else "COMPLETED"
    return MarlinValidationSampleResult(
        accession=sample.accession,
        sample_id=sample.sample_id,
        status=status,
        decision=report.decision,
        top_class=report.top_class,
        top_class_score=report.top_class_score,
        expected_class=sample.expected_class,
        concordant=concordant,
        feature_vector_sha256=report.feature_summary.feature_vector_sha256,
        raw_score_sha256=raw_score_sha256,
        observed_model_feature_count=report.feature_summary.observed_model_feature_count,
        observed_fraction=report.feature_summary.observed_fraction,
        runtime_seconds=runtime_seconds,
    )


def execute_marlin_validation(
    manifest: MarlinValidationManifest,
    *,
    classify_sample: ClassificationCallback,
    absolute_score_tolerance: float = 0.0,
) -> MarlinValidationCohortReport:
    """Run every checksummed sample twice and aggregate only deterministic outcomes."""
    if not math.isfinite(absolute_score_tolerance) or absolute_score_tolerance < 0:
        raise ValueError("MARLIN validation score tolerance must be finite and non-negative")
    _validate_registered_inputs(manifest)
    results: list[MarlinValidationSampleResult] = []
    deterministic_failures = 0

    for sample in manifest.samples:
        try:
            first, first_digest, first_seconds = classify_sample(sample)
            if sha256_file(sample.local_input_path) != sample.input_sha256:
                raise ValueError("input changed after first classification")
            second, second_digest, second_seconds = classify_sample(sample)
            if sha256_file(sample.local_input_path) != sample.input_sha256:
                raise ValueError("input changed after repeated classification")
            _validate_report_identity(sample, manifest, first)
            _validate_report_identity(sample, manifest, second)
            if not _reports_are_deterministic(
                first,
                second,
                first_digest,
                second_digest,
                absolute_score_tolerance=absolute_score_tolerance,
            ):
                deterministic_failures += 1
                results.append(
                    MarlinValidationSampleResult(
                        accession=sample.accession,
                        sample_id=sample.sample_id,
                        status="FAILED",
                        expected_class=sample.expected_class,
                        runtime_seconds=first_seconds + second_seconds,
                        reason=(
                            "Repeated classification was not deterministic within the frozen "
                            "tolerance."
                        ),
                    )
                )
                continue
            results.append(
                _sample_result(
                    sample,
                    first,
                    first_digest,
                    (first_seconds + second_seconds) / 2,
                )
            )
        except Exception as exc:
            results.append(
                MarlinValidationSampleResult(
                    accession=sample.accession,
                    sample_id=sample.sample_id,
                    status="FAILED",
                    expected_class=sample.expected_class,
                    runtime_seconds=0,
                    reason=f"Classification execution failed: {type(exc).__name__}",
                 )
            )

    high_confidence = sum(
        item.status == "COMPLETED"
        and item.decision is MarlinClassificationDecision.HIGH_CONFIDENCE
        for item in results
    )
    unknown = sum(
        item.status == "COMPLETED"
        and item.decision is MarlinClassificationDecision.UNKNOWN
        for item in results
    )
    no_call = sum(item.status == "NO_CALL" for item in results)
    failed = sum(item.status == "FAILED" for item in results)
    concordant = sum(item.concordant is True for item in results)
    discordant = sum(item.concordant is False for item in results)
    return MarlinValidationCohortReport(
        validation_id=manifest.validation_id,
        artifact_lock_id=manifest.artifact_lock_id,
        runtime_profile_id=manifest.runtime_profile_id,
        sample_count=len(results),
        high_confidence_count=high_confidence,
        unknown_count=unknown,
        no_call_count=no_call,
        failed_count=failed,
        concordant_count=concordant,
        discordant_count=discordant,
        deterministic_rerun_failures=deterministic_failures,
        results=results,
        generated_at=datetime.now(UTC),
    )


def _write_report_atomic(report: MarlinValidationCohortReport, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, staged_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    staged = Path(staged_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(report.model_dump_json(indent=2) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(staged, path)
    except BaseException:
        staged.unlink(missing_ok=True)
        raise
    return path


def run_validation_manifest(
    manifest_path: Path,
    *,
    artifact_lock_path: Path,
    runtime_profile_path: Path,
    artifact_paths: MarlinArtifactPaths,
    inference_script: Path,
    output_path: Path,
    rscript_path: str = "Rscript",
    work_dir: Path | None = None,
) -> Path:
    """Execute a local locked validation manifest through the canonical sample runner."""
    manifest = load_model(manifest_path, MarlinValidationManifest)
    lock = load_model(artifact_lock_path, MarlinArtifactLock)
    profile = load_model(runtime_profile_path, MarlinRuntimeCompatibilityProfile)
    if manifest.artifact_lock_id != lock.lock_id:
        raise ValueError("MARLIN validation manifest artifact lock differs from supplied lock")
    if manifest.runtime_profile_id != profile.profile_id:
        raise ValueError("MARLIN validation manifest runtime profile differs from supplied profile")
    if profile.reference_runtime_lock_id != lock.lock_id:
        raise ValueError("MARLIN runtime profile belongs to a different artifact/runtime lock")
    if profile.execution_backend != lock.execution_backend:
        raise ValueError("MARLIN runtime profile backend differs from artifact lock")

    resources = MarlinRunResources(artifact_paths, inference_script)
    base_work_dir = work_dir or (Path(output_path).parent / ".marlin-validation-runtime")
    runner = SubprocessRunner()

    def classify_sample(
        sample: MarlinValidationSample,
    ) -> tuple[MarlinPredictionReport, str | None, float]:
        started = time.perf_counter()
        report, runtime_result = run_precomputed_marlin_classification(
            sample_id=sample.sample_id,
            input_path=sample.local_input_path,
            genome_build=sample.genome_build,
            lock=lock,
            runtime_profile=profile,
            resources=resources,
            runner=runner,
            work_dir=base_work_dir / sample.sample_id,
            rscript_path=rscript_path,
        )
        elapsed = time.perf_counter() - started
        raw_digest = None if runtime_result is None else runtime_result.raw_score_sha256
        return report, raw_digest, elapsed

    report = execute_marlin_validation(
        manifest,
        classify_sample=classify_sample,
        absolute_score_tolerance=profile.absolute_score_tolerance,
    )
    return _write_report_atomic(report, output_path)

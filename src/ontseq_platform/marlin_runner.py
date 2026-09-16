from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .execution import CommandRunner
from .marlin_artifacts import (
    MarlinArtifactPaths,
    load_exported_feature_ids,
    verify_marlin_artifact_lock,
)
from .marlin_classification import classify_marlin_scores, load_marlin_class_annotations
from .marlin_contracts import (
    MarlinArtifactLock,
    MarlinPredictionReport,
    MarlinRuntimeCompatibilityProfile,
)
from .marlin_features import build_marlin_feature_vector
from .marlin_input import parse_marlin_probe_bed
from .marlin_runtime import MarlinRuntimeResult, run_marlin_inference
from .marlin_runtime_freeze import verify_marlin_v1_runtime_profile_policy
from .models import GenomeBuild


@dataclass(frozen=True)
class MarlinRunResources:
    artifact_paths: MarlinArtifactPaths
    inference_script: Path


def verify_runtime_profile_identity(
    profile: MarlinRuntimeCompatibilityProfile,
    lock: MarlinArtifactLock,
) -> None:
    if profile.reference_runtime_lock_id != lock.lock_id:
        raise ValueError("MARLIN runtime profile belongs to a different artifact/runtime lock")
    if profile.execution_backend != lock.execution_backend:
        raise ValueError("MARLIN runtime profile execution backend differs from artifact lock")
    verify_marlin_v1_runtime_profile_policy(profile)


def run_precomputed_marlin_classification(
    *,
    sample_id: str,
    input_path: Path,
    genome_build: GenomeBuild,
    lock: MarlinArtifactLock,
    runtime_profile: MarlinRuntimeCompatibilityProfile,
    resources: MarlinRunResources,
    runner: CommandRunner,
    work_dir: Path,
    rscript_path: str = "Rscript",
) -> tuple[MarlinPredictionReport, MarlinRuntimeResult | None]:
    """Canonical precomputed-CpG MARLIN path shared by CLI and validation."""
    if genome_build is not GenomeBuild.GRCH37:
        raise ValueError("MARLIN v1 runner currently requires GRCh37/hg19")
    verify_runtime_profile_identity(runtime_profile, lock)
    verify_marlin_artifact_lock(lock, resources.artifact_paths)

    feature_ids = load_exported_feature_ids(resources.artifact_paths.canonical_feature_list)
    source = parse_marlin_probe_bed(input_path, genome_build=genome_build)
    vector = build_marlin_feature_vector(
        source,
        feature_ids,
        feature_artifact_sha256=lock.canonical_feature_list_sha256,
    )
    annotations = load_marlin_class_annotations(resources.artifact_paths.class_annotations)

    runtime_result: MarlinRuntimeResult | None = None
    if vector.summary.observed_model_feature_count > 0:
        runtime_result = run_marlin_inference(
            vector,
            lock,
            model_path=resources.artifact_paths.model,
            class_labels=tuple(f"model-unit-{index}" for index in range(1, 43)),
            inference_script=resources.inference_script,
            runner=runner,
            work_dir=work_dir,
            rscript_path=rscript_path,
        )

    report = classify_marlin_scores(
        sample_id=sample_id,
        runtime_result=runtime_result,
        feature_summary=vector.summary,
        input_fingerprint=source.input_fingerprint,
        source_kind=source.source_kind,
        genome_build=source.genome_build,
        artifact_lock_id=lock.lock_id,
        runtime_profile_id=runtime_profile.profile_id,
        annotations=annotations,
    )
    # Rehash all controlled resources after execution so a mid-run mutation cannot be
    # accepted merely because preflight passed before the file changed.
    verify_marlin_artifact_lock(lock, resources.artifact_paths)
    return report, runtime_result

from __future__ import annotations

import argparse
import os
import tempfile
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

from pydantic import BaseModel

from .execution import SubprocessRunner
from .io import load_model
from .marlin_artifacts import (
    MarlinArtifactPaths,
    create_marlin_artifact_lock,
    load_exported_feature_ids,
)
from .marlin_contracts import MarlinArtifactLock, MarlinRuntimeCompatibilityProfile
from .marlin_features import build_marlin_feature_vector
from .marlin_input import parse_marlin_probe_bed
from .marlin_runner import MarlinRunResources, run_precomputed_marlin_classification
from .models import GenomeBuild
from .reference import sha256_file

COMMANDS = frozenset({"marlin-lock", "marlin-features", "marlin-classify", "marlin-validate"})


def _write_text_atomic(path: Path, text: str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, staged_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    staged = Path(staged_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(staged, path)
    except BaseException:
        staged.unlink(missing_ok=True)
        raise
    return path


def _write_model_atomic(model: BaseModel, path: Path) -> Path:
    return _write_text_atomic(path, model.model_dump_json(indent=2) + "\n")


def _add_common_artifact_paths(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--feature-rdata", type=Path, required=True)
    parser.add_argument("--feature-list", type=Path, required=True)
    parser.add_argument("--class-annotations", type=Path, required=True)
    parser.add_argument("--probe-bed", type=Path, required=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ontseq", description="ONTSeq MARLIN commands")
    subparsers = parser.add_subparsers(dest="command", required=True)

    lock = subparsers.add_parser("marlin-lock", help="Lock local MARLIN v1 artifacts")
    _add_common_artifact_paths(lock)
    lock.add_argument("--lock-id", required=True)
    lock.add_argument("--code-version", required=True)
    lock.add_argument("--code-manifest-sha256", required=True)
    lock.add_argument("--genome-build", choices=[GenomeBuild.GRCH37.value], required=True)
    lock.add_argument("--r-version", required=True)
    lock.add_argument("--keras-version", required=True)
    lock.add_argument("--tensorflow-version", required=True)
    lock.add_argument("--python-version-if-used")
    lock.add_argument("--execution-backend", required=True)
    lock.add_argument("--output", type=Path, required=True)

    features = subparsers.add_parser(
        "marlin-features", help="Build the locked MARLIN v1 feature vector"
    )
    features.add_argument("--input", type=Path, required=True)
    features.add_argument("--genome-build", choices=[GenomeBuild.GRCH37.value], required=True)
    features.add_argument("--artifact-lock", type=Path, required=True)
    features.add_argument("--feature-list", type=Path, required=True)
    features.add_argument("--output-vector", type=Path, required=True)
    features.add_argument("--output-summary", type=Path, required=True)

    classify = subparsers.add_parser("marlin-classify", help="Run locked MARLIN v1 classification")
    classify.add_argument("--sample-id", required=True)
    classify.add_argument("--input", type=Path, required=True)
    classify.add_argument("--genome-build", choices=[GenomeBuild.GRCH37.value], required=True)
    classify.add_argument("--artifact-lock", type=Path, required=True)
    classify.add_argument("--runtime-profile", type=Path, required=True)
    _add_common_artifact_paths(classify)
    classify.add_argument("--inference-script", type=Path, required=True)
    classify.add_argument("--rscript", default="Rscript")
    classify.add_argument("--work-dir", type=Path)
    classify.add_argument("--output", type=Path, required=True)

    validate = subparsers.add_parser(
        "marlin-validate", help="Validate locked MARLIN classifications against a manifest"
    )
    validate.add_argument("--manifest", type=Path, required=True)
    validate.add_argument("--artifact-lock", type=Path, required=True)
    validate.add_argument("--runtime-profile", type=Path, required=True)
    _add_common_artifact_paths(validate)
    validate.add_argument("--inference-script", type=Path, required=True)
    validate.add_argument("--rscript", default="Rscript")
    validate.add_argument("--work-dir", type=Path)
    validate.add_argument("--output", type=Path, required=True)
    return parser


def _artifact_paths(args: argparse.Namespace) -> MarlinArtifactPaths:
    return MarlinArtifactPaths(
        model=args.model,
        feature_rdata=args.feature_rdata,
        canonical_feature_list=args.feature_list,
        class_annotations=args.class_annotations,
        probe_bed=args.probe_bed,
    )


def _run_lock(args: argparse.Namespace) -> Path:
    lock = create_marlin_artifact_lock(
        _artifact_paths(args),
        lock_id=args.lock_id,
        code_version_or_commit=args.code_version,
        code_manifest_sha256=args.code_manifest_sha256,
        genome_build=GenomeBuild(args.genome_build),
        R_version=args.r_version,
        keras_version=args.keras_version,
        tensorflow_version=args.tensorflow_version,
        python_version_if_used=args.python_version_if_used,
        execution_backend=args.execution_backend,
        created_at=datetime.now(UTC),
    )
    return _write_model_atomic(lock, args.output)


def _load_locked_features(feature_list: Path, lock: MarlinArtifactLock) -> tuple[str, ...]:
    observed_hash = sha256_file(feature_list)
    if observed_hash != lock.canonical_feature_list_sha256:
        raise ValueError("MARLIN canonical feature-list SHA-256 differs from artifact lock")
    return load_exported_feature_ids(feature_list)


def _run_features(args: argparse.Namespace) -> tuple[Path, Path]:
    lock = load_model(args.artifact_lock, MarlinArtifactLock)
    feature_ids = _load_locked_features(args.feature_list, lock)
    source = parse_marlin_probe_bed(args.input, genome_build=GenomeBuild(args.genome_build))
    vector = build_marlin_feature_vector(
        source,
        feature_ids,
        feature_artifact_sha256=lock.canonical_feature_list_sha256,
    )
    vector_path = _write_text_atomic(
        args.output_vector, "".join(f"{int(value)}\n" for value in vector.values)
    )
    summary_path = _write_model_atomic(vector.summary, args.output_summary)
    return vector_path, summary_path


def _check_runtime_profile_identity(
    profile: MarlinRuntimeCompatibilityProfile,
    lock: MarlinArtifactLock,
) -> None:
    if profile.reference_runtime_lock_id != lock.lock_id:
        raise ValueError("MARLIN runtime profile belongs to a different artifact/runtime lock")
    if profile.execution_backend != lock.execution_backend:
        raise ValueError("MARLIN runtime profile execution backend differs from artifact lock")


def _run_classify(args: argparse.Namespace) -> Path:
    lock = load_model(args.artifact_lock, MarlinArtifactLock)
    profile = load_model(args.runtime_profile, MarlinRuntimeCompatibilityProfile)
    resources = MarlinRunResources(_artifact_paths(args), args.inference_script)
    work_dir = args.work_dir or (args.output.parent / ".marlin-runtime")
    report, _runtime_result = run_precomputed_marlin_classification(
        sample_id=args.sample_id,
        input_path=args.input,
        genome_build=GenomeBuild(args.genome_build),
        lock=lock,
        runtime_profile=profile,
        resources=resources,
        runner=SubprocessRunner(),
        work_dir=work_dir,
        rscript_path=args.rscript,
    )
    return _write_model_atomic(report, args.output)


def run_command(args: argparse.Namespace) -> None:
    if args.command == "marlin-lock":
        print(_run_lock(args))
    elif args.command == "marlin-features":
        for path in _run_features(args):
            print(path)
    elif args.command == "marlin-classify":
        print(_run_classify(args))
    elif args.command == "marlin-validate":
        from .marlin_validation import run_validation_manifest

        print(
            run_validation_manifest(
                args.manifest,
                artifact_lock_path=args.artifact_lock,
                runtime_profile_path=args.runtime_profile,
                artifact_paths=_artifact_paths(args),
                inference_script=args.inference_script,
                output_path=args.output,
                rscript_path=args.rscript,
                work_dir=args.work_dir,
            )
        )
    else:
        raise ValueError(f"Unsupported MARLIN command: {args.command}")


def main(argv: Sequence[str] | None = None) -> None:
    args = _parser().parse_args(argv)
    run_command(args)

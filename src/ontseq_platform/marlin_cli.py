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
from .marlin_contracts import (
    MarlinArtifactLock,
    MarlinRuntimeCompatibilityProfile,
    MarlinRuntimeIdentity,
)
from .marlin_features import build_marlin_feature_vector
from .marlin_input import parse_marlin_probe_bed
from .marlin_runner import MarlinRunResources, run_precomputed_marlin_classification
from .marlin_runtime import probe_marlin_runtime
from .marlin_runtime_compare import (
    execute_marlin_dual_runtime_comparison,
    runtime_identity_from_probe,
)
from .marlin_runtime_freeze import (
    freeze_runtime_compatibility,
    render_frozen_runtime_fixture,
    verify_runtime_probe_matches_lock,
)
from .models import GenomeBuild
from .reference import sha256_file

COMMANDS = frozenset(
    {
        "marlin-lock",
        "marlin-runtime-probe",
        "marlin-freeze-runtime",
        "marlin-compare-runtimes",
        "marlin-features",
        "marlin-classify",
        "marlin-validate",
    }
)


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


def _require_new_output(path: Path, *, label: str) -> Path:
    candidate = Path(path)
    if candidate.exists():
        raise FileExistsError(f"{label} output already exists: {candidate}")
    return candidate


def _require_new_freeze_outputs(
    fixture_path: Path,
    profile_path: Path,
    runtime_identity_path: Path | None = None,
) -> tuple[Path, ...]:
    paths = [Path(fixture_path), Path(profile_path)]
    if runtime_identity_path is not None:
        paths.append(Path(runtime_identity_path))
    resolved = [path.resolve(strict=False) for path in paths]
    if len(resolved) != len(set(resolved)):
        raise ValueError("MARLIN runtime freeze outputs must use different paths")
    for path in paths:
        if path.exists():
            raise FileExistsError(f"MARLIN runtime freeze output already exists: {path}")
    return tuple(paths)


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
    lock.add_argument("--runtime-probe-script", type=Path, required=True)
    lock.add_argument("--rscript", default="Rscript")
    lock.add_argument("--output", type=Path, required=True)

    probe = subparsers.add_parser(
        "marlin-runtime-probe", help="Probe the live MARLIN R/Keras/TensorFlow runtime"
    )
    probe.add_argument("--probe-script", type=Path, required=True)
    probe.add_argument("--rscript", default="Rscript")
    probe.add_argument("--output", type=Path, required=True)

    freeze = subparsers.add_parser(
        "marlin-freeze-runtime",
        help="Freeze a non-biological MARLIN runtime compatibility fixture",
    )
    freeze.add_argument("--artifact-lock", type=Path, required=True)
    freeze.add_argument("--model", type=Path, required=True)
    freeze.add_argument("--inference-script", type=Path, required=True)
    freeze.add_argument("--runtime-probe-script", type=Path, required=True)
    freeze.add_argument("--profile-id", required=True)
    freeze.add_argument("--rscript", default="Rscript")
    freeze.add_argument("--work-dir", type=Path)
    freeze.add_argument("--output-fixture", type=Path, required=True)
    freeze.add_argument("--output-profile", type=Path, required=True)
    freeze.add_argument("--output-runtime-identity", type=Path)

    compare = subparsers.add_parser(
        "marlin-compare-runtimes",
        help="Compare a live MARLIN runtime against a frozen reference",
    )
    compare.add_argument("--reference-artifact-lock", type=Path, required=True)
    compare.add_argument("--candidate-artifact-lock", type=Path, required=True)
    compare.add_argument("--reference-runtime-identity", type=Path, required=True)
    compare.add_argument("--reference-profile", type=Path, required=True)
    compare.add_argument("--runtime-fixture", type=Path, required=True)
    compare.add_argument("--model", type=Path, required=True)
    compare.add_argument("--inference-script", type=Path, required=True)
    compare.add_argument("--candidate-runtime-probe-script", type=Path, required=True)
    compare.add_argument("--candidate-rscript", required=True)
    compare.add_argument("--comparison-id", required=True)
    compare.add_argument("--output", type=Path, required=True)

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
    validate.add_argument("--runtime-fixture", type=Path, required=True)
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


def _run_runtime_probe(args: argparse.Namespace) -> Path:
    report = probe_marlin_runtime(
        probe_script=args.probe_script,
        runner=SubprocessRunner(),
        rscript_path=args.rscript,
        timeout_seconds=120,
    )
    return _write_model_atomic(report, args.output)


def _run_lock(args: argparse.Namespace) -> Path:
    runtime = probe_marlin_runtime(
        probe_script=args.runtime_probe_script,
        runner=SubprocessRunner(),
        rscript_path=args.rscript,
        timeout_seconds=120,
    )
    lock = create_marlin_artifact_lock(
        _artifact_paths(args),
        lock_id=args.lock_id,
        code_version_or_commit=args.code_version,
        code_manifest_sha256=args.code_manifest_sha256,
        genome_build=GenomeBuild(args.genome_build),
        R_version=runtime.R_version,
        keras_version=runtime.keras_version,
        tensorflow_version=runtime.tensorflow_version,
        python_version_if_used=runtime.python_version,
        execution_backend=runtime.execution_backend,
        created_at=datetime.now(UTC),
    )
    return _write_model_atomic(lock, args.output)


def _run_freeze(args: argparse.Namespace) -> tuple[Path, ...]:
    output_paths = _require_new_freeze_outputs(
        args.output_fixture,
        args.output_profile,
        args.output_runtime_identity,
    )
    fixture_path = output_paths[0]
    profile_path = output_paths[1]
    identity_path = output_paths[2] if len(output_paths) == 3 else None
    lock = load_model(args.artifact_lock, MarlinArtifactLock)
    runner = SubprocessRunner()
    live_runtime = probe_marlin_runtime(
        probe_script=args.runtime_probe_script,
        runner=runner,
        rscript_path=args.rscript,
        timeout_seconds=120,
    )
    verify_runtime_probe_matches_lock(lock, live_runtime)
    frozen_at = datetime.now(UTC)
    work_dir = args.work_dir or (profile_path.parent / ".marlin-freeze-runtime")
    vector, _runtime_result, profile = freeze_runtime_compatibility(
        lock,
        model_path=args.model,
        inference_script=args.inference_script,
        runner=runner,
        work_dir=work_dir,
        profile_id=args.profile_id,
        created_at=frozen_at,
        rscript_path=args.rscript,
    )
    runtime_identity = (
        runtime_identity_from_probe(
            lock,
            live_runtime,
            runtime_id=f"{lock.lock_id}:runtime",
            created_at=frozen_at,
        )
        if identity_path is not None
        else None
    )

    written: list[Path] = []
    try:
        written.append(
            _write_text_atomic(
                fixture_path,
                render_frozen_runtime_fixture(vector),
            )
        )
        written.append(_write_model_atomic(profile, profile_path))
        if identity_path is not None and runtime_identity is not None:
            written.append(_write_model_atomic(runtime_identity, identity_path))
    except BaseException:
        for path in output_paths:
            path.unlink(missing_ok=True)
        raise
    return tuple(written)


def _run_compare_runtimes(args: argparse.Namespace) -> Path:
    output = _require_new_output(args.output, label="MARLIN dual-runtime comparison")
    reference_lock = load_model(args.reference_artifact_lock, MarlinArtifactLock)
    candidate_lock = load_model(args.candidate_artifact_lock, MarlinArtifactLock)
    reference_identity = load_model(args.reference_runtime_identity, MarlinRuntimeIdentity)
    reference_profile = load_model(args.reference_profile, MarlinRuntimeCompatibilityProfile)
    report = execute_marlin_dual_runtime_comparison(
        comparison_id=args.comparison_id,
        reference_lock=reference_lock,
        candidate_lock=candidate_lock,
        reference_runtime_identity=reference_identity,
        reference_profile=reference_profile,
        runtime_fixture_path=args.runtime_fixture,
        model_path=args.model,
        inference_script=args.inference_script,
        candidate_probe_script=args.candidate_runtime_probe_script,
        runner=SubprocessRunner(),
        work_dir=output.parent / ".marlin-runtime-compare",
        candidate_rscript_path=args.candidate_rscript,
        created_at=datetime.now(UTC),
    )
    return _write_model_atomic(report, output)


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
    elif args.command == "marlin-runtime-probe":
        print(_run_runtime_probe(args))
    elif args.command == "marlin-freeze-runtime":
        for path in _run_freeze(args):
            print(path)
    elif args.command == "marlin-compare-runtimes":
        print(_run_compare_runtimes(args))
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
                runtime_fixture_path=args.runtime_fixture,
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

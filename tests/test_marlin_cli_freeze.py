from __future__ import annotations

import hashlib
import sys
from datetime import UTC, datetime
from pathlib import Path

import pytest

from ontseq_platform import entrypoint, marlin_cli
from ontseq_platform.marlin_cli import (
    COMMANDS,
    _parser,
    _require_new_freeze_outputs,
    _run_freeze,
)
from ontseq_platform.marlin_contracts import (
    MarlinArtifactLock,
    MarlinRuntimeCompatibilityProfile,
    MarlinRuntimeIdentity,
)
from ontseq_platform.marlin_runtime import MarlinRuntimeProbeReport
from ontseq_platform.marlin_runtime_freeze import generate_frozen_runtime_fixture
from ontseq_platform.models import GenomeBuild


def _lock(model_sha256: str = "0" * 64) -> MarlinArtifactLock:
    return MarlinArtifactLock(
        lock_id="MARLIN_REFERENCE_LOCK",
        model_version="1.0.0",
        model_source_uri="https://doi.org/10.5281/zenodo.15565404",
        model_sha256=model_sha256,
        code_source_uri="https://github.com/hovestadt/MARLIN",
        code_version_or_commit="442aa603415a54f62e7367794f9a31c6bc20fc2d",
        code_manifest_sha256="1" * 64,
        feature_source_uri="marlin_v1.features.RData",
        feature_sha256="2" * 64,
        canonical_feature_list_sha256="3" * 64,
        class_annotation_source_uri="marlin_v1.class_annotations.xlsx",
        class_annotation_sha256="4" * 64,
        probe_resource_source_uri="marlin_v1.probes_hg19.bed.gz",
        probe_resource_sha256="5" * 64,
        genome_build=GenomeBuild.GRCH37,
        R_version="4.1.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version_if_used="3.10.21",
        execution_backend="cpu",
        created_at=datetime(2026, 9, 16, 4, 0, tzinfo=UTC),
    )


def _probe() -> MarlinRuntimeProbeReport:
    return MarlinRuntimeProbeReport(
        R_version="4.1.3",
        keras_version="2.13.0",
        tensorflow_version="2.13.0",
        python_version="3.10.21",
        execution_backend="cpu",
    )


def _freeze_args(tmp_path: Path, *, include_identity: bool) -> list[str]:
    args = [
        "marlin-freeze-runtime",
        "--artifact-lock",
        str(tmp_path / "lock.json"),
        "--model",
        str(tmp_path / "model.hdf5"),
        "--inference-script",
        str(tmp_path / "infer.R"),
        "--runtime-probe-script",
        str(tmp_path / "probe.R"),
        "--profile-id",
        "MARLIN_V1_CPU_FROZEN",
        "--output-fixture",
        str(tmp_path / "runtime-fixture.txt"),
        "--output-profile",
        str(tmp_path / "runtime-profile.json"),
    ]
    if include_identity:
        args += [
            "--output-runtime-identity",
            str(tmp_path / "reference-runtime-identity.json"),
        ]
    return args


def _install_freeze_fakes(
    monkeypatch: pytest.MonkeyPatch,
    *,
    lock: MarlinArtifactLock,
    probe_calls: list[int],
) -> None:
    monkeypatch.setattr(marlin_cli, "load_model", lambda *_args, **_kwargs: lock)

    def fake_probe(**_kwargs: object) -> MarlinRuntimeProbeReport:
        probe_calls.append(1)
        return _probe()

    def fake_freeze(
        lock_arg: MarlinArtifactLock,
        **kwargs: object,
    ):
        created_at = kwargs["created_at"]
        assert isinstance(created_at, datetime)
        vector = generate_frozen_runtime_fixture(lock_arg)
        profile = MarlinRuntimeCompatibilityProfile(
            profile_id=str(kwargs["profile_id"]),
            reference_runtime_lock_id=lock_arg.lock_id,
            feature_vector_sha256=vector.summary.feature_vector_sha256,
            reference_scores=[1.0] + [0.0] * 41,
            absolute_score_tolerance=1e-7,
            score_sum_tolerance=1e-5,
            top_model_unit_index=0,
            execution_backend=lock_arg.execution_backend,
            created_at=created_at,
        )
        return vector, object(), profile

    monkeypatch.setattr(marlin_cli, "probe_marlin_runtime", fake_probe)
    monkeypatch.setattr(marlin_cli, "freeze_runtime_compatibility", fake_freeze)


def test_global_help_lists_runtime_freeze_command(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(sys, "argv", ["ontseq", "--help"])
    entrypoint.main()
    output = capsys.readouterr().out
    assert "marlin-freeze-runtime" in output
    assert "marlin-freeze-runtime" in COMMANDS


def test_runtime_freeze_requires_locked_runtime_inputs_and_outputs() -> None:
    args = _parser().parse_args(
        [
            "marlin-freeze-runtime",
            "--artifact-lock",
            "lock.json",
            "--model",
            "model.hdf5",
            "--inference-script",
            "infer.R",
            "--runtime-probe-script",
            "probe.R",
            "--profile-id",
            "MARLIN_V1_CPU_FROZEN",
            "--output-fixture",
            "runtime-fixture.txt",
            "--output-profile",
            "runtime-profile.json",
        ]
    )
    assert args.command == "marlin-freeze-runtime"
    assert args.artifact_lock == Path("lock.json")
    assert args.model == Path("model.hdf5")
    assert args.inference_script == Path("infer.R")
    assert args.runtime_probe_script == Path("probe.R")
    assert args.profile_id == "MARLIN_V1_CPU_FROZEN"
    assert args.output_fixture == Path("runtime-fixture.txt")
    assert args.output_profile == Path("runtime-profile.json")
    assert args.output_runtime_identity is None


def test_runtime_freeze_accepts_optional_reference_runtime_identity_output() -> None:
    args = _parser().parse_args(
        [
            "marlin-freeze-runtime",
            "--artifact-lock",
            "lock.json",
            "--model",
            "model.hdf5",
            "--inference-script",
            "infer.R",
            "--runtime-probe-script",
            "probe.R",
            "--profile-id",
            "MARLIN_V1_CPU_FROZEN",
            "--output-fixture",
            "runtime-fixture.txt",
            "--output-profile",
            "runtime-profile.json",
            "--output-runtime-identity",
            "reference-runtime-identity.json",
        ]
    )

    assert args.output_runtime_identity == Path("reference-runtime-identity.json")


def test_runtime_freeze_rejects_missing_live_probe_argument() -> None:
    with pytest.raises(SystemExit):
        _parser().parse_args(
            [
                "marlin-freeze-runtime",
                "--artifact-lock",
                "lock.json",
                "--model",
                "model.hdf5",
                "--inference-script",
                "infer.R",
                "--profile-id",
                "MARLIN_V1_CPU_FROZEN",
                "--output-fixture",
                "runtime-fixture.txt",
                "--output-profile",
                "runtime-profile.json",
            ]
        )


def test_runtime_freeze_refuses_existing_or_colliding_outputs(tmp_path: Path) -> None:
    fixture = tmp_path / "runtime-fixture.txt"
    profile = tmp_path / "runtime-profile.json"
    identity = tmp_path / "reference-runtime-identity.json"
    fixture.write_text("existing\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="already exists"):
        _require_new_freeze_outputs(fixture, profile, identity)

    fixture.unlink()
    with pytest.raises(ValueError, match="different paths"):
        _require_new_freeze_outputs(profile, profile, identity)
    with pytest.raises(ValueError, match="different paths"):
        _require_new_freeze_outputs(fixture, profile, profile)


def test_runtime_freeze_persists_identity_from_same_live_probe_and_timestamp(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lock = _lock(hashlib.sha256(b"model").hexdigest())
    probe_calls: list[int] = []
    _install_freeze_fakes(monkeypatch, lock=lock, probe_calls=probe_calls)
    args = _parser().parse_args(_freeze_args(tmp_path, include_identity=True))

    written = _run_freeze(args)

    assert probe_calls == [1]
    assert len(written) == 3
    profile = MarlinRuntimeCompatibilityProfile.model_validate_json(
        (tmp_path / "runtime-profile.json").read_text(encoding="utf-8")
    )
    identity = MarlinRuntimeIdentity.model_validate_json(
        (tmp_path / "reference-runtime-identity.json").read_text(encoding="utf-8")
    )
    assert identity.runtime_id == f"{lock.lock_id}:runtime"
    assert identity.R_version == lock.R_version
    assert identity.python_version == lock.python_version_if_used
    assert identity.execution_backend == lock.execution_backend
    assert identity.created_at == profile.created_at


def test_runtime_freeze_without_identity_remains_backward_compatible(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lock = _lock(hashlib.sha256(b"model").hexdigest())
    probe_calls: list[int] = []
    _install_freeze_fakes(monkeypatch, lock=lock, probe_calls=probe_calls)
    args = _parser().parse_args(_freeze_args(tmp_path, include_identity=False))

    written = _run_freeze(args)

    assert probe_calls == [1]
    assert len(written) == 2
    assert (tmp_path / "runtime-fixture.txt").is_file()
    assert (tmp_path / "runtime-profile.json").is_file()
    assert not (tmp_path / "reference-runtime-identity.json").exists()


def test_runtime_freeze_rolls_back_siblings_if_identity_write_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lock = _lock(hashlib.sha256(b"model").hexdigest())
    probe_calls: list[int] = []
    _install_freeze_fakes(monkeypatch, lock=lock, probe_calls=probe_calls)
    args = _parser().parse_args(_freeze_args(tmp_path, include_identity=True))
    identity_path = tmp_path / "reference-runtime-identity.json"
    original_write_model_atomic = marlin_cli._write_model_atomic

    def fail_identity(model: object, path: Path) -> Path:
        if path == identity_path:
            raise RuntimeError("synthetic identity write failure")
        return original_write_model_atomic(model, path)  # type: ignore[arg-type]

    monkeypatch.setattr(marlin_cli, "_write_model_atomic", fail_identity)

    with pytest.raises(RuntimeError, match="identity write failure"):
        _run_freeze(args)

    assert not (tmp_path / "runtime-fixture.txt").exists()
    assert not (tmp_path / "runtime-profile.json").exists()
    assert not identity_path.exists()

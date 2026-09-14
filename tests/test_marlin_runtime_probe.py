from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from types import SimpleNamespace

import pytest

from ontseq_platform import marlin_cli, marlin_runtime
from ontseq_platform.execution import CommandResult, ToolExecutionError


class ProbeRunner:
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
            stdout="probe-log",
            stderr=self.stderr,
        )


def _payload(
    *,
    r_version: str = "4.2.3",
    keras_version: str = "2.13.0",
    tensorflow_version: str = "2.13.0",
    python_version: str = "3.10.14",
    backend: str = "cpu",
) -> str:
    return (
        "key\tvalue\n"
        f"R_version\t{r_version}\n"
        f"keras_version\t{keras_version}\n"
        f"tensorflow_version\t{tensorflow_version}\n"
        f"python_version\t{python_version}\n"
        f"execution_backend\t{backend}\n"
    )


def _probe_function():
    probe = getattr(marlin_runtime, "probe_marlin_runtime", None)
    assert callable(probe), "MARLIN runtime probe is not implemented"
    return probe


def test_runtime_probe_reports_actual_versions_and_backend(tmp_path: Path) -> None:
    script = tmp_path / "probe.R"
    script.write_text("# fixture\n", encoding="utf-8")
    runner = ProbeRunner(_payload())

    report = _probe_function()(
        probe_script=script,
        runner=runner,
        rscript_path="Rscript",
        timeout_seconds=30,
    )

    assert report.R_version == "4.2.3"
    assert report.keras_version == "2.13.0"
    assert report.tensorflow_version == "2.13.0"
    assert report.python_version == "3.10.14"
    assert report.execution_backend == "cpu"
    assert runner.argv is not None
    assert runner.argv[0] == "Rscript"
    assert runner.argv[1] == str(script.resolve())


def test_runtime_probe_requires_patch_level_python_version(tmp_path: Path) -> None:
    script = tmp_path / "probe.R"
    script.write_text("# fixture\n", encoding="utf-8")

    with pytest.raises(ValueError, match="python_version"):
        _probe_function()(
            probe_script=script,
            runner=ProbeRunner(_payload(python_version="3.10")),
        )


def test_runtime_probe_rejects_missing_or_duplicate_fields(tmp_path: Path) -> None:
    script = tmp_path / "probe.R"
    script.write_text("# fixture\n", encoding="utf-8")
    missing = _payload().replace("python_version\t3.10.14\n", "")
    with pytest.raises(ValueError, match="runtime probe fields"):
        _probe_function()(
            probe_script=script,
            runner=ProbeRunner(missing),
        )

    duplicate = _payload() + "R_version\t4.2.3\n"
    with pytest.raises(ValueError, match="duplicate"):
        _probe_function()(
            probe_script=script,
            runner=ProbeRunner(duplicate),
        )


def test_runtime_probe_rejects_invalid_backend_and_process_failure(tmp_path: Path) -> None:
    script = tmp_path / "probe.R"
    script.write_text("# fixture\n", encoding="utf-8")
    with pytest.raises(ValueError, match="execution_backend"):
        _probe_function()(
            probe_script=script,
            runner=ProbeRunner(_payload(backend="tpu")),
        )

    with pytest.raises(ToolExecutionError, match="runtime probe failed"):
        _probe_function()(
            probe_script=script,
            runner=ProbeRunner(None, returncode=7, stderr="tensorflow import failed"),
        )


def _new_lock_argv() -> list[str]:
    return [
        "marlin-lock",
        "--model",
        "model.hdf5",
        "--feature-rdata",
        "features.RData",
        "--feature-list",
        "features.txt",
        "--class-annotations",
        "classes.xlsx",
        "--probe-bed",
        "probes.bed.gz",
        "--lock-id",
        "LOCK",
        "--code-version",
        "commit",
        "--code-manifest-sha256",
        "1" * 64,
        "--genome-build",
        "GRCh37",
        "--runtime-probe-script",
        "scripts/marlin_runtime_probe.R",
        "--rscript",
        "Rscript",
        "--output",
        "lock.json",
    ]


def test_marlin_cli_exposes_runtime_probe_and_lock_requires_live_probe() -> None:
    parser = marlin_cli._parser()
    args = parser.parse_args(
        [
            "marlin-runtime-probe",
            "--probe-script",
            "scripts/marlin_runtime_probe.R",
            "--output",
            "runtime.json",
        ]
    )
    assert args.command == "marlin-runtime-probe"
    assert args.probe_script == Path("scripts/marlin_runtime_probe.R")

    lock_args = parser.parse_args(_new_lock_argv())
    assert lock_args.runtime_probe_script == Path("scripts/marlin_runtime_probe.R")

    old_manual_args = [
        item
        for item in _new_lock_argv()
        if item not in {"--runtime-probe-script", "scripts/marlin_runtime_probe.R"}
    ]
    old_manual_args.extend(
        [
            "--r-version",
            "4.1.3",
            "--keras-version",
            "2.13.0",
            "--tensorflow-version",
            "2.13.0",
            "--python-version-if-used",
            "3.10.14",
            "--execution-backend",
            "cpu",
        ]
    )
    with pytest.raises(SystemExit):
        parser.parse_args(old_manual_args)


def test_marlin_lock_populates_runtime_identity_from_probe(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_probe(*, probe_script: Path, runner: object, rscript_path: str, timeout_seconds: int):
        captured["probe_script"] = probe_script
        captured["rscript_path"] = rscript_path
        return SimpleNamespace(
            R_version="4.2.3",
            keras_version="2.13.0",
            tensorflow_version="2.13.0",
            python_version="3.10.14",
            execution_backend="cpu",
        )

    def fake_create(_paths: object, **kwargs: object):
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(marlin_cli, "probe_marlin_runtime", fake_probe, raising=False)
    monkeypatch.setattr(marlin_cli, "create_marlin_artifact_lock", fake_create)
    monkeypatch.setattr(marlin_cli, "_write_model_atomic", lambda _model, path: path)

    args = marlin_cli._parser().parse_args(_new_lock_argv())
    args.output = tmp_path / "lock.json"
    output = marlin_cli._run_lock(args)

    assert output == tmp_path / "lock.json"
    assert captured["R_version"] == "4.2.3"
    assert captured["keras_version"] == "2.13.0"
    assert captured["tensorflow_version"] == "2.13.0"
    assert captured["python_version_if_used"] == "3.10.14"
    assert captured["execution_backend"] == "cpu"

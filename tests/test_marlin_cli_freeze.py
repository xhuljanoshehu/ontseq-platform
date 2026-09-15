from __future__ import annotations

import sys
from pathlib import Path

import pytest

from ontseq_platform import entrypoint
from ontseq_platform.marlin_cli import COMMANDS, _parser, _require_new_freeze_outputs


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
    assert args.profile_id == "MARLIN_V1_CPU_FROZEN"
    assert args.output_fixture == Path("runtime-fixture.txt")
    assert args.output_profile == Path("runtime-profile.json")


def test_runtime_freeze_refuses_existing_or_colliding_outputs(tmp_path: Path) -> None:
    fixture = tmp_path / "runtime-fixture.txt"
    profile = tmp_path / "runtime-profile.json"
    fixture.write_text("existing\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="already exists"):
        _require_new_freeze_outputs(fixture, profile)

    fixture.unlink()
    with pytest.raises(ValueError, match="different paths"):
        _require_new_freeze_outputs(profile, profile)

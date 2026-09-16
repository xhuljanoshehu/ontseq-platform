from __future__ import annotations

import sys
from pathlib import Path

import pytest

from ontseq_platform import entrypoint
from ontseq_platform.marlin_cli import COMMANDS, _parser


def _compare_args() -> list[str]:
    return [
        "marlin-compare-runtimes",
        "--reference-artifact-lock",
        "reference-lock.json",
        "--candidate-artifact-lock",
        "candidate-lock.json",
        "--reference-runtime-identity",
        "reference-runtime.json",
        "--reference-profile",
        "reference-profile.json",
        "--runtime-fixture",
        "runtime-fixture.txt",
        "--model",
        "model.hdf5",
        "--inference-script",
        "infer.R",
        "--candidate-runtime-probe-script",
        "probe.R",
        "--candidate-rscript",
        "candidate-Rscript",
        "--comparison-id",
        "MARLIN_R413_VS_R423",
        "--output",
        "runtime-comparison.json",
    ]


def test_global_help_lists_dual_runtime_comparison(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(sys, "argv", ["ontseq", "--help"])
    entrypoint.main()
    output = capsys.readouterr().out

    assert "marlin-compare-runtimes" in output
    assert "marlin-compare-runtimes" in COMMANDS


def test_dual_runtime_comparison_parser_requires_exact_inputs() -> None:
    args = _parser().parse_args(_compare_args())

    assert args.command == "marlin-compare-runtimes"
    assert args.reference_artifact_lock == Path("reference-lock.json")
    assert args.candidate_artifact_lock == Path("candidate-lock.json")
    assert args.reference_runtime_identity == Path("reference-runtime.json")
    assert args.reference_profile == Path("reference-profile.json")
    assert args.runtime_fixture == Path("runtime-fixture.txt")
    assert args.model == Path("model.hdf5")
    assert args.inference_script == Path("infer.R")
    assert args.candidate_runtime_probe_script == Path("probe.R")
    assert args.candidate_rscript == "candidate-Rscript"
    assert args.comparison_id == "MARLIN_R413_VS_R423"
    assert args.output == Path("runtime-comparison.json")


def test_dual_runtime_comparison_rejects_missing_required_argument() -> None:
    incomplete = _compare_args()
    index = incomplete.index("--reference-runtime-identity")
    del incomplete[index : index + 2]

    with pytest.raises(SystemExit):
        _parser().parse_args(incomplete)


def test_dual_runtime_comparison_refuses_existing_output(tmp_path: Path) -> None:
    from ontseq_platform.marlin_cli import _require_new_output

    output = tmp_path / "runtime-comparison.json"
    output.write_text("existing\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="already exists"):
        _require_new_output(output, label="MARLIN dual-runtime comparison")

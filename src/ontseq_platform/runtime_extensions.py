from __future__ import annotations

from pathlib import Path

from .io import load_model
from .target_coverage import TargetCoveragePolicy
from .target_coverage_extension import (
    TargetCoverageExtensionSettings,
    register_target_coverage_extension,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _target_coverage_policy() -> TargetCoveragePolicy:
    path = _repo_root() / "configs/qc/adaptive_target_coverage.technical.yaml"
    if path.is_file():
        return load_model(path, TargetCoveragePolicy)
    return TargetCoveragePolicy(
        profile_id="adaptive_target_coverage_technical_v1",
        status="technical_defaults_only",
        expected_version="0.3.14",
        thresholds=[1, 10, 20, 30],
        mapq=0,
        exclude_flags=1796,
        note=(
            "Technical descriptive defaults only. Coverage thresholds are not validated assay "
            "adequacy, fusion reportability, CNV reportability, or clinical no-call thresholds."
        ),
    )


def register_builtin_runtime_extensions() -> None:
    """Register execution extensions that are part of the packaged open-source runtime."""

    register_target_coverage_extension(
        TargetCoverageExtensionSettings(policy=_target_coverage_policy())
    )

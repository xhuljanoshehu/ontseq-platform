"""Resolve typed coverage artifacts from core and installed-runtime envelopes."""

from __future__ import annotations

from pathlib import Path

from .models import SampleManifest
from .target_coverage import TargetCoverageReport

COVERAGE_ARTIFACT_CONTRACT = "core-or-sample-coverage-v1"


def load_run_coverage(
    root: Path, manifest: SampleManifest, *, selection: bool = False
) -> TargetCoverageReport | None:
    """Accept either producer's exact filename, refusing inconsistent duplicates."""
    stem = "selection-coverage" if selection else "target-coverage"
    reports = []
    for name in (f"{stem}.json", f"{manifest.sample_id}.{stem}.json"):
        path = (root / "qc" / name).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError("Coverage artifact path escapes the run envelope")
        if not path.is_file():
            continue
        report = TargetCoverageReport.model_validate_json(path.read_text(encoding="utf-8"))
        if report.sample_id != manifest.sample_id:
            raise ValueError("Coverage artifact belongs to a different sample")
        if report.genome_build != manifest.assay.genome_build:
            raise ValueError("Coverage artifact uses a different genome build")
        reports.append(report)
    if not reports:
        return None
    if any(report != reports[0] for report in reports[1:]):
        raise ValueError("Conflicting coverage artifacts in the run envelope")
    return reports[0]

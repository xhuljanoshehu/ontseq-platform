"""Resolve typed coverage artifacts from run envelopes.

Inside a run, the SV and report stages read only the coverage the target-coverage stage
recorded for that run (:func:`ontseq_platform.pipeline.runner.load_current_coverage`).
:func:`load_run_coverage` is the reader for archived envelopes outside a run, where no stage
record is at hand; it also accepts the sample-named file the retired runtime extension of
ONTSeq <= 0.8.2 wrote, and refuses an envelope that holds two different answers.
"""

from __future__ import annotations

from pathlib import Path

from .models import SampleManifest
from .target_coverage import TargetCoverageReport

#: The in-run handoff: SV observability and reviewer reports consume only the current,
#: checksum-verified target-coverage stage artifacts. Recorded in their plans, so changing
#: the handoff changes their resume signatures.
COVERAGE_ARTIFACT_CONTRACT = "current-stage-coverage-v2"
#: The archived-envelope reader below, which accepts either historical producer name.
ARCHIVED_COVERAGE_READER = "core-or-sample-coverage-v1"


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

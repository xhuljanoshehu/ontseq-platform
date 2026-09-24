"""Reformat a concluded local run without writing to its immutable envelope."""

from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path

from ..cnv.qdnaseq import QDNAseqCallReport
from ..coverage_artifacts import load_run_coverage
from ..marlin_native_contracts import NativeMarlinReport
from ..methylation import MethylationReport
from ..models import PipelineResult
from ..pipeline.lock import LOCK_FILENAME
from ..pipeline.runner import METHYLATION_REPORT, QC_READ_LENGTH_HISTOGRAM, RUN_REPORT
from ..pipeline.state import RunReport
from ..qc import read_length_histogram_from_tsv
from ..report import render_html
from ..report_marlin import validate_marlin_identity
from ..report_methylation import validate_methylation_identity
from ..report_plots import ReadLengthBin
from .guard import resolve_within


class ArchivedEvidence:
    """Read only checksum-bound artifacts of the current, unlocked result assembly."""

    def __init__(self, envelope: Path, result: PipelineResult) -> None:
        self.envelope = envelope
        self.verified_bytes: dict[str, bytes] = {}
        if (envelope / LOCK_FILENAME).exists():
            raise ValueError("the run envelope is still locked")
        run_path = resolve_within(envelope / RUN_REPORT, [envelope])
        self.run = (
            RunReport.model_validate_json(self.read(run_path)) if run_path.is_file() else None
        )
        if self.run is None:
            return
        run = self.run
        if (run.run_id, run.sample_id, run.genome_build, run.manifest.analysis.profile) != (
            result.manifest.run_id,
            result.manifest.sample_id,
            result.manifest.assay.genome_build,
            result.manifest.analysis.profile,
        ):
            raise ValueError("run identity does not match the result")
        assembly = [stage for stage in run.stages if stage.stage.value == "assemble"]
        if len(assembly) != 1 or assembly[0].status.value != "COMPLETED":
            raise ValueError("the current run has no completed result assembly")
        relative = f"normalized/{result.manifest.sample_id}.result.json"
        self.checked(relative)
        if len([item for item in assembly[0].outputs if item.relative_path == relative]) != 1:
            raise ValueError("result evidence is missing from the current assembly")
        if PipelineResult.model_validate_json(self.verified_bytes[relative]) != result:
            raise ValueError("normalized result differs from its current assembly evidence")

    @staticmethod
    def read(path: Path) -> bytes:
        with path.open("rb") as stream:
            content = stream.read(64 * 1024 * 1024 + 1)
        if len(content) > 64 * 1024 * 1024:
            raise ValueError("report evidence exceeds the export size limit")
        return content

    def checked(self, relative: str) -> Path:
        path = resolve_within(self.envelope / relative, [self.envelope])
        if self.run is not None:
            records = [
                record
                for stage in self.run.stages
                for record in stage.outputs
                if record.relative_path == relative
            ]
            if records and not path.is_file():
                raise ValueError("recorded report evidence is missing")
            if not records and not path.is_file():
                return path
            content = self.read(path)
            if len(records) != 1 or (
                len(content) != records[0].size_bytes
                or hashlib.sha256(content).hexdigest() != records[0].sha256
            ):
                raise ValueError("report evidence does not match its recorded checksum")
            self.verified_bytes[relative] = content
        return path


def load_archived_marlin(
    envelope: Path,
    result: PipelineResult,
    *,
    evidence: ArchivedEvidence | None = None,
) -> NativeMarlinReport | None:
    """Load current typed MARLIN evidence; failed states never read leftover prediction files."""
    evidence = evidence if evidence is not None else ArchivedEvidence(envelope, result)
    run = evidence.run
    checked = evidence.checked
    sample = result.manifest.sample_id
    marlin_relative = f"evidence/marlin/{sample}.marlin.json"
    marlin_stages = [stage for stage in run.stages if stage.stage.value == "marlin"] if run else []
    marlin_outcomes = [outcome for outcome in result.modules if outcome.module.value == "marlin"]
    marlin = None
    marlin_requested = any(
        item.value in {"marlin", "methylation"} for item in result.manifest.analysis.modules
    )
    unrequested = (
        not marlin_requested
        and not marlin_outcomes
        and all(stage.status.value == "NOT_RUN" for stage in marlin_stages)
    )
    if unrequested and result.provenance.reference_checksums.get("marlin_report"):
        raise ValueError("MARLIN unrequested run cannot claim prediction evidence")
    if not unrequested and (
        marlin_stages
        or marlin_outcomes
        or result.provenance.reference_checksums.get("marlin_report")
    ):
        if len(marlin_stages) != 1 or len(marlin_outcomes) != 1:
            raise ValueError("MARLIN evidence has no unique recorded execution outcome")
        stage, outcome = marlin_stages[0], marlin_outcomes[0]
        stage_reason = stage.reason
        resume_prefix = "Resumed unchanged from a previous run. "
        while stage_reason.startswith(resume_prefix):
            stage_reason = stage_reason[len(resume_prefix) :]
        if (stage.status, stage_reason, stage.tools) != (
            outcome.status,
            outcome.reason,
            outcome.tools,
        ):
            raise ValueError("MARLIN execution does not match the normalized result")
        if stage.status.value in {"NOT_RUN", "FAILED"}:
            if result.provenance.reference_checksums.get("marlin_report"):
                raise ValueError("MARLIN blocked/failed stage cannot use stale prediction evidence")
            marlin = NativeMarlinReport(
                run_id=result.manifest.run_id,
                sample_id=sample,
                genome_build=result.manifest.assay.genome_build,
                status=stage.status,
                reason=stage.reason,
                tools=stage.tools,
                warnings=stage.warnings,
                limitations=stage.limitations,
            )
        else:
            marlin_path = checked(marlin_relative)
            records = [
                record for record in stage.outputs if record.relative_path == marlin_relative
            ]
            if len(records) != 1 or not marlin_path.is_file():
                raise ValueError("MARLIN concluded execution is missing recorded evidence")
            if result.provenance.reference_checksums.get("marlin_report") != records[0].sha256:
                raise ValueError("MARLIN evidence checksum differs from the normalized result")
            marlin = NativeMarlinReport.model_validate_json(
                evidence.verified_bytes[marlin_relative]
            )
    if marlin is not None:
        validate_marlin_identity(result, marlin)
    return marlin


def current_befund(envelope: Path, result: PipelineResult) -> bytes:
    evidence = ArchivedEvidence(envelope, result)
    checked = evidence.checked
    sample = result.manifest.sample_id
    checked(f"normalized/{sample}.result.json")
    cnv_path = checked(f"evidence/cnv/{sample}.qdnaseq.json")
    cnv = (
        QDNAseqCallReport.model_validate_json(cnv_path.read_bytes()) if cnv_path.is_file() else None
    )
    if cnv is not None:
        outcomes = [item for item in result.modules if item.module.value == "cnv"]
        if (
            len(outcomes) != 1
            or outcomes[0].status != cnv.status
            or outcomes[0].tools != cnv.tools
            or any(event not in result.events for event in cnv.events)
        ):
            raise ValueError("CNV evidence does not match the normalized result")
        for fit in cnv.fits:
            for name in (fit.copy_number_plot, fit.fit_plot):
                checked(f"evidence/cnv/qdnaseq/{name}")
    methyl_path = checked(METHYLATION_REPORT.format(sample=sample))
    methylation = (
        MethylationReport.model_validate_json(methyl_path.read_bytes())
        if methyl_path.is_file()
        else None
    )
    try:
        validate_methylation_identity(result, methylation)
    except ValueError as exc:
        raise ValueError("methylation evidence does not match the result") from exc
    marlin = load_archived_marlin(envelope, result, evidence=evidence)
    for name in (
        "target-coverage.json",
        "selection-coverage.json",
        f"{sample}.target-coverage.json",
        f"{sample}.selection-coverage.json",
    ):
        checked(f"qc/{name}")
    target = load_run_coverage(envelope, result.manifest)
    selection = load_run_coverage(envelope, result.manifest, selection=True)
    histogram_path = checked(QC_READ_LENGTH_HISTOGRAM)
    histogram = (
        [
            ReadLengthBin(start=a, end=b, count=c, bases=d)
            for a, b, c, d in read_length_histogram_from_tsv(histogram_path.read_text())
        ]
        if histogram_path.is_file()
        else None
    )
    with tempfile.TemporaryDirectory(prefix="ontseq-befund-") as directory:
        path = render_html(
            result,
            Path(directory) / "report.html",
            cnv_report=cnv,
            cnv_evidence_root=envelope / "evidence/cnv/qdnaseq",
            methylation_report=methylation,
            marlin_report=marlin,
            qc_histogram=histogram,
            target_coverage=target,
            selection_coverage=selection,
        )
        return path.read_bytes()

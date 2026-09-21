"""Reformat a concluded local run without writing to its immutable envelope."""

from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path

from ..cnv.qdnaseq import QDNAseqCallReport
from ..coverage_artifacts import load_run_coverage
from ..methylation import MethylationReport
from ..models import PipelineResult
from ..pipeline.lock import LOCK_FILENAME
from ..pipeline.runner import METHYLATION_REPORT, QC_READ_LENGTH_HISTOGRAM, RUN_REPORT
from ..pipeline.state import RunReport
from ..qc import read_length_histogram_from_tsv
from ..report import render_html
from ..report_plots import ReadLengthBin
from .guard import resolve_within


def current_befund(envelope: Path, result: PipelineResult) -> bytes:
    if (envelope / LOCK_FILENAME).exists():
        raise ValueError("the run envelope is still locked")
    run_path = resolve_within(envelope / RUN_REPORT, [envelope])
    run = RunReport.model_validate_json(run_path.read_bytes()) if run_path.is_file() else None
    if run is not None and (run.run_id, run.sample_id, run.genome_build) != (
        result.manifest.run_id,
        result.manifest.sample_id,
        result.manifest.assay.genome_build,
    ):
        raise ValueError("run identity does not match the result")

    def checked(relative: str) -> Path:
        path = resolve_within(envelope / relative, [envelope])
        if run is not None:
            records = [
                record
                for stage in run.stages
                for record in stage.outputs
                if record.relative_path == relative
            ]
            if records and not path.is_file():
                raise ValueError("recorded report evidence is missing")
            if not records and not path.is_file():
                return path
            if len(records) != 1 or (
                path.stat().st_size != records[0].size_bytes
                or hashlib.sha256(path.read_bytes()).hexdigest() != records[0].sha256
            ):
                raise ValueError("report evidence does not match its recorded checksum")
        return path

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
    if methylation is not None and (
        (methylation.sample_id, methylation.genome_build)
        != (sample, result.manifest.assay.genome_build)
        or result.provenance.reference_checksums.get("bedmethyl")
        != methylation.bedmethyl_fingerprint.sha256
        or not any(
            item.module.value == "methylation"
            and item.status == methylation.status
            and methylation.tool in item.tools
            for item in result.modules
        )
    ):
        raise ValueError("methylation evidence does not match the result")
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
            qc_histogram=histogram,
            target_coverage=target,
            selection_coverage=selection,
        )
        return path.read_bytes()

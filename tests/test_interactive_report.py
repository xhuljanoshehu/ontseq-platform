from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from ontseq_platform.demo import build_demo_result
from ontseq_platform.report import render_html


def payload(document: str) -> dict:
    match = re.search(
        r'<script id="ontseq-befund-data" type="application/json">(.*?)</script>', document, re.S
    )
    assert match, "The portable report must include its run-specific interactive contract"
    return json.loads(match.group(1))


def test_report_exports_actual_run_without_input_paths_or_script_injection(tmp_path: Path) -> None:
    result = build_demo_result()
    result.manifest.sample_id = "SYNTHETIC_INTERACTIVE"
    result.warnings += ['</script><script>alert("unsafe")</script>', "/private/SECRET.bam"]
    before = result.model_dump_json()
    document = render_html(result, tmp_path / "report.html").read_text()
    data = payload(document)
    assert data["view"]["sample_id"] == "SYNTHETIC_INTERACTIVE"
    assert data["view"]["release_status"] == result.release_status.value
    assert data["cnv"] is None
    assert data["requested_modules"] == [m.value for m in result.manifest.analysis.modules]
    assert "SECRET.bam" not in document
    assert "</script><script>alert" not in document
    assert result.model_dump_json() == before
    assert 'id="ontseq-befund-root"' in document
    assert 'id="ontseq-static-report"' in document
    assert "befund-interactive-v1" in document


def test_cnv_attachment_requires_same_sample_and_build(tmp_path: Path) -> None:
    from ontseq_platform.cnv.qdnaseq import CnvFit, QDNAseqCallReport

    result = build_demo_result()
    fit = CnvFit(
        bin_size_kbp=500,
        cellularity=0.5,
        ploidy=4,
        fit_error=0.01,
        candidate_count=2,
        segment_count=0,
        segment_file="segments.tsv",
        chromosome_file="chromosomes.tsv",
        fit_plot="fit.png",
        copy_number_plot="copy.png",
        rds_file="model.rds",
        alternatives=[dict(cellularity=0.8, ploidy=2, fit_error=0.02)],
    )
    cnv = QDNAseqCallReport(
        sample_id=result.manifest.sample_id,
        genome_build=result.manifest.assay.genome_build,
        status="NO_CALL",
        primary_fit=fit,
        fits=[fit],
        chromosome_consensus=[],
        events=[],
        tools=[],
        output_files=[],
    )
    document = render_html(result, tmp_path / "report.html", cnv_report=cnv).read_text()
    data = payload(document)
    assert data["cnv"]["primary_fit"]["alternatives"][0]["ploidy"] == 2
    cnv.sample_id = "SYNTHETIC_OTHER"
    with pytest.raises(ValueError, match="sample|identity"):
        render_html(result, tmp_path / "bad.html", cnv_report=cnv)


def test_report_without_module_outcomes_preserves_unknown(tmp_path: Path) -> None:
    result = build_demo_result()
    result.modules = []
    data = payload(render_html(result, tmp_path / "report.html").read_text())
    assert data["view"]["modules"] == []


def test_assets_reject_traversal_and_keep_missing_plot_explicit(tmp_path: Path) -> None:
    from ontseq_platform.report_interactive import read_plot

    assert read_plot(tmp_path, "missing.png") is None
    for name in ("../secret.png", "/etc/passwd", "nested/../../secret.png"):
        with pytest.raises(ValueError):
            read_plot(tmp_path, name)
    target = tmp_path / "fake.png"
    target.write_text("not a PNG")
    with pytest.raises(ValueError, match="PNG"):
        read_plot(tmp_path, "fake.png")


def test_current_export_rejects_missing_recorded_evidence(tmp_path: Path) -> None:
    import hashlib
    from datetime import UTC, datetime

    from ontseq_platform.pipeline.runner import QC_READ_LENGTH_HISTOGRAM, RUN_REPORT
    from ontseq_platform.pipeline.state import ArtifactRecord, RunReport, StageRecord
    from ontseq_platform.service.befund import current_befund

    result = build_demo_result()
    relative = f"normalized/{result.manifest.sample_id}.result.json"
    path = tmp_path / relative
    path.parent.mkdir()
    path.write_text(result.model_dump_json())
    records = [
        ArtifactRecord(
            relative_path=relative,
            size_bytes=path.stat().st_size,
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            exportable=True,
        ),
        ArtifactRecord(
            relative_path=QC_READ_LENGTH_HISTOGRAM, size_bytes=10, sha256="0" * 64, exportable=True
        ),
    ]
    now = datetime.now(UTC)
    run = RunReport(
        run_id=result.manifest.run_id,
        sample_id=result.manifest.sample_id,
        input_kind=result.manifest.input.kind,
        genome_build=result.manifest.assay.genome_build,
        manifest=result.manifest,
        passed=True,
        verdict_reason="Synthetic only",
        pipeline_version="0.8.2",
        git_commit="synthetic",
        started_at=now,
        finished_at=now,
        stages=[
            StageRecord(
                stage="assemble",
                title="Synthetic",
                status="COMPLETED",
                verification="verified_pure_python",
                required=True,
                reason="Synthetic",
                signature="0" * 64,
                outputs=records,
            )
        ],
    )
    (tmp_path / RUN_REPORT).parent.mkdir(parents=True, exist_ok=True)
    (tmp_path / RUN_REPORT).write_text(run.model_dump_json())
    with pytest.raises(ValueError, match="missing|checksum"):
        current_befund(tmp_path, result)

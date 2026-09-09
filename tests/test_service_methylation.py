"""Authenticated local methylation UI integration, with synthetic data only."""

from __future__ import annotations

import hashlib
import http.client
import json
from datetime import UTC, datetime
from pathlib import Path
from unittest.mock import patch

import pytest
from test_methylation import _policy, _row
from test_service_boundaries import _request, _service

from ontseq_platform.demo import build_demo_result
from ontseq_platform.execution import CommandResult, ToolExecutionError
from ontseq_platform.methylation import normalize_methylation
from ontseq_platform.methylation_probe import MethylationProbe, _stat_identity
from ontseq_platform.models import AnalysisModule, ModuleOutcome, ModuleRunStatus, ToolRecord
from ontseq_platform.pipeline.lock import LOCK_FILENAME
from ontseq_platform.pipeline.runner import METHYLATION_REPORT, RESULT_JSON, RUN_REPORT
from ontseq_platform.pipeline.stages import StageId, VerificationStatus
from ontseq_platform.pipeline.state import ArtifactRecord, RunReport, StageRecord
from ontseq_platform.service.app import _build_manifest, _verify_methylation_runtime


def test_probe_requires_token_origin_host_and_allowed_bam(tmp_path: Path) -> None:
    allowed = tmp_path / "allowed"
    allowed.mkdir()
    bam = allowed / "SYNTHETIC.bam"
    bam.write_bytes(b"synthetic placeholder")
    outside = tmp_path / "OUTSIDE.bam"
    outside.write_bytes(b"synthetic placeholder")
    with (
        _service(allowed, allowed / "results") as (config, port),
        patch("ontseq_platform.service.app.probe_bam_methylation") as probe,
    ):
        payload = {"bam_path": str(bam)}
        route = "/api/methylation/probe"
        assert _request(port, "POST", route, body=payload)[0] == 401
        assert (
            _request(port, "POST", route, token=config.token, body={"bam_path": str(outside)})[0]
            == 403
        )
        assert (
            _request(
                port, "POST", route, token=config.token, host=f"evil.invalid:{port}", body=payload
            )[0]
            == 403
        )
        connection = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
        connection.request(
            "POST",
            route,
            json.dumps(payload),
            {
                "X-ONTSeq-Token": config.token,
                "Origin": "https://evil.invalid",
            },
        )
        response = connection.getresponse()
        assert response.status == 403
        response.read()
        connection.close()
        probe.assert_not_called()
        probe.return_value = MethylationProbe(
            "detected", "Synthetic paired tags.", 1, bam_identity=_stat_identity(bam.stat())
        )
        status, body = _request(port, "POST", route, token=config.token, body=payload)
        assert status == 200
        assert json.loads(body)["status"] == "detected"
        assert json.loads(body)["bam_path"] == str(bam.resolve())
        assert "methylation_available" in json.loads(body)
        assert "methylation_unavailable_reason" in json.loads(body)
        probe.assert_called_once()
        assert probe.call_args.args == (bam.resolve(),)
        assert probe.call_args.kwargs["samtools"] == "samtools"
        assert probe.call_args.kwargs["confirmation_offset"] is None
        assert not probe.call_args.kwargs["cancel_event"].is_set()


@pytest.mark.parametrize("include", [None, False, True])
def test_manifest_requests_methylation_only_on_explicit_opt_in(
    tmp_path: Path, include: bool | None
) -> None:
    bam = tmp_path / "SYNTHETIC.bam"
    bam.write_bytes(b"synthetic placeholder")
    Path(f"{bam}.bai").write_bytes(b"synthetic index")
    payload = {
        "bam": str(bam),
        "sample_id": "SYNTHETIC_001",
        "run_id": "SYNTHETIC_RUN",
        "assay": "lcwgs",
        "genome_build": "GRCh38",
    }
    if include is not None:
        payload["include_methylation"] = include
    manifest = _build_manifest(payload, reference_id="SYNTHETIC_REF", allowed_roots=[tmp_path])
    assert (AnalysisModule.METHYLATION in manifest.analysis.modules) == (include is True)


@pytest.mark.parametrize("value", ["false", "true", 1, None, []])
def test_run_api_rejects_ambiguous_opt_in(tmp_path: Path, value: object) -> None:
    with _service(tmp_path, tmp_path / "results") as (config, port):
        status, body = _request(
            port, "POST", "/api/runs", token=config.token, body={"include_methylation": value}
        )
    assert status == 400
    assert b"explicit JSON boolean" in body


@pytest.mark.parametrize("failure", ["unavailable", "wrong_version"])
def test_missing_or_wrong_modkit_is_a_visible_prerequisite_failure(failure: str) -> None:
    with patch("ontseq_platform.service.app.SubprocessRunner.run") as run:
        if failure == "unavailable":
            run.side_effect = ToolExecutionError("synthetic private diagnostic")
        else:
            run.return_value = CommandResult(("modkit", "--version"), 0, "modkit 0.6.1", "")
        with pytest.raises(ValueError, match="Methylation requires") as error:
            _verify_methylation_runtime(_policy(), "modkit")
        assert "synthetic private diagnostic" not in str(error.value)


def _methylation_envelope(root: Path) -> tuple[str, str, Path]:
    result = build_demo_result()
    result.manifest.analysis.modules.append(AnalysisModule.METHYLATION)
    tool = ToolRecord(name="modkit", version="0.4.1")
    result.modules = [item for item in result.modules if item.module != AnalysisModule.METHYLATION]
    result.modules.append(
        ModuleOutcome(
            module=AnalysisModule.METHYLATION,
            status=ModuleRunStatus.COMPLETED,
            reason="Synthetic only",
            tools=[tool],
        )
    )
    run_id, sample_id = result.manifest.run_id, result.manifest.sample_id
    envelope = root / run_id / sample_id
    result_path = envelope / RESULT_JSON.format(sample=sample_id)
    result_path.parent.mkdir(parents=True)
    raw = envelope / "synthetic.bedmethyl"
    raw.write_text(_row("chr1", 100, "m", 20, 15) + "\n", encoding="utf-8")
    methylation = normalize_methylation(
        sample_id=sample_id,
        genome_build=result.manifest.assay.genome_build,
        bedmethyl_path=raw,
        policy=_policy(),
        tool=tool,
    )
    result.provenance.reference_checksums["bedmethyl"] = methylation.bedmethyl_fingerprint.sha256
    result_bytes = result.model_dump_json().encode("utf-8")
    result_path.write_bytes(result_bytes)
    relative = METHYLATION_REPORT.format(sample=sample_id)
    path = envelope / relative
    path.parent.mkdir(parents=True)
    encoded = methylation.model_dump_json().encode("utf-8")
    path.write_bytes(encoded)
    timestamp = datetime.now(UTC)
    run = RunReport(
        run_id=run_id,
        sample_id=sample_id,
        input_kind=result.manifest.input.kind,
        genome_build=result.manifest.assay.genome_build,
        manifest=result.manifest,
        passed=False,
        verdict_reason="Synthetic only",
        pipeline_version="0.7.0-test",
        git_commit="SYNTHETIC",
        started_at=timestamp,
        finished_at=timestamp,
        unverified_stages=[StageId.METHYLATION],
        stages=[
            StageRecord(
                stage=StageId.METHYLATION,
                title="Methylation",
                status=methylation.status,
                verification=VerificationStatus.UNVERIFIED_ADAPTER,
                required=False,
                reason="Synthetic only",
                signature="a" * 64,
                tools=[tool],
                outputs=[
                    ArtifactRecord(
                        relative_path=relative,
                        size_bytes=len(encoded),
                        sha256=hashlib.sha256(encoded).hexdigest(),
                        exportable=True,
                    )
                ],
            ),
            StageRecord(
                stage=StageId.ASSEMBLE,
                title="Assembly",
                status=ModuleRunStatus.COMPLETED,
                verification=VerificationStatus.VERIFIED_PURE_PYTHON,
                required=True,
                reason="Synthetic only",
                signature="b" * 64,
                outputs=[
                    ArtifactRecord(
                        relative_path=RESULT_JSON.format(sample=sample_id),
                        size_bytes=len(result_bytes),
                        sha256=hashlib.sha256(result_bytes).hexdigest(),
                        exportable=True,
                    )
                ],
            ),
        ],
    )
    run_path = envelope / RUN_REPORT
    run_path.parent.mkdir(parents=True)
    run_path.write_text(run.model_dump_json(), encoding="utf-8")
    return run_id, sample_id, path


def test_regional_values_require_identity_and_recorded_checksum(tmp_path: Path) -> None:
    output = tmp_path / "results"
    run_id, sample_id, path = _methylation_envelope(output)
    route = f"/api/methylation?run_id={run_id}&sample_id={sample_id}"
    with _service(tmp_path, output) as (config, port):
        assert _request(port, "GET", route)[0] == 401
        status, body = _request(port, "GET", route, token=config.token)
        assert status == 200, body
        payload = json.loads(body)
        assert payload["run_id"] == run_id and payload["sample_id"] == sample_id
        assert payload["regions"][0]["mean_modified_fraction"] == 0.75
        assert payload["regions"][0]["mean_valid_coverage"] == 20
        assert payload["research_only"] is True
        envelope = output / run_id / sample_id
        lock = envelope / LOCK_FILENAME
        lock.write_text("synthetic lock", encoding="utf-8")
        assert _request(port, "GET", route, token=config.token)[0] == 409
        lock.unlink()
        original = path.read_bytes()
        path.write_bytes(original + b" ")
        status, body = _request(port, "GET", route, token=config.token)
        assert status == 400 and b"checksum" in body
        assert (
            _request(
                port, "GET", f"/api/methylation?run_id=..&sample_id={sample_id}", token=config.token
            )[0]
            == 403
        )


@pytest.mark.parametrize("mismatch", ["sample", "outcome"])
def test_regional_report_rejects_mismatched_typed_provenance(tmp_path: Path, mismatch: str) -> None:
    output = tmp_path / "results"
    run_id, sample_id, path = _methylation_envelope(output)
    envelope = output / run_id / sample_id
    if mismatch == "sample":
        payload = json.loads(path.read_bytes())
        payload["sample_id"] = "OTHER_SYNTHETIC"
        data = json.dumps(payload).encode("utf-8")
        path.write_bytes(data)
        run_path = envelope / RUN_REPORT
        run = json.loads(run_path.read_bytes())
        record = run["stages"][0]["outputs"][0]
        record.update(size_bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
        run_path.write_text(json.dumps(run), encoding="utf-8")
    else:
        result_path = envelope / RESULT_JSON.format(sample=sample_id)
        payload = json.loads(result_path.read_bytes())
        outcome = next(item for item in payload["modules"] if item["module"] == "methylation")
        outcome["status"] = "NOT_RUN"
        data = json.dumps(payload).encode("utf-8")
        result_path.write_bytes(data)
        run_path = envelope / RUN_REPORT
        run = json.loads(run_path.read_bytes())
        record = run["stages"][1]["outputs"][0]
        record.update(size_bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
        run_path.write_text(json.dumps(run), encoding="utf-8")
    with _service(tmp_path, output) as (config, port):
        route = f"/api/methylation?run_id={run_id}&sample_id={sample_id}"
        status, body = _request(port, "GET", route, token=config.token)
    assert status == 400
    assert b"does not match" in body


@pytest.mark.parametrize("state", ["failed", "stale_checksum", "wrong_bedmethyl"])
def test_new_methylation_cannot_be_combined_with_stale_or_unbound_result(
    tmp_path: Path, state: str
) -> None:
    output = tmp_path / "results"
    run_id, sample_id, _ = _methylation_envelope(output)
    envelope = output / run_id / sample_id
    run_path = envelope / RUN_REPORT
    run = json.loads(run_path.read_bytes())
    assembly = run["stages"][1]
    if state == "failed":
        assembly.update(status="FAILED", outputs=[])
    elif state == "stale_checksum":
        assembly["outputs"][0]["sha256"] = "e" * 64
    else:
        result_path = envelope / RESULT_JSON.format(sample=sample_id)
        result = json.loads(result_path.read_bytes())
        result["provenance"]["reference_checksums"]["bedmethyl"] = "f" * 64
        data = json.dumps(result).encode("utf-8")
        result_path.write_bytes(data)
        assembly["outputs"][0].update(size_bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    run_path.write_text(json.dumps(run), encoding="utf-8")
    with _service(tmp_path, output) as (config, port):
        route = f"/api/methylation?run_id={run_id}&sample_id={sample_id}"
        status, body = _request(port, "GET", route, token=config.token)
    assert status == 400
    assert b"assembly" in body or b"source provenance" in body

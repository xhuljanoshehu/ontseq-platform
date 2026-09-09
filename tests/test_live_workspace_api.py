"""Real HTTP boundary tests for the local result workspace (synthetic data only)."""

from __future__ import annotations

import http.client
import json
from pathlib import Path
from unittest.mock import patch

import pytest
from test_service_boundaries import _request, _service

from ontseq_platform.demo import build_demo_result
from ontseq_platform.pipeline.runner import REPORT_HTML, RESULT_JSON
from ontseq_platform.service.app import WORKSPACE_PAGE


def _result(root: Path) -> tuple[str, str, Path]:
    result = build_demo_result()
    run_id = result.manifest.run_id
    sample_id = result.manifest.sample_id
    envelope = root / run_id / sample_id
    path = envelope / RESULT_JSON.format(sample=sample_id)
    path.parent.mkdir(parents=True)
    path.write_text(result.model_dump_json(), encoding="utf-8")
    report = envelope / REPORT_HTML.format(sample=sample_id)
    report.parent.mkdir()
    report.write_text("<!doctype html><p>SYNTHETIC TEST ONLY</p>", encoding="utf-8")
    return run_id, sample_id, envelope


def test_live_results_require_token_and_preserve_run_identity(tmp_path: Path) -> None:
    output = tmp_path / "results"
    run_id, sample_id, _ = _result(output)
    route = f"/api/results?run_id={run_id}&sample_id={sample_id}"
    with _service(tmp_path, output) as (config, port):
        assert _request(port, "GET", route)[0] == 401
        status, body = _request(port, "GET", route, token=config.token)
        assert status == 200
        payload = json.loads(body)
        assert payload["manifest"]["run_id"] == run_id
        assert payload["manifest"]["sample_id"] == sample_id


def test_config_reports_the_exact_desktop_launch_identity(tmp_path: Path) -> None:
    instance_id = "b" * 32
    with _service(tmp_path, tmp_path / "results") as (config, port):
        config.instance_id = instance_id
        status, body = _request(port, "GET", "/api/config", token=config.token)

    assert status == 200
    assert json.loads(body)["instance_id"] == instance_id


@pytest.mark.parametrize("kind", ["html", "json"])
def test_live_artifacts_are_fixed_named_outputs(tmp_path: Path, kind: str) -> None:
    output = tmp_path / "results"
    run_id, sample_id, _ = _result(output)
    with _service(tmp_path, output) as (config, port):
        route = f"/api/artifacts?run_id={run_id}&sample_id={sample_id}&kind={kind}"
        assert _request(port, "GET", route, token=config.token)[0] == 200
        assert _request(port, "GET", route + "-unknown", token=config.token)[0] == 400
        traversal = f"/api/results?run_id=..&sample_id={sample_id}"
        assert _request(port, "GET", traversal, token=config.token)[0] == 403


def test_live_results_reject_swapped_identity(tmp_path: Path) -> None:
    output = tmp_path / "results"
    run_id, sample_id, envelope = _result(output)
    path = envelope / RESULT_JSON.format(sample=sample_id)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["manifest"]["sample_id"] = "OTHER_SYNTHETIC"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with _service(tmp_path, output) as (config, port):
        route = f"/api/results?run_id={run_id}&sample_id={sample_id}"
        assert _request(port, "GET", route, token=config.token)[0] == 400


def test_workspace_token_is_not_exposed_to_rebinding_host(tmp_path: Path) -> None:
    with (
        _service(tmp_path, tmp_path / "results") as (config, port),
        patch.object(Path, "is_file", return_value=True),
        patch.object(Path, "read_text", return_value='<meta content="__ONTSEQ_TOKEN__">'),
    ):
        status, body = _request(port, "GET", "/workspace", host=f"evil.invalid:{port}")
        assert status == 403
        assert config.token.encode() not in body
        status, body = _request(port, "GET", "/workspace")
        assert status == 200
        assert config.token.encode() in body
    assert WORKSPACE_PAGE.name == "workspace.html"


@pytest.mark.parametrize("endpoint", ["results", "artifacts"])
@pytest.mark.parametrize("mismatch", ["build", "profile"])
def test_cross_build_result_cannot_be_exported(
    tmp_path: Path, endpoint: str, mismatch: str
) -> None:
    output = tmp_path / "results"
    run_id, sample_id, envelope = _result(output)
    path = envelope / RESULT_JSON.format(sample=sample_id)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["reference_context"] = {
        "profile_id": (
            payload["manifest"]["analysis"]["profile"] if mismatch == "build" else "OTHER_PROFILE"
        ),
        "profile_version": "test",
        "genome_build": "GRCh37" if mismatch == "build" else "GRCh38",
        "reference_bundle_id": "SYNTHETIC_GRCH37",
        "reference_bundle_version": "test",
        "knowledge_bundle_id": "SYNTHETIC_KNOWLEDGE",
        "knowledge_bundle_version": "test",
        "resource_root": "/synthetic-test-only",
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    with _service(tmp_path, output) as (config, port):
        route = f"/api/{endpoint}?run_id={run_id}&sample_id={sample_id}"
        status, body = _request(port, "GET", route, token=config.token)
        assert status == 400
        assert b"reference context" in body or b"ISCN proposal provenance" in body


def test_workspace_is_not_embeddable_and_is_not_cached(tmp_path: Path) -> None:
    with _service(tmp_path, tmp_path / "results") as (_, port):
        connection = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
        try:
            connection.request("GET", "/workspace")
            response = connection.getresponse()
            assert response.status == 200
            assert response.getheader("X-Frame-Options") == "DENY"
            assert "frame-ancestors 'none'" in response.getheader("Content-Security-Policy", "")
            assert response.getheader("Cache-Control") == "no-store"
            response.read()
        finally:
            connection.close()

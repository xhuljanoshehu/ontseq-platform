"""Transport checks for MARLIN readiness; never accepts model paths from the browser."""

import json
from pathlib import Path

from test_service_boundaries import _request, _service


def test_marlin_readiness_is_authenticated_and_reports_missing_installation(tmp_path: Path):
    with _service(tmp_path, tmp_path / "runs") as (config, port):
        code, _ = _request(port, "GET", "/api/marlin/readiness?genome_build=GRCh38")
        assert code == 401
        code, payload = _request(
            port, "GET", "/api/marlin/readiness?genome_build=GRCh38", token=config.token
        )
        assert code == 200
        data = json.loads(payload)
        assert data["ready"] is False
        assert data["reason"]
        assert data["genome_build"] == "GRCh38"
        assert data["validation_status"] == "UNVALIDATED_RESEARCH"


def test_marlin_readiness_rejects_ambiguous_build_or_external_configuration(tmp_path: Path):
    with _service(tmp_path, tmp_path / "runs") as (config, port):
        for query in (
            "",
            "genome_build=unknown",
            "genome_build=GRCh38&genome_build=GRCh37",
            "genome_build=GRCh38&installation=/tmp/foreign",
        ):
            code, _ = _request(port, "GET", "/api/marlin/readiness?" + query, token=config.token)
            assert code == 400


def test_configured_runtime_does_not_hide_missing_modkit(tmp_path: Path):
    from unittest.mock import patch

    from ontseq_platform.marlin_native_contracts import NativeMarlinReadiness

    with _service(tmp_path, tmp_path / "runs") as (config, port):
        with (
            patch(
                "ontseq_platform.marlin_native.check_native_marlin_readiness",
                return_value=NativeMarlinReadiness(ready=True, reason="Runtime configured"),
            ),
            patch(
                "ontseq_platform.modkit_build.identify_modkit_binary",
                side_effect=FileNotFoundError("Configured modkit missing"),
            ),
        ):
            code, payload = _request(
                port, "GET", "/api/marlin/readiness?genome_build=GRCh38", token=config.token
            )
        value = json.loads(payload)
        assert code == 200 and value["ready"] is False
        assert value["readiness_level"] == "UNAVAILABLE"
        assert "modkit" in value["reason"]

"""Real loopback API checks with only synthetic BAM files and controlled readers."""

from __future__ import annotations

import http.client
import json
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import pytest
from test_service_boundaries import _request, _service

from ontseq_platform.methylation_probe import MethylationProbe, _stat_identity
from ontseq_platform.service.app import Jobs, RunJob
from ontseq_platform.service.methylation_scans import (
    MAX_CACHED_PROBES,
    MAX_RETAINED_SCANS,
    AvailabilityCache,
    MethylationScans,
    ScanBusyError,
    ScanNotFoundError,
)

AVAILABLE = {
    "methylation_available": True,
    "methylation_unavailable_reason": None,
    "expected_modkit_version": "0.4.1",
}


def _bam(root: Path, name: str = "SYNTHETIC.bam") -> Path:
    path = root / name
    path.write_bytes(b"wholly synthetic placeholder")
    return path.resolve()


def _probe(
    bam: Path,
    *,
    status: str = "detected",
    mode: str = "quick",
    checked: int = 12,
    offset: int | None = 1234,
) -> MethylationProbe:
    return MethylationProbe(
        status,
        "Synthetic discovery result.",
        checked_reads=checked,
        complete=status == "not_detected",
        bam_identity=_stat_identity(bam.stat()),
        reason_code="detected_5mc" if status == "detected" else "complete_no_tags",
        scan_mode=mode,
        detected_modifications=("C+m:5mC",) if status == "detected" else (),
        confirmation_offset=offset if status == "detected" else None,
    )


def _json_request(
    port: int, method: str, route: str, token: str, body: dict[str, Any] | None = None
) -> tuple[int, dict[str, Any]]:
    status, data = _request(port, method, route, token=token, body=body)
    return status, json.loads(data)


def _finished(port: int, token: str, scan_id: str) -> dict[str, Any]:
    deadline = time.monotonic() + 3
    while time.monotonic() < deadline:
        status, snapshot = _json_request(port, "GET", f"/api/methylation/scans/{scan_id}", token)
        assert status == 200
        if snapshot["state"] != "running":
            return snapshot
        time.sleep(0.01)
    raise AssertionError("Synthetic scan did not finish")


@pytest.fixture(autouse=True)
def _available_runtime(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "ontseq_platform.service.app._methylation_availability", lambda _: dict(AVAILABLE)
    )


def test_scan_transport_authentication_path_scope_and_identifier_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    allowed = tmp_path / "allowed"
    allowed.mkdir()
    bam = _bam(allowed)
    outside = _bam(tmp_path, "OUTSIDE.bam")
    calls: list[Path] = []
    monkeypatch.setattr(
        "ontseq_platform.service.app.probe_bam_methylation",
        lambda path, **_: calls.append(path) or _probe(path),
    )
    with _service(allowed, allowed / "runs") as (config, port):
        create = "/api/methylation/scans"
        known_shape = create + "/" + "a" * 32
        for method, route, payload in [
            ("POST", create, {"bam_path": str(bam)}),
            ("GET", known_shape, None),
            ("POST", known_shape + "/cancel", {}),
        ]:
            assert _request(port, method, route, body=payload)[0] == 401
            assert (
                _request(
                    port,
                    method,
                    route,
                    token=config.token,
                    host=f"evil.invalid:{port}",
                    body=payload,
                )[0]
                == 403
            )
        assert (
            _json_request(port, "POST", create, config.token, {"bam_path": str(outside)})[0] == 403
        )
        for bad_id in ["..", "invalid", "A" * 32, "a" * 32 + "/extra"]:
            assert _json_request(port, "GET", create + "/" + bad_id, config.token)[0] == 400
            assert (
                _json_request(port, "POST", create + "/" + bad_id + "/cancel", config.token, {})[0]
                == 400
            )
        assert _json_request(port, "GET", known_shape, config.token)[0] == 404
        assert _json_request(port, "POST", known_shape + "/cancel", config.token, {})[0] == 404
        connection = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
        connection.request(
            "POST",
            create,
            json.dumps({"bam_path": str(bam)}),
            {"X-ONTSeq-Token": config.token, "Origin": "https://evil.invalid"},
        )
        response = connection.getresponse()
        assert response.status == 403
        response.read()
        connection.close()
        assert not calls


def test_background_scan_reports_progress_busy_and_cancel_without_caching_old_positive(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bam = _bam(tmp_path)
    started, stopped = threading.Event(), threading.Event()

    def reader(path: Path, **kwargs: Any) -> MethylationProbe:
        if kwargs.get("mode") == "thorough":
            kwargs["progress_callback"]({"checked_reads": 44, "elapsed_seconds": 0.1})
            started.set()
            assert kwargs["cancel_event"].wait(3)
            stopped.set()
            # A racing positive returned after cancellation must never enter the cache.
            return _probe(path, mode="thorough", checked=45)
        return _probe(path, status="not_detected")

    monkeypatch.setattr("ontseq_platform.service.app.probe_bam_methylation", reader)
    with _service(tmp_path, tmp_path / "runs") as (config, port):
        status, initial = _json_request(
            port, "POST", "/api/methylation/scans", config.token, {"bam_path": str(bam)}
        )
        assert status == 202 and initial["state"] == "running" and initial["result"] is None
        assert started.wait(1)
        status, progress = _json_request(
            port, "GET", "/api/methylation/scans/" + initial["scan_id"], config.token
        )
        assert (
            status == 200 and progress["checked_reads"] == 44 and progress["elapsed_seconds"] >= 0
        )
        status, busy = _json_request(
            port, "POST", "/api/methylation/scans", config.token, {"bam_path": str(bam)}
        )
        assert status == 409 and busy["reason_code"] == "busy"
        status, run_busy = _json_request(port, "POST", "/api/runs", config.token, {})
        assert status == 409 and run_busy["reason_code"] == "busy"
        status, busy = _json_request(
            port,
            "POST",
            "/api/methylation/probe",
            config.token,
            {"bam_path": str(bam), "force_refresh": True},
        )
        assert status == 200 and busy["status"] == "unknown" and busy["reason_code"] == "busy"
        status, cancelled = _json_request(
            port,
            "POST",
            "/api/methylation/scans/" + initial["scan_id"] + "/cancel",
            config.token,
            {},
        )
        assert status == 200 and cancelled["state"] in {"running", "cancelled"}
        assert stopped.wait(1)
        cancelled = _finished(port, config.token, initial["scan_id"])
        assert cancelled["state"] == "cancelled"
        assert cancelled["result"]["reason_code"] == "cancelled"
        deadline = time.monotonic() + 1
        while True:
            _, fresh = _json_request(
                port, "POST", "/api/methylation/probe", config.token, {"bam_path": str(bam)}
            )
            if fresh["reason_code"] != "busy":
                break
            assert time.monotonic() < deadline
        assert fresh["status"] == "not_detected"


def test_completed_thorough_result_has_private_witness_and_refresh_rechecks_it(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bam = _bam(tmp_path)
    calls: list[dict[str, Any]] = []

    def reader(path: Path, **kwargs: Any) -> MethylationProbe:
        calls.append(kwargs)
        return _probe(path, mode=kwargs.get("mode", "quick"), checked=20_000)

    monkeypatch.setattr("ontseq_platform.service.app.probe_bam_methylation", reader)
    with _service(tmp_path, tmp_path / "runs") as (config, port):
        _, initial = _json_request(
            port, "POST", "/api/methylation/scans", config.token, {"bam_path": str(bam)}
        )
        finished = _finished(port, config.token, initial["scan_id"])
        assert finished["state"] == "completed"
        result = finished["result"]
        assert result["status"] == "detected" and result["checked_reads"] == 20_000
        assert result["methylation_available"] is True and "confirmation_offset" not in result
        assert result["detected_modifications"] == ["C+m:5mC"]
        assert len(calls) == 1
        _, cached = _json_request(
            port, "POST", "/api/methylation/probe", config.token, {"bam_path": str(bam)}
        )
        assert cached["scan_mode"] == "thorough" and len(calls) == 1
        _, refreshed = _json_request(
            port,
            "POST",
            "/api/methylation/probe",
            config.token,
            {"bam_path": str(bam), "force_refresh": True},
        )
        assert refreshed["status"] == "detected" and len(calls) == 2
        assert calls[-1]["confirmation_offset"] == 1234
        assert "confirmation_offset" not in refreshed
        bam.write_bytes(b"changed synthetic content and size")
        _, changed_snapshot = _json_request(
            port, "GET", "/api/methylation/scans/" + initial["scan_id"], config.token
        )
        assert changed_snapshot["result"]["status"] == "unknown"
        assert changed_snapshot["result"]["reason_code"] == "file_changed"
        _, changed = _json_request(
            port,
            "POST",
            "/api/methylation/probe",
            config.token,
            {"bam_path": str(bam), "force_refresh": True},
        )
        assert len(calls) == 3 and calls[-1]["confirmation_offset"] is None
        assert changed["bam_identity"] != result["bam_identity"]


@pytest.mark.parametrize("value", [None, 0, "true", [], {}])
def test_quick_refresh_flag_requires_boolean(tmp_path: Path, value: object) -> None:
    bam = _bam(tmp_path)
    with _service(tmp_path, tmp_path / "runs") as (config, port):
        assert (
            _json_request(
                port,
                "POST",
                "/api/methylation/probe",
                config.token,
                {"bam_path": str(bam), "force_refresh": value},
            )[0]
            == 400
        )


def test_scan_worker_failure_and_service_shutdown_cleanup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bam = _bam(tmp_path)
    started, stopped = threading.Event(), threading.Event()

    def reader(_: Path, **kwargs: Any) -> MethylationProbe:
        if started.is_set():
            assert kwargs["cancel_event"].wait(3)
            stopped.set()
            raise RuntimeError("Private synthetic detail must not be returned")
        started.set()
        raise RuntimeError("Private synthetic detail must not be returned")

    monkeypatch.setattr("ontseq_platform.service.app.probe_bam_methylation", reader)
    with _service(tmp_path, tmp_path / "runs") as (config, port):
        _, initial = _json_request(
            port, "POST", "/api/methylation/scans", config.token, {"bam_path": str(bam)}
        )
        failed = _finished(port, config.token, initial["scan_id"])
        assert failed["state"] == "failed" and failed["result"]["reason_code"] == "worker_failed"
        assert "Private" not in json.dumps(failed)
        assert (
            _json_request(
                port, "POST", "/api/methylation/scans", config.token, {"bam_path": str(bam)}
            )[0]
            == 202
        )
    assert stopped.is_set(), "Service shutdown must cancel and join its active discovery worker"


def test_active_quick_probe_also_owns_the_thorough_slot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bam = _bam(tmp_path)
    started, release = threading.Event(), threading.Event()

    def reader(path: Path, **_: Any) -> MethylationProbe:
        started.set()
        assert release.wait(3)
        return _probe(path)

    monkeypatch.setattr("ontseq_platform.service.app.probe_bam_methylation", reader)
    with _service(tmp_path, tmp_path / "runs") as (config, port), ThreadPoolExecutor() as pool:
        pending = pool.submit(
            _json_request,
            port,
            "POST",
            "/api/methylation/probe",
            config.token,
            {"bam_path": str(bam)},
        )
        try:
            assert started.wait(1)
            status, response = _json_request(
                port, "POST", "/api/methylation/scans", config.token, {"bam_path": str(bam)}
            )
            assert status == 409 and response["reason_code"] == "busy"
        finally:
            release.set()
        assert pending.result(timeout=3)[1]["status"] == "detected"


def test_changed_file_during_scan_and_failed_witness_cannot_reuse_positive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bam = _bam(tmp_path)
    results = [_probe(bam, mode="thorough")]

    def reader(path: Path, **kwargs: Any) -> MethylationProbe:
        if results:
            return results.pop()
        assert kwargs["confirmation_offset"] == 1234
        return MethylationProbe(
            "unknown",
            "The current witness did not confirm 5mC; bounded sample was incomplete.",
            bam_identity=_stat_identity(path.stat()),
            reason_code="sample_incomplete",
        )

    monkeypatch.setattr("ontseq_platform.service.app.probe_bam_methylation", reader)
    with _service(tmp_path, tmp_path / "runs") as (config, port):
        _, initial = _json_request(
            port, "POST", "/api/methylation/scans", config.token, {"bam_path": str(bam)}
        )
        assert _finished(port, config.token, initial["scan_id"])["result"]["status"] == "detected"
        _, fresh = _json_request(
            port,
            "POST",
            "/api/methylation/probe",
            config.token,
            {"bam_path": str(bam), "force_refresh": True},
        )
        assert fresh["status"] == "unknown" and fresh["reason_code"] == "sample_incomplete"
        # This response is the current cache value, not the earlier positive job record.
        _, cached = _json_request(
            port, "POST", "/api/methylation/probe", config.token, {"bam_path": str(bam)}
        )
        assert cached["status"] == "unknown"

    before = _probe(bam, mode="thorough")

    def changed_reader(path: Path, **_: Any) -> MethylationProbe:
        path.write_bytes(b"modified synthetic input during worker execution")
        return before

    monkeypatch.setattr("ontseq_platform.service.app.probe_bam_methylation", changed_reader)
    with _service(tmp_path, tmp_path / "runs") as (config, port):
        _, initial = _json_request(
            port, "POST", "/api/methylation/scans", config.token, {"bam_path": str(bam)}
        )
        result = _finished(port, config.token, initial["scan_id"])["result"]
        assert result["status"] == "unknown" and result["reason_code"] == "file_changed"


def test_availability_cache_tracks_binary_identity_ttl_and_force_refresh() -> None:
    identity, now, checks = ["binary-a"], [0.0], []

    def check() -> dict[str, Any]:
        checks.append(identity[0])
        return {**AVAILABLE, "methylation_available": identity[0] == "binary-a"}

    cache = AvailabilityCache(check, lambda: identity[0], clock=lambda: now[0])
    assert cache.get()["methylation_available"]
    assert cache.get()["methylation_available"] and len(checks) == 1
    identity[0] = "binary-b"
    assert not cache.get()["methylation_available"] and len(checks) == 2
    now[0] = 11
    cache.get()
    cache.get(force_refresh=True)
    assert len(checks) == 4


def test_analysis_admission_waits_for_cancel_ack_and_scan_close_preserves_analysis(
    tmp_path: Path,
) -> None:
    bam = _bam(tmp_path)
    started, saw_cancel, acknowledge = threading.Event(), threading.Event(), threading.Event()
    jobs = Jobs()

    def reader(path: Path, **kwargs: Any) -> MethylationProbe:
        started.set()
        assert kwargs["cancel_event"].wait(2)
        saw_cancel.set()
        assert acknowledge.wait(2)
        return _probe(path, mode="thorough")

    manager = MethylationScans(
        reader,
        AvailabilityCache(lambda: dict(AVAILABLE), lambda: "synthetic"),
        samtools="unused",
        validate_path=lambda path: path,
        analysis_running=jobs.running,
    )
    job = RunJob("SYNTHETIC_RUN", "SYNTHETIC_SAMPLE", tmp_path / "run")
    try:
        initial = manager.start(bam)
        assert started.wait(1)
        snapshot = manager.cancel(initial["scan_id"])
        assert saw_cancel.wait(1) and snapshot["state"] == "running"
        with pytest.raises(ScanBusyError):
            manager.claim_analysis(lambda: jobs.claim(job))
        assert not jobs.running()
        acknowledge.set()
        deadline = time.monotonic() + 2
        while manager.get(initial["scan_id"])["state"] == "running":
            assert time.monotonic() < deadline
            time.sleep(0.001)
        manager.claim_analysis(lambda: jobs.claim(job))
        assert jobs.running()
        with pytest.raises(ScanBusyError):
            manager.start(bam)
        assert manager.quick(bam, force_refresh=True)["reason_code"] == "busy"
    finally:
        acknowledge.set()
        manager.close()
    assert jobs.get(job.run_id) is job and job.state == "running"


def test_cache_and_job_retention_are_bounded_and_expired_witness_requires_new_read(
    tmp_path: Path,
) -> None:
    now, calls = [0.0], []
    availability = AvailabilityCache(
        lambda: dict(AVAILABLE), lambda: "synthetic", clock=lambda: now[0]
    )

    def reader(path: Path, **kwargs: Any) -> MethylationProbe:
        calls.append(kwargs)
        return _probe(path, mode=kwargs.get("mode", "quick"))

    manager = MethylationScans(
        reader,
        availability,
        samtools="unused",
        validate_path=lambda path: path,
        clock=lambda: now[0],
    )
    try:
        bam = _bam(tmp_path)
        initial = manager.start(bam)
        deadline = time.monotonic() + 2
        while manager.get(initial["scan_id"])["state"] == "running":
            assert time.monotonic() < deadline
            time.sleep(0.001)
        now[0] = 31
        manager.quick(bam, force_refresh=True)
        assert calls[-1]["confirmation_offset"] == 1234, (
            "Expired display cache retains a freshly checked job witness"
        )
        now[0] = 601
        with pytest.raises(ScanNotFoundError):
            manager.get(initial["scan_id"])
        manager.quick(bam, force_refresh=True)
        assert calls[-1]["confirmation_offset"] is None
        for index in range(MAX_CACHED_PROBES + 3):
            manager.quick(_bam(tmp_path, f"SYNTHETIC_{index}.bam"))
        assert len(manager._cache) <= MAX_CACHED_PROBES
        for _ in range(MAX_RETAINED_SCANS + 3):
            snapshot = manager.start(bam)
            deadline = time.monotonic() + 2
            while manager.get(snapshot["scan_id"])["state"] == "running":
                assert time.monotonic() < deadline
                time.sleep(0.001)
        assert len(manager._scans) == MAX_RETAINED_SCANS
    finally:
        manager.close()

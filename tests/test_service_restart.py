"""Restart cannot race a browser analysis or target another service instance."""

import threading
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest
from test_service_boundaries import _request

from ontseq_platform.service.app import (
    JobRejected,
    Jobs,
    RunJob,
    ServiceConfig,
    close_handler_resources,
    make_handler,
)


@contextmanager
def service(root):
    config = ServiceConfig(
        reference_lock=None,
        qc_policy=root / "qc",
        sniffles_policy=root / "sv",
        target_coverage_policy=root / "coverage",
        output_dir=root / "runs",
        allowed_roots=[root],
        port=0,
        instance_id="a" * 32,
    )
    jobs = Jobs()
    handler = make_handler(config, jobs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    config.port = server.server_address[1]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield config, jobs, thread
    finally:
        server.shutdown()
        close_handler_resources(handler)
        server.server_close()
        thread.join(5)


def test_restart_requires_token_and_matching_identity(tmp_path: Path):
    with service(tmp_path) as (config, jobs, thread):
        for token, instance, expected in [
            (None, config.instance_id, 401),
            (config.token, "b" * 32, 409),
        ]:
            status, _ = _request(
                config.port,
                "POST",
                "/api/session/stop",
                token=token,
                body={"instance_id": instance},
            )
            assert status == expected
            assert thread.is_alive()
        jobs.claim(RunJob("SYN_RUN", "SYN_SAMPLE", Path("/synthetic")))


def test_restart_refuses_active_analysis_and_keeps_service_alive(tmp_path: Path):
    with service(tmp_path) as (config, jobs, thread):
        job = RunJob("SYN_RUN", "SYN_SAMPLE", Path("/synthetic"))
        jobs.claim(job)
        status, _ = _request(
            config.port,
            "POST",
            "/api/session/stop",
            token=config.token,
            body={"instance_id": config.instance_id},
        )
        assert status == 409
        assert jobs.running()
        assert thread.is_alive()
        job.state = "completed"
        status, _ = _request(
            config.port,
            "POST",
            "/api/session/stop",
            token=config.token,
            body={"instance_id": config.instance_id},
        )
        assert status == 200
        thread.join(3)
        assert not thread.is_alive()
        with pytest.raises(JobRejected):
            jobs.claim(RunJob("LATE_RUN", "LATE_SAMPLE", Path("/synthetic")))


def test_shutdown_and_run_claim_have_exactly_one_winner():
    for _ in range(40):
        jobs = Jobs()
        barrier = threading.Barrier(2)

        def claim(jobs=jobs, barrier=barrier):
            barrier.wait()
            try:
                jobs.claim(RunJob("SYN_RUN", "SYN_SAMPLE", Path("/synthetic")))
                return True
            except JobRejected:
                return False

        def stop(jobs=jobs, barrier=barrier):
            barrier.wait()
            try:
                jobs.begin_shutdown()
                return True
            except JobRejected:
                return False

        with ThreadPoolExecutor(2) as pool:
            a, b = pool.submit(claim), pool.submit(stop)
            assert a.result() != b.result()


def test_restart_waits_for_probe_cleanup(tmp_path: Path, monkeypatch):
    from test_service_methylation_scans import _probe

    started = threading.Event()
    released = threading.Event()

    def probe(bam, **kwargs):
        started.set()
        assert kwargs["cancel_event"].wait(3)
        released.set()
        return _probe(bam)

    monkeypatch.setattr("ontseq_platform.service.app.probe_bam_methylation", probe)
    monkeypatch.setattr("ontseq_platform.service.app._methylation_availability", lambda _: {})
    bam = tmp_path / "SYNTHETIC.bam"
    bam.write_bytes(b"synthetic")
    with service(tmp_path) as (config, jobs, thread):
        status, _ = _request(
            config.port,
            "POST",
            "/api/methylation/scans",
            token=config.token,
            body={"bam_path": str(bam)},
        )
        assert status == 202
        assert started.wait(3)
        status, _ = _request(
            config.port,
            "POST",
            "/api/session/stop",
            token=config.token,
            body={"instance_id": config.instance_id},
        )
        assert status == 200
        assert released.is_set()
        thread.join(3)
        assert not thread.is_alive()


def test_restart_refuses_unacknowledged_probe_cleanup(tmp_path, monkeypatch):
    # Simulate the real manager exhausting its bounded close wait. It must not
    # acknowledge a safe stop while a worker still owns input handles.
    monkeypatch.setattr("ontseq_platform.service.app.MethylationScans.close", lambda _: None)
    monkeypatch.setattr(
        "ontseq_platform.service.app.MethylationScans.discovery_active", lambda _: True
    )
    with service(tmp_path) as (config, jobs, thread):
        status, _ = _request(
            config.port,
            "POST",
            "/api/session/stop",
            token=config.token,
            body={"instance_id": config.instance_id},
        )
        assert status == 409
        assert thread.is_alive()
        with pytest.raises(JobRejected):
            jobs.claim(RunJob("LATE", "SYN_SAMPLE", tmp_path))

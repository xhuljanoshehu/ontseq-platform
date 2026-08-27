"""Boundaries the local service has to hold at the transport layer.

The guard functions are unit tested in ``test_service_guard``. What is checked here is that
the handler actually calls them — a correct check that a route forgets to make protects
nothing. Each test drives a real server over a real socket, because the things being checked
(the ``Host`` header, two requests arriving at once, a directory entry that cannot be
stat'ed) do not exist at the level of a function call.
"""

from __future__ import annotations

import http.client
import json
import os
import tempfile
import threading
import unittest
from collections.abc import Iterator
from contextlib import contextmanager
from http.server import ThreadingHTTPServer
from pathlib import Path

from ontseq_platform.service.app import JobRejected, Jobs, RunJob, ServiceConfig, make_handler
from ontseq_platform.service.guard import TOKEN_HEADER


@contextmanager
def _service(root: Path, output_dir: Path) -> Iterator[tuple[ServiceConfig, int]]:
    """Start the handler on an ephemeral port and stop it again."""
    config = ServiceConfig(
        # Only the routes exercised here are reached, and none of them loads a policy.
        reference_lock=root / "reference.lock.json",
        output_dir=output_dir,
        allowed_roots=[root],
        qc_policy=root / "qc.yaml",
        sniffles_policy=root / "sv.yaml",
        target_coverage_policy=root / "coverage.yaml",
        port=0,
    )
    server = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(config, Jobs()))
    port = server.server_address[1]
    assert isinstance(port, int)
    config.port = port
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield config, port
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def _request(
    port: int,
    method: str,
    path: str,
    *,
    token: str | None = None,
    host: str | None = None,
    body: dict[str, object] | None = None,
) -> tuple[int, bytes]:
    connection = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    headers: dict[str, str] = {"Host": host if host is not None else f"127.0.0.1:{port}"}
    if token is not None:
        headers[TOKEN_HEADER] = token
    payload = None
    if body is not None:
        payload = json.dumps(body).encode("utf-8")
        headers["Content-Type"] = "application/json"
    try:
        connection.request(method, path, body=payload, headers=headers)
        response = connection.getresponse()
        return response.status, response.read()
    finally:
        connection.close()


class ReviewRouteTests(unittest.TestCase):
    """``/api/review/<run>/<sample>`` builds a filesystem path out of the request URL."""

    def test_dot_dot_in_the_route_cannot_reach_a_sibling_envelope(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            output_dir = base / "runs"
            output_dir.mkdir()
            # A complete, reviewable envelope that the service was never pointed at.
            outside = base / "SAMPLE_001"
            (outside / "release").mkdir(parents=True)
            (outside / "release" / "release.json").write_text("{}", encoding="utf-8")

            with _service(base, output_dir) as (config, port):
                status, body = _request(
                    port,
                    "POST",
                    "/api/review/../SAMPLE_001",
                    token=config.token,
                    body={"decision": "accepted", "reviewer": "someone"},
                )

            self.assertEqual(status, 403, body)
            self.assertFalse(
                (outside / "review").exists(),
                "a review was written into an envelope outside the output directory",
            )

    def test_a_well_formed_route_for_a_missing_envelope_is_a_plain_bad_request(self) -> None:
        """Refusing traversal must not turn every unknown envelope into a 403."""
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            output_dir = base / "runs"
            output_dir.mkdir()
            with _service(base, output_dir) as (config, port):
                status, body = _request(
                    port,
                    "POST",
                    "/api/review/RUN_404/SAMPLE_404",
                    token=config.token,
                    body={"decision": "accepted", "reviewer": "someone"},
                )
            self.assertEqual(status, 400, body)


class PageRouteTests(unittest.TestCase):
    def test_the_page_is_served_to_a_loopback_host(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            with _service(base, base / "runs") as (config, port):
                status, body = _request(port, "GET", "/")
            self.assertEqual(status, 200)
            self.assertIn(config.token.encode("utf-8"), body)

    def test_a_rebound_hostname_is_not_handed_the_session_token(self) -> None:
        """The page carries the token, so it gets the rebinding check the API gets."""
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            with _service(base, base / "runs") as (config, port):
                status, body = _request(port, "GET", "/", host=f"evil.example:{port}")
            self.assertEqual(status, 403)
            self.assertNotIn(config.token.encode("utf-8"), body)


class BrowseRouteTests(unittest.TestCase):
    def test_a_dangling_symlink_does_not_break_the_listing(self) -> None:
        """One unreadable entry must not cost the operator the whole directory."""
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            (base / "runs").mkdir()
            real = base / "real.bam"
            real.write_bytes(b"BAM")
            Path(f"{real}.bai").write_bytes(b"BAI")
            os.symlink(base / "absent.bam", base / "dangling.bam")

            with _service(base, base / "runs") as (config, port):
                status, body = _request(port, "GET", "/api/browse", token=config.token)

            self.assertEqual(status, 200, body)
            listing = json.loads(body)
            names = [item["name"] for item in listing["entries"]]
            self.assertIn("real.bam", names)
            self.assertNotIn("dangling.bam", names)
            self.assertEqual(listing["unreadable_entries"], 1)

    def test_a_clean_directory_reports_nothing_unreadable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            (base / "runs").mkdir()
            (base / "real.bam").write_bytes(b"BAM")
            with _service(base, base / "runs") as (config, port):
                status, body = _request(port, "GET", "/api/browse", token=config.token)
            self.assertEqual(status, 200, body)
            self.assertEqual(json.loads(body)["unreadable_entries"], 0)


class LocateRouteTests(unittest.TestCase):
    """A name lookup must not follow a BAM symlink outside an allowed root."""

    def test_a_regular_bam_inside_the_root_is_returned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            allowed = base / "allowed"
            allowed.mkdir()
            (allowed / "real.bam").write_bytes(b"BAM")

            with _service(allowed, allowed / "runs") as (config, port):
                status, body = _request(
                    port,
                    "GET",
                    "/api/locate?name=real.bam",
                    token=config.token,
                )

            self.assertEqual(status, 200, body)
            payload = json.loads(body)
            self.assertEqual(len(payload["matches"]), 1)
            self.assertEqual(payload["matches"][0]["posix"], str(allowed / "real.bam"))

    def test_a_bam_symlink_outside_the_root_is_not_returned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            allowed = base / "allowed"
            allowed.mkdir()
            outside = base / "outside.bam"
            outside.write_bytes(b"BAM outside the service boundary")
            os.symlink(outside, allowed / "escaped.bam")

            with _service(allowed, allowed / "runs") as (config, port):
                status, body = _request(
                    port,
                    "GET",
                    "/api/locate?name=escaped.bam",
                    token=config.token,
                )

            self.assertEqual(status, 200, body)
            self.assertEqual(json.loads(body)["matches"], [])


class JobClaimTests(unittest.TestCase):
    """One analysis at a time, decided under the lock rather than between two calls."""

    @staticmethod
    def _job(run_id: str) -> RunJob:
        return RunJob(run_id=run_id, sample_id="SAMPLE_001", envelope=Path("envelope"))

    def test_a_second_run_is_refused_while_one_is_running(self) -> None:
        jobs = Jobs()
        jobs.claim(self._job("RUN_001"))
        with self.assertRaises(JobRejected):
            jobs.claim(self._job("RUN_002"))

    def test_a_run_may_start_once_the_previous_one_has_finished(self) -> None:
        jobs = Jobs()
        first = self._job("RUN_001")
        jobs.claim(first)
        first.state = "passed"
        second = self._job("RUN_002")
        jobs.claim(second)
        self.assertIs(jobs.get("RUN_002"), second)
        self.assertIs(jobs.get("RUN_001"), first)

    def test_reusing_a_finished_run_id_is_refused_rather_than_losing_the_record(self) -> None:
        jobs = Jobs()
        first = self._job("RUN_001")
        jobs.claim(first)
        first.state = "passed"
        with self.assertRaises(JobRejected):
            jobs.claim(self._job("RUN_001"))
        self.assertIs(jobs.get("RUN_001"), first)

    def test_only_one_of_many_simultaneous_claims_wins(self) -> None:
        """The check and the registration are one step, so a tie cannot start two runs.

        Asking ``running()`` and then calling ``add()`` passes every single-threaded test
        and still lets two POSTs that arrive together both start a pipeline.
        """
        jobs = Jobs()
        attempts = 16
        start = threading.Barrier(attempts)
        accepted: list[str] = []
        refused: list[str] = []
        guard = threading.Lock()

        def attempt(index: int) -> None:
            start.wait(timeout=10)
            try:
                jobs.claim(self._job(f"RUN_{index:03d}"))
            except JobRejected:
                with guard:
                    refused.append(f"RUN_{index:03d}")
            else:
                with guard:
                    accepted.append(f"RUN_{index:03d}")

        threads = [threading.Thread(target=attempt, args=(index,)) for index in range(attempts)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=10)

        self.assertEqual(len(accepted), 1, f"accepted {accepted}")
        self.assertEqual(len(refused), attempts - 1)


if __name__ == "__main__":
    unittest.main()

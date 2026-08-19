from __future__ import annotations

import http.client
import json
import threading
import unittest
from collections.abc import Iterator
from contextlib import contextmanager
from http import HTTPStatus
from http.server import ThreadingHTTPServer
from pathlib import Path
from typing import Any

from ontseq_platform.service.app import Jobs, ServiceConfig, make_handler
from ontseq_platform.service.guard import TOKEN_HEADER


@contextmanager
def recording_service() -> Iterator[
    tuple[ServiceConfig, ThreadingHTTPServer, list[dict[str, Any]]]
]:
    """Run the real HTTP handler while replacing only the analytical run launch."""
    config = ServiceConfig(
        reference_lock=Path("synthetic-reference-lock.json"),
        output_dir=Path("synthetic-output"),
        allowed_roots=[Path("synthetic-input")],
        qc_policy=Path("synthetic-qc-policy.yaml"),
        sniffles_policy=Path("synthetic-sniffles-policy.yaml"),
        port=0,
    )
    received: list[dict[str, Any]] = []
    handler = make_handler(config, Jobs())

    def record_start(self: Any, payload: dict[str, Any]) -> None:
        if not str(payload.get("bam", "")).strip():
            self._refuse(HTTPStatus.FORBIDDEN, "no BAM was selected")
            return
        received.append(
            {
                "payload": payload,
                "content_length": self.headers.get("Content-Length"),
                "transfer_encoding": self.headers.get("Transfer-Encoding"),
            }
        )
        self._json(HTTPStatus.ACCEPTED, {"state": "recorded"})

    handler._start_run = record_start
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    config.port = server.server_port
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield config, server, received
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


class ServiceRequestBodyTests(unittest.TestCase):
    payload = {
        "bam": r"C:\synthetic\sample.bam",
        "sample_id": "synthetic-sample",
        "run_id": "synthetic-run",
        "genome_build": "GRCh38",
        "assay": "lcwgs",
    }

    def send_chunked(
        self,
        config: ServiceConfig,
        server: ThreadingHTTPServer,
    ) -> tuple[int, bytes]:
        encoded = json.dumps(self.payload).encode("utf-8")
        headers = {
            TOKEN_HEADER: config.token,
            "Content-Type": "application/json",
        }
        connection = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=5)
        try:
            connection.request(
                "POST",
                "/api/runs",
                body=iter((encoded[:11], encoded[11:])),
                headers=headers,
                encode_chunked=True,
            )
            response = connection.getresponse()
            return response.status, response.read()
        finally:
            connection.close()

    def send_chunked_headers(
        self,
        server: ThreadingHTTPServer,
        *,
        token: str,
        origin: str | None = None,
    ) -> tuple[int, bytes]:
        """Send chunked headers only; an auth refusal must not wait for the body."""
        connection = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=5)
        try:
            connection.putrequest("POST", "/api/runs")
            connection.putheader("Transfer-Encoding", "chunked")
            connection.putheader("Content-Type", "application/json")
            connection.putheader(TOKEN_HEADER, token)
            if origin is not None:
                connection.putheader("Origin", origin)
            connection.endheaders()
            response = connection.getresponse()
            return response.status, response.read()
        finally:
            connection.close()

    def test_dotnet_style_chunked_json_reaches_the_run_handler(self) -> None:
        with recording_service() as (config, server, received):
            # An iterable body makes http.client use the same HTTP/1.1 chunked framing
            # observed from .NET PostAsJsonAsync when content length is not known.
            status, response_body = self.send_chunked(config, server)

        self.assertEqual(status, HTTPStatus.ACCEPTED, response_body)
        self.assertEqual(received[0]["payload"], self.payload)
        self.assertIsNone(received[0]["content_length"])
        self.assertEqual(received[0]["transfer_encoding"], "chunked")

    def test_content_length_json_remains_supported(self) -> None:
        encoded = json.dumps(self.payload).encode("utf-8")
        with recording_service() as (config, server, received):
            connection = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=5)
            try:
                connection.request(
                    "POST",
                    "/api/runs",
                    body=encoded,
                    headers={TOKEN_HEADER: config.token, "Content-Type": "application/json"},
                )
                response = connection.getresponse()
                response_body = response.read()
            finally:
                connection.close()

        self.assertEqual(response.status, HTTPStatus.ACCEPTED, response_body)
        self.assertEqual(received[0]["payload"], self.payload)
        self.assertEqual(received[0]["content_length"], str(len(encoded)))
        self.assertIsNone(received[0]["transfer_encoding"])

    def test_chunked_body_does_not_bypass_token_check(self) -> None:
        with recording_service() as (config, server, received):
            status, _body = self.send_chunked_headers(server, token="wrong")

        self.assertEqual(status, HTTPStatus.UNAUTHORIZED)
        self.assertEqual(received, [])

    def test_chunked_body_does_not_bypass_origin_check(self) -> None:
        with recording_service() as (config, server, received):
            status, _body = self.send_chunked_headers(
                server,
                token=config.token,
                origin="https://evil.example",
            )

        self.assertEqual(status, HTTPStatus.FORBIDDEN)
        self.assertEqual(received, [])


if __name__ == "__main__":
    unittest.main()

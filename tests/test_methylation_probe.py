"""Cross-platform worker lifecycle, private protocol and fail-closed result tests."""

from __future__ import annotations

import io
import json
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from ontseq_platform.methylation_probe import (
    MAX_WORKER_MESSAGE_BYTES,
    MethylationProbe,
    probe_bam_methylation,
)


class _Process:
    def __init__(self, output: bytes = b"", returncode: int = 0) -> None:
        self.stdout = io.BytesIO(output)
        self.stdin = io.BytesIO()
        self.returncode = returncode
        self.killed = False
        self.waited = False

    def wait(self, timeout: float | None = None) -> int:
        self.waited = True
        return self.returncode

    def poll(self) -> int | None:
        return self.returncode if self.waited or self.killed else None

    def kill(self) -> None:
        self.killed = True


def _message(code: str, *, complete: bool = False, offset: int | None = None) -> bytes:
    return (
        json.dumps(
            {
                "type": "result",
                "reason_code": code,
                "checked_reads": 2,
                "complete": complete,
                "confirmation_offset": offset,
            }
        ).encode()
        + b"\n"
    )


@pytest.fixture
def bam(tmp_path: Path) -> Path:
    path = tmp_path / "SYNTHETIC.bam"
    path.write_bytes(b"synthetic placeholder; worker is mocked")
    return path


def test_old_positional_result_contract_remains_valid() -> None:
    probe = MethylationProbe("unknown", "Synthetic", 2, False)
    assert probe.probe_version == "mm-ml-presence-v2"
    assert probe.reader == "pysam" and probe.confirmation_offset is None


@pytest.mark.parametrize(
    ("code", "complete", "offset", "status"),
    [
        ("detected_5mc", False, 1234, "detected"),
        ("complete_no_tags", True, None, "not_detected"),
        ("complete_no_tags", False, None, "unknown"),
        ("sample_incomplete", False, None, "unknown"),
        ("reader_unavailable", False, None, "unknown"),
        ("unsupported_modification", True, None, "unknown"),
        ("invalid_tags", True, None, "unknown"),
    ],
)
def test_protocol_results_are_typed_private_and_cleaned_up(
    bam: Path,
    code: str,
    complete: bool,
    offset: int | None,
    status: str,
) -> None:
    process = _Process(_message(code, complete=complete, offset=offset))
    with patch("ontseq_platform.methylation_probe.subprocess.Popen", return_value=process) as start:
        result = probe_bam_methylation(bam, samtools="deliberately unavailable legacy tool")
    assert result.status == status and result.complete == complete
    assert result.reason_code == code
    assert result.bam_identity and len(result.bam_identity) == 64
    assert result.confirmation_offset == offset
    assert str(bam) not in repr(result)
    assert str(bam) not in start.call_args.args[0]
    assert "samtools" not in start.call_args.args[0]
    assert start.call_args.kwargs["shell"] is False
    assert start.call_args.kwargs["stderr"] == subprocess.DEVNULL
    assert process.waited and process.stdout.closed and process.stdin.closed


@pytest.mark.parametrize(
    "output",
    [
        b"SYNTHETIC_PRIVATE_READ\n",
        b"x" * (MAX_WORKER_MESSAGE_BYTES + 1) + b"\n",
        _message("SYNTHETIC_PRIVATE_READ"),
        _message("detected_5mc"),
        b"",
    ],
)
def test_invalid_worker_output_never_leaks_raw_diagnostics(bam: Path, output: bytes) -> None:
    process = _Process(output)
    with patch("ontseq_platform.methylation_probe.subprocess.Popen", return_value=process):
        result = probe_bam_methylation(bam)
    assert result.status == "unknown" and result.reason_code == "worker_failed"
    assert "SYNTHETIC_PRIVATE_READ" not in repr(result)
    assert process.waited and process.stdout.closed


def test_changed_file_invalidates_positive_witness(bam: Path) -> None:
    process = _Process(_message("detected_5mc", offset=1234))
    original_wait = process.wait

    def wait(timeout: float | None = None) -> int:
        bam.write_bytes(b"changed synthetic placeholder")
        return original_wait(timeout)

    with (
        patch("ontseq_platform.methylation_probe.subprocess.Popen", return_value=process),
        patch.object(process, "wait", side_effect=wait),
    ):
        result = probe_bam_methylation(bam)
    assert result.reason_code == "file_changed"
    assert result.bam_identity is None and result.confirmation_offset is None


@pytest.mark.parametrize("action", ["timeout", "cancelled"])
def test_blocking_reader_is_killed_on_deadline_or_cancel(bam: Path, action: str) -> None:
    stopped = threading.Event()
    cancelled = threading.Event()
    process = _Process()

    class BlockingStream(io.BytesIO):
        def readline(self, size: int | None = -1) -> bytes:
            assert stopped.wait(3)
            return b""

    process.stdout = BlockingStream()
    original_kill = process.kill

    def kill() -> None:
        original_kill()
        stopped.set()

    timer = threading.Timer(0.02, cancelled.set)
    if action == "cancelled":
        timer.start()
    try:
        with (
            patch("ontseq_platform.methylation_probe.subprocess.Popen", return_value=process),
            patch.object(process, "kill", side_effect=kill),
        ):
            result = probe_bam_methylation(
                bam,
                timeout_seconds=0.03 if action == "timeout" else 2,
                cancel_event=cancelled,
            )
    finally:
        timer.cancel()
    assert result.reason_code == action and result.elapsed_seconds < 1
    assert process.killed and process.waited and process.stdout.closed


def test_progress_contains_only_safe_counts_and_elapsed(bam: Path) -> None:
    output = b'{"type":"progress","checked_reads":1,"secret":"PRIVATE"}\n'
    process = _Process(output + _message("complete_no_tags", complete=True))
    progress: list[dict[str, object]] = []
    with patch("ontseq_platform.methylation_probe.subprocess.Popen", return_value=process):
        probe_bam_methylation(bam, mode="thorough", progress_callback=progress.append)
    assert set(progress[0]) == {"checked_reads", "elapsed_seconds", "scan_mode", "reader"}
    assert "PRIVATE" not in repr(progress)


def test_unavailable_worker_and_precancel_are_typed(bam: Path) -> None:
    with patch(
        "ontseq_platform.methylation_probe.subprocess.Popen", side_effect=OSError("PRIVATE")
    ):
        assert probe_bam_methylation(bam).reason_code == "reader_unavailable"
    cancelled = threading.Event()
    cancelled.set()
    with patch("ontseq_platform.methylation_probe.subprocess.Popen") as start:
        assert probe_bam_methylation(bam, cancel_event=cancelled).reason_code == "cancelled"
        start.assert_not_called()


@pytest.mark.parametrize("action", ["timeout", "cancelled"])
def test_real_isolated_process_is_reaped_on_deadline_or_cancel(bam: Path, action: str) -> None:
    original_start = subprocess.Popen
    children: list[subprocess.Popen[bytes]] = []
    cancelled = threading.Event()

    def start(_command: list[str], **kwargs: Any) -> subprocess.Popen[bytes]:
        child = original_start([sys.executable, "-c", "import time; time.sleep(30)"], **kwargs)
        children.append(child)
        return child

    timer = threading.Timer(0.05, cancelled.set)
    if action == "cancelled":
        timer.start()
    try:
        with patch("ontseq_platform.methylation_probe.subprocess.Popen", side_effect=start):
            result = probe_bam_methylation(
                bam,
                timeout_seconds=0.05 if action == "timeout" else 2,
                cancel_event=cancelled,
            )
    finally:
        timer.cancel()
        for child in children:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=2)
    assert result.reason_code == action and result.elapsed_seconds < 1
    assert len(children) == 1 and children[0].poll() is not None
    assert children[0].stdout is not None and children[0].stdout.closed


def test_native_memory_guard_can_be_installed_in_an_isolated_child() -> None:
    child = subprocess.run(
        [
            sys.executable,
            "-c",
            "from ontseq_platform._methylation_worker import _memory_limit; assert _memory_limit()",
        ],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        timeout=10,
        check=False,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
    )
    assert child.returncode == 0, "isolated memory limit could not be installed"


def test_failing_progress_callback_is_typed_and_cleans_up(bam: Path) -> None:
    process = _Process(b'{"type":"progress","checked_reads":1}\n')

    def failure(_progress: dict[str, Any]) -> None:
        raise RuntimeError("SYNTHETIC_PRIVATE_ERROR")

    with patch("ontseq_platform.methylation_probe.subprocess.Popen", return_value=process):
        result = probe_bam_methylation(bam, progress_callback=failure)
    assert result.reason_code == "worker_failed" and "PRIVATE" not in repr(result)
    assert process.waited and process.stdout.closed

"""Cancellable, isolated BAM modification discovery with aggregate-only diagnostics."""

from __future__ import annotations

import hashlib
import json
import os
import queue
import subprocess
import sys
import threading
import time
from collections.abc import Callable
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

PROBE_VERSION = "mm-ml-presence-v2"
MAX_PROBE_READS = 10_000
MAX_ALIGNMENT_LINE_BYTES = 262_144  # Compatibility only: v2 never reads SAM lines.
PROBE_TIMEOUT_SECONDS = 3.0
MAX_WORKER_MESSAGE_BYTES = 4096

_REASONS = {
    "detected_5mc": "Valid paired MM/ML annotations contain supported cytosine 5mC calls.",
    "complete_no_tags": "The complete BAM scan contained no modified-base tags.",
    "sample_incomplete": "The bounded sample contained no supported 5mC calls; scan incomplete.",
    "timeout": "Tag discovery reached its time limit.",
    "cancelled": "Tag discovery was cancelled.",
    "reader_unavailable": "The direct BAM reader pysam is unavailable in this runtime.",
    "input_unavailable": "The selected BAM is unavailable for tag discovery.",
    "file_changed": "The BAM changed during tag discovery.",
    "record_limit": "A BAM record exceeded the bounded reader memory or record limits.",
    "invalid_tags": "Modified-base tags were incomplete or invalid.",
    "unsupported_modification": "Modified-base tags contain no supported cytosine 5mC calls.",
    "read_error": "The direct BAM reader could not complete the scan.",
    "worker_failed": "The isolated BAM reader could not complete the scan.",
}


@dataclass(frozen=True)
class MethylationProbe:
    status: Literal["detected", "not_detected", "unknown"]
    reason: str
    checked_reads: int = 0
    complete: bool = False
    probe_version: str = PROBE_VERSION
    bam_identity: str | None = None
    bam_identity_kind: Literal["stat-fingerprint-v1"] = "stat-fingerprint-v1"
    reason_code: str = ""
    elapsed_seconds: float = 0.0
    scan_mode: str = "quick"
    reader: str = "pysam"
    detected_modifications: tuple[str, ...] = ()
    # Private cache witness. API serializers MUST omit this field.
    confirmation_offset: int | None = None


def _stat_identity(info: os.stat_result) -> str:
    """Change detection only; this is not a BAM content checksum."""
    fields = (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns)
    return hashlib.sha256(b"stat-fingerprint-v1|" + "|".join(map(str, fields)).encode()).hexdigest()


def _read_messages(process: subprocess.Popen[bytes], messages: queue.Queue[dict[str, Any]]) -> None:
    """Drain a bounded, trusted worker protocol; never relay arbitrary worker strings."""
    assert process.stdout is not None
    try:
        while line := process.stdout.readline(MAX_WORKER_MESSAGE_BYTES + 1):
            if len(line) > MAX_WORKER_MESSAGE_BYTES or not line.endswith(b"\n"):
                raise ValueError("invalid worker message")
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError("invalid worker message")
            messages.put(value, timeout=1)
    except (OSError, ValueError, queue.Full):
        with suppress(queue.Full):
            messages.put_nowait({"type": "invalid"})
    finally:
        with suppress(queue.Full):
            messages.put_nowait({"type": "closed"})


def probe_bam_methylation(
    bam: Path,
    *,
    samtools: str = "samtools",
    maximum_reads: int = MAX_PROBE_READS,
    maximum_line_bytes: int = MAX_ALIGNMENT_LINE_BYTES,
    timeout_seconds: float = PROBE_TIMEOUT_SECONDS,
    mode: Literal["quick", "thorough"] = "quick",
    cancel_event: threading.Event | None = None,
    progress_callback: Callable[[dict[str, Any]], None] | None = None,
    confirmation_offset: int | None = None,
) -> MethylationProbe:
    """Discover supported MM/ML annotations without exposing alignment contents.

    The caller must apply its allowed-input-path gate first. Quick mode samples indexed
    BAM positions, falling back to a bounded sequential scan. Thorough mode ignores the
    read-count limit and streams until EOF, cancellation or the caller's deadline.
    Legacy samtools/line-size arguments remain accepted but are no longer used.
    Progress callbacks receive only counts, elapsed time, mode and reader identity.
    """
    del samtools
    if mode not in {"quick", "thorough"}:
        raise ValueError("Unknown methylation scan mode")
    if maximum_reads < 1 or maximum_line_bytes < 1 or not 0 < timeout_seconds <= 1800:
        raise ValueError("Methylation probe limits must be positive and at most 1800 seconds")
    if confirmation_offset is not None and (
        type(confirmation_offset) is not int or not 0 <= confirmation_offset < 2**63
    ):
        raise ValueError("Invalid BAM confirmation offset")
    started = time.monotonic()
    checked = 0
    identity: str | None = None

    def result(code: str, *, complete: bool = False, offset: int | None = None) -> MethylationProbe:
        nonlocal identity
        if identity is not None:
            try:
                current = _stat_identity(bam.stat())
            except OSError:
                current = None
            if current != identity:
                code, complete, offset, identity = "file_changed", False, None, None
        status: Literal["detected", "not_detected", "unknown"] = "unknown"
        if code == "detected_5mc":
            status = "detected"
        elif code == "complete_no_tags" and complete:
            status = "not_detected"
        return MethylationProbe(
            status,
            _REASONS[code],
            checked,
            complete,
            bam_identity=identity,
            reason_code=code,
            elapsed_seconds=round(time.monotonic() - started, 3),
            scan_mode=mode,
            detected_modifications=("C+m:5mC",) if status == "detected" else (),
            confirmation_offset=offset if status == "detected" else None,
        )

    if cancel_event is not None and cancel_event.is_set():
        return result("cancelled")
    try:
        if not bam.is_file():
            return result("input_unavailable")
        identity = _stat_identity(bam.stat())
        process = subprocess.Popen(
            [sys.executable, str(Path(__file__).with_name("_methylation_worker.py"))],
            shell=False,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0,
        )
    except (OSError, ValueError):
        return result("reader_unavailable")
    messages: queue.Queue[dict[str, Any]] = queue.Queue(maxsize=64)
    receiver = threading.Thread(target=_read_messages, args=(process, messages), daemon=True)
    receiver.start()
    try:
        assert process.stdin is not None
        process.stdin.write(
            json.dumps(
                {
                    "bam": str(bam),
                    "mode": mode,
                    "maximum_reads": maximum_reads,
                    "confirmation_offset": confirmation_offset,
                }
            ).encode("utf-8")
            + b"\n"
        )
        process.stdin.close()
        while True:
            if cancel_event is not None and cancel_event.is_set():
                return result("cancelled")
            if time.monotonic() - started >= timeout_seconds:
                return result("timeout")
            try:
                message = messages.get(timeout=min(0.05, timeout_seconds))
            except queue.Empty:
                if process.poll() is not None and not receiver.is_alive():
                    return result("worker_failed")
                continue
            count = message.get("checked_reads", checked)
            if type(count) is not int or count < checked:
                return result("worker_failed")
            checked = count
            kind = message.get("type")
            if kind == "progress":
                if progress_callback is not None:
                    try:
                        progress_callback(
                            {
                                "checked_reads": checked,
                                "elapsed_seconds": time.monotonic() - started,
                                "scan_mode": mode,
                                "reader": "pysam",
                            }
                        )
                    except Exception:
                        return result("worker_failed")
            elif kind == "result":
                code = message.get("reason_code")
                complete = message.get("complete")
                offset = message.get("confirmation_offset")
                if (
                    not isinstance(code, str)
                    or code not in _REASONS
                    or type(complete) is not bool
                    or (offset is not None and (type(offset) is not int or not 0 <= offset < 2**63))
                ):
                    return result("worker_failed")
                if code == "detected_5mc" and offset is None:
                    return result("worker_failed")
                remaining = max(0.001, timeout_seconds - (time.monotonic() - started))
                try:
                    if process.wait(timeout=min(remaining, 1.0)) != 0:
                        return result("worker_failed")
                except subprocess.TimeoutExpired:
                    return result("timeout")
                return result(str(code), complete=complete, offset=offset)
            else:
                return result("worker_failed")
    except (OSError, ValueError):
        return result("worker_failed")
    finally:
        if process.poll() is None:
            with suppress(OSError):
                process.kill()
        with suppress(subprocess.TimeoutExpired):
            process.wait(timeout=2)
        receiver.join(timeout=1)
        if process.stdout is not None:
            process.stdout.close()
        if process.stdin is not None and not process.stdin.closed:
            process.stdin.close()

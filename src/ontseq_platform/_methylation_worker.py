"""Private isolated pysam worker. Its stdout protocol contains no alignment contents."""

from __future__ import annotations

import gzip
import importlib
import json
import os
import struct
import sys
import time
import zlib
from collections.abc import Callable
from pathlib import Path
from typing import Any, BinaryIO, cast

MEMORY_LIMIT_BYTES = 768 * 1024 * 1024
MAX_QUERY_BASES = 4_000_000
MAX_MODIFICATION_CALLS = 1_000_000
MAX_TAG_BYTES = 16 * 1024 * 1024
MAX_INDEX_BYTES = 64 * 1024 * 1024
MAX_INDEX_ENTRIES = 2_000_000
INDEX_WINDOWS = 32


def _memory_limit() -> bool:
    if os.name != "nt":
        resource = importlib.import_module("resource")
        resource.setrlimit(resource.RLIMIT_AS, (MEMORY_LIMIT_BYTES, MEMORY_LIMIT_BYTES))
        return True
    # A Job Object bounds native allocations as well as Python allocations on Windows.
    import ctypes
    from ctypes import wintypes

    class Basic(ctypes.Structure):
        _fields_ = [
            ("process_time", ctypes.c_longlong),
            ("job_time", ctypes.c_longlong),
            ("flags", wintypes.DWORD),
            ("minimum_working_set", ctypes.c_size_t),
            ("maximum_working_set", ctypes.c_size_t),
            ("active_processes", wintypes.DWORD),
            ("affinity", ctypes.c_size_t),
            ("priority", wintypes.DWORD),
            ("scheduling", wintypes.DWORD),
        ]

    class IoCounters(ctypes.Structure):
        _fields_ = [(f"counter_{number}", ctypes.c_ulonglong) for number in range(6)]

    class Extended(ctypes.Structure):
        _fields_ = [
            ("basic", Basic),
            ("io", IoCounters),
            ("process_memory", ctypes.c_size_t),
            ("job_memory", ctypes.c_size_t),
            ("peak_process_memory", ctypes.c_size_t),
            ("peak_job_memory", ctypes.c_size_t),
        ]

    # ctypes exports these names only on Windows; resolve them at this FFI boundary.
    kernel = vars(ctypes)["WinDLL"]("kernel32", use_last_error=True)
    kernel.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
    kernel.CreateJobObjectW.restype = wintypes.HANDLE
    kernel.SetInformationJobObject.argtypes = [
        wintypes.HANDLE,
        ctypes.c_int,
        ctypes.c_void_p,
        wintypes.DWORD,
    ]
    kernel.SetInformationJobObject.restype = wintypes.BOOL
    kernel.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
    kernel.AssignProcessToJobObject.restype = wintypes.BOOL
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    job = kernel.CreateJobObjectW(None, None)
    if not job:
        return False
    limit = Extended()
    limit.basic.flags = 0x100  # JOB_OBJECT_LIMIT_PROCESS_MEMORY
    limit.process_memory = MEMORY_LIMIT_BYTES
    if not kernel.SetInformationJobObject(job, 9, ctypes.byref(limit), ctypes.sizeof(limit)):
        kernel.CloseHandle(job)
        return False
    if not kernel.AssignProcessToJobObject(job, kernel.GetCurrentProcess()):
        kernel.CloseHandle(job)
        return False
    # Keep the handle open until this short-lived worker exits.
    return True


class _IndexReader:
    def __init__(self, stream: BinaryIO) -> None:
        self.stream = stream
        self.consumed = 0

    def read(self, size: int) -> bytes:
        if size < 0 or self.consumed + size > MAX_INDEX_BYTES:
            raise ValueError("index limit")
        value = self.stream.read(size)
        self.consumed += len(value)
        if len(value) != size:
            raise ValueError("incomplete index")
        return value

    def count(self) -> int:
        value = struct.unpack("<i", self.read(4))[0]
        if not 0 <= value <= MAX_INDEX_ENTRIES:
            raise ValueError("index count limit")
        return int(value)


def _indexed_offsets(bam: Path) -> list[int]:
    """Read only bounded BAI/CSI offset metadata, never identifiers or sequence data.

    Index offsets are BAM record boundaries, unlike tell() immediately before a region
    iterator advances: the latter can silently seek or skip records. Sampling these exact
    boundaries permits a fresh, independently seekable positive witness.
    """
    candidates = [Path(str(bam) + suffix) for suffix in (".bai", ".csi")]
    candidates.extend(bam.with_suffix(suffix) for suffix in (".bai", ".csi"))
    for path in candidates:
        if path.is_symlink() or not path.is_file() or path.stat().st_size > MAX_INDEX_BYTES:
            continue
        try:
            with path.open("rb") as raw:
                compressed = raw.read(2) == b"\x1f\x8b"
                raw.seek(0)
                stream = cast(BinaryIO, gzip.GzipFile(fileobj=raw) if compressed else raw)
                reader = _IndexReader(stream)
                magic = reader.read(4)
                if magic == b"CSI\x01":
                    min_shift, depth, auxiliary_size = struct.unpack("<iii", reader.read(12))
                    if not 0 <= min_shift <= 31 or not 0 <= depth <= 10:
                        raise ValueError("unsupported index geometry")
                    reader.read(auxiliary_size)
                    normal_bin_limit = ((1 << (3 * (depth + 1))) - 1) // 7
                elif magic == b"BAI\x01":
                    normal_bin_limit = 37449
                else:
                    continue
                offsets: set[int] = set()
                for _ in range(reader.count()):
                    for _ in range(reader.count()):
                        bin_id = struct.unpack("<I", reader.read(4))[0]
                        if magic == b"CSI\x01":
                            reader.read(8)  # loffset is not an extra BAM record.
                        chunks = reader.count()
                        for _ in range(chunks):
                            start, end = struct.unpack("<QQ", reader.read(16))
                            if bin_id < normal_bin_limit and 0 < start < end < 2**63:
                                offsets.add(start)
                    if magic == b"BAI\x01":
                        for _ in range(reader.count()):
                            offset = struct.unpack("<Q", reader.read(8))[0]
                            if 0 < offset < 2**63:
                                offsets.add(offset)
                ordered = sorted(offsets)
                if len(ordered) <= INDEX_WINDOWS:
                    return ordered
                return [
                    ordered[round(i * (len(ordered) - 1) / (INDEX_WINDOWS - 1))]
                    for i in range(INDEX_WINDOWS)
                ]
        except (OSError, ValueError, EOFError, struct.error, zlib.error):
            continue
    return []


def _modification_state(record: Any) -> str:
    """Validate tag pairing, type, cardinality and positions through pysam/HTSlib."""
    if record.query_length > MAX_QUERY_BASES:
        return "record_limit"
    tags = record.get_tags(with_value_type=True)
    mm = [(value, kind) for name, value, kind in tags if name in {"MM", "Mm"}]
    ml = [(value, kind) for name, value, kind in tags if name in {"ML", "Ml"}]
    if not mm and not ml:
        return "invalid_tags" if any(name == "MN" for name, _, _ in tags) else "no_tags"
    if len(mm) != 1 or len(ml) != 1:
        return "invalid_tags"
    mm_value, mm_type = mm[0]
    ml_value, ml_type = ml[0]
    if mm_type != "Z" or not isinstance(mm_value, str) or not mm_value:
        return "invalid_tags"
    if ml_type != "B" or getattr(ml_value, "typecode", None) != "B":
        return "invalid_tags"
    if len(mm_value) > MAX_TAG_BYTES or len(ml_value) > MAX_MODIFICATION_CALLS:
        return "record_limit"
    mn = [(value, kind) for name, value, kind in tags if name == "MN"]
    if len(mn) > 1 or (
        mn
        and (
            type(mn[0][0]) is not int or mn[0][1] not in "cCsSiI" or mn[0][0] != record.query_length
        )
    ):
        return "invalid_tags"
    modifications = record.modified_bases
    if modifications is None:
        return "invalid_tags"
    calls = sum(len(values) for values in modifications.values())
    if calls != len(ml_value) or calls == 0:
        return "invalid_tags"
    for values in modifications.values():
        if any(
            not 0 <= pos < record.query_length or not 0 <= quality <= 255 for pos, quality in values
        ):
            return "invalid_tags"
    if any(
        base == "C" and modification == "m" and values
        for (base, _strand, modification), values in modifications.items()
    ):
        return "detected_5mc"
    return "unsupported_modification"


def _scan(config: dict[str, Any], emit: Callable[[dict[str, Any]], None]) -> None:
    checked = 0
    last_progress = 0.0
    problem: str | None = None

    def finish(code: str, *, complete: bool = False, offset: int | None = None) -> None:
        emit(
            {
                "type": "result",
                "reason_code": code,
                "checked_reads": checked,
                "complete": complete,
                "confirmation_offset": offset,
            }
        )

    def inspect_record(record: Any, offset: int) -> bool:
        nonlocal checked, last_progress, problem
        checked += 1
        now = time.monotonic()
        if checked == 1 or now - last_progress >= 0.2:
            emit({"type": "progress", "checked_reads": checked})
            last_progress = now
        state = _modification_state(record)
        if state == "detected_5mc":
            finish(state, offset=offset)
            return True
        if state == "record_limit":
            finish(state)
            return True
        if state == "invalid_tags" or (state != "no_tags" and problem is None):
            problem = state
        return False

    try:
        pysam = importlib.import_module("pysam")
    except (ImportError, OSError):
        finish("reader_unavailable")
        return
    try:
        if not _memory_limit():
            finish("worker_failed")
            return
        pysam.set_verbosity(0)
        path = Path(config["bam"])
        with pysam.AlignmentFile(str(path), "rb", check_sq=False, threads=1) as bam:
            if not bam.is_bam:
                finish("read_error")
                return
            data_start = bam.tell()
            witness = config.get("confirmation_offset")
            if witness is not None:
                try:
                    bam.seek(witness)
                    offset = bam.tell()
                    if inspect_record(next(bam), offset):
                        return
                except (OSError, ValueError, StopIteration):
                    pass
            thorough = config["mode"] == "thorough"
            maximum_reads = int(config["maximum_reads"])
            offsets = [] if thorough else _indexed_offsets(path)
            if offsets:
                # Start away from the prefix, then alternate across the index's spread.
                offsets = offsets[len(offsets) // 2 :] + offsets[: len(offsets) // 2]
                budget = max(1, maximum_reads // len(offsets))
                for start in offsets:
                    bam.seek(start)
                    for _ in range(budget):
                        if checked >= maximum_reads:
                            break
                        offset = bam.tell()
                        try:
                            record = next(bam)
                        except StopIteration:
                            break
                        if inspect_record(record, offset):
                            return
                finish(problem or "sample_incomplete")
                return
            bam.seek(data_start)
            while thorough or checked < maximum_reads:
                offset = bam.tell()
                try:
                    record = next(bam)
                except StopIteration:
                    finish(problem or "complete_no_tags", complete=True)
                    return
                if inspect_record(record, offset):
                    return
            finish(problem or "sample_incomplete")
    except MemoryError:
        finish("record_limit")
    except (OSError, ValueError, TypeError, OverflowError):
        finish("read_error")


def main() -> None:
    def emit(message: dict[str, Any]) -> None:
        print(json.dumps(message, separators=(",", ":")), flush=True)

    try:
        line = sys.stdin.buffer.readline(16385)
        if len(line) > 16384 or not line.endswith(b"\n"):
            raise ValueError("invalid input")
        config = json.loads(line)
        _scan(config, emit)
    except Exception:
        # Never print exceptions: HTSlib and Python errors can contain private fields.
        emit(
            {
                "type": "result",
                "reason_code": "worker_failed",
                "checked_reads": 0,
                "complete": False,
                "confirmation_offset": None,
            }
        )


if __name__ == "__main__":
    main()

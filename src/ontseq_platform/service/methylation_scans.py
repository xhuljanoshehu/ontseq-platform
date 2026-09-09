"""Short-lived, identity-bound BAM discovery state for one local service instance.

The cache is an interactive optimization, never analytical evidence or user consent.
Before-start refreshes always invoke the reader again. A previous positive may supply
only a private BAM virtual offset so the reader can freshly inspect its witness record.
"""

from __future__ import annotations

import threading
import time
import uuid
from collections import OrderedDict
from collections.abc import Callable, Hashable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from ..methylation_probe import PROBE_VERSION, MethylationProbe, _stat_identity
from .guard import GuardError

PROBE_CACHE_TTL_SECONDS = 30.0
MAX_CACHED_PROBES = 128
SCAN_RETENTION_SECONDS = 600.0
MAX_RETAINED_SCANS = 32
THOROUGH_TIMEOUT_SECONDS = 300.0
AVAILABILITY_CACHE_TTL_SECONDS = 10.0


class ScanBusyError(RuntimeError):
    pass


class ScanNotFoundError(KeyError):
    pass


class AvailabilityCache:
    """Cache a bounded runtime capability check independently from BAM evidence."""

    def __init__(
        self,
        check: Callable[[], dict[str, Any]],
        identity: Callable[[], Hashable],
        *,
        ttl_seconds: float = AVAILABILITY_CACHE_TTL_SECONDS,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        self._check = check
        self._identity = identity
        self._ttl = ttl_seconds
        self._clock = clock
        self._lock = threading.Lock()
        self._cached: tuple[Hashable, float, dict[str, Any]] | None = None

    def get(self, *, force_refresh: bool = False) -> dict[str, Any]:
        with self._lock:
            identity = self._identity()
            now = self._clock()
            cached = self._cached
            if not force_refresh and cached and cached[0] == identity and cached[1] > now:
                return dict(cached[2])
            result = self._check()
            # A replacement during --version must not seed the new binary's cache.
            if self._identity() == identity:
                self._cached = (identity, self._clock() + self._ttl, dict(result))
            else:
                self._cached = None
                result = {
                    **result,
                    "methylation_available": False,
                    "methylation_unavailable_reason": (
                        "The methylation runtime changed during its check."
                    ),
                }
            return result


@dataclass
class _Scan:
    scan_id: str
    bam: Path
    bam_identity: str
    started_at: float
    cancel_event: threading.Event
    state: Literal["running", "completed", "cancelled", "failed"] = "running"
    checked_reads: int = 0
    finished_at: float | None = None
    result: MethylationProbe | None = None
    thread: threading.Thread | None = None


class MethylationScans:
    """Own one discovery slot, bounded display cache and cancellable thorough jobs."""

    def __init__(
        self,
        probe: Callable[..., MethylationProbe],
        availability: AvailabilityCache,
        *,
        samtools: str,
        validate_path: Callable[[Path], Path],
        analysis_running: Callable[[], bool] = lambda: False,
        cache_ttl_seconds: float = PROBE_CACHE_TTL_SECONDS,
        retention_seconds: float = SCAN_RETENTION_SECONDS,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        self._probe = probe
        self._availability = availability
        self._samtools = samtools
        self._validate_path = validate_path
        self._analysis_running = analysis_running
        self._cache_ttl = cache_ttl_seconds
        self._retention = retention_seconds
        self._clock = clock
        self._lock = threading.RLock()
        self._cache: OrderedDict[tuple[str, str, str], tuple[float, MethylationProbe]] = (
            OrderedDict()
        )
        self._scans: OrderedDict[str, _Scan] = OrderedDict()
        self._active_cancel: threading.Event | None = None
        self._closed = False

    def availability(self, *, force_refresh: bool = False) -> dict[str, Any]:
        return self._availability.get(force_refresh=force_refresh)

    def _identity(self, bam: Path) -> str:
        if self._validate_path(bam) != bam or not bam.is_file():
            raise ValueError("The selected BAM changed its resolved location.")
        return _stat_identity(bam.stat())

    def _prune_locked(self) -> None:
        now = self._clock()
        for key, (expires, _) in list(self._cache.items()):
            if expires <= now:
                del self._cache[key]
        for scan_id, scan in list(self._scans.items()):
            if (
                scan.finished_at is not None
                and now - scan.finished_at >= self._retention
                and (scan.thread is None or not scan.thread.is_alive())
            ):
                del self._scans[scan_id]

    def _claim_locked(self) -> threading.Event:
        if self._closed or self._active_cancel is not None or self._analysis_running():
            raise ScanBusyError("Another BAM discovery scan is active or the service is stopping.")
        cancel = threading.Event()
        self._active_cancel = cancel
        return cancel

    def claim_analysis(self, claim: Callable[[], None]) -> None:
        """Serialize pipeline admission with discovery admission, including cancellation."""
        with self._lock:
            if self._closed or self._active_cancel is not None:
                raise ScanBusyError("BAM discovery must finish or acknowledge cancellation first.")
            claim()

    def discovery_active(self) -> bool:
        with self._lock:
            return self._active_cancel is not None

    def _unknown(
        self,
        reason_code: str,
        reason: str,
        *,
        checked: int = 0,
        elapsed: float = 0.0,
        mode: str = "quick",
        identity: str | None = None,
    ) -> MethylationProbe:
        return MethylationProbe(
            "unknown",
            reason,
            checked_reads=checked,
            reason_code=reason_code,
            elapsed_seconds=elapsed,
            scan_mode=mode,
            bam_identity=identity,
        )

    def _witness_locked(self, bam: Path, identity: str) -> int | None:
        cached = self._cache.get((str(bam), identity, PROBE_VERSION))
        candidates = [cached[1]] if cached else []
        # Display cache expiry must not lose a thorough witness: its record still has
        # to be read again, and retention is bounded independently from the UI cache.
        candidates.extend(
            scan.result
            for scan in reversed(self._scans.values())
            if scan.bam == bam
            and scan.bam_identity == identity
            and scan.state == "completed"
            and scan.result is not None
        )
        for result in candidates:
            if result.status == "detected" and result.probe_version == PROBE_VERSION:
                return result.confirmation_offset
        return None

    def _cache_result(self, bam: Path, identity: str, result: MethylationProbe) -> MethylationProbe:
        try:
            current = self._identity(bam)
        except (OSError, ValueError, GuardError):
            current = None
        if current != identity or (
            result.bam_identity is not None and result.bam_identity != identity
        ):
            result = self._unknown(
                "file_changed",
                "The BAM changed during discovery.",
                checked=result.checked_reads,
                elapsed=result.elapsed_seconds,
                mode=result.scan_mode,
            )
        with self._lock:
            for key in list(self._cache):
                if key[0] == str(bam) and key[1] != current:
                    del self._cache[key]
            if (
                not self._closed
                and current == identity
                and result.bam_identity == identity
                and result.probe_version == PROBE_VERSION
                and result.reason_code != "cancelled"
            ):
                key = (str(bam), identity, PROBE_VERSION)
                self._cache[key] = (self._clock() + self._cache_ttl, result)
                self._cache.move_to_end(key)
                while len(self._cache) > MAX_CACHED_PROBES:
                    self._cache.popitem(last=False)
        return result

    def _public_result(
        self, bam: Path, result: MethylationProbe, *, force_availability: bool = False
    ) -> dict[str, Any]:
        available = self.availability(force_refresh=force_availability)
        try:
            current_identity = self._identity(bam)
        except (OSError, ValueError, GuardError):
            current_identity = None
        if result.bam_identity is not None and current_identity != result.bam_identity:
            result = self._unknown(
                "file_changed",
                "The BAM changed since this discovery result.",
                checked=result.checked_reads,
                elapsed=result.elapsed_seconds,
                mode=result.scan_mode,
            )
        elif result.status in {"detected", "not_detected"} and result.bam_identity is None:
            result = self._unknown(
                "read_error",
                "The discovery result is not bound to the selected BAM.",
                checked=result.checked_reads,
                elapsed=result.elapsed_seconds,
                mode=result.scan_mode,
            )
        payload = asdict(result)
        payload.pop("confirmation_offset", None)
        return {
            **payload,
            "bam_path": str(bam),
            **available,
        }

    def quick(self, bam: Path, *, force_refresh: bool = False) -> dict[str, Any]:
        identity = self._identity(bam)
        cached_result = None
        with self._lock:
            self._prune_locked()
            cached = self._cache.get((str(bam), identity, PROBE_VERSION))
            if cached and not force_refresh:
                cached_result = cached[1]
            witness = self._witness_locked(bam, identity)
            try:
                cancel = None if cached_result else self._claim_locked()
            except ScanBusyError:
                cached_result = self._unknown(
                    "busy", "Another BAM discovery scan is in progress.", identity=identity
                )
                cancel = None
        if cached_result is not None:
            return self._public_result(bam, cached_result, force_availability=force_refresh)
        assert cancel is not None
        try:
            self._identity(bam)
            result = self._probe(
                bam, samtools=self._samtools, cancel_event=cancel, confirmation_offset=witness
            )
            if cancel.is_set():
                result = self._unknown(
                    "cancelled",
                    "BAM discovery was cancelled.",
                    checked=result.checked_reads,
                    elapsed=result.elapsed_seconds,
                )
            result = self._cache_result(bam, identity, result)
        except Exception:
            result = self._unknown("worker_failed", "The BAM discovery worker could not complete.")
        finally:
            with self._lock:
                self._active_cancel = None
        return self._public_result(bam, result, force_availability=force_refresh)

    def start(self, bam: Path) -> dict[str, Any]:
        identity = self._identity(bam)
        with self._lock:
            self._prune_locked()
            cancel = self._claim_locked()
            while len(self._scans) >= MAX_RETAINED_SCANS:
                self._scans.popitem(last=False)
            scan = _Scan(uuid.uuid4().hex, bam, identity, self._clock(), cancel)
            self._scans[scan.scan_id] = scan
            scan.thread = threading.Thread(
                target=self._work, args=(scan,), daemon=True, name="ontseq-methylation-scan"
            )
            try:
                scan.thread.start()
            except RuntimeError:
                self._active_cancel = None
                del self._scans[scan.scan_id]
                raise
            snapshot = self._snapshot_locked(scan)
        return snapshot

    def _work(self, scan: _Scan) -> None:
        def progress(payload: dict[str, Any]) -> None:
            with self._lock:
                count = payload.get("checked_reads")
                if (
                    scan.state == "running"
                    and isinstance(count, int)
                    and not isinstance(count, bool)
                ):
                    scan.checked_reads = max(scan.checked_reads, count)

        try:
            if self._identity(scan.bam) != scan.bam_identity:
                result = self._unknown(
                    "file_changed", "The BAM changed before discovery started.", mode="thorough"
                )
            else:
                result = self._probe(
                    scan.bam,
                    samtools=self._samtools,
                    mode="thorough",
                    timeout_seconds=THOROUGH_TIMEOUT_SECONDS,
                    cancel_event=scan.cancel_event,
                    progress_callback=progress,
                )
        except Exception:
            result = self._unknown(
                "worker_failed",
                "The BAM discovery worker could not complete.",
                checked=scan.checked_reads,
                mode="thorough",
            )
        with self._lock:
            scan.checked_reads = max(scan.checked_reads, result.checked_reads)
            if scan.cancel_event.is_set():
                scan.state = "cancelled"
                result = self._unknown(
                    "cancelled",
                    "The thorough BAM scan was cancelled.",
                    checked=scan.checked_reads,
                    mode="thorough",
                    elapsed=max(0.0, self._clock() - scan.started_at),
                )
            else:
                result = self._cache_result(scan.bam, scan.bam_identity, result)
                scan.state = (
                    "failed"
                    if result.reason_code
                    in {"worker_failed", "reader_unavailable", "input_unavailable", "read_error"}
                    else "completed"
                )
            scan.result = result
            scan.finished_at = self._clock()
            self._active_cancel = None

    def _snapshot_locked(self, scan: _Scan) -> dict[str, Any]:
        return {
            "scan_id": scan.scan_id,
            "bam_path": str(scan.bam),
            "state": scan.state,
            "checked_reads": scan.checked_reads,
            "elapsed_seconds": max(0.0, (scan.finished_at or self._clock()) - scan.started_at),
            "result": None,
        }

    def get(self, scan_id: str) -> dict[str, Any]:
        with self._lock:
            self._prune_locked()
            scan = self._scans.get(scan_id)
            if scan is None:
                raise ScanNotFoundError(scan_id)
            snapshot = self._snapshot_locked(scan)
            result = scan.result
        if result is not None:
            snapshot["result"] = self._public_result(scan.bam, result)
        return snapshot

    def cancel(self, scan_id: str) -> dict[str, Any]:
        with self._lock:
            self._prune_locked()
            scan = self._scans.get(scan_id)
            if scan is None:
                raise ScanNotFoundError(scan_id)
            if scan.state == "running":
                scan.cancel_event.set()
                # Remain running until the reader has acknowledged cancellation and
                # released its process/resources; terminal means the slot is reusable.
        return self.get(scan_id)

    def close(self) -> None:
        with self._lock:
            self._closed = True
            if self._active_cancel is not None:
                self._active_cancel.set()
            threads = [scan.thread for scan in self._scans.values() if scan.thread is not None]
            self._cache.clear()
        deadline = time.monotonic() + 5.0
        for thread in threads:
            if thread is not threading.current_thread():
                thread.join(timeout=max(0.0, deadline - time.monotonic()))

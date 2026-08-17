from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.pipeline.lock import (
    LOCK_FILENAME,
    LockHolder,
    RunAlreadyRunning,
    read_holder,
    run_lock,
)


def _holder(**overrides: object) -> LockHolder:
    values: dict[str, object] = {
        "pid": 4242,
        "hostname": socket.gethostname(),
        "acquired_at": "2026-08-17T09:00:00+00:00",
        "run_id": "RUN_001",
        "sample_id": "SAMPLE_001",
        "pipeline_version": "0.0.0-test",
    }
    values.update(overrides)
    return LockHolder(**values)  # type: ignore[arg-type]


class LockCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.root = Path(self._temporary.name) / "envelope"
        self.lock_path = self.root / LOCK_FILENAME

    def _lock(self, **overrides: object):
        values: dict[str, object] = {
            "run_id": "RUN_001",
            "sample_id": "SAMPLE_001",
            "pipeline_version": "0.0.0-test",
        }
        values.update(overrides)
        return run_lock(self.root, **values)  # type: ignore[arg-type]

    def _write_lock(self, holder: LockHolder) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        self.lock_path.write_text(json.dumps(holder.__dict__), encoding="utf-8")


class AcquisitionTests(LockCase):
    def test_the_lock_file_exists_while_held_and_is_gone_afterwards(self) -> None:
        with self._lock():
            self.assertTrue(self.lock_path.is_file())
        self.assertFalse(self.lock_path.exists())

    def test_a_clean_acquisition_warns_about_nothing(self) -> None:
        with self._lock() as warnings:
            self.assertEqual(warnings, [])

    def test_the_holder_records_who_took_it(self) -> None:
        """Without this, a held lock is a dead end rather than a thing to investigate."""
        with self._lock():
            holder = read_holder(self.lock_path)
        self.assertIsNotNone(holder)
        self.assertEqual(holder.pid, os.getpid())
        self.assertEqual(holder.hostname, socket.gethostname())
        self.assertEqual(holder.run_id, "RUN_001")
        self.assertEqual(holder.sample_id, "SAMPLE_001")

    def test_the_envelope_directory_is_created_if_missing(self) -> None:
        self.assertFalse(self.root.exists())
        with self._lock():
            self.assertTrue(self.root.is_dir())

    def test_the_lock_is_released_when_the_block_raises(self) -> None:
        """A failed run must not leave the envelope locked against the retry."""
        with self.assertRaises(RuntimeError):
            with self._lock():
                raise RuntimeError("stage exploded")
        self.assertFalse(self.lock_path.exists())

    def test_the_lock_can_be_taken_again_after_release(self) -> None:
        with self._lock():
            pass
        with self._lock() as warnings:
            self.assertEqual(warnings, [])


class ContentionTests(LockCase):
    def test_a_live_local_holder_blocks(self) -> None:
        self._write_lock(_holder(pid=os.getpid()))
        with self.assertRaises(RunAlreadyRunning) as raised:
            with self._lock():
                pass
        self.assertIn("already in use", str(raised.exception))

    def test_the_error_names_the_holder_and_the_file_to_remove(self) -> None:
        self._write_lock(_holder(pid=os.getpid(), run_id="OTHER_RUN"))
        with self.assertRaises(RunAlreadyRunning) as raised:
            with self._lock():
                pass
        self.assertEqual(raised.exception.holder.run_id, "OTHER_RUN")
        self.assertIn(LOCK_FILENAME, str(raised.exception))

    def test_a_blocked_acquisition_leaves_the_existing_lock_intact(self) -> None:
        self._write_lock(_holder(pid=os.getpid()))
        before = self.lock_path.read_text(encoding="utf-8")
        with self.assertRaises(RunAlreadyRunning):
            with self._lock():
                pass
        self.assertEqual(self.lock_path.read_text(encoding="utf-8"), before)

    def test_nesting_the_same_envelope_is_refused(self) -> None:
        """Two runs in one process are still two runs."""
        with self._lock():
            with self.assertRaises(RunAlreadyRunning):
                with self._lock():
                    pass

    def test_a_different_envelope_is_unaffected(self) -> None:
        other = Path(self._temporary.name) / "other-envelope"
        with self._lock():
            with run_lock(
                other, run_id="RUN_002", sample_id="SAMPLE_002", pipeline_version="0.0.0-test"
            ) as warnings:
                self.assertEqual(warnings, [])


class StaleTests(LockCase):
    def _dead_pid(self) -> int:
        """A PID that has certainly exited, obtained by starting a process and reaping it."""
        process = subprocess.Popen([sys.executable, "-c", ""])
        process.wait()
        return process.pid

    def test_a_dead_local_holder_is_reclaimed(self) -> None:
        self._write_lock(_holder(pid=self._dead_pid()))
        with self._lock() as warnings:
            self.assertEqual(len(warnings), 1)
            self.assertIn("Reclaimed a stale run lock", warnings[0])

    def test_the_reclaim_warning_names_the_previous_holder(self) -> None:
        """The reclaim goes into the run report, so it has to say what was stepped over."""
        self._write_lock(_holder(pid=self._dead_pid(), run_id="CRASHED_RUN"))
        with self._lock() as warnings:
            self.assertIn("CRASHED_RUN", warnings[0])

    def test_an_empty_lock_file_is_reclaimed(self) -> None:
        """Left by a process that died between creating the file and writing to it."""
        self.root.mkdir(parents=True, exist_ok=True)
        self.lock_path.write_text("", encoding="utf-8")
        with self._lock() as warnings:
            self.assertIn("unreadable or incomplete", warnings[0])

    def test_a_corrupt_lock_file_is_reclaimed(self) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        self.lock_path.write_text("{not json", encoding="utf-8")
        with self._lock() as warnings:
            self.assertEqual(len(warnings), 1)

    def test_a_lock_from_another_host_is_never_reclaimed(self) -> None:
        """On shared storage a crashed remote run and a live one look identical."""
        self._write_lock(_holder(pid=self._dead_pid(), hostname="some-other-machine"))
        with self.assertRaises(RunAlreadyRunning) as raised:
            with self._lock():
                pass
        self.assertIn("another host", str(raised.exception))

    def test_reclaiming_can_be_switched_off(self) -> None:
        self._write_lock(_holder(pid=self._dead_pid()))
        with self.assertRaises(RunAlreadyRunning):
            with self._lock(reclaim_stale=False):
                pass

    def test_a_reclaimed_lock_belongs_to_the_new_holder(self) -> None:
        self._write_lock(_holder(pid=self._dead_pid()))
        with self._lock():
            holder = read_holder(self.lock_path)
        self.assertEqual(holder.pid, os.getpid())


class ReleaseTests(LockCase):
    def test_a_lock_taken_over_by_someone_else_is_not_deleted_on_exit(self) -> None:
        """If we were wrongly judged dead, the lock is theirs now — leaving is not our call."""
        with self._lock():
            self._write_lock(_holder(pid=os.getpid() + 1, run_id="TOOK_OVER"))
        self.assertTrue(self.lock_path.is_file())
        self.assertEqual(read_holder(self.lock_path).run_id, "TOOK_OVER")

    def test_a_manually_removed_lock_does_not_break_release(self) -> None:
        with self._lock():
            self.lock_path.unlink()
        self.assertFalse(self.lock_path.exists())


class RealConcurrencyTests(LockCase):
    """Two actual processes, because the whole point is a race no single process can show."""

    #: Loads lock.py by path rather than importing the package, which would pull in
    #: pydantic. The module is stdlib-only by design, and this keeps the test that way too.
    SCRIPT = """
import sys, time, pathlib, importlib.util
spec = importlib.util.spec_from_file_location("ontseq_lock", sys.argv[1])
lock = importlib.util.module_from_spec(spec)
# Registered before execution: dataclass() resolves the postponed annotations through
# sys.modules, and an unregistered module makes that lookup fail.
sys.modules[spec.name] = lock
spec.loader.exec_module(lock)
run_lock, RunAlreadyRunning = lock.run_lock, lock.RunAlreadyRunning
try:
    with run_lock(
        pathlib.Path(sys.argv[2]),
        run_id="RUN_001",
        sample_id="SAMPLE_001",
        pipeline_version="0.0.0-test",
    ):
        print("ACQUIRED", flush=True)
        time.sleep(float(sys.argv[3]))
except RunAlreadyRunning:
    print("BLOCKED", flush=True)
"""

    def _spawn(self, hold_seconds: float) -> subprocess.Popen[str]:
        module = (
            Path(__file__).resolve().parent.parent
            / "src/ontseq_platform/pipeline/lock.py"
        )
        return subprocess.Popen(
            [sys.executable, "-c", self.SCRIPT, str(module), str(self.root), str(hold_seconds)],
            stdout=subprocess.PIPE,
            text=True,
        )

    def test_a_second_process_is_refused_while_the_first_holds(self) -> None:
        first = self._spawn(2.0)
        self.addCleanup(first.wait)
        self.assertEqual(first.stdout.readline().strip(), "ACQUIRED")

        second = self._spawn(0.0)
        self.assertEqual(second.communicate()[0].strip(), "BLOCKED")

    def test_the_lock_is_free_once_the_first_process_exits(self) -> None:
        first = self._spawn(0.0)
        self.assertEqual(first.communicate()[0].strip(), "ACQUIRED")
        second = self._spawn(0.0)
        self.assertEqual(second.communicate()[0].strip(), "ACQUIRED")


if __name__ == "__main__":
    unittest.main()

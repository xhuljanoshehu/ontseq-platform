from __future__ import annotations

import json
import os
import tempfile
import time
import unittest
from pathlib import Path

from ontseq_platform.pipeline.watch import (
    LEDGER_FILENAME,
    Attempt,
    Candidate,
    Ledger,
    Outcome,
    Readiness,
    discover,
    inspect_directory,
    sample_id_from_directory,
    should_attempt,
)

NOW = 1_000_000.0


class WatchCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.root = Path(self._temporary.name)
        self.watch = self.root / "drop"
        self.watch.mkdir()

    def _drop(self, name: str, *, age_seconds: float = 10_000, files: int = 2) -> Path:
        directory = self.watch / name
        directory.mkdir(parents=True, exist_ok=True)
        for index in range(files):
            item = directory / f"reads_{index}.bam"
            item.write_bytes(b"BAM")
            os.utime(item, (NOW - age_seconds, NOW - age_seconds))
        return directory

    def _discover(self, **overrides: object) -> list[Candidate]:
        values: dict[str, object] = {"quiet_seconds": 300.0, "now": NOW}
        values.update(overrides)
        return discover(self.watch, **values)  # type: ignore[arg-type]


class DiscoveryTests(WatchCase):
    def test_a_missing_watch_folder_is_an_error_not_an_empty_result(self) -> None:
        """Silently finding nothing is how a mistyped path looks like an idle system."""
        with self.assertRaises(NotADirectoryError):
            discover(self.root / "nowhere")

    def test_an_empty_watch_folder_yields_nothing(self) -> None:
        self.assertEqual(self._discover(), [])

    def test_candidates_are_returned_sorted_by_name(self) -> None:
        for name in ("SAMPLE_C", "SAMPLE_A", "SAMPLE_B"):
            self._drop(name)
        self.assertEqual(
            [item.name for item in self._discover()], ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"]
        )

    def test_loose_files_are_not_candidates(self) -> None:
        (self.watch / "stray.bam").write_bytes(b"BAM")
        self.assertEqual(self._discover(), [])

    def test_dot_and_underscore_directories_are_ignored(self) -> None:
        self._drop(".hidden")
        self._drop("_scratch")
        self._drop("SAMPLE_A")
        self.assertEqual([item.name for item in self._discover()], ["SAMPLE_A"])

    def test_files_are_counted_across_subdirectories(self) -> None:
        directory = self._drop("SAMPLE_A", files=1)
        nested = directory / "pod5"
        nested.mkdir()
        (nested / "batch.pod5").write_bytes(b"signal")
        os.utime(nested / "batch.pod5", (NOW - 10_000, NOW - 10_000))
        self.assertEqual(self._discover()[0].file_count, 2)


class QuiescenceTests(WatchCase):
    def test_a_directory_untouched_beyond_the_window_is_ready(self) -> None:
        self._drop("SAMPLE_A", age_seconds=1_000)
        candidate = self._discover(quiet_seconds=300.0)[0]
        self.assertIs(candidate.readiness, Readiness.READY)

    def test_a_directory_still_being_written_is_not_ready(self) -> None:
        self._drop("SAMPLE_A", age_seconds=10)
        candidate = self._discover(quiet_seconds=300.0)[0]
        self.assertIs(candidate.readiness, Readiness.STILL_CHANGING)
        self.assertIn("quiet window", candidate.detail)

    def test_the_newest_file_decides_not_the_oldest(self) -> None:
        """One late-arriving file means the run is still being written."""
        directory = self._drop("SAMPLE_A", age_seconds=10_000)
        late = directory / "final.bam"
        late.write_bytes(b"BAM")
        os.utime(late, (NOW - 5, NOW - 5))
        self.assertIs(self._discover()[0].readiness, Readiness.STILL_CHANGING)

    def test_the_readiness_reason_admits_that_quiescence_is_a_heuristic(self) -> None:
        self._drop("SAMPLE_A", age_seconds=1_000)
        self.assertIn("heuristic", self._discover()[0].detail)

    def test_an_empty_directory_is_never_ready(self) -> None:
        (self.watch / "SAMPLE_A").mkdir()
        candidate = self._discover()[0]
        self.assertIs(candidate.readiness, Readiness.EMPTY)

    def test_a_directory_with_only_subdirectories_is_empty(self) -> None:
        (self.watch / "SAMPLE_A" / "pod5").mkdir(parents=True)
        self.assertIs(self._discover()[0].readiness, Readiness.EMPTY)

    def test_a_file_dated_in_the_future_is_not_ready(self) -> None:
        """Clock skew must not make an in-progress run look long finished."""
        directory = self._drop("SAMPLE_A", age_seconds=10_000)
        skewed = directory / "skewed.bam"
        skewed.write_bytes(b"BAM")
        os.utime(skewed, (NOW + 5_000, NOW + 5_000))
        self.assertIs(self._discover()[0].readiness, Readiness.STILL_CHANGING)


class MarkerTests(WatchCase):
    def test_a_marker_makes_a_directory_ready_regardless_of_age(self) -> None:
        """An explicit signal from the producer beats guessing from timestamps."""
        directory = self._drop("SAMPLE_A", age_seconds=1)
        (directory / ".ready").write_text("", encoding="utf-8")
        candidate = self._discover(ready_marker=".ready")[0]
        self.assertIs(candidate.readiness, Readiness.READY)
        self.assertIn(".ready", candidate.detail)

    def test_without_the_marker_a_quiet_directory_is_still_not_ready(self) -> None:
        """Configuring a marker means the marker is the contract, not a second opinion."""
        self._drop("SAMPLE_A", age_seconds=10_000)
        candidate = self._discover(ready_marker=".ready")[0]
        self.assertIs(candidate.readiness, Readiness.MARKER_MISSING)

    def test_an_empty_directory_with_a_marker_is_still_empty(self) -> None:
        directory = self.watch / "SAMPLE_A"
        directory.mkdir()
        (directory / ".ready").write_text("", encoding="utf-8")
        # The marker is itself a file, so the directory is not empty — but it carries no
        # payload. Counting it is honest: emptiness is about files, not about meaning.
        self.assertIs(self._discover(ready_marker=".ready")[0].readiness, Readiness.READY)


class SampleIdTests(unittest.TestCase):
    def test_a_conforming_directory_name_is_used_as_is(self) -> None:
        self.assertEqual(sample_id_from_directory("AML_0031"), "AML_0031")

    def test_a_nonconforming_name_is_refused_rather_than_cleaned(self) -> None:
        """Inventing an identifier puts a name nobody chose onto a reviewer artifact."""
        for name in ("has space", "ä-umlaut", "-leading-dash", "ab", "x" * 65, ""):
            self.assertIsNone(sample_id_from_directory(name), name)

    def test_permitted_punctuation_is_accepted(self) -> None:
        self.assertEqual(sample_id_from_directory("260611_RAD114_AS_S700"), "260611_RAD114_AS_S700")


def _candidate(readiness: Readiness = Readiness.READY, detail: str = "ok") -> Candidate:
    return Candidate(path=Path("/drop/S"), name="S", readiness=readiness, detail=detail)


def _attempt(outcome: Outcome, detail: str = "detail") -> Attempt:
    return Attempt(name="S", outcome=outcome, detail=detail, attempted_at="2026-08-17T09:00:00Z")


class ShouldAttemptTests(unittest.TestCase):
    def test_a_ready_unseen_candidate_is_attempted(self) -> None:
        run, why = should_attempt(_candidate(), None)
        self.assertTrue(run)
        self.assertIn("not attempted", why)

    def test_a_not_ready_candidate_is_skipped_with_its_own_reason(self) -> None:
        run, why = should_attempt(_candidate(Readiness.STILL_CHANGING, "still writing"), None)
        self.assertFalse(run)
        self.assertEqual(why, "still writing")

    def test_a_completed_candidate_is_not_repeated(self) -> None:
        run, why = should_attempt(_candidate(), _attempt(Outcome.COMPLETED))
        self.assertFalse(run)
        self.assertIn("already completed", why)

    def test_a_failure_is_not_retried_automatically(self) -> None:
        """A deterministic failure repeated every poll hides the sample that needs a human."""
        run, why = should_attempt(_candidate(), _attempt(Outcome.FAILED, "qc gate failed"))
        self.assertFalse(run)
        self.assertIn("qc gate failed", why)
        self.assertIn("--retry-failed", why)

    def test_a_failure_is_retried_when_asked(self) -> None:
        run, _ = should_attempt(_candidate(), _attempt(Outcome.FAILED), retry_failed=True)
        self.assertTrue(run)

    def test_a_lock_conflict_is_always_retried(self) -> None:
        """Being blocked says something about another process, not about this sample."""
        run, why = should_attempt(_candidate(), _attempt(Outcome.LOCKED))
        self.assertTrue(run)
        self.assertIn("retrying", why)

    def test_a_rejected_candidate_stays_rejected(self) -> None:
        run, why = should_attempt(_candidate(), _attempt(Outcome.REJECTED, "unusable name"))
        self.assertFalse(run)
        self.assertIn("unusable name", why)

    def test_a_rejected_candidate_is_not_revived_by_retry_failed(self) -> None:
        run, _ = should_attempt(
            _candidate(), _attempt(Outcome.REJECTED, "unusable name"), retry_failed=True
        )
        self.assertFalse(run)


class LedgerTests(WatchCase):
    def setUp(self) -> None:
        super().setUp()
        self.ledger_path = self.root / "out" / LEDGER_FILENAME

    def test_a_missing_ledger_loads_empty(self) -> None:
        self.assertEqual(list(Ledger.load(self.ledger_path).attempts()), [])

    def test_a_recorded_attempt_survives_a_reload(self) -> None:
        ledger = Ledger.load(self.ledger_path)
        ledger.record(_attempt(Outcome.COMPLETED, "all stages passed"))
        reloaded = Ledger.load(self.ledger_path)
        self.assertEqual(reloaded.get("S").outcome, Outcome.COMPLETED)
        self.assertEqual(reloaded.get("S").detail, "all stages passed")

    def test_recording_the_same_name_twice_keeps_the_newer_outcome(self) -> None:
        ledger = Ledger.load(self.ledger_path)
        ledger.record(_attempt(Outcome.FAILED))
        ledger.record(_attempt(Outcome.COMPLETED))
        self.assertEqual(Ledger.load(self.ledger_path).get("S").outcome, Outcome.COMPLETED)

    def test_a_corrupt_ledger_loads_empty_rather_than_stopping_the_watcher(self) -> None:
        """One bad file must not turn into a stopped watcher; resume makes the redo cheap."""
        self.ledger_path.parent.mkdir(parents=True, exist_ok=True)
        self.ledger_path.write_text("{not json", encoding="utf-8")
        self.assertEqual(list(Ledger.load(self.ledger_path).attempts()), [])

    def test_an_unparseable_entry_is_dropped_and_the_rest_survive(self) -> None:
        self.ledger_path.parent.mkdir(parents=True, exist_ok=True)
        self.ledger_path.write_text(
            json.dumps(
                {
                    "attempts": {
                        "GOOD": {
                            "outcome": "completed",
                            "detail": "d",
                            "attempted_at": "t",
                        },
                        "BAD": {"outcome": "not-an-outcome"},
                    }
                }
            ),
            encoding="utf-8",
        )
        ledger = Ledger.load(self.ledger_path)
        self.assertIsNotNone(ledger.get("GOOD"))
        self.assertIsNone(ledger.get("BAD"))

    def test_the_ledger_is_written_atomically(self) -> None:
        """A watcher killed mid-write must not come back to a half-written ledger."""
        ledger = Ledger.load(self.ledger_path)
        ledger.record(_attempt(Outcome.COMPLETED))
        leftovers = [item for item in self.ledger_path.parent.iterdir() if item.suffix == ".tmp"]
        self.assertEqual(leftovers, [])
        json.loads(self.ledger_path.read_text(encoding="utf-8"))

    def test_attempts_are_returned_in_a_stable_order(self) -> None:
        ledger = Ledger.load(self.ledger_path)
        for name in ("S_C", "S_A", "S_B"):
            ledger.record(Attempt(name, Outcome.COMPLETED, "d", "t"))
        self.assertEqual([item.name for item in ledger.attempts()], ["S_A", "S_B", "S_C"])


if __name__ == "__main__":
    unittest.main()


class GridionMarkerTests(unittest.TestCase):
    """MinKNOW's completion signal is a glob one level down, not a fixed top-level name.

    A GridION writes ``final_summary_<flowcell>_<run>_<hash>.txt`` into the run directory
    when sequencing finishes. Matching only a literal name at the top level would leave the
    one authoritative signal the instrument emits unusable, and force every real run onto
    the quiescence heuristic the design explicitly calls a heuristic.
    """

    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.sample = Path(self._temporary.name) / "AML_00123"
        self.run = self.sample / "20260818_1030_X1_FAV12345_abcdef12"
        (self.run / "pod5_pass").mkdir(parents=True)
        (self.run / "pod5_pass" / "reads.pod5").write_bytes(b"POD5")

    def _inspect(self, marker: str | None):
        return inspect_directory(
            self.sample, ready_marker=marker, quiet_seconds=0.0, now=time.time() + 10_000
        )

    def _final_summary(self) -> None:
        (self.run / "final_summary_FAV12345_abcdef12_0123abcd.txt").write_text("x")

    def test_a_run_still_sequencing_is_not_ready(self) -> None:
        readiness, detail, _, _ = self._inspect("final_summary_*.txt")
        self.assertIs(readiness, Readiness.MARKER_MISSING)
        self.assertIn("final_summary_*.txt", detail)

    def test_a_finished_run_is_ready(self) -> None:
        self._final_summary()
        readiness, detail, _, _ = self._inspect("final_summary_*.txt")
        self.assertIs(readiness, Readiness.READY)
        self.assertIn("final_summary_FAV12345", detail)

    def test_the_marker_is_found_below_the_sample_directory(self) -> None:
        """MinKNOW writes it into the run directory, not the sample directory above it."""
        self._final_summary()
        self.assertFalse((self.sample / "final_summary_FAV12345_abcdef12_0123abcd.txt").exists())
        self.assertIs(self._inspect("final_summary_*.txt")[0], Readiness.READY)

    def test_a_literal_marker_name_still_works(self) -> None:
        (self.sample / "READY").write_text("")
        self.assertIs(self._inspect("READY")[0], Readiness.READY)

    def test_a_marker_matching_nothing_leaves_the_run_waiting(self) -> None:
        self._final_summary()
        self.assertIs(self._inspect("sequencing_finished.txt")[0], Readiness.MARKER_MISSING)

    def test_a_directory_matching_the_glob_is_not_a_marker(self) -> None:
        """Only a file signals completion; a directory of that name is not the producer's word."""
        (self.run / "final_summary_dir.txt").mkdir()
        self.assertIs(self._inspect("final_summary_*.txt")[0], Readiness.MARKER_MISSING)

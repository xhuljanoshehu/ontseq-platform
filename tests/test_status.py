from __future__ import annotations

import hashlib
import json
import os
import socket
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.pipeline.lock import LOCK_FILENAME
from ontseq_platform.pipeline.review import REVIEW_LOG, Decision, ReviewState, append_entry
from ontseq_platform.status import (
    RunState,
    exit_code,
    render_json,
    render_ledger,
    render_text,
    scan,
)


def _run_report(*, run_id: str, sample_id: str, passed: bool, unverified: list[str]) -> dict:
    """A minimal but valid RunReport payload.

    ``RunReport`` validates that ``unverified_stages`` names exactly the concluded stages
    whose adapter is unexecuted, so the basecall record is present only when the caller
    asks for it. Listing an unverified stage that is not in ``stages`` — or omitting one
    that is — makes the payload unloadable, and every test then sees UNREADABLE.
    """
    stages: list[dict] = [
        {
            "stage": "intake",
            "title": "Aligned-BAM integrity gate",
            "status": "COMPLETED",
            "verification": "verified_with_real_tool",
            "required": True,
            "reason": "gate returned PASS",
            "signature": "b" * 64,
        }
    ]
    if "basecall" in unverified:
        stages.append(
            {
                "stage": "basecall",
                "title": "Dorado basecalling",
                "status": "COMPLETED",
                "verification": "unverified_adapter",
                "required": True,
                "reason": "basecalled",
                "signature": "c" * 64,
            }
        )
    return {
        "schema_version": "0.1.0",
        "run_id": run_id,
        "sample_id": sample_id,
        "input_kind": "aligned_bam",
        "genome_build": "GRCh38",
        "manifest": {
            "schema_version": "0.1.0",
            "sample_id": sample_id,
            "run_id": run_id,
            "input": {
                "kind": "aligned_bam",
                "path": "/nowhere/in.bam",
                "index_path": "/nowhere/in.bam.bai",
            },
            "assay": {"mode": "lcwgs", "genome_build": "GRCh38", "reference_id": "REF_V1"},
            "analysis": {"profile": "lcwgs", "modules": ["qc"]},
        },
        "passed": passed,
        "verdict_reason": "PASS - everything ran" if passed else "FAIL - qc broke",
        "stages": stages,
        "pipeline_version": "0.0.0-test",
        "git_commit": "0" * 40,
        "started_at": "2026-08-17T09:00:00Z",
        "finished_at": "2026-08-17T09:10:00Z",
        "unverified_stages": unverified,
    }


class StatusCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.output = Path(self._temporary.name) / "runs"
        self.output.mkdir(parents=True)

    def _envelope(self, run_id: str, sample_id: str) -> Path:
        root = self.output / run_id / sample_id
        (root / "manifest").mkdir(parents=True, exist_ok=True)
        (root / "provenance").mkdir(parents=True, exist_ok=True)
        return root

    def _with_report(
        self,
        run_id: str = "RUN_001",
        sample_id: str = "SAMPLE_A",
        *,
        passed: bool = True,
        unverified: list[str] | None = None,
    ) -> Path:
        root = self._envelope(run_id, sample_id)
        (root / "provenance" / "run.json").write_text(
            json.dumps(
                _run_report(
                    run_id=run_id,
                    sample_id=sample_id,
                    passed=passed,
                    unverified=unverified or [],
                )
            ),
            encoding="utf-8",
        )
        return root

    def _lock(self, root: Path, *, pid: int, hostname: str | None = None) -> None:
        (root / LOCK_FILENAME).write_text(
            json.dumps(
                {
                    "pid": pid,
                    "hostname": hostname or socket.gethostname(),
                    "acquired_at": "2026-08-17T09:00:00+00:00",
                    "run_id": root.parent.name,
                    "sample_id": root.name,
                    "pipeline_version": "0.0.0-test",
                }
            ),
            encoding="utf-8",
        )

    @staticmethod
    def _dead_pid() -> int:
        process = subprocess.Popen([sys.executable, "-c", ""])
        process.wait()
        return process.pid


class ScanTests(StatusCase):
    def test_a_missing_output_directory_is_an_error(self) -> None:
        with self.assertRaises(NotADirectoryError):
            scan(self.output.parent / "nowhere")

    def test_an_empty_output_directory_yields_nothing(self) -> None:
        self.assertEqual(scan(self.output), [])

    def test_a_finished_passing_run_is_reported_as_passed(self) -> None:
        self._with_report()
        status = scan(self.output)[0]
        self.assertIs(status.state, RunState.PASSED)
        self.assertEqual(status.run_id, "RUN_001")
        self.assertEqual(status.sample_id, "SAMPLE_A")

    def test_a_failing_run_is_reported_as_failed(self) -> None:
        self._with_report(passed=False)
        self.assertIs(scan(self.output)[0].state, RunState.FAILED)

    def test_envelopes_are_returned_in_a_stable_order(self) -> None:
        self._with_report("RUN_B", "SAMPLE_2")
        self._with_report("RUN_A", "SAMPLE_1")
        self.assertEqual(
            [(item.run_id, item.sample_id) for item in scan(self.output)],
            [("RUN_A", "SAMPLE_1"), ("RUN_B", "SAMPLE_2")],
        )

    def test_a_single_run_can_be_selected(self) -> None:
        self._with_report("RUN_A", "SAMPLE_1")
        self._with_report("RUN_B", "SAMPLE_2")
        self.assertEqual([item.run_id for item in scan(self.output, run_id="RUN_B")], ["RUN_B"])

    def test_directories_that_are_not_envelopes_are_ignored(self) -> None:
        (self.output / "RUN_A" / "not-an-envelope").mkdir(parents=True)
        self.assertEqual(scan(self.output), [])

    def test_an_unreadable_report_is_reported_rather_than_skipped(self) -> None:
        """Silently omitting it would make a broken run look like a run that never existed."""
        root = self._envelope("RUN_001", "SAMPLE_A")
        (root / "provenance" / "run.json").write_text("{not json", encoding="utf-8")
        status = scan(self.output)[0]
        self.assertIs(status.state, RunState.UNREADABLE)
        self.assertIn("could not be read", status.detail)


class LockStateTests(StatusCase):
    def test_a_live_lock_is_reported_as_running(self) -> None:
        root = self._with_report()
        self._lock(root, pid=os.getpid())
        status = scan(self.output)[0]
        self.assertIs(status.state, RunState.RUNNING)
        self.assertTrue(status.holder_alive)

    def test_a_dead_lock_is_reported_as_interrupted(self) -> None:
        """Where a run died. Worth knowing before somebody deletes the directory."""
        root = self._with_report()
        self._lock(root, pid=self._dead_pid())
        status = scan(self.output)[0]
        self.assertIs(status.state, RunState.INTERRUPTED)
        self.assertFalse(status.holder_alive)
        self.assertIn("reclaims the lock", status.detail)

    def test_a_lock_from_another_host_is_running_with_unknown_liveness(self) -> None:
        """Unknown is the honest answer, and must not be collapsed into "probably fine"."""
        root = self._with_report()
        self._lock(root, pid=self._dead_pid(), hostname="some-other-machine")
        status = scan(self.output)[0]
        self.assertIs(status.state, RunState.RUNNING)
        self.assertIsNone(status.holder_alive)
        self.assertIn("another host", status.detail)

    def test_a_held_lock_outranks_the_report_on_disk(self) -> None:
        """The report is a snapshot from partway through the run that is still going."""
        root = self._with_report(passed=True)
        self._lock(root, pid=os.getpid())
        self.assertIs(scan(self.output)[0].state, RunState.RUNNING)

    def test_an_envelope_with_neither_lock_nor_report_is_unfinished(self) -> None:
        self._envelope("RUN_001", "SAMPLE_A")
        status = scan(self.output)[0]
        self.assertIs(status.state, RunState.UNFINISHED)
        self.assertIn("never reached a verdict", status.detail)


class RenderTests(StatusCase):
    def test_an_empty_scan_says_so(self) -> None:
        self.assertIn("no run envelopes", render_text([]))

    def test_the_summary_counts_every_state(self) -> None:
        self._with_report("RUN_A", "SAMPLE_1", passed=True)
        self._with_report("RUN_B", "SAMPLE_2", passed=False)
        rendered = render_text(scan(self.output))
        self.assertIn("1 passed", rendered)
        self.assertIn("1 failed", rendered)
        self.assertIn("total: 2 envelope(s)", rendered)

    def test_unverified_adapters_are_surfaced_in_the_summary(self) -> None:
        """A completed run resting on an unexecuted adapter must not read as simply fine."""
        self._with_report(unverified=["basecall"])
        self.assertIn("UNVERIFIED ADAPTERS COMPLETED: basecall", render_text(scan(self.output)))

    def test_the_lock_holder_is_shown_when_one_exists(self) -> None:
        root = self._with_report()
        self._lock(root, pid=os.getpid())
        self.assertIn(f"pid {os.getpid()}", render_text(scan(self.output)))

    def test_stages_appear_only_when_verbose(self) -> None:
        self._with_report()
        self.assertNotIn("intake", render_text(scan(self.output)))
        self.assertIn("intake", render_text(scan(self.output), verbose=True))

    def test_json_output_is_parseable_and_carries_the_state(self) -> None:
        self._with_report(passed=False, unverified=["basecall"])
        payload = json.loads(render_json(scan(self.output)))
        self.assertEqual(payload[0]["state"], "failed")
        self.assertEqual(payload[0]["run_id"], "RUN_001")
        self.assertEqual(len(payload[0]["stages"]), 2)

    def test_json_records_unknown_liveness_as_null_not_false(self) -> None:
        root = self._with_report()
        self._lock(root, pid=self._dead_pid(), hostname="elsewhere")
        payload = json.loads(render_json(scan(self.output)))
        self.assertIsNone(payload[0]["holder"]["alive"])


class LedgerRenderTests(StatusCase):
    def test_no_ledger_renders_nothing(self) -> None:
        self.assertEqual(render_ledger(self.output), "")

    def test_a_ledger_is_summarised(self) -> None:
        (self.output / ".ontseq-watch.json").write_text(
            json.dumps(
                {
                    "attempts": {
                        "SAMPLE_A": {
                            "outcome": "completed",
                            "detail": "PASS",
                            "attempted_at": "t",
                        }
                    }
                }
            ),
            encoding="utf-8",
        )
        rendered = render_ledger(self.output)
        self.assertIn("SAMPLE_A", rendered)
        self.assertIn("COMPLETED", rendered)


class ExitCodeTests(StatusCase):
    def test_all_passing_is_zero(self) -> None:
        self._with_report()
        self.assertEqual(exit_code(scan(self.output)), 0)

    def test_a_failure_is_two(self) -> None:
        self._with_report(passed=False)
        self.assertEqual(exit_code(scan(self.output)), 2)

    def test_an_interrupted_run_is_six(self) -> None:
        root = self._with_report()
        self._lock(root, pid=self._dead_pid())
        self.assertEqual(exit_code(scan(self.output)), 6)

    def test_a_run_in_progress_is_not_a_problem(self) -> None:
        """A check that fires while the pipeline works trains people to ignore it."""
        root = self._with_report()
        self._lock(root, pid=os.getpid())
        self.assertEqual(exit_code(scan(self.output)), 0)

    def test_a_failure_outranks_an_interruption(self) -> None:
        self._with_report("RUN_A", "SAMPLE_1", passed=False)
        root = self._with_report("RUN_B", "SAMPLE_2")
        self._lock(root, pid=self._dead_pid())
        self.assertEqual(exit_code(scan(self.output)), 2)

    def test_an_empty_output_directory_is_not_an_alert(self) -> None:
        self.assertEqual(exit_code([]), 0)


class ReviewVisibilityTests(StatusCase):
    """An operator scanning many envelopes should not need a second command per envelope."""

    def _release(self, root: Path) -> str:
        (root / "release").mkdir(parents=True, exist_ok=True)
        path = root / "release" / "release.json"
        path.write_text('{"synthetic": true}', encoding="utf-8")
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def _sign_off(self, root: Path, *, digest: str, decision: str = "accepted") -> None:
        append_entry(
            root / REVIEW_LOG,
            decision=Decision(decision),
            reviewer="dr.mueller",
            run_id=root.parent.name,
            sample_id=root.name,
            release_sha256=digest,
        )

    def test_an_envelope_nobody_reviewed_reports_no_review(self) -> None:
        self._with_report()
        self.assertIsNone(scan(self.output)[0].review)

    def test_an_accepted_envelope_is_shown(self) -> None:
        root = self._with_report()
        self._sign_off(root, digest=self._release(root))
        status = scan(self.output)[0]
        self.assertIs(status.review, ReviewState.ACCEPTED)
        self.assertIn("review: ACCEPTED", render_text([status]))

    def test_a_changed_release_makes_the_review_stale_here_too(self) -> None:
        """The two surfaces must agree; a stale review shown as accepted would be worse."""
        root = self._with_report()
        self._sign_off(root, digest=self._release(root))
        (root / "release" / "release.json").write_text('{"synthetic": false}', encoding="utf-8")
        self.assertIs(scan(self.output)[0].review, ReviewState.STALE)

    def test_the_review_state_reaches_the_json(self) -> None:
        root = self._with_report()
        self._sign_off(root, digest=self._release(root))
        payload = json.loads(render_json(scan(self.output)))
        self.assertEqual(payload[0]["review"], "accepted")

    def test_an_unreviewed_run_is_still_exit_code_zero(self) -> None:
        """`status` answers whether the runs worked, not whether they may be released."""
        self._with_report()
        self.assertEqual(exit_code(scan(self.output)), 0)

    def test_a_rejected_review_does_not_change_the_run_exit_code_either(self) -> None:
        root = self._with_report()
        self._sign_off(root, digest=self._release(root), decision="rejected")
        self.assertIs(scan(self.output)[0].review, ReviewState.REJECTED)
        self.assertEqual(exit_code(scan(self.output)), 0)


if __name__ == "__main__":
    unittest.main()

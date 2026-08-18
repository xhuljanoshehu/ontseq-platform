from __future__ import annotations

import json
import tempfile
import unittest
from datetime import UTC, datetime
from pathlib import Path

from ontseq_platform.pipeline.review import (
    Decision,
    ReviewError,
    ReviewState,
    accepted_reviewers,
    append_entry,
    current_state,
    entry_digest,
    exit_code,
    read_log,
    verify_chain,
)

BUNDLE = "a" * 64
OTHER_BUNDLE = "b" * 64


class ReviewCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.log = Path(self._temporary.name) / "review" / "review.log.jsonl"

    def record(
        self,
        decision: Decision = Decision.ACCEPTED,
        *,
        reviewer: str = "dr.mueller",
        release: str = BUNDLE,
        note: str = "",
        when: str = "2026-08-18T09:00:00+00:00",
    ):
        return append_entry(
            self.log,
            decision=decision,
            reviewer=reviewer,
            run_id="RUN_001",
            sample_id="AML_00123",
            release_sha256=release,
            note=note,
            now=datetime.fromisoformat(when).astimezone(UTC),
        )

    def lines(self) -> list[dict]:
        return [json.loads(line) for line in self.log.read_text().splitlines() if line.strip()]

    def rewrite(self, payloads: list[dict]) -> None:
        lines = [json.dumps(item, sort_keys=True, separators=(",", ":")) for item in payloads]
        self.log.write_text("\n".join(lines) + "\n")


class AppendTests(ReviewCase):
    def test_a_missing_log_is_an_empty_history_not_an_error(self) -> None:
        self.assertEqual(read_log(self.log), [])

    def test_the_first_entry_has_no_predecessor(self) -> None:
        entry = self.record()
        self.assertEqual(entry.sequence, 1)
        self.assertIsNone(entry.previous_entry_sha256)

    def test_each_entry_names_the_one_before_it(self) -> None:
        first = self.record()
        second = self.record(Decision.REJECTED, reviewer="dr.schmidt", note="repeat")
        self.assertEqual(second.previous_entry_sha256, first.entry_sha256)
        self.assertEqual(second.sequence, 2)

    def test_a_later_review_does_not_replace_an_earlier_one(self) -> None:
        """An audit trail whose history can be overwritten is not an audit trail."""
        self.record()
        self.record(Decision.REJECTED, reviewer="dr.schmidt")
        self.assertEqual(len(read_log(self.log)), 2)

    def test_an_anonymous_review_is_refused(self) -> None:
        with self.assertRaises(ReviewError):
            self.record(reviewer="   ")

    def test_the_identity_is_recorded_as_asserted_not_authenticated(self) -> None:
        """Nothing here authenticates anybody, so the record must not imply otherwise."""
        self.assertEqual(self.record().identity_source, "asserted")

    def test_the_review_binds_to_what_was_reviewed(self) -> None:
        self.assertEqual(self.record(release=BUNDLE).release_sha256, BUNDLE)


class ChainTests(ReviewCase):
    def test_an_untouched_log_verifies(self) -> None:
        self.record()
        self.record(reviewer="dr.schmidt")
        intact, _ = verify_chain(read_log(self.log))
        self.assertTrue(intact)

    def test_an_empty_log_verifies(self) -> None:
        self.assertTrue(verify_chain([])[0])

    def test_editing_an_entry_breaks_the_chain(self) -> None:
        self.record(note="looks fine")
        payloads = self.lines()
        payloads[0]["note"] = "looks terrible"
        self.rewrite(payloads)
        intact, detail = verify_chain(read_log(self.log))
        self.assertFalse(intact)
        self.assertIn("edited", detail)

    def test_removing_an_entry_breaks_the_chain(self) -> None:
        self.record()
        self.record(reviewer="dr.schmidt")
        self.record(reviewer="dr.weber")
        payloads = self.lines()
        self.rewrite([payloads[0], payloads[2]])
        self.assertFalse(verify_chain(read_log(self.log))[0])

    def test_reordering_entries_breaks_the_chain(self) -> None:
        self.record()
        self.record(reviewer="dr.schmidt")
        payloads = self.lines()
        self.rewrite([payloads[1], payloads[0]])
        self.assertFalse(verify_chain(read_log(self.log))[0])

    def test_appending_to_a_broken_log_is_refused(self) -> None:
        """A record that looks continuous and is not would be worse than refusing."""
        self.record(note="original")
        payloads = self.lines()
        payloads[0]["note"] = "tampered"
        self.rewrite(payloads)
        with self.assertRaises(ReviewError):
            self.record(reviewer="dr.schmidt")

    def test_a_truncated_final_line_is_reported_not_ignored(self) -> None:
        self.record()
        self.log.write_text(self.log.read_text()[:-20])
        with self.assertRaises(ReviewError):
            read_log(self.log)

    def test_the_digest_ignores_only_itself(self) -> None:
        entry = self.record()
        payload = self.lines()[0]
        self.assertEqual(entry_digest(payload), entry.entry_sha256)
        payload["reviewer"] = "somebody.else"
        self.assertNotEqual(entry_digest(payload), entry.entry_sha256)


class StateTests(ReviewCase):
    def test_no_review_is_pending(self) -> None:
        state, detail = current_state([], BUNDLE)
        self.assertIs(state, ReviewState.PENDING)
        self.assertIn("no review", detail)

    def test_an_acceptance_of_the_current_content_is_accepted(self) -> None:
        self.record()
        self.assertIs(current_state(read_log(self.log), BUNDLE)[0], ReviewState.ACCEPTED)

    def test_a_rejection_of_the_current_content_is_rejected(self) -> None:
        self.record(Decision.REJECTED, note="coverage too low")
        state, detail = current_state(read_log(self.log), BUNDLE)
        self.assertIs(state, ReviewState.REJECTED)
        self.assertIn("coverage too low", detail)

    def test_changed_content_makes_an_acceptance_stale_not_valid(self) -> None:
        """The judgement stands for what it saw; it must not carry over to something else."""
        self.record()
        state, detail = current_state(read_log(self.log), OTHER_BUNDLE)
        self.assertIs(state, ReviewState.STALE)
        self.assertIn("says nothing about the content here now", detail)

    def test_a_vanished_release_bundle_is_stale(self) -> None:
        self.record()
        self.assertIs(current_state(read_log(self.log), None)[0], ReviewState.STALE)

    def test_the_latest_decision_is_the_current_one(self) -> None:
        self.record(Decision.ACCEPTED)
        self.record(Decision.REJECTED, reviewer="dr.schmidt", note="withdrawn")
        self.assertIs(current_state(read_log(self.log), BUNDLE)[0], ReviewState.REJECTED)

    def test_a_broken_chain_outranks_the_decision_in_it(self) -> None:
        self.record()
        payloads = self.lines()
        payloads[0]["decision"] = "rejected"
        self.rewrite(payloads)
        self.assertIs(current_state(read_log(self.log), BUNDLE)[0], ReviewState.BROKEN)


class FourEyesTests(ReviewCase):
    def test_two_reviewers_of_the_same_content_both_count(self) -> None:
        self.record(reviewer="dr.mueller")
        self.record(reviewer="dr.schmidt")
        self.assertEqual(len(accepted_reviewers(read_log(self.log), BUNDLE)), 2)

    def test_one_person_accepting_twice_is_still_one_reviewer(self) -> None:
        """Otherwise a four-eyes rule is satisfied by running the same command again."""
        self.record(reviewer="dr.mueller")
        self.record(reviewer="dr.mueller")
        self.assertEqual(accepted_reviewers(read_log(self.log), BUNDLE), ("dr.mueller",))

    def test_an_acceptance_of_other_content_does_not_count(self) -> None:
        self.record(reviewer="dr.mueller", release=OTHER_BUNDLE)
        self.record(reviewer="dr.schmidt", release=BUNDLE)
        self.assertEqual(accepted_reviewers(read_log(self.log), BUNDLE), ("dr.schmidt",))

    def test_a_rejection_is_not_an_acceptance(self) -> None:
        self.record(Decision.REJECTED, reviewer="dr.mueller")
        self.assertEqual(accepted_reviewers(read_log(self.log), BUNDLE), ())


class ExitCodeTests(unittest.TestCase):
    def test_accepted_is_zero(self) -> None:
        self.assertEqual(exit_code(ReviewState.ACCEPTED), 0)

    def test_rejected_is_two(self) -> None:
        self.assertEqual(exit_code(ReviewState.REJECTED), 2)

    def test_a_broken_trail_is_two(self) -> None:
        self.assertEqual(exit_code(ReviewState.BROKEN), 2)

    def test_pending_is_six(self) -> None:
        self.assertEqual(exit_code(ReviewState.PENDING), 6)

    def test_stale_is_six(self) -> None:
        self.assertEqual(exit_code(ReviewState.STALE), 6)

    def test_too_few_reviewers_blocks_an_otherwise_accepted_run(self) -> None:
        self.assertEqual(exit_code(ReviewState.ACCEPTED, reviewers=1, required_reviewers=2), 6)

    def test_enough_reviewers_passes(self) -> None:
        self.assertEqual(exit_code(ReviewState.ACCEPTED, reviewers=2, required_reviewers=2), 0)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import json
import unittest

from ontseq_platform.pipeline.checks import (
    Check,
    CheckList,
    CheckStatus,
    exit_code,
    render_json,
    render_text,
    required_tools,
    worst,
)
from ontseq_platform.pipeline.stages import SPEC_BY_STAGE, InputKindName, StageId


class RequiredToolTests(unittest.TestCase):
    def _named(self, kind: InputKindName) -> dict[str, bool]:
        return {item.name: item.required for item in required_tools(kind)}

    def test_an_aligned_bam_run_needs_no_basecaller_or_aligner(self) -> None:
        tools = self._named(InputKindName.ALIGNED_BAM)
        self.assertNotIn("dorado", tools)
        self.assertNotIn("minimap2", tools)
        self.assertIn("samtools", tools)
        self.assertIn("cramino", tools)

    def test_an_unaligned_bam_run_needs_the_aligner_but_not_the_basecaller(self) -> None:
        tools = self._named(InputKindName.UNALIGNED_BAM)
        self.assertNotIn("dorado", tools)
        self.assertTrue(tools["minimap2"])

    def test_a_pod5_run_needs_the_basecaller(self) -> None:
        self.assertTrue(self._named(InputKindName.POD5)["dorado"])

    def test_a_tool_serving_only_an_optional_stage_is_not_required(self) -> None:
        """Sniffles serves SV alone, and SV is optional, so its absence must not block."""
        self.assertFalse(SPEC_BY_STAGE[StageId.SV].required)
        self.assertFalse(self._named(InputKindName.ALIGNED_BAM)["sniffles"])

    def test_a_tool_shared_by_a_required_and_an_optional_stage_is_required(self) -> None:
        tools = {item.name: item for item in required_tools(InputKindName.UNALIGNED_BAM)}
        self.assertIn(StageId.ALIGN, tools["samtools"].stages)
        self.assertIn(StageId.INTAKE, tools["samtools"].stages)
        self.assertTrue(tools["samtools"].required)

    def test_requirements_are_returned_in_a_stable_order(self) -> None:
        names = [item.name for item in required_tools(InputKindName.POD5)]
        self.assertEqual(names, sorted(names))

    def test_every_tool_names_at_least_one_planned_stage(self) -> None:
        for kind in InputKindName:
            for requirement in required_tools(kind):
                self.assertTrue(requirement.stages, requirement.name)

    def test_the_expected_version_can_be_attached_without_mutating(self) -> None:
        original = required_tools(InputKindName.ALIGNED_BAM)[0]
        updated = original.with_expected_version("1.2.3")
        self.assertIsNone(original.expected_version)
        self.assertEqual(updated.expected_version, "1.2.3")
        self.assertEqual(updated.name, original.name)


class WorstTests(unittest.TestCase):
    @staticmethod
    def _of(*statuses: CheckStatus) -> list[Check]:
        return [
            Check(name=str(index), status=status, detail="")
            for index, status in enumerate(statuses)
        ]

    def test_failed_outranks_everything(self) -> None:
        self.assertIs(
            worst(self._of(CheckStatus.OK, CheckStatus.WARNING, CheckStatus.FAILED)),
            CheckStatus.FAILED,
        )

    def test_warning_outranks_unknown(self) -> None:
        self.assertIs(
            worst(self._of(CheckStatus.UNKNOWN, CheckStatus.WARNING)), CheckStatus.WARNING
        )

    def test_unknown_outranks_ok(self) -> None:
        """Not having looked must never be reported as having looked and found it fine."""
        self.assertIs(worst(self._of(CheckStatus.OK, CheckStatus.UNKNOWN)), CheckStatus.UNKNOWN)

    def test_only_skipped_checks_report_skipped(self) -> None:
        self.assertIs(worst(self._of(CheckStatus.SKIPPED)), CheckStatus.SKIPPED)

    def test_no_checks_at_all_is_skipped(self) -> None:
        self.assertIs(worst([]), CheckStatus.SKIPPED)


class ExitCodeTests(unittest.TestCase):
    @staticmethod
    def _of(*statuses: CheckStatus) -> list[Check]:
        return [
            Check(name=str(index), status=status, detail="")
            for index, status in enumerate(statuses)
        ]

    def test_a_failure_is_two(self) -> None:
        self.assertEqual(exit_code(self._of(CheckStatus.OK, CheckStatus.FAILED)), 2)

    def test_a_warning_does_not_block(self) -> None:
        self.assertEqual(exit_code(self._of(CheckStatus.WARNING)), 0)

    def test_an_unknown_does_not_block(self) -> None:
        """Refusing to start because a GPU is unknowable from here makes this unusable."""
        self.assertEqual(exit_code(self._of(CheckStatus.UNKNOWN)), 0)

    def test_nothing_checked_is_not_a_failure(self) -> None:
        self.assertEqual(exit_code([]), 0)


class CheckListTests(unittest.TestCase):
    def test_checks_keep_the_order_they_were_made_in(self) -> None:
        checks = CheckList()
        checks.ok("first", "")
        checks.failed("second", "")
        checks.warning("third", "")
        self.assertEqual([item.name for item in checks.checks], ["first", "second", "third"])

    def test_only_failed_checks_are_blocking(self) -> None:
        checks = CheckList()
        self.assertTrue(checks.failed("a", "").blocking)
        self.assertFalse(checks.warning("b", "").blocking)
        self.assertFalse(checks.unknown("c", "").blocking)
        self.assertFalse(checks.skipped("d", "").blocking)
        self.assertFalse(checks.ok("e", "").blocking)

    def test_a_stage_can_be_attached_to_a_check(self) -> None:
        checks = CheckList()
        self.assertIs(checks.failed("a", "", stage=StageId.ALIGN).stage, StageId.ALIGN)


class RenderTests(unittest.TestCase):
    def _sample(self) -> list[Check]:
        checks = CheckList()
        checks.ok("input.bam", "1.20 GiB")
        checks.failed("tool.samtools", "not on PATH", remedy="install samtools")
        checks.skipped("basecall.model", "this run does not basecall")
        checks.unknown("disk.free", "42.0 GiB free")
        return checks.checks

    def test_skipped_checks_are_hidden_by_default(self) -> None:
        rendered = render_text(self._sample())
        self.assertNotIn("basecall.model", rendered)
        self.assertIn("basecall.model", render_text(self._sample(), verbose=True))

    def test_passing_checks_are_shown(self) -> None:
        """Printing only problems leaves a reader unable to tell clean from incomplete."""
        self.assertIn("input.bam", render_text(self._sample()))

    def test_a_remedy_is_shown_under_its_check(self) -> None:
        self.assertIn("install samtools", render_text(self._sample()))

    def test_the_summary_reports_the_worst_status(self) -> None:
        self.assertIn("preflight: FAILED", render_text(self._sample()))

    def test_the_summary_counts_skipped_checks_it_did_not_print(self) -> None:
        self.assertIn("1 skipped", render_text(self._sample()))

    def test_nothing_to_check_says_so(self) -> None:
        self.assertEqual(render_text([]), "no checks apply")

    def test_json_names_every_blocking_check(self) -> None:
        payload = json.loads(render_json(self._sample()))
        self.assertEqual(payload["blocking"], ["tool.samtools"])
        self.assertEqual(payload["verdict"], "failed")

    def test_json_keeps_skipped_checks_a_reader_would_not_see(self) -> None:
        payload = json.loads(render_json(self._sample()))
        self.assertEqual(len(payload["checks"]), 4)

    def test_json_records_a_missing_stage_as_null(self) -> None:
        payload = json.loads(render_json(self._sample()))
        self.assertIsNone(payload["checks"][0]["stage"])


if __name__ == "__main__":
    unittest.main()

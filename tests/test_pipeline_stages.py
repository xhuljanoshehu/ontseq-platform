from __future__ import annotations

import unittest

from ontseq_platform.pipeline.stages import (
    EXECUTION_ORDER,
    SPEC_BY_STAGE,
    STAGE_SPECS,
    InputKindName,
    StageId,
    StageOutcome,
    VerificationStatus,
    applies_to,
    blocking_dependency,
    effective_dependencies,
    planned_stages,
    summarize,
    unverified_specs,
)


class GraphIntegrityTests(unittest.TestCase):
    def test_execution_order_covers_every_stage_exactly_once(self) -> None:
        self.assertEqual(len(EXECUTION_ORDER), len(STAGE_SPECS))
        self.assertEqual(set(EXECUTION_ORDER), set(SPEC_BY_STAGE))

    def test_dependencies_precede_dependents(self) -> None:
        position = {stage: index for index, stage in enumerate(EXECUTION_ORDER)}
        for spec in STAGE_SPECS:
            for dependency in spec.depends_on:
                self.assertLess(
                    position[dependency],
                    position[spec.stage],
                    f"{dependency} must precede {spec.stage}",
                )

    def test_every_stage_declares_a_purpose_and_verification(self) -> None:
        for spec in STAGE_SPECS:
            self.assertTrue(spec.purpose, f"{spec.stage} has no purpose")
            self.assertIsInstance(spec.verification, VerificationStatus)

    def test_non_applicable_stages_explain_themselves(self) -> None:
        for spec in STAGE_SPECS:
            if spec.applicable_for != frozenset(InputKindName):
                self.assertTrue(
                    spec.not_applicable_reason,
                    f"{spec.stage} can be skipped but does not say why",
                )


class ApplicabilityTests(unittest.TestCase):
    def test_aligned_bam_run_skips_basecalling_and_alignment(self) -> None:
        planned = planned_stages(InputKindName.ALIGNED_BAM)
        self.assertNotIn(StageId.BASECALL, planned)
        self.assertNotIn(StageId.ALIGN, planned)
        self.assertIn(StageId.INTAKE, planned)

    def test_unaligned_bam_run_aligns_but_does_not_basecall(self) -> None:
        planned = planned_stages(InputKindName.UNALIGNED_BAM)
        self.assertNotIn(StageId.BASECALL, planned)
        self.assertIn(StageId.ALIGN, planned)

    def test_pod5_run_includes_the_whole_chain(self) -> None:
        planned = planned_stages(InputKindName.POD5)
        self.assertIn(StageId.BASECALL, planned)
        self.assertIn(StageId.ALIGN, planned)
        self.assertIn(StageId.RELEASE, planned)

    def test_planned_stages_are_returned_in_execution_order(self) -> None:
        for kind in InputKindName:
            planned = planned_stages(kind)
            self.assertEqual(list(planned), [s for s in EXECUTION_ORDER if s in set(planned)])


class DependencyResolutionTests(unittest.TestCase):
    """Skipped stages must be bridged, not waited on."""

    def test_intake_waits_for_alignment_when_alignment_applies(self) -> None:
        self.assertEqual(
            effective_dependencies(StageId.INTAKE, InputKindName.UNALIGNED_BAM),
            (StageId.ALIGN,),
        )

    def test_intake_has_no_dependency_for_an_aligned_bam_run(self) -> None:
        self.assertEqual(effective_dependencies(StageId.INTAKE, InputKindName.ALIGNED_BAM), ())

    def test_alignment_waits_for_basecalling_only_from_pod5(self) -> None:
        self.assertEqual(
            effective_dependencies(StageId.ALIGN, InputKindName.POD5), (StageId.BASECALL,)
        )
        self.assertEqual(effective_dependencies(StageId.ALIGN, InputKindName.UNALIGNED_BAM), ())

    def test_resolution_is_transitive_through_two_skipped_stages(self) -> None:
        # intake -> align -> basecall, and neither applies to an aligned-BAM run.
        self.assertEqual(effective_dependencies(StageId.INTAKE, InputKindName.ALIGNED_BAM), ())

    def test_applies_to_matches_the_spec(self) -> None:
        self.assertTrue(applies_to(StageId.BASECALL, InputKindName.POD5))
        self.assertFalse(applies_to(StageId.BASECALL, InputKindName.ALIGNED_BAM))


class FailurePropagationTests(unittest.TestCase):
    """A stage downstream of a problem is NOT_RUN, never FAILED."""

    def test_stage_with_satisfied_dependencies_may_run(self) -> None:
        outcomes = {StageId.INTAKE: StageOutcome.COMPLETED}
        self.assertIsNone(blocking_dependency(StageId.QC, InputKindName.ALIGNED_BAM, outcomes))

    def test_failed_dependency_blocks_and_names_the_cause(self) -> None:
        outcomes = {StageId.INTAKE: StageOutcome.FAILED}
        blocked = blocking_dependency(StageId.QC, InputKindName.ALIGNED_BAM, outcomes)
        self.assertIsNotNone(blocked)
        self.assertEqual(blocked.blocking_stage, StageId.INTAKE)
        self.assertIn("never started", blocked.describe())

    def test_missing_dependency_outcome_blocks(self) -> None:
        blocked = blocking_dependency(StageId.QC, InputKindName.ALIGNED_BAM, {})
        self.assertIsNotNone(blocked)
        self.assertEqual(blocked.blocking_outcome, StageOutcome.NOT_RUN)

    def test_no_call_dependency_blocks_downstream(self) -> None:
        """NO_CALL means the stage produced no artifact, so a dependent cannot proceed."""
        outcomes = {StageId.INTAKE: StageOutcome.NO_CALL}
        self.assertIsNotNone(blocking_dependency(StageId.QC, InputKindName.ALIGNED_BAM, outcomes))


class VerdictTests(unittest.TestCase):
    def _complete_aligned_run(self) -> dict[StageId, StageOutcome]:
        return {
            stage: StageOutcome.COMPLETED for stage in planned_stages(InputKindName.ALIGNED_BAM)
        }

    def test_complete_run_passes(self) -> None:
        verdict = summarize(InputKindName.ALIGNED_BAM, self._complete_aligned_run())
        self.assertTrue(verdict.passed)
        self.assertEqual(verdict.failed_stages, ())
        self.assertIn("every applicable stage", verdict.describe())

    def test_failed_required_stage_fails_the_run(self) -> None:
        outcomes = self._complete_aligned_run()
        outcomes[StageId.QC] = StageOutcome.FAILED
        verdict = summarize(InputKindName.ALIGNED_BAM, outcomes)
        self.assertFalse(verdict.passed)
        self.assertEqual(verdict.failed_stages, (StageId.QC,))
        self.assertIn("FAILED", verdict.describe())

    def test_incomplete_required_stage_fails_the_run(self) -> None:
        outcomes = self._complete_aligned_run()
        outcomes[StageId.REPORT] = StageOutcome.NOT_RUN
        verdict = summarize(InputKindName.ALIGNED_BAM, outcomes)
        self.assertFalse(verdict.passed)
        self.assertIn(StageId.REPORT, verdict.incomplete_required_stages)

    def test_skipped_optional_stage_does_not_fail_the_run_but_is_reported(self) -> None:
        outcomes = self._complete_aligned_run()
        outcomes[StageId.CNV] = StageOutcome.NOT_RUN
        verdict = summarize(InputKindName.ALIGNED_BAM, outcomes)
        self.assertTrue(verdict.passed)
        self.assertIn(StageId.CNV, verdict.skipped_optional_stages)
        self.assertIn("not a negative biological finding", verdict.describe())

    def test_no_call_counts_as_a_completed_stage(self) -> None:
        """A caller that looked and declined has done its job."""
        outcomes = self._complete_aligned_run()
        outcomes[StageId.SV] = StageOutcome.NO_CALL
        verdict = summarize(InputKindName.ALIGNED_BAM, outcomes)
        self.assertTrue(verdict.passed)
        self.assertNotIn(StageId.SV, verdict.skipped_optional_stages)

    def test_unverified_stages_are_surfaced_even_on_a_passing_run(self) -> None:
        outcomes = {stage: StageOutcome.COMPLETED for stage in planned_stages(InputKindName.POD5)}
        verdict = summarize(InputKindName.POD5, outcomes)
        self.assertTrue(verdict.passed)
        # Basecalling completed, but its adapter has never met a real Dorado binary.
        self.assertIn(StageId.BASECALL, verdict.unverified_stages)

    def test_aligned_run_does_not_report_basecalling_as_incomplete(self) -> None:
        verdict = summarize(InputKindName.ALIGNED_BAM, self._complete_aligned_run())
        self.assertNotIn(StageId.BASECALL, verdict.incomplete_required_stages)
        self.assertNotIn(StageId.BASECALL, verdict.skipped_optional_stages)


class VerificationTests(unittest.TestCase):
    def test_basecalling_is_declared_unverified(self) -> None:
        self.assertEqual(
            SPEC_BY_STAGE[StageId.BASECALL].verification,
            VerificationStatus.UNVERIFIED_ADAPTER,
        )

    def test_real_tool_stages_are_declared_as_such(self) -> None:
        """Verification is claimed only where CI actually runs the binary."""
        for stage in (StageId.ALIGN, StageId.INTAKE, StageId.QC, StageId.SV):
            self.assertEqual(
                SPEC_BY_STAGE[stage].verification,
                VerificationStatus.VERIFIED_WITH_REAL_TOOL,
            )

    def test_unverified_specs_filters_correctly(self) -> None:
        specs = unverified_specs(planned_stages(InputKindName.POD5))
        stages = {spec.stage for spec in specs}
        self.assertIn(StageId.BASECALL, stages)
        self.assertNotIn(StageId.INTAKE, stages)
        self.assertNotIn(StageId.ALIGN, stages)

    def test_an_unaligned_bam_run_flags_only_the_unwired_stages(self) -> None:
        """Below POD5 every wired adapter is exercised by CI; only stubs remain."""
        specs = unverified_specs(planned_stages(InputKindName.UNALIGNED_BAM))
        self.assertEqual(
            {spec.stage for spec in specs}, {StageId.TARGET_COVERAGE, StageId.CNV}
        )


if __name__ == "__main__":
    unittest.main()

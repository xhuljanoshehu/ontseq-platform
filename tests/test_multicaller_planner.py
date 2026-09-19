from __future__ import annotations

import hashlib
import unittest

from ontseq_platform.multicaller_contracts import (
    CallerAssayRegime,
    CallerInputArtifact,
    CallerInputRole,
    CallerLaneRequest,
    CallerMode,
    CallerPlanningDecision,
)
from ontseq_platform.multicaller_planner import plan_caller_lane, seal_multicaller_plan


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _input(role: CallerInputRole, suffix: str | None = None) -> CallerInputArtifact:
    label = suffix or role.value
    return CallerInputArtifact(
        artifact_id=f"synthetic-{label}",
        role=role,
        sha256=_sha(label),
    )


def _request(
    mode: CallerMode,
    regime: CallerAssayRegime,
    *inputs: CallerInputArtifact,
    requested: bool = True,
    parents: tuple[str, ...] = (),
) -> CallerLaneRequest:
    return CallerLaneRequest(
        lane_id=f"lane-{mode.value.replace(':', '-')}",
        mode_id=mode,
        assay_regime=regime,
        inputs=list(inputs),
        parent_lane_ids=list(parents),
        requested=requested,
    )


class MultiCallerEligibilityPlannerTests(unittest.TestCase):
    def test_savana_paired_missing_normal_is_ineligible_without_mode_fallback(self) -> None:
        request = _request(
            CallerMode.SAVANA_PAIRED,
            CallerAssayRegime.TUMOR_NORMAL_LONG_READ,
            _input(CallerInputRole.TUMOR_BAM),
            _input(CallerInputRole.SNP_VCF),
        )

        planned = plan_caller_lane(request, runtime_available=True)

        self.assertEqual(planned.mode_id, CallerMode.SAVANA_PAIRED)
        self.assertEqual(planned.decision, CallerPlanningDecision.INELIGIBLE)
        self.assertTrue(any("matched normal" in reason.lower() for reason in planned.reasons))
        self.assertNotIn("tumor_only", " ".join(planned.reasons).lower())

    def test_missing_runtime_is_explicit_not_a_no_call(self) -> None:
        request = _request(
            CallerMode.QDNASEQ_ACE_MULTIBIN,
            CallerAssayRegime.LCWGS,
            _input(CallerInputRole.ALIGNED_BAM),
        )

        planned = plan_caller_lane(request, runtime_available=False)

        self.assertEqual(planned.decision, CallerPlanningDecision.RUNTIME_UNAVAILABLE)
        self.assertTrue(any("runtime" in reason.lower() for reason in planned.reasons))

    def test_adaptive_sampling_does_not_satisfy_ichorcna_ulp_wgs_regime(self) -> None:
        request = _request(
            CallerMode.ICHORCNA_ULP_WGS,
            CallerAssayRegime.ADAPTIVE_SAMPLING,
            _input(CallerInputRole.COVERAGE_PROFILE),
        )

        planned = plan_caller_lane(request, runtime_available=True)

        self.assertEqual(planned.decision, CallerPlanningDecision.INELIGIBLE)
        self.assertTrue(any("regime" in reason.lower() for reason in planned.reasons))

    def test_wakhan_requires_explicit_phased_input(self) -> None:
        request = _request(
            CallerMode.WAKHAN_PHASED_CNA,
            CallerAssayRegime.TUMOR_ONLY_LONG_READ,
            _input(CallerInputRole.TUMOR_BAM),
        )

        planned = plan_caller_lane(request, runtime_available=True)

        self.assertEqual(planned.decision, CallerPlanningDecision.INELIGIBLE)
        self.assertTrue(any("phased" in reason.lower() for reason in planned.reasons))

    def test_unrequested_lane_stays_not_requested_even_if_runtime_exists(self) -> None:
        request = _request(
            CallerMode.QDNASEQ_ACE_MULTIBIN,
            CallerAssayRegime.LCWGS,
            _input(CallerInputRole.ALIGNED_BAM),
            requested=False,
        )

        planned = plan_caller_lane(request, runtime_available=True)

        self.assertEqual(planned.decision, CallerPlanningDecision.NOT_REQUESTED)
        self.assertTrue(any("not requested" in reason.lower() for reason in planned.reasons))

    def test_lane_request_rejects_ambiguous_duplicate_input_roles(self) -> None:
        with self.assertRaises(ValueError):
            CallerLaneRequest(
                lane_id="lane-duplicate-input",
                mode_id=CallerMode.QDNASEQ_ACE_MULTIBIN,
                assay_regime=CallerAssayRegime.LCWGS,
                inputs=[
                    _input(CallerInputRole.ALIGNED_BAM, "one"),
                    _input(CallerInputRole.ALIGNED_BAM, "two"),
                ],
            )


class MultiCallerPlanSealingTests(unittest.TestCase):
    def _runtime_map(self, *modes: CallerMode) -> dict[CallerMode, bool]:
        return {mode: True for mode in modes}

    def test_dependency_cycle_is_rejected_before_execution(self) -> None:
        sniffles = _request(
            CallerMode.SNIFFLES2_STANDARD,
            CallerAssayRegime.ADAPTIVE_SAMPLING,
            _input(CallerInputRole.ALIGNED_BAM, "sniffles-bam"),
            parents=("lane-spectre-sniffles-supported",),
        ).model_copy(update={"lane_id": "lane-sniffles"})
        spectre = _request(
            CallerMode.SPECTRE_SNIFFLES_SUPPORTED,
            CallerAssayRegime.LCWGS,
            _input(CallerInputRole.COVERAGE_PROFILE),
            parents=("lane-sniffles",),
        ).model_copy(update={"lane_id": "lane-spectre-sniffles-supported"})

        with self.assertRaisesRegex(ValueError, "cycle"):
            seal_multicaller_plan(
                plan_id="synthetic-cycle",
                requests=[sniffles, spectre],
                runtime_availability=self._runtime_map(
                    CallerMode.SNIFFLES2_STANDARD,
                    CallerMode.SPECTRE_SNIFFLES_SUPPORTED,
                ),
            )

    def test_spectre_depth_only_rejects_hidden_sniffles_parent(self) -> None:
        sniffles = _request(
            CallerMode.SNIFFLES2_STANDARD,
            CallerAssayRegime.ADAPTIVE_SAMPLING,
            _input(CallerInputRole.ALIGNED_BAM, "sniffles"),
        ).model_copy(update={"lane_id": "lane-sniffles"})
        depth_only = _request(
            CallerMode.SPECTRE_DEPTH_ONLY,
            CallerAssayRegime.LCWGS,
            _input(CallerInputRole.COVERAGE_PROFILE),
            parents=("lane-sniffles",),
        )

        with self.assertRaisesRegex(ValueError, "does not accept parent"):
            seal_multicaller_plan(
                plan_id="synthetic-hidden-dependency",
                requests=[sniffles, depth_only],
                runtime_availability=self._runtime_map(
                    CallerMode.SNIFFLES2_STANDARD,
                    CallerMode.SPECTRE_DEPTH_ONLY,
                ),
            )

    def test_parent_content_rekeys_dependent_but_not_unrelated_lane(self) -> None:
        def build(sniffles_suffix: str):
            sniffles = _request(
                CallerMode.SNIFFLES2_STANDARD,
                CallerAssayRegime.ADAPTIVE_SAMPLING,
                _input(CallerInputRole.ALIGNED_BAM, sniffles_suffix),
            ).model_copy(update={"lane_id": "lane-sniffles"})
            spectre = _request(
                CallerMode.SPECTRE_SNIFFLES_SUPPORTED,
                CallerAssayRegime.LCWGS,
                _input(CallerInputRole.COVERAGE_PROFILE),
                parents=("lane-sniffles",),
            ).model_copy(update={"lane_id": "lane-spectre-sniffles-supported"})
            qdnaseq = _request(
                CallerMode.QDNASEQ_ACE_MULTIBIN,
                CallerAssayRegime.LCWGS,
                _input(CallerInputRole.ALIGNED_BAM, "qdnaseq"),
            )
            return seal_multicaller_plan(
                plan_id="synthetic-parent-lock",
                requests=[sniffles, spectre, qdnaseq],
                runtime_availability=self._runtime_map(
                    CallerMode.SNIFFLES2_STANDARD,
                    CallerMode.SPECTRE_SNIFFLES_SUPPORTED,
                    CallerMode.QDNASEQ_ACE_MULTIBIN,
                ),
            )

        first = build("sniffles-v1")
        second = build("sniffles-v2")
        first_lanes = {item.lane_id: item for item in first.lanes}
        second_lanes = {item.lane_id: item for item in second.lanes}

        self.assertNotEqual(first.plan_sha256, second.plan_sha256)
        self.assertNotEqual(
            first_lanes["lane-sniffles"].lane_sha256,
            second_lanes["lane-sniffles"].lane_sha256,
        )
        self.assertNotEqual(
            first_lanes["lane-spectre-sniffles-supported"].lane_sha256,
            second_lanes["lane-spectre-sniffles-supported"].lane_sha256,
        )
        self.assertEqual(
            first_lanes["lane-qdnaseq_ace-multibin"].lane_sha256,
            second_lanes["lane-qdnaseq_ace-multibin"].lane_sha256,
        )

    def test_plan_lock_is_deterministic_across_request_order(self) -> None:
        first_request = _request(
            CallerMode.QDNASEQ_ACE_MULTIBIN,
            CallerAssayRegime.LCWGS,
            _input(CallerInputRole.ALIGNED_BAM, "one"),
        )
        second_request = _request(
            CallerMode.CUTESV_STANDARD,
            CallerAssayRegime.ADAPTIVE_SAMPLING,
            _input(CallerInputRole.ALIGNED_BAM, "two"),
        )

        first = seal_multicaller_plan(
            plan_id="synthetic-deterministic",
            requests=[first_request, second_request],
            runtime_availability=self._runtime_map(
                CallerMode.QDNASEQ_ACE_MULTIBIN,
                CallerMode.CUTESV_STANDARD,
            ),
        )
        second = seal_multicaller_plan(
            plan_id="synthetic-deterministic",
            requests=[second_request, first_request],
            runtime_availability=self._runtime_map(
                CallerMode.QDNASEQ_ACE_MULTIBIN,
                CallerMode.CUTESV_STANDARD,
            ),
        )

        self.assertEqual(first.plan_sha256, second.plan_sha256)
        self.assertEqual(first.topological_order, second.topological_order)

    def test_duplicate_lane_ids_are_rejected(self) -> None:
        first = _request(
            CallerMode.QDNASEQ_ACE_MULTIBIN,
            CallerAssayRegime.LCWGS,
            _input(CallerInputRole.ALIGNED_BAM, "one"),
        ).model_copy(update={"lane_id": "duplicate-lane"})
        second = _request(
            CallerMode.CUTESV_STANDARD,
            CallerAssayRegime.ADAPTIVE_SAMPLING,
            _input(CallerInputRole.ALIGNED_BAM, "two"),
        ).model_copy(update={"lane_id": "duplicate-lane"})

        with self.assertRaisesRegex(ValueError, "lane IDs"):
            seal_multicaller_plan(
                plan_id="synthetic-duplicates",
                requests=[first, second],
                runtime_availability=self._runtime_map(
                    CallerMode.QDNASEQ_ACE_MULTIBIN,
                    CallerMode.CUTESV_STANDARD,
                ),
            )


if __name__ == "__main__":
    unittest.main()

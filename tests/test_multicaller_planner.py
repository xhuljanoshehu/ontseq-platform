from __future__ import annotations

import hashlib
import unittest

from ontseq_platform.multicaller_planner import plan_caller_lane

from ontseq_platform.multicaller_contracts import (
    CallerAssayRegime,
    CallerInputArtifact,
    CallerInputRole,
    CallerLaneRequest,
    CallerMode,
    CallerPlanningDecision,
)


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
) -> CallerLaneRequest:
    return CallerLaneRequest(
        lane_id=f"lane-{mode.value.replace(':', '-')}",
        mode_id=mode,
        assay_regime=regime,
        inputs=list(inputs),
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


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import hashlib
import unittest

from ontseq_platform.multicaller_report import (
    CallerDependencyStatus,
    CallerExecutionOutcome,
    CallerLaneExecutionRecord,
    ComparisonGroupState,
    MultiCallerComparisonGroup,
    MultiCallerComparisonMember,
    MultiCallerComparisonReport,
    MultiCallerEvidenceLink,
    build_multicaller_comparison_report,
)
from pydantic import ValidationError

from ontseq_platform.multicaller_contracts import (
    CallerAssayRegime,
    CallerInputArtifact,
    CallerInputRole,
    CallerMode,
    CallerPlanningDecision,
)
from ontseq_platform.multicaller_identity import (
    CallerEvidenceSourceKind,
    build_multicaller_collection_address,
)
from ontseq_platform.multicaller_planner import seal_multicaller_plan


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _input(role: CallerInputRole, label: str) -> CallerInputArtifact:
    return CallerInputArtifact(
        artifact_id=label,
        role=role,
        sha256=_sha(label),
    )


def _dependent_plan():
    from ontseq_platform.multicaller_contracts import CallerLaneRequest

    severus = CallerLaneRequest(
        lane_id="lane-severus",
        mode_id=CallerMode.SEVERUS_SINGLE_SAMPLE,
        assay_regime=CallerAssayRegime.TUMOR_ONLY_LONG_READ,
        inputs=[_input(CallerInputRole.TUMOR_BAM, "tumor-severus")],
    )
    wakhan = CallerLaneRequest(
        lane_id="lane-wakhan",
        mode_id=CallerMode.WAKHAN_PHASED_CNA,
        assay_regime=CallerAssayRegime.TUMOR_ONLY_LONG_READ,
        inputs=[
            _input(CallerInputRole.TUMOR_BAM, "tumor-wakhan"),
            _input(CallerInputRole.PHASED_VARIANTS, "phased-vcf"),
        ],
        parent_lane_ids=["lane-severus"],
    )
    return seal_multicaller_plan(
        plan_id="plan-dependent",
        requests=[severus, wakhan],
        runtime_availability={
            CallerMode.SEVERUS_SINGLE_SAMPLE: True,
            CallerMode.WAKHAN_PHASED_CNA: True,
        },
    )


def _lane(plan, lane_id: str):
    return next(item for item in plan.lanes if item.lane_id == lane_id)


class MultiCallerComparisonReportTests(unittest.TestCase):
    def test_direct_dependency_is_visible_and_not_called_independent_confirmation(self) -> None:
        plan = _dependent_plan()
        severus_lane = _lane(plan, "lane-severus")
        wakhan_lane = _lane(plan, "lane-wakhan")

        severus_address = build_multicaller_collection_address(
            lane_sha256=severus_lane.lane_sha256,
            source_kind=CallerEvidenceSourceKind.NATIVE_ARTIFACT,
            source_id="severus-native-vcf",
        )
        wakhan_address = build_multicaller_collection_address(
            lane_sha256=wakhan_lane.lane_sha256,
            source_kind=CallerEvidenceSourceKind.NATIVE_ARTIFACT,
            source_id="wakhan-native-cna",
        )
        report = build_multicaller_comparison_report(
            plan=plan,
            executions=[
                CallerLaneExecutionRecord(
                    lane_id="lane-severus",
                    outcome=CallerExecutionOutcome.COMPLETED,
                    evidence_links=[
                        MultiCallerEvidenceLink(
                            address=severus_address,
                            label="Severus native somatic SV evidence",
                            native_artifact_ids=["severus-native-vcf"],
                        )
                    ],
                ),
                CallerLaneExecutionRecord(
                    lane_id="lane-wakhan",
                    outcome=CallerExecutionOutcome.COMPLETED,
                    evidence_links=[
                        MultiCallerEvidenceLink(
                            address=wakhan_address,
                            label="Wakhan native CNA evidence",
                            native_artifact_ids=["wakhan-native-cna"],
                        )
                    ],
                ),
            ],
        )

        lanes = {item.lane_id: item for item in report.lanes}
        self.assertEqual(
            lanes["lane-severus"].dependency_status,
            CallerDependencyStatus.NO_DECLARED_CALLER_DEPENDENCY,
        )
        self.assertEqual(
            lanes["lane-wakhan"].dependency_status,
            CallerDependencyStatus.DIRECT_CALLER_DEPENDENCY,
        )
        self.assertEqual(lanes["lane-wakhan"].parent_lane_ids, ["lane-severus"])
        self.assertEqual(
            lanes["lane-wakhan"].parent_lane_sha256s,
            {"lane-severus": severus_lane.lane_sha256},
        )
        self.assertNotIn("independent_confirmation", report.model_dump_json())

    def test_failed_unavailable_ineligible_and_no_call_lanes_remain_visible(self) -> None:
        from ontseq_platform.multicaller_contracts import CallerLaneRequest

        failed = CallerLaneRequest(
            lane_id="lane-failed",
            mode_id=CallerMode.QDNASEQ_ACE_MULTIBIN,
            assay_regime=CallerAssayRegime.LCWGS,
            inputs=[_input(CallerInputRole.ALIGNED_BAM, "failed-bam")],
        )
        no_call = CallerLaneRequest(
            lane_id="lane-nocall",
            mode_id=CallerMode.QDNASEQ_ACE_MULTIBIN,
            assay_regime=CallerAssayRegime.LCWGS,
            inputs=[_input(CallerInputRole.ALIGNED_BAM, "nocall-bam")],
        )
        unavailable = CallerLaneRequest(
            lane_id="lane-unavailable",
            mode_id=CallerMode.ICHORCNA_ULP_WGS,
            assay_regime=CallerAssayRegime.CFDNA_ULP_WGS,
            inputs=[_input(CallerInputRole.COVERAGE_PROFILE, "ichor-coverage")],
        )
        ineligible = CallerLaneRequest(
            lane_id="lane-ineligible",
            mode_id=CallerMode.SAVANA_PAIRED,
            assay_regime=CallerAssayRegime.TUMOR_NORMAL_LONG_READ,
            inputs=[
                _input(CallerInputRole.TUMOR_BAM, "savana-tumor"),
                _input(CallerInputRole.SNP_VCF, "savana-snps"),
            ],
        )
        plan = seal_multicaller_plan(
            plan_id="plan-visible-states",
            requests=[failed, no_call, unavailable, ineligible],
            runtime_availability={
                CallerMode.QDNASEQ_ACE_MULTIBIN: True,
                CallerMode.ICHORCNA_ULP_WGS: False,
                CallerMode.SAVANA_PAIRED: True,
            },
        )

        report = build_multicaller_comparison_report(
            plan=plan,
            executions=[
                CallerLaneExecutionRecord(
                    lane_id="lane-failed",
                    outcome=CallerExecutionOutcome.FAILED,
                    reasons=["Synthetic caller failure."],
                ),
                CallerLaneExecutionRecord(
                    lane_id="lane-nocall",
                    outcome=CallerExecutionOutcome.NO_CALL,
                    reasons=["No event passed the registered technical policy."],
                ),
                CallerLaneExecutionRecord(
                    lane_id="lane-unavailable",
                    outcome=CallerExecutionOutcome.NOT_RUN,
                    reasons=["Pinned runtime is unavailable."],
                ),
                CallerLaneExecutionRecord(
                    lane_id="lane-ineligible",
                    outcome=CallerExecutionOutcome.NOT_RUN,
                    reasons=["Matched normal input is missing."],
                ),
            ],
        )

        lanes = {item.lane_id: item for item in report.lanes}
        self.assertEqual(set(lanes), set(plan.topological_order))
        self.assertEqual(lanes["lane-failed"].execution_outcome, CallerExecutionOutcome.FAILED)
        self.assertEqual(lanes["lane-nocall"].execution_outcome, CallerExecutionOutcome.NO_CALL)
        self.assertEqual(
            lanes["lane-unavailable"].planning_decision,
            CallerPlanningDecision.RUNTIME_UNAVAILABLE,
        )
        self.assertEqual(
            lanes["lane-ineligible"].planning_decision,
            CallerPlanningDecision.INELIGIBLE,
        )
        self.assertEqual(
            lanes["lane-unavailable"].execution_outcome,
            CallerExecutionOutcome.NOT_RUN,
        )
        self.assertEqual(
            lanes["lane-ineligible"].execution_outcome,
            CallerExecutionOutcome.NOT_RUN,
        )

    def test_native_evidence_links_remain_reachable(self) -> None:
        plan = _dependent_plan()
        lane = _lane(plan, "lane-severus")
        address = build_multicaller_collection_address(
            lane_sha256=lane.lane_sha256,
            source_kind=CallerEvidenceSourceKind.NATIVE_ARTIFACT,
            source_id="severus-native-vcf",
        )
        report = build_multicaller_comparison_report(
            plan=plan,
            executions=[
                CallerLaneExecutionRecord(
                    lane_id="lane-severus",
                    outcome=CallerExecutionOutcome.COMPLETED,
                    evidence_links=[
                        MultiCallerEvidenceLink(
                            address=address,
                            label="Native VCF",
                            native_artifact_ids=["severus-native-vcf"],
                        )
                    ],
                ),
                CallerLaneExecutionRecord(
                    lane_id="lane-wakhan",
                    outcome=CallerExecutionOutcome.NOT_ASSESSABLE,
                    reasons=["Phased CNA output unavailable for this synthetic comparison."],
                ),
            ],
        )

        link = next(
            item
            for item in report.lanes
            if item.lane_id == "lane-severus"
        ).evidence_links[0]
        self.assertEqual(link.address.collection_id, address.collection_id)
        self.assertEqual(link.address.source_id, "severus-native-vcf")
        self.assertEqual(link.native_artifact_ids, ["severus-native-vcf"])

    def test_conflicting_numeric_observations_are_preserved_without_aggregation(self) -> None:
        plan = _dependent_plan()
        severus_lane = _lane(plan, "lane-severus")
        wakhan_lane = _lane(plan, "lane-wakhan")
        severus_address = build_multicaller_collection_address(
            lane_sha256=severus_lane.lane_sha256,
            source_kind=CallerEvidenceSourceKind.FULL_EVIDENCE,
            source_id="severus-copy-observation",
        )
        wakhan_address = build_multicaller_collection_address(
            lane_sha256=wakhan_lane.lane_sha256,
            source_kind=CallerEvidenceSourceKind.FULL_EVIDENCE,
            source_id="wakhan-copy-observation",
        )
        group = MultiCallerComparisonGroup(
            group_id="group-chr7-copy-number",
            subject="chr7:1000000-5000000 copy number",
            state=ComparisonGroupState.CONFLICTING,
            members=[
                MultiCallerComparisonMember(
                    lane_id="lane-severus",
                    evidence_collection_id=severus_address.collection_id,
                    observation="copy_number",
                    numeric_value=1.0,
                    unit="copies",
                ),
                MultiCallerComparisonMember(
                    lane_id="lane-wakhan",
                    evidence_collection_id=wakhan_address.collection_id,
                    observation="copy_number",
                    numeric_value=3.0,
                    unit="copies",
                ),
            ],
        )
        report = build_multicaller_comparison_report(
            plan=plan,
            executions=[
                CallerLaneExecutionRecord(
                    lane_id="lane-severus",
                    outcome=CallerExecutionOutcome.COMPLETED,
                    evidence_links=[
                        MultiCallerEvidenceLink(
                            address=severus_address,
                            label="Severus copy-number observation",
                        )
                    ],
                ),
                CallerLaneExecutionRecord(
                    lane_id="lane-wakhan",
                    outcome=CallerExecutionOutcome.COMPLETED,
                    evidence_links=[
                        MultiCallerEvidenceLink(
                            address=wakhan_address,
                            label="Wakhan copy-number observation",
                        )
                    ],
                ),
            ],
            groups=[group],
        )

        values = [member.numeric_value for member in report.groups[0].members]
        self.assertEqual(values, [1.0, 3.0])
        payload = report.model_dump()
        self.assertNotIn("aggregate_value", payload["groups"][0])
        self.assertNotIn("average", payload["groups"][0])
        self.assertNotIn("winner", payload)
        self.assertNotIn("majority_vote", payload)

    def test_report_contract_rejects_winner_or_majority_fields(self) -> None:
        plan = _dependent_plan()
        report = build_multicaller_comparison_report(
            plan=plan,
            executions=[
                CallerLaneExecutionRecord(
                    lane_id="lane-severus",
                    outcome=CallerExecutionOutcome.NO_CALL,
                    reasons=["Synthetic no-call."],
                ),
                CallerLaneExecutionRecord(
                    lane_id="lane-wakhan",
                    outcome=CallerExecutionOutcome.NOT_ASSESSABLE,
                    reasons=["Synthetic not-assessable."],
                ),
            ],
        )
        payload = report.model_dump(mode="json")

        with self.assertRaises(ValidationError):
            MultiCallerComparisonReport.model_validate(
                {**payload, "winner": "lane-severus"}
            )
        with self.assertRaises(ValidationError):
            MultiCallerComparisonReport.model_validate(
                {**payload, "majority_vote": "supported"}
            )

    def test_evidence_link_must_match_the_lane_hash(self) -> None:
        plan = _dependent_plan()
        wrong_address = build_multicaller_collection_address(
            lane_sha256=_sha("wrong-lane"),
            source_kind=CallerEvidenceSourceKind.NATIVE_ARTIFACT,
            source_id="wrong-native-evidence",
        )

        with self.assertRaisesRegex(ValueError, "lane SHA"):
            build_multicaller_comparison_report(
                plan=plan,
                executions=[
                    CallerLaneExecutionRecord(
                        lane_id="lane-severus",
                        outcome=CallerExecutionOutcome.COMPLETED,
                        evidence_links=[
                            MultiCallerEvidenceLink(
                                address=wrong_address,
                                label="Mismatched native evidence",
                            )
                        ],
                    ),
                    CallerLaneExecutionRecord(
                        lane_id="lane-wakhan",
                        outcome=CallerExecutionOutcome.NOT_ASSESSABLE,
                        reasons=["Synthetic not-assessable."],
                    ),
                ],
            )


if __name__ == "__main__":
    unittest.main()

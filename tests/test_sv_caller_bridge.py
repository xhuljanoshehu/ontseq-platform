from __future__ import annotations

import unittest

from ontseq_platform.cutesv import CuteSVCallReport, CuteSVPolicy
from ontseq_platform.models import (
    EventType,
    Evidence,
    FileFingerprint,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    SnifflesCallReport,
    SnifflesPolicy,
    ToolRecord,
)
from ontseq_platform.sv_caller_bridge import compare_sniffles_and_cutesv
from ontseq_platform.sv_concordance import SVConcordancePolicy, SVConcordanceStatus


def _event(event_id: str, caller: str, version: str) -> GenomicEvent:
    return GenomicEvent(
        event_id=event_id,
        event_type=EventType.TRANSLOCATION,
        primary=Locus(chromosome="chr1", start=999, end=1000),
        secondary=Locus(chromosome="chr2", start=4999, end=5000),
        evidence=[
            Evidence(
                caller=caller,
                caller_version=version,
                support_reads=8,
                filters=["PASS"],
            )
        ],
        confidence="unclassified",
        reportable=False,
    )


def _sniffles_report(*, sample_id: str = "SYNTHETIC_001") -> SnifflesCallReport:
    event = _event("SNIFFLES2-000001", "Sniffles2", "2.8.0")
    return SnifflesCallReport(
        sample_id=sample_id,
        genome_build=GenomeBuild.GRCH38,
        status=ModuleRunStatus.COMPLETED,
        policy=SnifflesPolicy(
            profile_id="synthetic",
            status="technical_defaults_only",
            min_support=5,
            min_sv_length=50,
            mapq=20,
            note="Synthetic software test only.",
        ),
        events=[event],
        raw_record_count=1,
        accepted_record_count=1,
        rejected_record_count=0,
        rejection_counts={},
        tool=ToolRecord(name="Sniffles2", version="2.8.0"),
        vcf_fingerprint=FileFingerprint(size_bytes=1),
    )


def _cutesv_report(
    *,
    sample_id: str = "SYNTHETIC_001",
    no_call: bool = False,
) -> CuteSVCallReport:
    events = [] if no_call else [_event("CUTESV-000001", "cuteSV", "2.1.4")]
    return CuteSVCallReport(
        sample_id=sample_id,
        genome_build=GenomeBuild.GRCH38,
        status=ModuleRunStatus.NO_CALL if no_call else ModuleRunStatus.COMPLETED,
        policy=CuteSVPolicy(min_support=5, min_sv_length=50),
        events=events,
        raw_record_count=0 if no_call else 1,
        accepted_record_count=0 if no_call else 1,
        rejected_record_count=0,
        rejection_counts={},
        tool=ToolRecord(name="cuteSV", version="2.1.4"),
        vcf_fingerprint=FileFingerprint(size_bytes=1),
    )


def _policy() -> SVConcordancePolicy:
    return SVConcordancePolicy(
        maximum_breakpoint_distance_bp=50,
        note="Synthetic software comparison only.",
    )


class SVCallerBridgeTests(unittest.TestCase):
    def test_exact_multicaller_agreement_remains_support_only(self) -> None:
        report = compare_sniffles_and_cutesv(
            _sniffles_report(),
            _cutesv_report(),
            _policy(),
        )

        self.assertEqual(len(report.pairs), 1)
        pair = report.pairs[0]
        self.assertEqual(pair.status, SVConcordanceStatus.EXACT_MATCH)
        self.assertEqual(pair.evidence_semantics, "support_only_not_truth")
        self.assertTrue(pair.research_only)
        self.assertFalse(pair.reportable)
        self.assertEqual(report.conclusion_semantics, "caller_concordance_is_not_truth")

    def test_cutesv_no_call_does_not_become_negative_evidence(self) -> None:
        report = compare_sniffles_and_cutesv(
            _sniffles_report(),
            _cutesv_report(no_call=True),
            _policy(),
        )

        self.assertEqual(report.pairs, [])
        self.assertEqual(report.unmatched_left_observation_ids, ["SNIFFLES2-000001"])
        self.assertEqual(report.unmatched_right_observation_ids, [])
        self.assertTrue(any("NO_CALL" in warning for warning in report.warnings))
        self.assertTrue(any("not a validated negative" in warning for warning in report.warnings))

    def test_sample_mismatch_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "same sample"):
            compare_sniffles_and_cutesv(
                _sniffles_report(),
                _cutesv_report(sample_id="SYNTHETIC_002"),
                _policy(),
            )


if __name__ == "__main__":
    unittest.main()

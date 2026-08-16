from __future__ import annotations

import unittest

from ontseq_platform.cnv.evaluate import build_case_mask, compare_methods, evaluate_case
from ontseq_platform.cnv.models import (
    CnvBenchmarkCase,
    CnvCallSet,
    CnvDataBasis,
    CnvEvaluationOptions,
    CnvSegment,
    CnvStrata,
    CnvTruthSet,
    CnvTruthSource,
    GenomicRegion,
    segment_to_genomic_event,
)
from ontseq_platform.cnv.states import CopyNumberState
from ontseq_platform.models import EventType, GenomeBuild, ModuleRunStatus

CONTIGS = {"chr5": 181_538_259, "chr8": 145_138_636}


def _truth(segments, *, background=CopyNumberState.NEUTRAL, resolution=100_000, **kwargs):
    return CnvTruthSet(
        truth_id="TRUTH_001",
        sample_id="SYNTHETIC_AML_001",
        genome_build=GenomeBuild.GRCH38,
        source=CnvTruthSource.SIMULATED,
        source_version="v1",
        background_state=background,
        resolution_bp=resolution,
        segments=segments,
        **kwargs,
    )


def _call_set(segments, *, call_set_id="CALL_001", method="baseline", **kwargs):
    return CnvCallSet(
        call_set_id=call_set_id,
        sample_id="SYNTHETIC_AML_001",
        genome_build=GenomeBuild.GRCH38,
        method=method,
        method_version="0.1.0",
        data_basis=CnvDataBasis.SIMULATED,
        background_state=kwargs.pop("background_state", CopyNumberState.NEUTRAL),
        status=ModuleRunStatus.COMPLETED if segments else ModuleRunStatus.NO_CALL,
        segments=segments,
        **kwargs,
    )


def _case(truth, call_set, **kwargs):
    return CnvBenchmarkCase(
        case_id="CASE_001",
        genome_build=GenomeBuild.GRCH38,
        contig_lengths=CONTIGS,
        truth=truth,
        call_set=call_set,
        **kwargs,
    )


LOSS = CnvSegment(
    contig="chr5", start=70_000_000, end=160_000_000, state=CopyNumberState.LOSS, copy_number=1.0
)
GAIN = CnvSegment(
    contig="chr8", start=0, end=145_138_636, state=CopyNumberState.GAIN, copy_number=3.0
)


class ContractTests(unittest.TestCase):
    def test_round_trip_produces_a_serializable_report(self) -> None:
        report = evaluate_case(_case(_truth([LOSS]), _call_set([LOSS])))
        payload = report.model_dump_json()
        self.assertIn("detection_rate", payload)
        self.assertIs(report.research_only, True)
        self.assertTrue(report.limitations)
        self.assertEqual(report.detection_rate.point, 1.0)

    def test_case_rejects_mismatched_genome_builds(self) -> None:
        truth = _truth([LOSS])
        call_set = _call_set([LOSS]).model_copy(update={"genome_build": GenomeBuild.GRCH37})
        with self.assertRaises(ValueError):
            _case(truth, call_set)

    def test_case_rejects_mismatched_samples(self) -> None:
        call_set = _call_set([LOSS]).model_copy(update={"sample_id": "OTHER_SAMPLE"})
        with self.assertRaises(ValueError):
            _case(_truth([LOSS]), call_set)

    def test_completed_call_set_must_contain_segments(self) -> None:
        with self.assertRaises(ValueError):
            CnvCallSet(
                call_set_id="BAD",
                sample_id="SYNTHETIC_AML_001",
                genome_build=GenomeBuild.GRCH38,
                method="x",
                method_version="1",
                data_basis=CnvDataBasis.SIMULATED,
                background_state=CopyNumberState.NEUTRAL,
                status=ModuleRunStatus.COMPLETED,
                segments=[],
            )

    def test_no_call_set_must_not_contain_segments(self) -> None:
        with self.assertRaises(ValueError):
            CnvCallSet(
                call_set_id="BAD",
                sample_id="SYNTHETIC_AML_001",
                genome_build=GenomeBuild.GRCH38,
                method="x",
                method_version="1",
                data_basis=CnvDataBasis.SIMULATED,
                background_state=CopyNumberState.NEUTRAL,
                status=ModuleRunStatus.NO_CALL,
                segments=[LOSS],
            )


class MaskIntegrationTests(unittest.TestCase):
    def test_analysis_scope_restricts_the_evaluable_genome(self) -> None:
        case = _case(
            _truth([LOSS]),
            _call_set([LOSS]),
            analysis_scope=[GenomicRegion(contig="chr5", start=0, end=181_538_259)],
        )
        mask = build_case_mask(case)
        self.assertEqual(mask.evaluable_bases, 181_538_259)
        self.assertEqual(mask.excluded_bases_by_reason["outside_analysis_scope"], 145_138_636)

    def test_caller_no_call_regions_become_unassessable(self) -> None:
        case = _case(
            _truth([LOSS]),
            _call_set(
                [],
                no_call_regions=[
                    GenomicRegion(contig="chr5", start=0, end=181_538_259),
                    GenomicRegion(contig="chr8", start=0, end=145_138_636),
                ],
            ),
        )
        report = evaluate_case(case)
        self.assertEqual(report.truth_events[0].outcome, "NOT_ASSESSABLE")
        self.assertIsNone(report.detection_rate.point)
        self.assertGreater(report.partition.excluded_bases_by_reason["caller_no_call"], 0)

    def test_truth_uninformative_regions_are_excluded(self) -> None:
        case = _case(
            _truth(
                [LOSS],
                uninformative_regions=[GenomicRegion(contig="chr5", start=0, end=181_538_259)],
            ),
            _call_set([LOSS]),
        )
        report = evaluate_case(case)
        self.assertIn("truth_not_informative", report.partition.excluded_bases_by_reason)

    def test_excluded_regions_are_attributed(self) -> None:
        case = _case(
            _truth([LOSS]),
            _call_set([LOSS]),
            excluded_regions=[GenomicRegion(contig="chr5", start=0, end=50_000_000)],
        )
        report = evaluate_case(case)
        self.assertEqual(report.partition.excluded_bases_by_reason["blacklist"], 50_000_000)


class WarningTests(unittest.TestCase):
    def test_calls_below_truth_resolution_are_flagged_not_counted_as_errors(self) -> None:
        small = CnvSegment(
            contig="chr5", start=0, end=50_000, state=CopyNumberState.GAIN, copy_number=3.0
        )
        case = _case(
            _truth([LOSS], resolution=10_000_000),
            _call_set([small, LOSS]),
        )
        report = evaluate_case(case)
        self.assertTrue(
            any("smaller than the truth set's declared resolution" in w for w in report.warnings)
        )


class MethodComparisonTests(unittest.TestCase):
    def test_shared_mask_prevents_rewarding_a_caller_for_its_own_blind_spots(self) -> None:
        """A method that declines to call must not shrink only its own denominator."""
        confident = _call_set([LOSS, GAIN], call_set_id="CONFIDENT", method="confident")
        cautious = _call_set(
            [LOSS],
            call_set_id="CAUTIOUS",
            method="cautious",
            no_call_regions=[GenomicRegion(contig="chr8", start=0, end=145_138_636)],
        )
        template = _case(_truth([LOSS, GAIN]), confident)
        reports = compare_methods(template, [confident, cautious])

        self.assertEqual(len(reports), 2)
        # Both methods are scored on the same evaluable genome.
        self.assertEqual(reports[0].partition.evaluable_bases, reports[1].partition.evaluable_bases)
        # chr8 is unassessable for both, so neither is credited nor blamed there.
        for report in reports:
            self.assertEqual(report.detection_rate.total, 1)


class ProjectionTests(unittest.TestCase):
    def test_whole_chromosome_segment_becomes_a_chromosome_event(self) -> None:
        event = segment_to_genomic_event(GAIN, event_id="cnv-1", contig_length=145_138_636)
        self.assertIsNotNone(event)
        self.assertEqual(event.event_type, EventType.CHROMOSOME_GAIN)
        self.assertIs(event.reportable, False)

    def test_subchromosomal_segment_becomes_a_deletion(self) -> None:
        event = segment_to_genomic_event(LOSS, event_id="cnv-2", contig_length=181_538_259)
        self.assertEqual(event.event_type, EventType.DELETION)

    def test_copy_neutral_loh_is_not_projected_onto_a_dosage_event(self) -> None:
        loh = CnvSegment(
            contig="chr5",
            start=0,
            end=1_000_000,
            state=CopyNumberState.COPY_NEUTRAL_LOH,
            copy_number=2.0,
        )
        self.assertIsNone(segment_to_genomic_event(loh, event_id="cnv-3"))


class OptionEchoTests(unittest.TestCase):
    def test_options_are_echoed_into_the_report(self) -> None:
        options = CnvEvaluationOptions(detection_overlap_fraction=0.9, copy_number_tolerance=0.1)
        report = evaluate_case(
            _case(
                _truth([LOSS]),
                _call_set([LOSS]),
                options=options,
                strata=CnvStrata(tumor_fraction=0.25, mean_coverage_x=3.0),
            )
        )
        self.assertEqual(report.options.detection_overlap_fraction, 0.9)
        self.assertEqual(report.options.copy_number_tolerance, 0.1)
        self.assertEqual(report.strata.tumor_fraction, 0.25)


if __name__ == "__main__":
    unittest.main()

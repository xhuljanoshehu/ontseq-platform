from __future__ import annotations

import unittest

from ontseq_platform.cnv.core import StateSegment, evaluate
from ontseq_platform.cnv.segment import (
    DepthBin,
    SegmentationParameters,
    call_segments,
    neutral_background_segments,
)
from ontseq_platform.cnv.simulate import (
    SimulationParameters,
    closed_world_truth,
    simulate_bins,
    simulate_dilution_series,
    truth_copy_number_at,
)
from ontseq_platform.cnv.states import CopyNumberState, expected_mixture_copy_number

CONTIGS = {"5": 181_538_259, "7": 159_345_973, "8": 145_138_636, "2": 242_193_529}
EVENTS = [
    ("5", 70_000_000, 160_000_000, 1.0),
    ("7", 65_000_000, 159_345_973, 1.0),
    ("8", 0, 145_138_636, 3.0),
]


def _truth() -> list[StateSegment]:
    return closed_world_truth(EVENTS, CONTIGS)


class MixtureModelTests(unittest.TestCase):
    def test_pure_tumor_returns_the_tumor_copy_number(self) -> None:
        self.assertEqual(expected_mixture_copy_number(1.0, tumor_fraction=1.0), 1.0)

    def test_pure_normal_returns_the_normal_copy_number(self) -> None:
        self.assertEqual(expected_mixture_copy_number(1.0, tumor_fraction=0.0), 2.0)

    def test_dilution_compresses_towards_baseline(self) -> None:
        self.assertAlmostEqual(expected_mixture_copy_number(1.0, tumor_fraction=0.2), 1.8)
        self.assertAlmostEqual(expected_mixture_copy_number(3.0, tumor_fraction=0.2), 2.2)

    def test_truth_copy_number_is_length_weighted_across_a_window(self) -> None:
        segments = [StateSegment("5", 0, 500, CopyNumberState.LOSS, 1.0)]
        value = truth_copy_number_at(segments, "5", 0, 1000, default_copy_number=2.0)
        self.assertAlmostEqual(value, 1.5)


class DeterminismTests(unittest.TestCase):
    def test_same_seed_reproduces_identical_counts(self) -> None:
        parameters = SimulationParameters(seed=42, tumor_fraction=0.5)
        first = simulate_bins(
            contig_lengths=CONTIGS, truth_segments=_truth(), parameters=parameters
        )
        second = simulate_bins(
            contig_lengths=CONTIGS, truth_segments=_truth(), parameters=parameters
        )
        self.assertEqual([b.count for b in first], [b.count for b in second])

    def test_different_seeds_diverge(self) -> None:
        first = simulate_bins(
            contig_lengths=CONTIGS,
            truth_segments=_truth(),
            parameters=SimulationParameters(seed=1),
        )
        second = simulate_bins(
            contig_lengths=CONTIGS,
            truth_segments=_truth(),
            parameters=SimulationParameters(seed=2),
        )
        self.assertNotEqual([b.count for b in first], [b.count for b in second])

    def test_dilution_series_is_reproducible_and_covers_every_level(self) -> None:
        fractions = [1.0, 0.5, 0.2]
        levels = simulate_dilution_series(
            contig_lengths=CONTIGS,
            truth_segments=_truth(),
            tumor_fractions=fractions,
            replicates=2,
            base_parameters=SimulationParameters(seed=11),
        )
        self.assertEqual(len(levels), 6)
        self.assertEqual(sorted({level.tumor_fraction for level in levels}), sorted(fractions))
        repeated = simulate_dilution_series(
            contig_lengths=CONTIGS,
            truth_segments=_truth(),
            tumor_fractions=fractions,
            replicates=2,
            base_parameters=SimulationParameters(seed=11),
        )
        self.assertEqual(
            [b.count for level in levels for b in level.bins],
            [b.count for level in repeated for b in level.bins],
        )

    def test_replicates_at_one_level_are_independent_draws(self) -> None:
        levels = simulate_dilution_series(
            contig_lengths=CONTIGS,
            truth_segments=_truth(),
            tumor_fractions=[0.5],
            replicates=2,
            base_parameters=SimulationParameters(seed=5),
        )
        self.assertNotEqual(
            [b.count for b in levels[0].bins], [b.count for b in levels[1].bins]
        )


class ClosedWorldTruthTests(unittest.TestCase):
    def test_truth_profile_is_gapless(self) -> None:
        truth = _truth()
        by_contig: dict[str, list[StateSegment]] = {}
        for segment in truth:
            by_contig.setdefault(segment.contig, []).append(segment)
        for contig, segments in by_contig.items():
            segments.sort(key=lambda item: item.start)
            self.assertEqual(segments[0].start, 0)
            self.assertEqual(segments[-1].end, CONTIGS[contig])
            for previous, current in zip(segments, segments[1:], strict=False):
                self.assertEqual(previous.end, current.start)


class EndToEndTests(unittest.TestCase):
    """Simulate, call and evaluate must close the loop without external tools."""

    def _run(self, tumor_fraction: float, seed: int = 7):
        parameters = SimulationParameters(
            bin_size_bp=1_000_000,
            mean_coverage_x=3.0,
            tumor_fraction=tumor_fraction,
            seed=seed,
        )
        bins = simulate_bins(
            contig_lengths=CONTIGS, truth_segments=_truth(), parameters=parameters
        )
        called = call_segments(
            [DepthBin(b.contig, b.start, b.end, b.count) for b in bins],
            SegmentationParameters(
                tumor_fraction=tumor_fraction if tumor_fraction < 1.0 else None
            ),
        )
        query = neutral_background_segments(called.segments, CONTIGS)
        return evaluate(
            truth_segments=_truth(),
            query_segments=query,
            evaluable={contig: [(0, length)] for contig, length in CONTIGS.items()},
            reference_bases=sum(CONTIGS.values()),
            truth_background=CopyNumberState.NEUTRAL,
            query_background=CopyNumberState.NEUTRAL,
        )

    def test_pure_tumor_recovers_every_event_with_high_base_concordance(self) -> None:
        report = self._run(1.0)
        self.assertEqual(report.detection_rate.successes, 3)
        self.assertEqual(report.detection_rate.total, 3)
        self.assertGreater(report.base_level.concordance, 0.99)
        self.assertEqual(
            sum(1 for item in report.query_events if item.outcome.value == "UNCONFIRMED"), 0
        )

    def test_moderate_dilution_still_recovers_every_event(self) -> None:
        report = self._run(0.5)
        self.assertEqual(report.detection_rate.successes, 3)
        self.assertGreater(report.base_level.concordance, 0.95)

    def test_base_level_concordance_degrades_before_event_detection_does(self) -> None:
        """The reason base-level metrics are reported alongside event detection.

        At 5% tumor fraction the naive baseline caller still 'detects' all three events,
        because it calls whole chromosomes. Event-level detection alone would report a
        perfect result; base-level concordance exposes the loss of specificity.
        """
        strong = self._run(1.0)
        weak = self._run(0.05)
        self.assertEqual(weak.detection_rate.successes, 3)
        self.assertLess(weak.base_level.concordance, strong.base_level.concordance)
        self.assertLess(weak.base_level.concordance, 0.8)


if __name__ == "__main__":
    unittest.main()

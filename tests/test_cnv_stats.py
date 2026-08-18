from __future__ import annotations

import unittest

from ontseq_platform.cnv.stats import (
    fit_logistic,
    inverse_normal_cdf,
    logistic_threshold,
    mcnemar_exact,
    mean_absolute_error,
    minimum_attainable_p_value,
    root_mean_square_error,
    wilson_interval,
)


class NormalQuantileTests(unittest.TestCase):
    def test_known_quantiles(self) -> None:
        self.assertAlmostEqual(inverse_normal_cdf(0.5), 0.0, places=9)
        self.assertAlmostEqual(inverse_normal_cdf(0.975), 1.959964, places=5)
        self.assertAlmostEqual(inverse_normal_cdf(0.995), 2.575829, places=5)

    def test_symmetry(self) -> None:
        self.assertAlmostEqual(inverse_normal_cdf(0.1), -inverse_normal_cdf(0.9), places=9)

    def test_out_of_range_is_rejected(self) -> None:
        for value in (0.0, 1.0, -0.1, 1.5):
            with self.assertRaises(ValueError):
                inverse_normal_cdf(value)


class WilsonIntervalTests(unittest.TestCase):
    def test_zero_denominator_is_undefined_not_zero(self) -> None:
        estimate = wilson_interval(0, 0)
        self.assertIsNone(estimate.point)
        self.assertIsNone(estimate.lower)
        self.assertIsNone(estimate.upper)

    def test_perfect_small_sample_is_not_reported_as_certain(self) -> None:
        estimate = wilson_interval(3, 3)
        self.assertEqual(estimate.point, 1.0)
        self.assertLess(estimate.lower, 0.5)
        self.assertLessEqual(estimate.upper, 1.0)

    def test_zero_successes_keeps_a_non_zero_upper_bound(self) -> None:
        estimate = wilson_interval(0, 5)
        self.assertEqual(estimate.point, 0.0)
        self.assertEqual(estimate.lower, 0.0)
        self.assertGreater(estimate.upper, 0.0)

    def test_interval_narrows_as_the_sample_grows(self) -> None:
        small = wilson_interval(8, 10)
        large = wilson_interval(800, 1000)
        self.assertLess(large.upper - large.lower, small.upper - small.lower)

    def test_invalid_counts_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            wilson_interval(5, 3)


class ErrorMetricTests(unittest.TestCase):
    def test_weighted_mean_absolute_error(self) -> None:
        value = mean_absolute_error([(1.0, 1.5), (2.0, 2.0)], [100, 300])
        self.assertAlmostEqual(value, 0.125)

    def test_rmse_penalises_large_errors_more(self) -> None:
        mae = mean_absolute_error([(0.0, 0.0), (0.0, 2.0)], [1, 1])
        rmse = root_mean_square_error([(0.0, 0.0), (0.0, 2.0)], [1, 1])
        self.assertGreater(rmse, mae)

    def test_zero_weight_is_undefined(self) -> None:
        self.assertIsNone(mean_absolute_error([], []))
        self.assertIsNone(root_mean_square_error([], []))


class McNemarTests(unittest.TestCase):
    def test_no_discordant_pairs_is_undefined_not_one(self) -> None:
        # Agreement on everything is not evidence of equivalence.
        self.assertIsNone(mcnemar_exact(0, 0))

    def test_symmetry(self) -> None:
        self.assertEqual(mcnemar_exact(3, 8), mcnemar_exact(8, 3))

    def test_known_exact_values(self) -> None:
        # n=1: two-sided p = 2 * 0.5 = 1.0
        self.assertAlmostEqual(mcnemar_exact(1, 0), 1.0)
        # n=2, smaller=0: 2 * (1 * 0.25) = 0.5
        self.assertAlmostEqual(mcnemar_exact(2, 0), 0.5)
        # n=6, smaller=0: 2 * (1/64) = 0.03125
        self.assertAlmostEqual(mcnemar_exact(6, 0), 0.03125)
        # n=10, smaller=1: 2 * (1 + 10) / 1024 = 0.021484375
        self.assertAlmostEqual(mcnemar_exact(9, 1), 0.021484375)

    def test_balanced_discordance_is_never_significant(self) -> None:
        self.assertAlmostEqual(mcnemar_exact(5, 5), 1.0)

    def test_p_value_never_exceeds_one(self) -> None:
        for first in range(6):
            for second in range(6):
                value = mcnemar_exact(first, second)
                if value is not None:
                    self.assertLessEqual(value, 1.0)
                    self.assertGreaterEqual(value, 0.0)

    def test_negative_counts_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            mcnemar_exact(-1, 2)


class LogisticFitTests(unittest.TestCase):
    def test_single_level_cannot_be_fitted(self) -> None:
        self.assertIsNone(fit_logistic([0.5], [3], [5]))

    def test_perfectly_separated_design_is_refused(self) -> None:
        # Every level is all-or-nothing, so the maximum likelihood slope is infinite.
        self.assertIsNone(fit_logistic([0.1, 0.9], [0, 5], [5, 5]))

    def test_graded_design_converges_with_a_positive_slope(self) -> None:
        fit = fit_logistic([0.05, 0.1, 0.2, 0.5], [1, 4, 8, 10], [10, 10, 10, 10])
        self.assertIsNotNone(fit)
        self.assertTrue(fit.converged)
        self.assertGreater(fit.slope, 0)

    def test_threshold_recovers_a_plausible_detection_limit(self) -> None:
        fit = fit_logistic([0.05, 0.1, 0.2, 0.5], [1, 4, 8, 10], [10, 10, 10, 10])
        limit = logistic_threshold(fit, 0.95)
        self.assertIsNotNone(limit)
        self.assertGreater(limit, 0.05)

    def test_threshold_is_withheld_for_a_negative_slope(self) -> None:
        fit = fit_logistic([0.05, 0.1, 0.2, 0.5], [10, 8, 4, 1], [10, 10, 10, 10])
        self.assertIsNotNone(fit)
        self.assertIsNone(logistic_threshold(fit, 0.95))


if __name__ == "__main__":
    unittest.main()


class MinimumAttainablePValueTests(unittest.TestCase):
    """A comparison that could never have reached significance is a design fault.

    Reported next to the p-value because the two readings are entirely different: "we
    compared them and found no difference" is a result, "this comparison could not have
    detected one" is not, and a bare p=0.125 looks identical to both.
    """

    def test_no_discordant_pair_has_no_attainable_value(self) -> None:
        self.assertIsNone(minimum_attainable_p_value(0))

    def test_four_pairs_cannot_reach_five_percent(self) -> None:
        self.assertAlmostEqual(minimum_attainable_p_value(4), 0.125)

    def test_six_pairs_is_the_smallest_count_that_can(self) -> None:
        self.assertGreater(minimum_attainable_p_value(5), 0.05)
        self.assertLess(minimum_attainable_p_value(6), 0.05)

    def test_the_floor_matches_the_most_extreme_actual_split(self) -> None:
        """The floor is a claim about the test, so it must equal what the test returns."""
        for pairs in range(1, 12):
            with self.subTest(pairs=pairs):
                self.assertAlmostEqual(minimum_attainable_p_value(pairs), mcnemar_exact(pairs, 0))

    def test_a_negative_count_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            minimum_attainable_p_value(-1)

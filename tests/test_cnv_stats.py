from __future__ import annotations

import unittest

from ontseq_platform.cnv.stats import (
    fit_logistic,
    inverse_normal_cdf,
    logistic_threshold,
    mean_absolute_error,
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

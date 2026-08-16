"""Small, dependency-free statistics helpers for CNV benchmarking.

Benchmark cohorts in this domain are small: a validation series may contain a handful of
specimens per dilution level. Reporting a bare point estimate such as "recall 0.8" from
five events invites over-reading, so every proportion this package reports is
accompanied by an interval and by the sample size the interval was computed from.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class ProportionEstimate:
    """A proportion with a Wilson score interval and the count it came from."""

    successes: int
    total: int
    point: float | None
    lower: float | None
    upper: float | None
    confidence_level: float


def inverse_normal_cdf(probability: float) -> float:
    """Return the standard normal quantile using Acklam's rational approximation.

    The absolute error is below 1.15e-9 across the open unit interval, which is far
    tighter than anything a benchmark interval needs, and it avoids a scipy dependency.
    """
    if not 0.0 < probability < 1.0:
        raise ValueError("probability must lie strictly between 0 and 1")
    a = (
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    )
    b = (
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    )
    c = (
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    )
    d = (
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e00,
        3.754408661907416e00,
    )
    lower_break = 0.02425
    upper_break = 1.0 - lower_break
    if probability < lower_break:
        q = math.sqrt(-2.0 * math.log(probability))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )
    if probability > upper_break:
        q = math.sqrt(-2.0 * math.log(1.0 - probability))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )
    q = probability - 0.5
    r = q * q
    return (
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5])
        * q
        / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    )


def wilson_interval(
    successes: int, total: int, *, confidence_level: float = 0.95
) -> ProportionEstimate:
    """Return a Wilson score interval for a binomial proportion.

    The Wilson interval is used instead of the normal approximation because benchmark
    proportions routinely sit at 0 or 1 with small denominators, where the normal
    interval collapses to zero width and asserts certainty that does not exist.

    A zero denominator yields ``None`` for every field rather than 0.0. An undefined
    proportion must stay undefined; substituting a number would let an unevaluated
    stratum look like a perfect or failed one.
    """
    if successes < 0 or total < 0 or successes > total:
        raise ValueError("successes must lie between 0 and total")
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence level must lie strictly between 0 and 1")
    if total == 0:
        return ProportionEstimate(
            successes=0,
            total=0,
            point=None,
            lower=None,
            upper=None,
            confidence_level=confidence_level,
        )
    z = inverse_normal_cdf(1.0 - (1.0 - confidence_level) / 2.0)
    proportion = successes / total
    denominator = 1.0 + (z * z) / total
    center = (proportion + (z * z) / (2 * total)) / denominator
    spread = (
        z
        / denominator
        * math.sqrt(proportion * (1.0 - proportion) / total + (z * z) / (4 * total * total))
    )
    return ProportionEstimate(
        successes=successes,
        total=total,
        point=proportion,
        lower=max(0.0, center - spread),
        upper=min(1.0, center + spread),
        confidence_level=confidence_level,
    )


def mean_absolute_error(
    pairs: Sequence[tuple[float, float]], weights: Sequence[int]
) -> float | None:
    """Return the weight-weighted mean absolute error, or ``None`` when unweighted."""
    total_weight = sum(weights)
    if total_weight <= 0:
        return None
    return sum(abs(a - b) * w for (a, b), w in zip(pairs, weights, strict=True)) / total_weight


def root_mean_square_error(
    pairs: Sequence[tuple[float, float]], weights: Sequence[int]
) -> float | None:
    """Return the weight-weighted RMSE, or ``None`` when unweighted."""
    total_weight = sum(weights)
    if total_weight <= 0:
        return None
    return math.sqrt(
        sum(((a - b) ** 2) * w for (a, b), w in zip(pairs, weights, strict=True)) / total_weight
    )


def mcnemar_exact(only_first: int, only_second: int) -> float | None:
    """Return the two-sided exact p-value of McNemar's test on discordant pairs.

    Two CNV methods scored on the same truth events produce paired binary outcomes, and
    a paired test is the only honest way to compare them: an unpaired comparison of two
    detection rates ignores that the same events drive both numbers and therefore
    overstates the uncertainty of the *difference*.

    Only discordant pairs carry information. Under the null hypothesis that the methods
    are equivalent, each discordant pair is a fair coin, so the exact test is a two-sided
    binomial test with ``p = 0.5`` on ``n = only_first + only_second`` trials.

    Returns ``None`` when there are no discordant pairs at all, because a test with no
    informative observations has no p-value. Reporting 1.0 there would suggest the
    methods were shown to be equivalent, when in fact nothing was measured.
    """
    if only_first < 0 or only_second < 0:
        raise ValueError("discordant counts must not be negative")
    total = only_first + only_second
    if total == 0:
        return None
    smaller = min(only_first, only_second)
    tail = sum(math.comb(total, k) for k in range(smaller + 1)) * (0.5**total)
    return min(1.0, 2.0 * tail)


@dataclass(frozen=True)
class LogisticFit:
    """A one-predictor logistic fit with an explicit convergence flag."""

    intercept: float
    slope: float
    converged: bool
    iterations: int


def fit_logistic(
    predictors: Sequence[float],
    successes: Sequence[int],
    totals: Sequence[int],
    *,
    maximum_iterations: int = 100,
    tolerance: float = 1e-9,
) -> LogisticFit | None:
    """Fit ``logit(p) = intercept + slope * x`` by Newton-Raphson on grouped counts.

    Returns ``None`` when the design cannot support a fit at all: fewer than two
    distinct predictor levels, no observations, or a perfectly separated design where
    the maximum-likelihood slope is infinite. Reporting a finite slope for a separated
    design would manufacture a limit of detection out of data that cannot support one.
    """
    if not (len(predictors) == len(successes) == len(totals)):
        raise ValueError("predictors, successes and totals must have equal length")
    levels = [(x, s, n) for x, s, n in zip(predictors, successes, totals, strict=True) if n > 0]
    if len({x for x, _, _ in levels}) < 2:
        return None
    observed = [s / n for _, s, n in levels]
    if all(p in (0.0, 1.0) for p in observed):
        # Complete separation: every level is all-or-nothing, so the likelihood has no
        # finite maximum.
        return None

    intercept = 0.0
    slope = 0.0
    for iteration in range(1, maximum_iterations + 1):
        gradient = [0.0, 0.0]
        hessian = [[0.0, 0.0], [0.0, 0.0]]
        for x, s, n in levels:
            eta = intercept + slope * x
            # Clamp to avoid overflow in exp for extreme intermediate steps.
            eta = max(-500.0, min(500.0, eta))
            p = 1.0 / (1.0 + math.exp(-eta))
            residual = s - n * p
            weight = n * p * (1.0 - p)
            gradient[0] += residual
            gradient[1] += residual * x
            hessian[0][0] += weight
            hessian[0][1] += weight * x
            hessian[1][0] += weight * x
            hessian[1][1] += weight * x * x
        determinant = hessian[0][0] * hessian[1][1] - hessian[0][1] * hessian[1][0]
        if abs(determinant) < 1e-12:
            return None
        step_intercept = (hessian[1][1] * gradient[0] - hessian[0][1] * gradient[1]) / determinant
        step_slope = (-hessian[1][0] * gradient[0] + hessian[0][0] * gradient[1]) / determinant
        intercept += step_intercept
        slope += step_slope
        if abs(step_intercept) < tolerance and abs(step_slope) < tolerance:
            return LogisticFit(
                intercept=intercept, slope=slope, converged=True, iterations=iteration
            )
    return LogisticFit(
        intercept=intercept, slope=slope, converged=False, iterations=maximum_iterations
    )


def logistic_threshold(fit: LogisticFit, target_probability: float) -> float | None:
    """Return the predictor value at which the fit reaches ``target_probability``.

    This is the model-based limit of detection, for example ``LoD95`` at
    ``target_probability = 0.95``. Returns ``None`` for a non-converged fit or a
    non-positive slope, because extrapolating a detection limit from a model that did
    not converge, or that predicts detection falling as signal rises, is meaningless.
    """
    if not fit.converged or fit.slope <= 0:
        return None
    if not 0.0 < target_probability < 1.0:
        raise ValueError("target probability must lie strictly between 0 and 1")
    return (math.log(target_probability / (1.0 - target_probability)) - fit.intercept) / fit.slope

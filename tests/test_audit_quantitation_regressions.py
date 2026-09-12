"""Numerical edge cases; these are statistical-model checks, not assay validation."""

from __future__ import annotations

import math

import pytest

from ontseq_platform.quantitation import (
    QuantitationError,
    _binomial_upper_tail,
    minimum_detectable_cancer_cell_fraction,
    minimum_detectable_vaf,
    minimum_variant_reads,
)


def test_large_symmetric_binomial_does_not_overflow() -> None:
    # Independent identity: P(X>=n/2) = 1/2 + P(X=n/2)/2 for p=1/2.
    expected = 0.5 + (math.comb(2000, 1000) / (2**2000)) / 2
    assert _binomial_upper_tail(1000, 2000, 0.5) == pytest.approx(expected, abs=1e-10)


def test_large_depth_rare_event_uses_exact_at_least_one_identity() -> None:
    expected = -math.expm1(10000 * math.log1p(-0.00001))
    assert _binomial_upper_tail(1, 10000, 0.00001) == pytest.approx(expected, abs=1e-10)


@pytest.mark.parametrize("power", [0.0, 1.0, -0.1, 1.1, math.nan, math.inf])
def test_invalid_requested_power_is_rejected(power: float) -> None:
    with pytest.raises(QuantitationError, match="power"):
        minimum_detectable_vaf(2, power=power)


def test_invalid_purity_is_not_hidden_by_an_unresolvable_depth() -> None:
    with pytest.raises(QuantitationError, match="tumour_fraction"):
        minimum_detectable_cancer_cell_fraction(1, tumour_fraction=-1, error_rate=0.9)


def test_high_depth_detection_floor_is_bounded_and_bracketed() -> None:
    needed = minimum_variant_reads(10000, error_rate=0.01, alpha=0.01)
    assert 100 < needed < 150
    assert _binomial_upper_tail(needed - 1, 10000, 0.01) > 0.01
    assert _binomial_upper_tail(needed, 10000, 0.01) <= 0.01
    vaf = minimum_detectable_vaf(10000, error_rate=0.01, alpha=0.01, power=0.95)
    assert vaf is not None and 0.01 < vaf < 0.1
    assert _binomial_upper_tail(needed, 10000, vaf) == pytest.approx(0.95, abs=1e-8)


@pytest.mark.parametrize("n", [1, 2, 5, 12, 40])
@pytest.mark.parametrize("p", [0.0001, 0.1, 0.5, 0.9, 0.9999])
def test_binomial_tail_agrees_with_independent_small_exact_sum(n: int, p: float) -> None:
    for k in range(n + 2):
        expected = math.fsum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k, n + 1))
        assert _binomial_upper_tail(k, n, p) == pytest.approx(expected, abs=1e-12)

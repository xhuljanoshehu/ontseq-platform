"""Copy-number state semantics.

The CNV subsystem compares *states per base*, not tool-specific numbers, because
different callers report copy number on incompatible scales (absolute integer copy
number, ploidy-corrected copy number, log2 ratio, or a categorical call). Reducing every
representation to one explicit state vocabulary is what makes methods comparable.

This module is dependency-free on purpose so that the state vocabulary can be unit
tested without the contract layer.
"""

from __future__ import annotations

from enum import StrEnum


class CopyNumberState(StrEnum):
    """Per-base copy-number state.

    ``NO_CALL`` is a first-class state and never a synonym for ``NEUTRAL``. A region a
    method could not assess must stay distinguishable from a region it assessed and
    found unchanged, otherwise sensitivity and specificity are both unmeasurable.
    """

    HOMOZYGOUS_LOSS = "homozygous_loss"
    LOSS = "loss"
    NEUTRAL = "neutral"
    GAIN = "gain"
    HIGH_AMPLIFICATION = "high_amplification"
    COPY_NEUTRAL_LOH = "copy_neutral_loh"
    NO_CALL = "no_call"


class StateDirection(StrEnum):
    """Coarse direction used when exact state agreement is too strict."""

    LOSS = "loss"
    NEUTRAL = "neutral"
    GAIN = "gain"
    UNDETERMINED = "undetermined"


class ConcordanceMode(StrEnum):
    """How two states are compared.

    ``DIRECTIONAL`` is the default for detection because clinical karyotyping asks
    whether a region is lost or gained, not whether two methods agreed on an exact
    integer copy number at low coverage. ``STRICT`` is retained for method development
    where exact state agreement is the question being studied.
    """

    STRICT = "strict"
    DIRECTIONAL = "directional"


#: States that assert something about copy number. ``NO_CALL`` is deliberately absent.
ASSERTED_STATES: frozenset[CopyNumberState] = frozenset(
    {
        CopyNumberState.HOMOZYGOUS_LOSS,
        CopyNumberState.LOSS,
        CopyNumberState.NEUTRAL,
        CopyNumberState.GAIN,
        CopyNumberState.HIGH_AMPLIFICATION,
        CopyNumberState.COPY_NEUTRAL_LOH,
    }
)

#: States that represent a copy-number change away from the sample baseline.
ALTERED_STATES: frozenset[CopyNumberState] = frozenset(
    {
        CopyNumberState.HOMOZYGOUS_LOSS,
        CopyNumberState.LOSS,
        CopyNumberState.GAIN,
        CopyNumberState.HIGH_AMPLIFICATION,
        CopyNumberState.COPY_NEUTRAL_LOH,
    }
)

_DIRECTIONS: dict[CopyNumberState, StateDirection] = {
    CopyNumberState.HOMOZYGOUS_LOSS: StateDirection.LOSS,
    CopyNumberState.LOSS: StateDirection.LOSS,
    CopyNumberState.NEUTRAL: StateDirection.NEUTRAL,
    # Copy-neutral LOH carries no dosage change, so its direction is neutral. It stays a
    # separate state because it is a distinct biological finding that a dosage-only
    # method cannot detect at all.
    CopyNumberState.COPY_NEUTRAL_LOH: StateDirection.NEUTRAL,
    CopyNumberState.GAIN: StateDirection.GAIN,
    CopyNumberState.HIGH_AMPLIFICATION: StateDirection.GAIN,
    CopyNumberState.NO_CALL: StateDirection.UNDETERMINED,
}


def direction(state: CopyNumberState) -> StateDirection:
    """Return the coarse direction of a state."""
    return _DIRECTIONS[state]


def is_asserted(state: CopyNumberState) -> bool:
    """Return whether the state makes a claim about copy number."""
    return state in ASSERTED_STATES


def is_altered(state: CopyNumberState) -> bool:
    """Return whether the state differs from an unchanged baseline."""
    return state in ALTERED_STATES


def concordant(
    truth: CopyNumberState,
    query: CopyNumberState,
    mode: ConcordanceMode = ConcordanceMode.DIRECTIONAL,
) -> bool:
    """Return whether a called state agrees with a truth state under ``mode``.

    ``NO_CALL`` is never concordant with anything, including itself. Bases where either
    side is ``NO_CALL`` must be removed from the evaluable genome before scoring rather
    than being counted as agreement.
    """
    if truth == CopyNumberState.NO_CALL or query == CopyNumberState.NO_CALL:
        return False
    if mode == ConcordanceMode.STRICT:
        return truth == query
    if truth == CopyNumberState.COPY_NEUTRAL_LOH or query == CopyNumberState.COPY_NEUTRAL_LOH:
        # Copy-neutral LOH is only ever concordant with itself: collapsing it into
        # "neutral" would let a dosage-only caller claim credit for detecting it.
        return truth == query
    return direction(truth) == direction(query)


def state_from_copy_number(
    copy_number: float,
    *,
    baseline_ploidy: float = 2.0,
    loss_threshold: float = 0.5,
    gain_threshold: float = 0.5,
    homozygous_loss_maximum: float = 0.5,
    amplification_copy_number: float = 6.0,
) -> CopyNumberState:
    """Derive a state from an absolute copy number.

    ``loss_threshold`` and ``gain_threshold`` are absolute copy-number distances from
    ``baseline_ploidy``. The defaults place the neutral band at ``ploidy +/- 0.5`` so
    that a diploid sample calls ``CN < 1.5`` a loss and ``CN > 2.5`` a gain, which is the
    conventional rounding boundary for integer copy number.

    The thresholds are engineering defaults for benchmarking. They are not validated
    reportability thresholds and must not be reused as such.
    """
    if copy_number < 0:
        raise ValueError("copy number must not be negative")
    if baseline_ploidy <= 0:
        raise ValueError("baseline ploidy must be positive")
    if loss_threshold <= 0 or gain_threshold <= 0:
        raise ValueError("loss and gain thresholds must be positive")
    if copy_number <= homozygous_loss_maximum:
        return CopyNumberState.HOMOZYGOUS_LOSS
    if copy_number >= amplification_copy_number:
        return CopyNumberState.HIGH_AMPLIFICATION
    if copy_number < baseline_ploidy - loss_threshold:
        return CopyNumberState.LOSS
    if copy_number > baseline_ploidy + gain_threshold:
        return CopyNumberState.GAIN
    return CopyNumberState.NEUTRAL


def expected_mixture_copy_number(
    tumor_copy_number: float,
    *,
    tumor_fraction: float,
    normal_copy_number: float = 2.0,
) -> float:
    """Return the copy number a bulk sample presents for a given tumor fraction.

    A specimen with blast fraction ``f`` and tumor copy number ``CN_t`` behaves like a
    mixture ``f * CN_t + (1 - f) * CN_n``. This is the reason a low-blast specimen
    compresses every alteration towards the neutral baseline, and it is the model the
    dilution-series simulator and the limit-of-detection analysis both use.
    """
    if not 0.0 <= tumor_fraction <= 1.0:
        raise ValueError("tumor fraction must be between 0 and 1")
    if tumor_copy_number < 0 or normal_copy_number < 0:
        raise ValueError("copy numbers must not be negative")
    return tumor_fraction * tumor_copy_number + (1.0 - tumor_fraction) * normal_copy_number


def tumor_copy_number_from_mixture(
    observed_copy_number: float,
    *,
    tumor_fraction: float,
    normal_copy_number: float = 2.0,
) -> float:
    """Invert :func:`expected_mixture_copy_number` to recover the tumor copy number.

    Raises when the tumor fraction is zero because the observed value then carries no
    information about the tumor compartment. Returning a plausible-looking number in
    that situation would be a silent fabrication.
    """
    if not 0.0 < tumor_fraction <= 1.0:
        raise ValueError("tumor fraction must be greater than 0 to invert the mixture")
    if observed_copy_number < 0 or normal_copy_number < 0:
        raise ValueError("copy numbers must not be negative")
    recovered = (
        observed_copy_number - (1.0 - tumor_fraction) * normal_copy_number
    ) / tumor_fraction
    return max(0.0, recovered)

"""Tumour fraction, copy-number baseline and subclone resolution from allele fractions.

A sample is a mixture. Every measurement this pipeline makes -- a log2 copy-number ratio, a
methylation beta value, a variant allele fraction -- is diluted by the normal cells in the
tube, and none of them says by how much. Reporting a copy-number call without the tumour
fraction is reporting a ratio and calling it a count.

The estimator has to be independent of copy number, or it is circular: the tumour fraction is
what converts an observed ratio into an integer copy number, so it cannot be derived from
those ratios. The way out is the allele fraction of a clonal, heterozygous somatic SNV in a
copy-neutral diploid region. Half the alleles in a tumour cell carry it and none in a normal
cell do, so

    VAF = f / 2          and therefore        f = 2 x VAF

for tumour fraction ``f``. A VAF of 0.25 means the tumour is half the sample. Once ``f`` is
known, an observed copy ratio can be resolved to an integer copy number, a methylation beta
value can be corrected for normal-cell dilution, and a second mutation's allele fraction
becomes a statement about what proportion of the *tumour* carries it -- the cancer cell
fraction, which is how subclones are read.

The identity holds only where the region really is copy-neutral and diploid. Under a
deletion, an amplification or copy-neutral LOH the allele fraction shifts, and a tumour
fraction taken from there is wrong in a direction nobody can see afterwards. So the
copy-number state is a required argument, and every function here refuses rather than
guesses. Copy-number calls select the regions; they never supply the number.

What this module does NOT do: call variants. It consumes allele counts a caller produced, and
ONTSeq has no small-variant caller wired in, so nothing in the pipeline calls it today. It is
written now because the arithmetic and its limits are decidable now -- in particular the
detection limits below, which say what a run of a given depth can and cannot resolve before
any dilution experiment is run.

Nothing here is validated for clinical use. No threshold is a clinical threshold, and no
default is an assay adequacy criterion. Research use only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

COPY_NEUTRAL_DIPLOID = "copy_neutral_diploid"
"""The only copy-number state in which the VAF-to-tumour-fraction identity holds."""

DETERMINABLE = "determinable"
NOT_COPY_NEUTRAL = "region_is_not_copy_neutral_diploid"
DEPTH_TOO_LOW = "depth_below_the_requested_minimum"
VAF_EXCEEDS_HETEROZYGOUS_MODEL = "allele_fraction_exceeds_a_clonal_heterozygous_model"
NO_VARIANT_READS = "no_variant_reads_observed"

TECHNICAL_MIN_DEPTH = 30
"""Technical default only. Not an assay adequacy threshold and not clinically validated.

Chosen so a Wilson interval on the allele fraction is narrower than roughly +/-0.2, which is
the point below which a tumour-fraction estimate stops constraining anything. A laboratory
that validates this assay sets its own value.
"""

TECHNICAL_ERROR_RATE = 0.01
"""Technical default only: an order-of-magnitude per-base error rate for modern ONT reads.

Used to place the detection floor. The real value depends on basecaller, chemistry, context
and the caller's own error model, and must be measured per assay rather than assumed.
"""


class QuantitationError(ValueError):
    """Raised when an input cannot describe an observation at all."""


@dataclass(frozen=True)
class Interval:
    """A closed interval. Always reported alongside a point estimate, never instead of one."""

    low: float
    high: float

    def __post_init__(self) -> None:
        if self.low > self.high:
            raise QuantitationError(f"interval is reversed: [{self.low}, {self.high}]")

    @property
    def width(self) -> float:
        return self.high - self.low


@dataclass(frozen=True)
class AlleleObservation:
    """Read counts supporting a variant at one site."""

    variant_reads: int
    total_reads: int

    def __post_init__(self) -> None:
        if self.total_reads <= 0:
            raise QuantitationError("total_reads must be positive")
        if self.variant_reads < 0:
            raise QuantitationError("variant_reads cannot be negative")
        if self.variant_reads > self.total_reads:
            raise QuantitationError(
                f"variant_reads ({self.variant_reads}) exceeds total_reads ({self.total_reads})"
            )

    @property
    def vaf(self) -> float:
        return self.variant_reads / self.total_reads


def wilson_interval(observation: AlleleObservation, *, confidence: float = 0.95) -> Interval:
    """Score interval for a binomial proportion.

    Chosen over the textbook normal interval because allele fractions are routinely near 0
    and at the depths this assay reaches; the normal interval misbehaves there, returning
    bounds below zero and a zero-width interval when no variant read is seen.
    """
    if not 0.0 < confidence < 1.0:
        raise QuantitationError(f"confidence must be in (0, 1), got {confidence}")
    z = _normal_quantile(1.0 - (1.0 - confidence) / 2.0)
    n = observation.total_reads
    k = observation.variant_reads
    denominator = n + z * z
    centre = (k + z * z / 2.0) / denominator
    spread = z / denominator * math.sqrt(k * (n - k) / n + z * z / 4.0)
    return Interval(low=max(0.0, centre - spread), high=min(1.0, centre + spread))


@dataclass(frozen=True)
class TumourFractionEstimate:
    """What proportion of the sample is tumour, or why that cannot be said."""

    status: str
    point: float | None
    interval: Interval | None
    observation: AlleleObservation
    copy_number_state: str

    @property
    def determinable(self) -> bool:
        return self.status == DETERMINABLE

    def reason(self) -> str:
        if self.status == DETERMINABLE:
            assert self.interval is not None and self.point is not None
            return (
                f"Tumour fraction {self.point:.3f} "
                f"(95% CI {self.interval.low:.3f}-{self.interval.high:.3f}) from "
                f"{self.observation.variant_reads}/{self.observation.total_reads} reads in a "
                "copy-neutral diploid region."
            )
        if self.status == NOT_COPY_NEUTRAL:
            return (
                f"Region is {self.copy_number_state!r}, not {COPY_NEUTRAL_DIPLOID!r}. The "
                "VAF-to-tumour-fraction identity does not hold under altered copy number or "
                "LOH, and the error would not be visible in the result."
            )
        if self.status == DEPTH_TOO_LOW:
            return (
                f"{self.observation.total_reads} reads is below the requested minimum depth. "
                "An allele fraction from too few reads carries an interval wide enough to "
                "admit almost any tumour fraction."
            )
        if self.status == NO_VARIANT_READS:
            return (
                "No variant reads observed. This bounds the tumour fraction from above but "
                "does not estimate it, and it is not evidence that the tumour is absent."
            )
        return (
            f"Allele fraction {self.observation.vaf:.3f} implies a tumour fraction above 1. "
            "The site is not a clonal heterozygous variant in a diploid region: it may be "
            "homozygous, in LOH, in an amplification, or germline."
        )


def tumour_fraction_from_clonal_snv(
    observation: AlleleObservation,
    *,
    copy_number_state: str,
    min_depth: int = TECHNICAL_MIN_DEPTH,
    confidence: float = 0.95,
) -> TumourFractionEstimate:
    """Estimate tumour fraction as ``2 x VAF``, or refuse and say why.

    ``copy_number_state`` has no default. The identity is only valid in a copy-neutral diploid
    region, and a default would let a caller obtain a confident, wrong number by omission --
    which is precisely the failure this module exists to prevent.

    The caller is responsible for the variant being somatic and clonal. Neither is checkable
    from read counts alone: a germline heterozygous SNP sits at VAF 0.5 regardless of tumour
    content and would report a tumour fraction of 1.
    """

    def refuse(status: str) -> TumourFractionEstimate:
        return TumourFractionEstimate(
            status=status,
            point=None,
            interval=None,
            observation=observation,
            copy_number_state=copy_number_state,
        )

    if copy_number_state != COPY_NEUTRAL_DIPLOID:
        return refuse(NOT_COPY_NEUTRAL)
    if observation.total_reads < min_depth:
        return refuse(DEPTH_TOO_LOW)
    if observation.variant_reads == 0:
        return refuse(NO_VARIANT_READS)
    if 2.0 * observation.vaf > 1.0:
        return refuse(VAF_EXCEEDS_HETEROZYGOUS_MODEL)

    band = wilson_interval(observation, confidence=confidence)
    return TumourFractionEstimate(
        status=DETERMINABLE,
        point=2.0 * observation.vaf,
        interval=Interval(low=min(1.0, 2.0 * band.low), high=min(1.0, 2.0 * band.high)),
        observation=observation,
        copy_number_state=copy_number_state,
    )


@dataclass(frozen=True)
class CancerCellFraction:
    """What proportion of the *tumour* carries a variant. 1.0 means clonal."""

    status: str
    point: float | None
    interval: Interval | None

    @property
    def determinable(self) -> bool:
        return self.status == DETERMINABLE

    def subclonal_at(self, threshold: float) -> bool | None:
        """Whether the variant is subclonal at a caller-supplied threshold.

        Returns ``None`` when the fraction is not determinable, so that "cannot tell" cannot
        be silently read as "clonal".
        """
        if self.point is None:
            return None
        return self.point < threshold


def cancer_cell_fraction(
    observation: AlleleObservation,
    *,
    tumour_fraction: float,
    copy_number_state: str,
    min_depth: int = TECHNICAL_MIN_DEPTH,
    confidence: float = 0.95,
) -> CancerCellFraction:
    """Resolve an allele fraction into a cancer cell fraction: ``CCF = 2 x VAF / f``.

    This is how subclones are read. A variant at half the allele fraction of the clonal ones
    sits in roughly half the tumour cells. The interval matters more than the point here --
    at the depths this assay reaches, two subclones can be statistically indistinguishable.
    """
    if not 0.0 < tumour_fraction <= 1.0:
        raise QuantitationError(f"tumour_fraction must be in (0, 1], got {tumour_fraction}")
    if copy_number_state != COPY_NEUTRAL_DIPLOID:
        return CancerCellFraction(status=NOT_COPY_NEUTRAL, point=None, interval=None)
    if observation.total_reads < min_depth:
        return CancerCellFraction(status=DEPTH_TOO_LOW, point=None, interval=None)

    band = wilson_interval(observation, confidence=confidence)
    scale = 2.0 / tumour_fraction
    return CancerCellFraction(
        status=DETERMINABLE,
        point=min(1.0, observation.vaf * scale),
        interval=Interval(low=min(1.0, band.low * scale), high=min(1.0, band.high * scale)),
    )


def _binomial_upper_tail(k: int, n: int, p: float) -> float:
    """``P(X >= k)`` for ``X ~ Binomial(n, p)``, computed exactly.

    Exact rather than a normal approximation because the whole question here lives in the
    tail, at small ``k`` and small ``p``, where the approximation is worst.
    """
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    return math.fsum(math.comb(n, i) * p**i * (1.0 - p) ** (n - i) for i in range(k, n + 1))


def minimum_variant_reads(
    depth: int, *, error_rate: float = TECHNICAL_ERROR_RATE, alpha: float = 0.01
) -> int:
    """Fewest variant reads that sequencing error alone would produce with probability <= alpha.

    This is the floor a real variant has to clear. It is a property of depth and error rate,
    not of any caller: a caller with a better error model does better than this, never worse.
    """
    if depth <= 0:
        raise QuantitationError("depth must be positive")
    if not 0.0 <= error_rate < 1.0:
        raise QuantitationError(f"error_rate must be in [0, 1), got {error_rate}")
    if not 0.0 < alpha < 1.0:
        raise QuantitationError(f"alpha must be in (0, 1), got {alpha}")
    for k in range(1, depth + 1):
        if _binomial_upper_tail(k, depth, error_rate) <= alpha:
            return k
    return depth + 1


def minimum_detectable_vaf(
    depth: int,
    *,
    error_rate: float = TECHNICAL_ERROR_RATE,
    alpha: float = 0.01,
    power: float = 0.95,
) -> float | None:
    """Smallest allele fraction detectable at ``depth`` with the given power, or ``None``.

    ``None`` means no allele fraction is detectable at that depth -- the read count needed to
    clear the error floor exceeds the depth available. That is the honest answer for shallow
    coverage, and it is the reason this assay's off-target fraction cannot carry somatic
    variant calling however good the caller is.
    """
    needed = minimum_variant_reads(depth, error_rate=error_rate, alpha=alpha)
    if needed > depth:
        return None
    low, high = 0.0, 1.0
    for _ in range(60):
        middle = (low + high) / 2.0
        if _binomial_upper_tail(needed, depth, middle) >= power:
            high = middle
        else:
            low = middle
    return high


def minimum_detectable_cancer_cell_fraction(
    depth: int,
    *,
    tumour_fraction: float,
    error_rate: float = TECHNICAL_ERROR_RATE,
    alpha: float = 0.01,
    power: float = 0.95,
) -> float | None:
    """Smallest subclone resolvable at a given depth and tumour fraction, or ``None``.

    The question a haematologist actually asks -- "how small a subclone can this see?" -- and
    it has two answers stacked: the allele fraction has to clear the error floor, and the
    tumour has to be a large enough part of the sample for that allele fraction to correspond
    to a small subclone. A value above 1 means no subclone is resolvable: even a variant in
    every tumour cell would sit below the detection floor.
    """
    floor = minimum_detectable_vaf(depth, error_rate=error_rate, alpha=alpha, power=power)
    if floor is None:
        return None
    if not 0.0 < tumour_fraction <= 1.0:
        raise QuantitationError(f"tumour_fraction must be in (0, 1], got {tumour_fraction}")
    resolvable = 2.0 * floor / tumour_fraction
    return resolvable if resolvable <= 1.0 else None


def expected_vaf(*, tumour_fraction: float, cancer_cell_fraction: float) -> float:
    """Forward model: the allele fraction a variant should show. Used to check the inverses.

    Also the generator for an in-silico dilution series -- pick a tumour fraction, pick a
    subclone, and this is the allele fraction the mixed reads should carry.
    """
    if not 0.0 <= tumour_fraction <= 1.0:
        raise QuantitationError(f"tumour_fraction must be in [0, 1], got {tumour_fraction}")
    if not 0.0 <= cancer_cell_fraction <= 1.0:
        raise QuantitationError(
            f"cancer_cell_fraction must be in [0, 1], got {cancer_cell_fraction}"
        )
    return tumour_fraction * cancer_cell_fraction / 2.0


def copy_number_from_ratio(*, ratio: float, tumour_fraction: float) -> float:
    """Resolve an observed copy ratio into a tumour copy number, given the tumour fraction.

    ``ratio = (f x CN + (1 - f) x 2) / 2``, inverted for ``CN``. This is what "setting the
    baseline to 2" means once dilution is accounted for: at a tumour fraction of 0.5 a true
    one-copy deletion shows a ratio of 0.75, not 0.5, and a fixed ratio threshold would miss
    it. The returned value is continuous on purpose -- rounding it to an integer is a
    modelling decision with its own uncertainty, and it is not made here.
    """
    if not 0.0 < tumour_fraction <= 1.0:
        raise QuantitationError(f"tumour_fraction must be in (0, 1], got {tumour_fraction}")
    if ratio < 0.0:
        raise QuantitationError(f"ratio cannot be negative, got {ratio}")
    return (2.0 * ratio - 2.0 * (1.0 - tumour_fraction)) / tumour_fraction


def _normal_quantile(probability: float) -> float:
    """Inverse standard normal CDF (Acklam's rational approximation).

    Carried here rather than imported so this module keeps no dependencies; accurate to about
    1e-9, far beyond what read counts in the tens or hundreds can justify.
    """
    if not 0.0 < probability < 1.0:
        raise QuantitationError(f"probability must be in (0, 1), got {probability}")
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
    break_low, break_high = 0.02425, 1.0 - 0.02425
    if probability < break_low:
        q = math.sqrt(-2.0 * math.log(probability))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )
    if probability > break_high:
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


ON_TARGET_DEPTH = 80
"""Measured mean on-target depth of run 260611_RAD114_AS_S700: 1.37 Gb over a 17.03 Mb panel."""

OFF_TARGET_DEPTH = 9
"""Measured genome-wide off-target depth of the same run, 8.8x, rounded to whole reads."""

_DEPTH_LADDER = (OFF_TARGET_DEPTH, 20, 30, 50, ON_TARGET_DEPTH, 100, 200, 500)
_FRACTION_LADDER = (0.2, 0.3, 0.5, 0.7, 0.9)


def _shallowest_resolvable(depth: int) -> float | None:
    """Best cancer cell fraction any tumour fraction in the ladder reaches at this depth."""
    reachable = [
        value
        for fraction in _FRACTION_LADDER
        if (value := minimum_detectable_cancer_cell_fraction(depth, tumour_fraction=fraction))
        is not None
    ]
    return min(reachable) if reachable else None


def _percent(value: float | None) -> str:
    return "not resolvable" if value is None else f"{value * 100:.1f}%"


def format_report() -> str:
    """Render the quantitative model and its limits as Markdown, computed rather than asserted."""
    lines = [
        "# Quantitative model: tumour fraction, copy-number baseline and subclones",
        "",
        "GENERATED FILE - do not edit by hand. Regenerate with:",
        "",
        "    python -m ontseq_platform.quantitation",
        "",
        "Every number below is computed by `ontseq_platform.quantitation` at generation time,",
        "so the document cannot drift from the code. Nothing here is validated for clinical",
        "use, and no value is an assay adequacy threshold. Research use only.",
        "",
        "## The problem",
        "",
        "A sample is a mixture of tumour and normal cells. A copy-number ratio, a methylation",
        "beta value and an allele fraction are all diluted by the normal fraction, and none of",
        "them states by how much. A copy-number call reported without the tumour fraction is a",
        "ratio presented as a count.",
        "",
        "## Why the estimator cannot come from copy number",
        "",
        "The tumour fraction is what converts an observed ratio into an integer copy number, so",
        "deriving it from those same ratios is circular. The anchor has to be independent of",
        "copy number: the allele fraction of a clonal, heterozygous somatic SNV in a",
        "copy-neutral diploid region. Half the alleles in a tumour cell carry it and none in a",
        "normal cell do, so `VAF = f / 2`, and therefore `f = 2 x VAF`.",
        "",
        "The identity holds *only* where the region is genuinely copy-neutral and diploid. Under",
        "a deletion, an amplification or copy-neutral LOH the allele fraction shifts and the",
        "resulting tumour fraction is wrong invisibly. Copy-number calls select the regions the",
        "estimate may be taken from; they never supply the number. `copy_number_state` is",
        "therefore a required argument with no default.",
        "",
        "## What follows once the tumour fraction is known",
        "",
        "1. **Copy-number baseline.** An observed ratio resolves to a tumour copy number.",
        "2. **Methylation.** Beta values can be corrected for normal-cell dilution.",
        "3. **Subclones.** A second variant's allele fraction becomes a cancer cell fraction,",
        "   `CCF = 2 x VAF / f` - the proportion of the tumour carrying it.",
        "",
        "## Why a fixed copy-ratio threshold is unsafe",
        "",
        "The same true copy number produces a different observed ratio at every tumour fraction.",
        "A threshold chosen for pure tumour misses real events in a diluted sample.",
        "",
        "| True copy number | " + " | ".join(f"f = {f:.0%}" for f in _FRACTION_LADDER) + " |",
        "|---|" + "---|" * len(_FRACTION_LADDER),
    ]
    for copy_number in (0, 1, 2, 3, 4):
        cells = [
            f"{(fraction * copy_number + (1.0 - fraction) * 2.0) / 2.0:.2f}"
            for fraction in _FRACTION_LADDER
        ]
        lines.append(f"| {copy_number} | " + " | ".join(cells) + " |")
    lines += [
        "",
        "A one-copy deletion shows a ratio of 0.75 at half tumour content, not 0.50. A filter",
        "set at a fixed distance from 1.0 is a filter whose meaning changes with every sample.",
        "",
        "## Detection limits by depth",
        "",
        f"Per-base error rate {TECHNICAL_ERROR_RATE:.0%} (technical default), alpha 0.01,",
        "power 0.95. `Min reads` is the fewest variant reads that sequencing error",
        "alone would produce with probability at most alpha - the floor any real variant must",
        "clear. It is a property of depth and error rate, not of a caller: a caller with a",
        "better error model does better than this, never worse.",
        "",
        "| Depth | Min variant reads | Smallest detectable VAF |",
        "|---|---|---|",
    ]
    for depth in _DEPTH_LADDER:
        marker = ""
        if depth == ON_TARGET_DEPTH:
            marker = " (measured on-target)"
        elif depth == OFF_TARGET_DEPTH:
            marker = " (measured off-target)"
        lines.append(
            f"| {depth}x{marker} | {minimum_variant_reads(depth)} | "
            f"{_percent(minimum_detectable_vaf(depth))} |"
        )

    lines += [
        "",
        "## Smallest resolvable subclone",
        "",
        "The question a haematologist asks. Two limits stack: the allele fraction has to clear",
        "the error floor, and the tumour has to be a large enough part of the sample for that",
        "allele fraction to correspond to a small subclone. Values are cancer cell fractions -",
        "the proportion of *tumour* cells carrying the variant.",
        "",
        "| Depth | " + " | ".join(f"f = {f:.0%}" for f in _FRACTION_LADDER) + " |",
        "|---|" + "---|" * len(_FRACTION_LADDER),
    ]
    for depth in _DEPTH_LADDER:
        cells = [
            _percent(minimum_detectable_cancer_cell_fraction(depth, tumour_fraction=fraction))
            for fraction in _FRACTION_LADDER
        ]
        lines.append(f"| {depth}x | " + " | ".join(cells) + " |")

    on_target_vaf = minimum_detectable_vaf(ON_TARGET_DEPTH)
    on_target_ccf = minimum_detectable_cancer_cell_fraction(ON_TARGET_DEPTH, tumour_fraction=0.5)
    lines += [
        "",
        "Read against this assay's measured depths:",
        "",
        f"* **On-target ({ON_TARGET_DEPTH}x)** resolves allele fractions down to about "
        f"{_percent(on_target_vaf)}. At a tumour fraction of 50% that is a subclone of about "
        f"{_percent(on_target_ccf)} of tumour cells - major subclones only, not minor ones.",
        f"* **Off-target ({OFF_TARGET_DEPTH}x)** resolves nothing below "
        f"{_percent(_shallowest_resolvable(OFF_TARGET_DEPTH))}, and only at the highest tumour",
        "  fraction in this table; below that it resolves no subclone at all. That is a",
        "  property of the depth, not of the caller, and no choice of software changes it.",
        "",
        "So the entire quantitative model lives inside the panel. Genome-wide, this assay can",
        "carry copy number and methylation, and cannot carry allele-fraction quantitation.",
        "",
        "## What is not implemented",
        "",
        "**All of the above consumes allele counts, and ONTSeq has no small-variant caller",
        "wired in.** There is no stage, no pinned caller, and no validation. This module is the",
        "arithmetic and its limits, written now because both are decidable now; nothing in the",
        "pipeline calls it. The dependency is the same one that blocks seven of the",
        "twenty-four drafted guideline criteria.",
        "",
        "Also not addressed here, and each a separate problem:",
        "",
        "* Distinguishing somatic from germline without a matched normal. A germline",
        "  heterozygous SNP sits at VAF 0.5 whatever the tumour content and would report a",
        "  tumour fraction of 1. The module refuses that value rather than returning it, which",
        "  is a guard, not a solution.",
        "* Clonality. The module cannot tell a clonal variant from a subclonal one without",
        "  being told the tumour fraction, and the tumour fraction is what a clonal variant is",
        "  used to establish. Breaking that requires either a known-clonal marker or joint",
        "  estimation over many sites.",
        "* Multi-site inference. Real purity and subclone calling fits a mixture model over",
        "  many variants at once. These functions are per-site and deliberately so.",
        "",
        "## In-silico dilution",
        "",
        "`expected_vaf` is the forward model, and it is the generator for a dilution series:",
        "choose a tumour fraction and a subclone size, and it gives the allele fraction the",
        "mixed reads should carry. Mixing reads from a known-positive sample with reads from a",
        "normal at defined ratios measures the sensitivity curve against a known truth, needs",
        "no new patient material and no new sequencing run, and validates the table above",
        "rather than trusting it. The table is a ceiling derived from counting statistics; a",
        "real caller will do worse, and the gap between them is the number worth having.",
        "",
    ]
    return "\n".join(lines)


def _repo_root() -> Path:
    candidate = Path(__file__).resolve().parents[2]
    if (candidate / "docs").is_dir():
        return candidate
    return Path.cwd()


def main() -> int:
    root = _repo_root()
    destination = root / "docs" / "QUANTITATIVE_MODEL.md"
    destination.write_text(format_report() + "\n", encoding="utf-8")
    print(f"wrote {destination.relative_to(root)}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

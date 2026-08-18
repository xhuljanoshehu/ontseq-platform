"""Cross-run aggregation, dilution-series analysis and limit-of-detection estimation.

A single evaluation answers "how did this method do on this sample". The questions the
project actually needs answered are different in kind:

- how does sensitivity change as coverage falls;
- how does it change as blast fraction falls;
- at what blast fraction does detection drop below an acceptable rate;
- does method A beat method B on the same samples, or only on easier ones.

All four require aggregation across many evaluations grouped by stratum, which is why
:class:`ontseq_platform.cnv.models.CnvStrata` is typed rather than free-form.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from datetime import UTC, datetime
from typing import Literal

from pydantic import Field, model_validator

from ..models import StrictModel
from .models import CnvEvaluationReport, ProportionResult
from .stats import (
    ProportionEstimate,
    fit_logistic,
    logistic_threshold,
    mcnemar_exact,
    minimum_attainable_p_value,
    wilson_interval,
)


class StratumSummary(StrictModel):
    """Pooled outcome for one stratum across several evaluations."""

    label: str
    evaluations: int = Field(ge=0)
    samples: int = Field(ge=0)
    detected_events: int = Field(ge=0)
    missed_events: int = Field(ge=0)
    not_assessable_events: int = Field(ge=0)
    unconfirmed_events: int = Field(ge=0)
    detection_rate: ProportionResult
    #: Base-level concordance averaged over evaluations, weighted by evaluable bases.
    weighted_base_concordance: float | None = Field(default=None, ge=0, le=1)
    evaluable_bases: int = Field(ge=0)
    excluded_bases: int = Field(ge=0)


class LimitOfDetection(StrictModel):
    """A detection limit with an explicit basis and an explicit failure mode.

    ``model_based_value`` is ``None`` whenever the logistic fit did not converge, the
    design was perfectly separated, or fewer than two dilution levels were supplied.
    Reporting a number in those cases would invent a detection limit the data cannot
    support, which is precisely the kind of claim this repository is built to avoid.
    """

    target_detection_rate: float = Field(gt=0, lt=1)
    predictor: Literal["tumor_fraction", "mean_coverage_x"]
    #: Lowest observed level whose *lower* confidence bound still meets the target.
    empirical_value: float | None = None
    model_based_value: float | None = None
    model_converged: bool = False
    levels_used: int = Field(ge=0)
    note: str


class SpecimenClustering(StrictModel):
    """How the scored events are distributed over specimens.

    Every interval in this report is computed over *events*, and several events routinely
    come from the same specimen. Events within one specimen share its purity, its library,
    its coverage and its artefacts, so they are not independent observations. Treating
    them as independent narrows every confidence interval: the interval describes a
    hypothetical population of independent events that does not exist.

    The correction is a specimen-aware endpoint, which is a study-design decision and not
    something this module can make on the caller's behalf. What it can do is refuse to
    hide the problem: the numbers below let a reader see the clustering, and
    ``intervals_are_anticonservative`` states in one field that the intervals beside them
    are narrower than the data support.
    """

    specimens: int = Field(ge=0)
    events: int = Field(ge=0)
    #: Largest number of events contributed by any single specimen.
    largest_specimen_events: int = Field(ge=0)
    mean_events_per_specimen: float | None = Field(default=None, ge=0)
    #: True whenever any specimen contributed more than one event.
    intervals_are_anticonservative: bool = False


class CnvAggregateReport(StrictModel):
    """Aggregated benchmark evidence across many evaluations."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    aggregate_id: str = Field(min_length=1)
    method: str
    method_version: str
    evaluations: int = Field(ge=0)
    overall_detection_rate: ProportionResult
    #: Event-level intervals are reported alongside the clustering that qualifies them.
    clustering: SpecimenClustering | None = None
    #: Detection rate computed one specimen at a time, then averaged over specimens, so
    #: each specimen carries equal weight regardless of how many events it contributed.
    #: Reported next to the event-level rate rather than replacing it: the two answer
    #: different questions, and a specimen-heavy dataset can pull them far apart.
    specimen_level_detection_rate: ProportionResult | None = None
    by_tumor_fraction: list[StratumSummary] = Field(default_factory=list)
    by_coverage: list[StratumSummary] = Field(default_factory=list)
    by_size_class: list[StratumSummary] = Field(default_factory=list)
    by_data_basis: list[StratumSummary] = Field(default_factory=list)
    limits_of_detection: list[LimitOfDetection] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def evaluations_are_consistent(self) -> CnvAggregateReport:
        if self.evaluations == 0 and self.overall_detection_rate.total > 0:
            raise ValueError("an aggregate with no evaluations cannot contain outcomes")
        return self


def _proportion(estimate: ProportionEstimate) -> ProportionResult:
    return ProportionResult(
        successes=estimate.successes,
        total=estimate.total,
        point=estimate.point,
        lower=estimate.lower,
        upper=estimate.upper,
        confidence_level=estimate.confidence_level,
    )


def _summarize(
    label: str, reports: Sequence[CnvEvaluationReport], confidence_level: float
) -> StratumSummary:
    detected = sum(r.detection_rate.successes for r in reports)
    assessable = sum(r.detection_rate.total for r in reports)
    not_assessable = sum(
        sum(1 for e in r.truth_events if e.outcome == "NOT_ASSESSABLE") for r in reports
    )
    unconfirmed = sum(sum(1 for e in r.query_events if e.outcome == "UNCONFIRMED") for r in reports)
    evaluable = sum(r.base_level.evaluable_bases for r in reports)
    weighted = (
        sum((r.base_level.concordance or 0.0) * r.base_level.evaluable_bases for r in reports)
        / evaluable
        if evaluable
        else None
    )
    return StratumSummary(
        label=label,
        evaluations=len(reports),
        samples=len({r.sample_id for r in reports}),
        detected_events=detected,
        missed_events=assessable - detected,
        not_assessable_events=not_assessable,
        unconfirmed_events=unconfirmed,
        detection_rate=_proportion(
            wilson_interval(detected, assessable, confidence_level=confidence_level)
        ),
        weighted_base_concordance=weighted,
        evaluable_bases=evaluable,
        excluded_bases=sum(r.partition.excluded_bases for r in reports),
    )


def _group(
    reports: Sequence[CnvEvaluationReport],
    key: Callable[[CnvEvaluationReport], str | None],
    confidence_level: float,
) -> list[StratumSummary]:
    buckets: dict[str, list[CnvEvaluationReport]] = {}
    for report in reports:
        label = key(report)
        if label is not None:
            buckets.setdefault(label, []).append(report)
    return [_summarize(label, buckets[label], confidence_level) for label in sorted(buckets)]


def _size_class_summaries(
    reports: Sequence[CnvEvaluationReport], confidence_level: float
) -> list[StratumSummary]:
    """Pool per-size-class detection across evaluations.

    Size class lives inside each report rather than on the report itself, so this
    aggregation walks the nested strata instead of grouping whole reports.
    """
    detected: dict[str, int] = {}
    missed: dict[str, int] = {}
    not_assessable: dict[str, int] = {}
    counts: dict[str, int] = {}
    for report in reports:
        for stratum in report.detection_by_size_class:
            detected[stratum.label] = detected.get(stratum.label, 0) + stratum.detected
            missed[stratum.label] = missed.get(stratum.label, 0) + stratum.missed
            not_assessable[stratum.label] = (
                not_assessable.get(stratum.label, 0) + stratum.not_assessable
            )
            counts[stratum.label] = counts.get(stratum.label, 0) + 1
    summaries: list[StratumSummary] = []
    for label in sorted(counts):
        hits = detected.get(label, 0)
        total = hits + missed.get(label, 0)
        summaries.append(
            StratumSummary(
                label=label,
                evaluations=counts[label],
                samples=0,
                detected_events=hits,
                missed_events=missed.get(label, 0),
                not_assessable_events=not_assessable.get(label, 0),
                unconfirmed_events=0,
                detection_rate=_proportion(
                    wilson_interval(hits, total, confidence_level=confidence_level)
                ),
                weighted_base_concordance=None,
                evaluable_bases=0,
                excluded_bases=0,
            )
        )
    return summaries


def estimate_limit_of_detection(
    levels: Sequence[tuple[float, int, int]],
    *,
    predictor: Literal["tumor_fraction", "mean_coverage_x"],
    target_detection_rate: float = 0.95,
    confidence_level: float = 0.95,
) -> LimitOfDetection:
    """Estimate the level at which detection reaches ``target_detection_rate``.

    ``levels`` are ``(predictor_value, detected, assessable)`` triples.

    Two estimates are produced and both are reported, because they fail differently. The
    empirical value is the lowest tested level whose *lower* confidence bound still meets
    the target; it is conservative, cannot extrapolate, and is limited to levels actually
    tested. The model-based value comes from a logistic fit and can interpolate, but is
    only meaningful when the fit converged on a non-separated design.
    """
    usable = [(x, d, n) for x, d, n in levels if n > 0]
    empirical: float | None = None
    for value, detected, total in sorted(usable):
        interval = wilson_interval(detected, total, confidence_level=confidence_level)
        if interval.lower is not None and interval.lower >= target_detection_rate:
            empirical = value
            break

    fit = fit_logistic(
        [x for x, _, _ in usable],
        [d for _, d, _ in usable],
        [n for _, _, n in usable],
    )
    model_value = logistic_threshold(fit, target_detection_rate) if fit else None
    if model_value is not None and model_value < 0:
        model_value = None

    if not usable:
        note = "No level contributed an assessable event; no limit of detection is defined."
    elif fit is None:
        note = (
            "The logistic fit was not attempted: fewer than two distinct levels, or a "
            "perfectly separated design. Only the empirical value is meaningful."
        )
    elif not fit.converged:
        note = "The logistic fit did not converge; the model-based value is withheld."
    else:
        note = (
            "Model-based value interpolates between tested levels and must not be "
            "extrapolated beyond them."
        )

    return LimitOfDetection(
        target_detection_rate=target_detection_rate,
        predictor=predictor,
        empirical_value=empirical,
        model_based_value=model_value,
        model_converged=bool(fit and fit.converged),
        levels_used=len(usable),
        note=note,
    )


def aggregate(
    reports: Sequence[CnvEvaluationReport],
    *,
    aggregate_id: str,
    confidence_level: float = 0.95,
    target_detection_rate: float = 0.95,
) -> CnvAggregateReport:
    """Pool evaluations of one method into a stratified benchmark summary.

    Raises when the reports mix methods. Pooling several methods into one detection rate
    would produce a number that describes nothing; use one aggregate per method and
    compare the aggregates.
    """
    if not reports:
        raise ValueError("cannot aggregate an empty report list")
    methods = {(r.method, r.method_version) for r in reports}
    if len(methods) > 1:
        raise ValueError(
            "aggregate() requires a single method; pooling methods yields a rate that "
            f"describes no method: {sorted(methods)}"
        )
    method, method_version = methods.pop()

    warnings: list[str] = []
    detected = sum(r.detection_rate.successes for r in reports)
    assessable = sum(r.detection_rate.total for r in reports)
    if assessable == 0:
        warnings.append(
            "No truth event was assessable in any evaluation; the detection rate is "
            "undefined rather than zero."
        )

    # Events cluster inside specimens, and every interval in this report is computed over
    # events. Rather than silently presenting intervals that are too narrow, count the
    # clustering, say so, and put a specimen-weighted rate beside the event-weighted one.
    by_specimen: dict[str, tuple[int, int]] = {}
    for report in reports:
        hits, total = by_specimen.get(report.sample_id, (0, 0))
        by_specimen[report.sample_id] = (
            hits + report.detection_rate.successes,
            total + report.detection_rate.total,
        )
    specimen_counts = [total for _, total in by_specimen.values()]
    clustered = any(count > 1 for count in specimen_counts)
    clustering = SpecimenClustering(
        specimens=len(by_specimen),
        events=assessable,
        largest_specimen_events=max(specimen_counts, default=0),
        mean_events_per_specimen=(assessable / len(by_specimen) if by_specimen else None),
        intervals_are_anticonservative=clustered,
    )
    if clustered:
        warnings.append(
            f"{assessable} scored event(s) come from {len(by_specimen)} specimen(s), up to "
            f"{clustering.largest_specimen_events} from one. Events within a specimen share "
            "its purity, library, coverage and artefacts, so they are not independent. "
            "Every confidence interval in this report is computed over events and is "
            "therefore narrower than the data support. Analytical validation needs a "
            "specimen-level endpoint, which is a study-design decision rather than a "
            "post-hoc correction."
        )
    # Each specimen contributes at most one success and one trial, so no specimen can
    # dominate the interval by contributing many events. This is a deliberately crude
    # specimen-level endpoint, and it is labelled as one rather than presented as the
    # cluster-robust analysis a validation study would specify.
    specimen_successes = sum(1 for hits, total in by_specimen.values() if total and hits == total)
    specimen_trials = sum(1 for _, total in by_specimen.values() if total)
    specimen_rate = (
        _proportion(
            wilson_interval(specimen_successes, specimen_trials, confidence_level=confidence_level)
        )
        if specimen_trials
        else None
    )

    tumor_levels: dict[float, tuple[int, int]] = {}
    coverage_levels: dict[float, tuple[int, int]] = {}
    for report in reports:
        if report.strata.tumor_fraction is not None:
            hits, total = tumor_levels.get(report.strata.tumor_fraction, (0, 0))
            tumor_levels[report.strata.tumor_fraction] = (
                hits + report.detection_rate.successes,
                total + report.detection_rate.total,
            )
        if report.strata.mean_coverage_x is not None:
            hits, total = coverage_levels.get(report.strata.mean_coverage_x, (0, 0))
            coverage_levels[report.strata.mean_coverage_x] = (
                hits + report.detection_rate.successes,
                total + report.detection_rate.total,
            )

    limits: list[LimitOfDetection] = []
    if len(tumor_levels) >= 2:
        limits.append(
            estimate_limit_of_detection(
                [(x, d, n) for x, (d, n) in sorted(tumor_levels.items())],
                predictor="tumor_fraction",
                target_detection_rate=target_detection_rate,
                confidence_level=confidence_level,
            )
        )
    if len(coverage_levels) >= 2:
        limits.append(
            estimate_limit_of_detection(
                [(x, d, n) for x, (d, n) in sorted(coverage_levels.items())],
                predictor="mean_coverage_x",
                target_detection_rate=target_detection_rate,
                confidence_level=confidence_level,
            )
        )

    return CnvAggregateReport(
        aggregate_id=aggregate_id,
        method=method,
        method_version=method_version,
        evaluations=len(reports),
        overall_detection_rate=_proportion(
            wilson_interval(detected, assessable, confidence_level=confidence_level)
        ),
        clustering=clustering,
        specimen_level_detection_rate=specimen_rate,
        by_tumor_fraction=_group(
            reports,
            lambda r: None if r.strata.tumor_fraction is None else f"{r.strata.tumor_fraction:g}",
            confidence_level,
        ),
        by_coverage=_group(
            reports,
            lambda r: (
                None if r.strata.mean_coverage_x is None else f"{r.strata.mean_coverage_x:g}x"
            ),
            confidence_level,
        ),
        by_size_class=_size_class_summaries(reports, confidence_level),
        by_data_basis=_group(reports, lambda r: r.data_basis.value, confidence_level),
        limits_of_detection=limits,
        warnings=warnings,
        limitations=[
            "Aggregated synthetic or public-benchmark results do not establish "
            "intended-use performance on AML specimens.",
            "Detection rates pool events of different sizes and classes unless read "
            "through the per-stratum tables.",
            "A limit of detection derived from one truth profile does not generalise to "
            "other event classes or sizes.",
            "Confidence intervals are event-level. Where a specimen contributed several "
            "events they are not independent observations, and the intervals are "
            "correspondingly optimistic; see the clustering block.",
            "The specimen-level rate counts a specimen as a success only when every "
            "assessable event in it was detected. It is a conservative screening summary, "
            "not the cluster-robust endpoint an analytical validation would pre-specify.",
        ],
    )


class PairedMethodComparison(StrictModel):
    """A paired comparison of two methods on the same truth events.

    The contingency counts are reported in full so the reader can judge the comparison
    rather than only its p-value.
    """

    schema_version: Literal["0.1.0"] = "0.1.0"
    method_a: str
    method_b: str
    paired_events: int = Field(ge=0)
    both_detected: int = Field(ge=0)
    only_a_detected: int = Field(ge=0)
    only_b_detected: int = Field(ge=0)
    neither_detected: int = Field(ge=0)
    #: Events dropped because they were not assessable under both methods.
    unpaired_events: int = Field(ge=0)
    p_value: float | None = Field(default=None, ge=0, le=1)
    #: Significance level the inferential claim is made at. Belongs in the report because
    #: a direction asserted at an alpha chosen after seeing the counts is not a result.
    alpha: float = Field(default=0.05, gt=0, lt=1)
    #: Which method the *test* supports. Only ever "a" or "b" when the paired test is
    #: significant at ``alpha``; otherwise "neither", however lopsided the counts look.
    favours: Literal["a", "b", "neither"] = "neither"
    #: Which way the discordant counts happen to lean. Descriptive only, and explicitly
    #: not a finding: with four discordant pairs the smallest attainable two-sided exact
    #: p-value is 0.125, so a 4-0 split can look decisive and prove nothing.
    observed_direction: Literal["a", "b", "neither"] = "neither"
    #: True when the discordant count is too small for any split to reach ``alpha``.
    #: Separates "we looked and found no difference" from "this test could never have
    #: found one", which a bare non-significant p-value cannot express.
    underpowered: bool = False
    #: Smallest two-sided exact p-value attainable at this discordant count.
    minimum_attainable_p_value: float | None = Field(default=None, ge=0, le=1)
    #: Distinct specimens contributing the discordant pairs. McNemar's exact test assumes
    #: each pair is an independent coin flip; pairs drawn from the same specimen are not,
    #: so a p-value computed over fewer specimens than pairs is too small.
    discordant_specimens: int = Field(default=0, ge=0)
    #: True when the discordant pairs come from fewer specimens than there are pairs.
    p_value_is_anticonservative: bool = False
    note: str
    research_only: Literal[True] = True

    @property
    def discordant(self) -> int:
        return self.only_a_detected + self.only_b_detected


def paired_detection_comparison(
    reports_a: Sequence[CnvEvaluationReport],
    reports_b: Sequence[CnvEvaluationReport],
    *,
    alpha: float = 0.05,
) -> PairedMethodComparison:
    """Compare two methods on the truth events both of them could assess.

    Events are paired by ``(sample_id, event_id)``. An event counts only when it is
    assessable under **both** methods: if one method could not look there, the pair
    carries no information about which method is better, and including it would let a
    method's blind spots influence the comparison.

    The p-value comes from McNemar's exact test on the discordant pairs. It is ``None``
    when no discordant pair exists, which means the comparison found no evidence either
    way rather than evidence of equivalence.

    ``favours`` names a method only when the test is significant at ``alpha``. The
    direction the counts happen to lean is reported separately as ``observed_direction``,
    because those are different claims and only one of them is a result: four discordant
    pairs split 4-0 look decisive and cannot reach any conventional threshold, since the
    smallest attainable two-sided p-value at that count is 0.125. ``underpowered`` marks
    exactly that situation, so a reader can tell a comparison that found no difference
    from one that could never have found one.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0, 1)")
    methods_a = {report.method for report in reports_a}
    methods_b = {report.method for report in reports_b}
    if len(methods_a) > 1 or len(methods_b) > 1:
        raise ValueError("each side of a paired comparison must contain exactly one method")
    if not methods_a or not methods_b:
        raise ValueError("both sides of a paired comparison must contain at least one report")

    def outcomes(reports: Sequence[CnvEvaluationReport]) -> dict[tuple[str, str], bool]:
        table: dict[tuple[str, str], bool] = {}
        for report in reports:
            for event in report.truth_events:
                if event.outcome in {"DETECTED", "MISSED"}:
                    table[(report.sample_id, event.event_id)] = event.outcome == "DETECTED"
        return table

    left = outcomes(reports_a)
    right = outcomes(reports_b)
    shared = sorted(set(left) & set(right))
    unpaired = len(set(left) ^ set(right))

    both = sum(1 for key in shared if left[key] and right[key])
    only_a = sum(1 for key in shared if left[key] and not right[key])
    only_b = sum(1 for key in shared if right[key] and not left[key])
    neither = sum(1 for key in shared if not left[key] and not right[key])
    p_value = mcnemar_exact(only_a, only_b)

    observed: Literal["a", "b", "neither"] = (
        "a" if only_a > only_b else "b" if only_b > only_a else "neither"
    )
    floor = minimum_attainable_p_value(only_a + only_b)
    underpowered = floor is None or floor > alpha
    discordant_specimens = len({key[0] for key in shared if left[key] != right[key]})
    clustered_pairs = discordant_specimens < (only_a + only_b)

    if not shared:
        note = (
            "No truth event was assessable under both methods, so the methods were not "
            "compared. This is not a tie."
        )
        favours: Literal["a", "b", "neither"] = "neither"
    elif p_value is None:
        note = (
            "The methods agreed on every paired event. With no discordant pair there is "
            "no evidence either way; this is not evidence of equivalence."
        )
        favours = "neither"
    elif underpowered:
        # The decisive case for reviewer trust: a 4-0 split reads as an obvious winner
        # and is not one, because no split of four pairs can reach 0.05. Naming a winner
        # here would be an inferential claim the data cannot support at any threshold.
        favours = "neither"
        note = (
            f"McNemar exact test on {only_a + only_b} discordant pair(s). The smallest "
            f"two-sided p-value attainable at this count is {floor:.4g}, above the "
            f"pre-specified alpha of {alpha:g}, so no observation could have reached "
            "significance. The comparison is underpowered by design, not inconclusive by "
            f"result; the counts lean towards {observed}, which is a description and not "
            "a finding."
        )
    elif p_value <= alpha:
        favours = observed
        note = (
            f"McNemar exact test on {only_a + only_b} discordant pair(s), p={p_value:.4g} "
            f"at a pre-specified alpha of {alpha:g}. The difference is significant in "
            f"favour of method {observed}."
        )
    else:
        favours = "neither"
        note = (
            f"McNemar exact test on {only_a + only_b} discordant pair(s), p={p_value:.4g} "
            f"at a pre-specified alpha of {alpha:g}. Not significant, so no method is "
            f"favoured; the counts lean towards {observed}, which is a description and "
            "not a finding. A non-significant result is not evidence of equivalence."
        )

    if clustered_pairs:
        note += (
            f" The {only_a + only_b} discordant pair(s) come from {discordant_specimens} "
            "specimen(s). McNemar's exact test treats each pair as an independent coin "
            "flip, and pairs from one specimen are not independent, so this p-value is "
            "smaller than the data support."
        )

    return PairedMethodComparison(
        method_a=methods_a.pop(),
        method_b=methods_b.pop(),
        paired_events=len(shared),
        both_detected=both,
        only_a_detected=only_a,
        only_b_detected=only_b,
        neither_detected=neither,
        unpaired_events=unpaired,
        p_value=p_value,
        alpha=alpha,
        favours=favours,
        observed_direction=observed,
        underpowered=underpowered,
        minimum_attainable_p_value=floor,
        discordant_specimens=discordant_specimens,
        p_value_is_anticonservative=clustered_pairs,
        note=note,
    )


def compare_aggregates(
    aggregates: Sequence[CnvAggregateReport],
) -> list[tuple[str, str, float | None, float | None]]:
    """Return a flat comparison table of ``(method, version, rate, lower_bound)``.

    Intentionally returns overlapping intervals rather than a ranking. With the cohort
    sizes available here, declaring a winner from point estimates whose intervals overlap
    would be a statistical error, so the caller is handed the intervals and must decide.
    """
    rows: list[tuple[str, str, float | None, float | None]] = []
    for item in aggregates:
        rows.append(
            (
                item.method,
                item.method_version,
                item.overall_detection_rate.point,
                item.overall_detection_rate.lower,
            )
        )
    return sorted(rows, key=lambda row: (row[3] is None, -(row[3] or 0.0)))

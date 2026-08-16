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


class CnvAggregateReport(StrictModel):
    """Aggregated benchmark evidence across many evaluations."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    aggregate_id: str = Field(min_length=1)
    method: str
    method_version: str
    evaluations: int = Field(ge=0)
    overall_detection_rate: ProportionResult
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
    unconfirmed = sum(
        sum(1 for e in r.query_events if e.outcome == "UNCONFIRMED") for r in reports
    )
    evaluable = sum(r.base_level.evaluable_bases for r in reports)
    weighted = (
        sum(
            (r.base_level.concordance or 0.0) * r.base_level.evaluable_bases for r in reports
        )
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
        by_data_basis=_group(
            reports, lambda r: r.data_basis.value, confidence_level
        ),
        limits_of_detection=limits,
        warnings=warnings,
        limitations=[
            "Aggregated synthetic or public-benchmark results do not establish "
            "intended-use performance on AML specimens.",
            "Detection rates pool events of different sizes and classes unless read "
            "through the per-stratum tables.",
            "A limit of detection derived from one truth profile does not generalise to "
            "other event classes or sizes.",
        ],
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

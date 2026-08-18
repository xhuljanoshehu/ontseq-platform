"""End-to-end synthetic CNV benchmark: simulate, call, evaluate, aggregate.

This module exists to make the CNV benchmark harness *executable evidence* rather than
documentation. It runs the entire path with no external binary, no reference genome and
no genomic data, which makes it safe for continuous integration and reproducible on any
machine.

It answers, for a deliberately simple baseline caller on deliberately simple synthetic
data, the questions the harness was built to answer: how detection behaves across blast
fraction and coverage, where specificity degrades before sensitivity does, and what a
limit of detection looks like when it is estimated honestly.

None of its numbers describe real assay performance.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ..io import write_json
from ..models import GenomeBuild, ModuleRunStatus
from .core import StateSegment
from .evaluate import evaluate_case
from .models import (
    CnvBenchmarkCase,
    CnvCallSet,
    CnvDataBasis,
    CnvEvaluationOptions,
    CnvEvaluationReport,
    CnvSegment,
    CnvStrata,
    CnvTruthSet,
    CnvTruthSource,
    GenomicRegion,
)
from .segment import (
    METHOD_NAME,
    METHOD_VERSION,
    DepthBin,
    SegmentationParameters,
    call_segments,
    neutral_background_segments,
)
from .simulate import SimulationParameters, closed_world_truth, simulate_bins
from .states import CopyNumberState
from .strata import (
    CnvAggregateReport,
    PairedMethodComparison,
    aggregate,
    paired_detection_comparison,
)

#: Synthetic genome. Contig lengths are GRCh38 values so that coordinates line up with
#: the other synthetic fixtures; a subset is used to keep the demo fast.
DEMO_CONTIGS: dict[str, int] = {
    "chr2": 242_193_529,
    "chr5": 181_538_259,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr17": 83_257_441,
    "chr20": 64_444_167,
}

#: A synthetic profile shaped like a complex AML karyotype: del(5q), -7q, +8, del(20q).
#: The coordinates are invented and carry no biological meaning.
DEMO_EVENTS: tuple[tuple[str, int, int, float], ...] = (
    ("chr5", 70_000_000, 160_000_000, 1.0),
    ("chr7", 65_000_000, 159_345_973, 1.0),
    ("chr8", 0, 145_138_636, 3.0),
    ("chr20", 39_000_000, 64_444_167, 1.0),
)

DEMO_TUMOR_FRACTIONS: tuple[float, ...] = (1.0, 0.6, 0.4, 0.25, 0.15, 0.1, 0.05)
DEMO_COVERAGES: tuple[float, ...] = (5.0, 3.0, 1.0)
DEMO_REPLICATES = 3


@dataclass(frozen=True)
class MethodVariant:
    """One caller configuration to benchmark.

    The demo runs two, on byte-identical simulated data, so that the paired method
    comparison is exercised on a real difference rather than only in unit tests. They
    differ solely in the segmentation split threshold, which makes the conservative
    variant less willing to declare a boundary.
    """

    key: str
    method: str
    parameters: SegmentationParameters
    rationale: str


DEMO_VARIANTS: tuple[MethodVariant, ...] = (
    MethodVariant(
        key="default",
        method=METHOD_NAME,
        parameters=SegmentationParameters(),
        rationale="Default split threshold of 4.0.",
    ),
    MethodVariant(
        key="conservative",
        method=f"{METHOD_NAME}-conservative",
        parameters=SegmentationParameters(split_threshold=8.0),
        rationale=(
            "Split threshold raised to 8.0, so a boundary must be twice as convincing. "
            "Expected to trade sensitivity for fewer spurious segments."
        ),
    ),
)


@dataclass(frozen=True)
class DemoOutputs:
    """Paths written by a demo run, plus the in-memory results."""

    truth_path: Path
    evaluation_paths: list[Path]
    aggregate_paths: list[Path]
    comparison_path: Path
    aggregates: list[CnvAggregateReport]
    comparison: PairedMethodComparison
    reports_by_variant: dict[str, list[CnvEvaluationReport]]

    @property
    def aggregate(self) -> CnvAggregateReport:
        """The default variant's aggregate, which the summary line describes."""
        return self.aggregates[0]


def demo_truth_segments() -> list[StateSegment]:
    """Return the gapless closed-world truth profile used by the demo."""
    return closed_world_truth(DEMO_EVENTS, DEMO_CONTIGS)


def build_demo_truth_set(sample_id: str = "SYNTHETIC_CNV_DEMO") -> CnvTruthSet:
    """Build the contract-level truth set for the demo profile."""
    segments = [
        CnvSegment(
            contig=contig,
            start=start,
            end=end,
            state=(CopyNumberState.LOSS if copy_number < 2.0 else CopyNumberState.GAIN),
            copy_number=copy_number,
            # Simulated truth knows its breakpoints exactly.
            start_uncertainty_bp=0,
            end_uncertainty_bp=0,
            notes=["Synthetic coordinates. No biological meaning."],
        )
        for contig, start, end, copy_number in DEMO_EVENTS
    ]
    return CnvTruthSet(
        truth_id="SYNTHETIC_CNV_DEMO_TRUTH",
        sample_id=sample_id,
        genome_build=GenomeBuild.GRCH38,
        source=CnvTruthSource.SIMULATED,
        source_version="ontseq-cnv-simulator-0.1.0",
        background_state=CopyNumberState.NEUTRAL,
        # The simulator asserts every base, but a bin-based method cannot resolve below
        # one bin, so the truth declares the bin size as its effective resolution.
        resolution_bp=1_000_000,
        segments=segments,
        baseline_ploidy=2.0,
        limitations=[
            "Fully synthetic. Models dosage and counting noise only.",
            "No GC bias, mappability, alignment error, or subclonal structure beyond a "
            "single mixture fraction.",
            "Establishes harness behaviour, never assay performance.",
        ],
    )


def _call_set_for(
    *,
    call_set_id: str,
    sample_id: str,
    tumor_fraction: float,
    mean_coverage_x: float,
    bin_size_bp: int,
    seed: int,
    variant: MethodVariant,
) -> tuple[CnvCallSet, list[GenomicRegion]]:
    parameters = SimulationParameters(
        bin_size_bp=bin_size_bp,
        mean_coverage_x=mean_coverage_x,
        tumor_fraction=tumor_fraction,
        seed=seed,
    )
    bins = simulate_bins(
        contig_lengths=DEMO_CONTIGS,
        truth_segments=demo_truth_segments(),
        parameters=parameters,
    )
    # The tumor fraction is deliberately NOT supplied to the caller. A benchmark that
    # hands the method the answer measures nothing, and inverting the mixture with a
    # known fraction amplifies noise by 1/f, which flatters or destroys the result
    # depending only on that one number.
    result = call_segments(
        [DepthBin(item.contig, item.start, item.end, item.count) for item in bins],
        variant.parameters,
    )
    filled = neutral_background_segments(result.segments, DEMO_CONTIGS)
    segments = [
        CnvSegment(
            contig=item.contig,
            start=item.start,
            end=item.end,
            state=item.state,
            copy_number=item.copy_number,
        )
        for item in filled
    ]
    no_call = [
        GenomicRegion(contig=contig, start=start, end=end, label="below_minimum_bin_count")
        for contig, start, end in result.no_call_regions
    ]
    call_set = CnvCallSet(
        call_set_id=call_set_id,
        sample_id=sample_id,
        genome_build=GenomeBuild.GRCH38,
        method=variant.method,
        method_version=METHOD_VERSION,
        data_basis=CnvDataBasis.SIMULATED,
        background_state=CopyNumberState.NEUTRAL,
        # The baseline segments the whole evaluable genome, so an empty result means it
        # looked everywhere and found nothing to report — a biological negative bounded by
        # the mask, not a failure to look. Reporting that as NO_CALL would be the exact
        # conflation the status vocabulary exists to prevent.
        status=ModuleRunStatus.COMPLETED,
        reports_biological_negative=not segments,
        segments=segments,
        no_call_regions=no_call,
        bin_size_bp=bin_size_bp,
        mean_coverage_x=mean_coverage_x,
        warnings=result.warnings,
    )
    return call_set, no_call


def run_demo_benchmark(
    *,
    tumor_fractions: Sequence[float] = DEMO_TUMOR_FRACTIONS,
    coverages: Sequence[float] = DEMO_COVERAGES,
    replicates: int = DEMO_REPLICATES,
    bin_size_bp: int = 1_000_000,
    seed: int = 20260816,
    variants: Sequence[MethodVariant] = DEMO_VARIANTS,
) -> tuple[CnvTruthSet, dict[str, list[CnvEvaluationReport]], list[CnvAggregateReport]]:
    """Run the full simulate-call-evaluate-aggregate loop for every variant.

    Each variant sees byte-identical simulated data for a given
    ``(coverage, fraction, replicate)`` cell, because the simulation seed depends only on
    that cell. That is what makes the downstream comparison genuinely paired.
    """
    truth = build_demo_truth_set()
    reports: dict[str, list[CnvEvaluationReport]] = {variant.key: [] for variant in variants}
    for coverage_index, coverage in enumerate(coverages):
        for fraction_index, fraction in enumerate(tumor_fractions):
            for replicate in range(replicates):
                derived_seed = (
                    seed + 100_000 * (coverage_index + 1) + 1_000 * (fraction_index + 1) + replicate
                )
                cell = f"{coverage:g}x-{fraction:g}-{replicate}".replace(".", "_")
                sample_id = f"SIM-{cell}"
                for variant in variants:
                    call_set, _ = _call_set_for(
                        call_set_id=f"CS-{variant.key}-{cell}",
                        sample_id=sample_id,
                        tumor_fraction=fraction,
                        mean_coverage_x=coverage,
                        bin_size_bp=bin_size_bp,
                        seed=derived_seed,
                        variant=variant,
                    )
                    case = CnvBenchmarkCase(
                        case_id=f"CASE-{variant.key}-{cell}",
                        genome_build=GenomeBuild.GRCH38,
                        contig_lengths=DEMO_CONTIGS,
                        truth=truth.model_copy(update={"sample_id": sample_id}),
                        call_set=call_set,
                        options=CnvEvaluationOptions(),
                        strata=CnvStrata(
                            assay_mode="simulated",
                            data_basis=CnvDataBasis.SIMULATED,
                            mean_coverage_x=coverage,
                            tumor_fraction=fraction,
                            bin_size_bp=bin_size_bp,
                            replicate=replicate,
                            sample_class="synthetic_positive",
                        ),
                    )
                    reports[variant.key].append(evaluate_case(case))
    aggregates = [
        aggregate(reports[variant.key], aggregate_id=f"SYNTHETIC_CNV_DEMO_{variant.key.upper()}")
        for variant in variants
    ]
    return truth, reports, aggregates


def write_demo_benchmark(
    output_dir: Path,
    *,
    tumor_fractions: Sequence[float] = DEMO_TUMOR_FRACTIONS,
    coverages: Sequence[float] = DEMO_COVERAGES,
    replicates: int = DEMO_REPLICATES,
    bin_size_bp: int = 1_000_000,
    seed: int = 20260816,
    variants: Sequence[MethodVariant] = DEMO_VARIANTS,
) -> DemoOutputs:
    """Run the demo and write every artifact as JSON.

    Requires at least two variants, because the paired comparison this writes is the
    point of the demo. A single-variant run should call :func:`run_demo_benchmark`.
    """
    if len(variants) < 2:
        raise ValueError("write_demo_benchmark needs at least two variants to compare")
    truth, reports, aggregates = run_demo_benchmark(
        tumor_fractions=tumor_fractions,
        coverages=coverages,
        replicates=replicates,
        bin_size_bp=bin_size_bp,
        seed=seed,
        variants=variants,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    truth_path = write_json(truth, output_dir / "cnv-demo.truth.json")

    evaluation_paths: list[Path] = []
    aggregate_paths: list[Path] = []
    for variant, summary in zip(variants, aggregates, strict=True):
        variant_dir = output_dir / "evaluations" / variant.key
        variant_dir.mkdir(parents=True, exist_ok=True)
        evaluation_paths.extend(
            write_json(report, variant_dir / f"{report.evaluation_id}.json")
            for report in reports[variant.key]
        )
        aggregate_paths.append(
            write_json(summary, output_dir / f"cnv-demo.aggregate.{variant.key}.json")
        )

    comparison = paired_detection_comparison(reports[variants[0].key], reports[variants[1].key])
    comparison_path = write_json(comparison, output_dir / "cnv-demo.comparison.json")
    return DemoOutputs(
        truth_path=truth_path,
        evaluation_paths=evaluation_paths,
        aggregate_paths=aggregate_paths,
        comparison_path=comparison_path,
        aggregates=aggregates,
        comparison=comparison,
        reports_by_variant=reports,
    )


def summarize_comparison(comparison: PairedMethodComparison) -> list[str]:
    """Render a short human-readable summary of a paired method comparison."""
    p_value = "undefined" if comparison.p_value is None else f"{comparison.p_value:.4g}"
    floor = (
        "n/a"
        if comparison.minimum_attainable_p_value is None
        else f"{comparison.minimum_attainable_p_value:.4g}"
    )
    lines = [
        f"paired comparison: {comparison.method_a} vs {comparison.method_b}",
        f"  paired events: {comparison.paired_events} (unpaired {comparison.unpaired_events})",
        f"  both={comparison.both_detected} onlyA={comparison.only_a_detected} "
        f"onlyB={comparison.only_b_detected} neither={comparison.neither_detected}",
        f"  McNemar exact p={p_value} at alpha={comparison.alpha:g}, favours={comparison.favours}",
        f"  observed direction={comparison.observed_direction} (descriptive), "
        f"smallest attainable p={floor}",
    ]
    if comparison.underpowered:
        lines.append("  UNDERPOWERED: no split of this many pairs could have reached alpha")
    if comparison.p_value_is_anticonservative:
        lines.append(
            f"  CLUSTERED: discordant pairs come from {comparison.discordant_specimens} "
            "specimen(s); the p-value is smaller than the data support"
        )
    lines.append(f"  {comparison.note}")
    return lines


def summarize_demo(summary: CnvAggregateReport) -> list[str]:
    """Render a short human-readable summary of an aggregate report."""
    lines = [
        f"method: {summary.method} {summary.method_version}",
        f"evaluations: {summary.evaluations}",
    ]
    overall = summary.overall_detection_rate
    if overall.point is None:
        lines.append("overall detection: undefined (no assessable truth event)")
    else:
        lines.append(
            f"overall detection: {overall.successes}/{overall.total} = {overall.point:.3f} "
            f"[{overall.lower:.3f}, {overall.upper:.3f}]"
        )
    for stratum in summary.by_tumor_fraction:
        rate = stratum.detection_rate
        concordance = (
            f" base_concordance={stratum.weighted_base_concordance:.3f}"
            if stratum.weighted_base_concordance is not None
            else ""
        )
        point = "undefined" if rate.point is None else f"{rate.point:.3f}"
        lines.append(
            f"  tumor_fraction={stratum.label:>5}: detection={point} "
            f"({rate.successes}/{rate.total}) unconfirmed={stratum.unconfirmed_events}"
            f"{concordance}"
        )
    for limit in summary.limits_of_detection:
        empirical = "none" if limit.empirical_value is None else f"{limit.empirical_value:g}"
        modelled = (
            "withheld" if limit.model_based_value is None else f"{limit.model_based_value:.4f}"
        )
        lines.append(
            f"  LoD{int(limit.target_detection_rate * 100)} by {limit.predictor}: "
            f"empirical={empirical} model={modelled}"
        )
    return lines

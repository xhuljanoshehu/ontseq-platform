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
from .strata import CnvAggregateReport, aggregate

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
class DemoOutputs:
    """Paths written by a demo run, plus the in-memory aggregate."""

    truth_path: Path
    evaluation_paths: list[Path]
    aggregate_path: Path
    aggregate: CnvAggregateReport
    reports: list[CnvEvaluationReport]


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
    result = call_segments(
        [DepthBin(item.contig, item.start, item.end, item.count) for item in bins],
        # The tumor fraction is deliberately NOT supplied to the caller. A benchmark that
        # hands the method the answer measures nothing, and inverting the mixture with a
        # known fraction amplifies noise by 1/f, which flatters or destroys the result
        # depending only on that one number.
        SegmentationParameters(),
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
        method=METHOD_NAME,
        method_version=METHOD_VERSION,
        data_basis=CnvDataBasis.SIMULATED,
        background_state=CopyNumberState.NEUTRAL,
        status=ModuleRunStatus.COMPLETED if segments else ModuleRunStatus.NO_CALL,
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
) -> tuple[CnvTruthSet, list[CnvEvaluationReport], CnvAggregateReport]:
    """Run the full simulate-call-evaluate-aggregate loop over a factorial design."""
    truth = build_demo_truth_set()
    reports: list[CnvEvaluationReport] = []
    for coverage_index, coverage in enumerate(coverages):
        for fraction_index, fraction in enumerate(tumor_fractions):
            for replicate in range(replicates):
                derived_seed = (
                    seed + 100_000 * (coverage_index + 1) + 1_000 * (fraction_index + 1) + replicate
                )
                sample_id = f"SIM-{coverage:g}x-{fraction:g}-{replicate}".replace(".", "_")
                call_set, _ = _call_set_for(
                    call_set_id=f"CS-{coverage:g}x-{fraction:g}-{replicate}".replace(".", "_"),
                    sample_id=sample_id,
                    tumor_fraction=fraction,
                    mean_coverage_x=coverage,
                    bin_size_bp=bin_size_bp,
                    seed=derived_seed,
                )
                case = CnvBenchmarkCase(
                    case_id=f"CASE-{coverage:g}x-{fraction:g}-{replicate}".replace(".", "_"),
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
                reports.append(evaluate_case(case))
    return truth, reports, aggregate(reports, aggregate_id="SYNTHETIC_CNV_DEMO")


def write_demo_benchmark(output_dir: Path, **kwargs: object) -> DemoOutputs:
    """Run the demo and write every artifact as JSON."""
    truth, reports, summary = run_demo_benchmark(**kwargs)  # type: ignore[arg-type]
    output_dir.mkdir(parents=True, exist_ok=True)
    evaluations_dir = output_dir / "evaluations"
    evaluations_dir.mkdir(parents=True, exist_ok=True)
    truth_path = write_json(truth, output_dir / "cnv-demo.truth.json")
    evaluation_paths = [
        write_json(report, evaluations_dir / f"{report.evaluation_id}.json") for report in reports
    ]
    aggregate_path = write_json(summary, output_dir / "cnv-demo.aggregate.json")
    return DemoOutputs(
        truth_path=truth_path,
        evaluation_paths=evaluation_paths,
        aggregate_path=aggregate_path,
        aggregate=summary,
        reports=reports,
    )


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

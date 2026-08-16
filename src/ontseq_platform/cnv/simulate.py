"""Deterministic synthetic read-depth generation for CNV benchmarking.

Purpose
-------

Two things are impossible without a simulator. First, proving that the comparison core
behaves correctly requires inputs whose truth is known exactly, including breakpoints.
Second, the questions this project needs answered - how sensitivity falls with coverage
and with blast fraction, and where the limit of detection sits - require dilution and
coverage series that no real cohort provides on demand.

Model
-----

For a locus with tumor copy number ``CN_t`` in a specimen of tumor fraction ``f``, the
bulk sample presents ``CN_mix = f * CN_t + (1 - f) * CN_n``. Read counts in a bin are
proportional to ``CN_mix`` relative to the genome-wide mean, and are drawn from a
negative binomial (a gamma-Poisson mixture) rather than a Poisson, because sequencing
depth is consistently overdispersed relative to Poisson. The dispersion parameter is
explicit and recorded.

What this is not
----------------

The simulator models dosage and counting noise only. It does not model GC bias,
mappability, alignment error, chimeric reads, adaptive-sampling enrichment kinetics, or
subclonal structure beyond a single mixture fraction. It therefore establishes that a
method's *arithmetic* behaves as intended across coverage and purity; it cannot
establish real-world performance, and no result derived from it may be presented as
analytical validation.
"""

from __future__ import annotations

import math
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from .core import StateSegment
from .intervals import canonical_contig
from .states import CopyNumberState, expected_mixture_copy_number, state_from_copy_number


@dataclass(frozen=True)
class Bin:
    """One fixed-width counting window."""

    contig: str
    start: int
    end: int
    count: int
    expected_copy_number: float

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class SimulationParameters:
    """Everything that determines a simulated run, echoed into provenance."""

    bin_size_bp: int = 1_000_000
    mean_coverage_x: float = 3.0
    read_length_bp: int = 10_000
    tumor_fraction: float = 1.0
    normal_copy_number: float = 2.0
    #: Negative-binomial dispersion. Larger means closer to Poisson. Typical short-read
    #: WGS bin counts sit around 50-200; long-read low-coverage data is noisier.
    dispersion: float = 50.0
    seed: int = 0

    def __post_init__(self) -> None:
        if self.bin_size_bp <= 0:
            raise ValueError("bin size must be positive")
        if self.mean_coverage_x < 0:
            raise ValueError("mean coverage must not be negative")
        if self.read_length_bp <= 0:
            raise ValueError("read length must be positive")
        if not 0.0 <= self.tumor_fraction <= 1.0:
            raise ValueError("tumor fraction must lie between 0 and 1")
        if self.dispersion <= 0:
            raise ValueError("dispersion must be positive")

    @property
    def expected_reads_per_bin(self) -> float:
        """Reads per bin at neutral copy number, from coverage and read length."""
        return self.mean_coverage_x * self.bin_size_bp / self.read_length_bp


def _poisson(rng: random.Random, mean: float) -> int:
    """Draw a Poisson variate.

    Knuth's product method is exact and fast for small means. Above 30 it needs too many
    iterations, so a normal approximation with continuity correction is used instead; at
    that magnitude the relative error is well below the noise the simulator is modelling.
    """
    if mean <= 0:
        return 0
    if mean < 30:
        limit = math.exp(-mean)
        count = 0
        product = rng.random()
        while product > limit:
            count += 1
            product *= rng.random()
        return count
    value = rng.gauss(mean, math.sqrt(mean))
    return max(0, int(round(value)))


def _negative_binomial(rng: random.Random, mean: float, dispersion: float) -> int:
    """Draw an overdispersed count as a gamma-Poisson mixture."""
    if mean <= 0:
        return 0
    rate = rng.gammavariate(dispersion, mean / dispersion)
    return _poisson(rng, rate)


def truth_copy_number_at(
    segments: Sequence[StateSegment],
    contig: str,
    start: int,
    end: int,
    *,
    default_copy_number: float,
) -> float:
    """Return the length-weighted tumor copy number across a window."""
    target = canonical_contig(contig)
    total = end - start
    if total <= 0:
        raise ValueError("window must be non-empty")
    covered = 0
    accumulated = 0.0
    for segment in segments:
        if canonical_contig(segment.contig) != target or segment.copy_number is None:
            continue
        overlap = min(end, segment.end) - max(start, segment.start)
        if overlap > 0:
            covered += overlap
            accumulated += segment.copy_number * overlap
    accumulated += default_copy_number * (total - covered)
    return accumulated / total


def simulate_bins(
    *,
    contig_lengths: Mapping[str, int],
    truth_segments: Sequence[StateSegment],
    parameters: SimulationParameters,
) -> list[Bin]:
    """Generate deterministic bin counts for a genome under a truth copy-number profile.

    The genome-wide mean mixture copy number is used as the normalisation baseline, which
    reproduces a real and important property of read-depth CNV analysis: depth data are
    inherently *relative*. A sample where every chromosome is duplicated is
    indistinguishable from a diploid one on depth alone.
    """
    rng = random.Random(parameters.seed)
    windows: list[tuple[str, int, int, float]] = []
    for contig in sorted(contig_lengths, key=lambda name: (len(name), name)):
        length = contig_lengths[contig]
        position = 0
        while position < length:
            end = min(position + parameters.bin_size_bp, length)
            if end - position >= parameters.bin_size_bp // 2:
                tumor_cn = truth_copy_number_at(
                    truth_segments,
                    contig,
                    position,
                    end,
                    default_copy_number=parameters.normal_copy_number,
                )
                mixture = expected_mixture_copy_number(
                    tumor_cn,
                    tumor_fraction=parameters.tumor_fraction,
                    normal_copy_number=parameters.normal_copy_number,
                )
                windows.append((canonical_contig(contig), position, end, mixture))
            position = end

    if not windows:
        return []
    baseline = sum(item[3] for item in windows) / len(windows)
    if baseline <= 0:
        raise ValueError("the simulated genome has zero mean copy number")

    bins: list[Bin] = []
    for contig, start, end, mixture in windows:
        scale = (end - start) / parameters.bin_size_bp
        mean = parameters.expected_reads_per_bin * scale * (mixture / baseline)
        bins.append(
            Bin(
                contig=contig,
                start=start,
                end=end,
                count=_negative_binomial(rng, mean, parameters.dispersion),
                expected_copy_number=mixture,
            )
        )
    return bins


def truth_profile(
    events: Sequence[tuple[str, int, int, float]],
    *,
    baseline_ploidy: float = 2.0,
) -> list[StateSegment]:
    """Build truth segments from ``(contig, start, end, tumor_copy_number)`` tuples."""
    segments: list[StateSegment] = []
    for contig, start, end, copy_number in events:
        segments.append(
            StateSegment(
                contig=canonical_contig(contig),
                start=start,
                end=end,
                state=state_from_copy_number(copy_number, baseline_ploidy=baseline_ploidy),
                copy_number=copy_number,
            )
        )
    return segments


def closed_world_truth(
    events: Sequence[tuple[str, int, int, float]],
    contig_lengths: Mapping[str, int],
    *,
    baseline_ploidy: float = 2.0,
) -> list[StateSegment]:
    """Expand event tuples into a gapless genome-wide truth profile.

    Every base not covered by an event becomes an explicit neutral segment. This is what
    a closed-world truth source asserts, and making it explicit rather than implicit means
    the evaluator never has to guess.
    """
    altered = truth_profile(events, baseline_ploidy=baseline_ploidy)
    by_contig: dict[str, list[StateSegment]] = {}
    for segment in altered:
        by_contig.setdefault(segment.contig, []).append(segment)

    complete: list[StateSegment] = []
    for contig, length in contig_lengths.items():
        key = canonical_contig(contig)
        items = sorted(by_contig.get(key, []), key=lambda item: item.start)
        cursor = 0
        for segment in items:
            if segment.start > cursor:
                complete.append(
                    StateSegment(
                        key, cursor, segment.start, CopyNumberState.NEUTRAL, baseline_ploidy
                    )
                )
            complete.append(segment)
            cursor = segment.end
        if cursor < length:
            complete.append(
                StateSegment(key, cursor, length, CopyNumberState.NEUTRAL, baseline_ploidy)
            )
    return complete


@dataclass(frozen=True)
class DilutionLevel:
    """One rung of a dilution series."""

    tumor_fraction: float
    replicate: int
    bins: list[Bin]
    parameters: SimulationParameters


def simulate_dilution_series(
    *,
    contig_lengths: Mapping[str, int],
    truth_segments: Sequence[StateSegment],
    tumor_fractions: Sequence[float],
    replicates: int,
    base_parameters: SimulationParameters,
) -> list[DilutionLevel]:
    """Generate a reproducible blast-fraction dilution series.

    Each ``(fraction, replicate)`` pair gets a distinct derived seed so the series is
    deterministic as a whole while replicates remain independent draws.
    """
    if replicates < 1:
        raise ValueError("replicates must be at least 1")
    levels: list[DilutionLevel] = []
    for fraction_index, fraction in enumerate(tumor_fractions):
        for replicate in range(replicates):
            seed = base_parameters.seed + 1000 * (fraction_index + 1) + replicate
            parameters = SimulationParameters(
                bin_size_bp=base_parameters.bin_size_bp,
                mean_coverage_x=base_parameters.mean_coverage_x,
                read_length_bp=base_parameters.read_length_bp,
                tumor_fraction=fraction,
                normal_copy_number=base_parameters.normal_copy_number,
                dispersion=base_parameters.dispersion,
                seed=seed,
            )
            levels.append(
                DilutionLevel(
                    tumor_fraction=fraction,
                    replicate=replicate,
                    bins=simulate_bins(
                        contig_lengths=contig_lengths,
                        truth_segments=truth_segments,
                        parameters=parameters,
                    ),
                    parameters=parameters,
                )
            )
    return levels

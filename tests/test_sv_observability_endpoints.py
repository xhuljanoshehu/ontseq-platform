"""Regression tests for breakpoint-level Adaptive Sampling SV observability."""

from __future__ import annotations

import pytest

from ontseq_platform.models import (
    AssayMode,
    EventType,
    FileFingerprint,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    SvObservability,
    TargetBedRole,
    ToolRecord,
)
from ontseq_platform.sv_observability import apply_sv_observability
from ontseq_platform.target_coverage import (
    TargetCoveragePolicy,
    TargetCoverageRegion,
    TargetCoverageReport,
)


def _region(*, start: int, end: int, mean_depth: float, region_id: str) -> TargetCoverageRegion:
    length = end - start
    thresholds = (1, 10, 20, 30)
    bases = {f"{threshold}x": length if mean_depth >= threshold else 0 for threshold in thresholds}
    return TargetCoverageRegion(
        chromosome="chr1",
        start=start,
        end=end,
        region_id=region_id,
        mean_depth=mean_depth,
        bases_at_threshold=bases,
        fraction_at_threshold={key: value / length for key, value in bases.items()},
    )


def _coverage(*regions: TargetCoverageRegion) -> TargetCoverageReport:
    return TargetCoverageReport(
        sample_id="SYNTHETIC_SV_OBS",
        genome_build=GenomeBuild.GRCH38,
        target_bed_version="synthetic-endpoints-v1",
        target_bed_role=TargetBedRole.ANALYSIS_ROI_UNBUFFERED,
        status=ModuleRunStatus.COMPLETED,
        policy=TargetCoveragePolicy(
            profile_id="synthetic-endpoints",
            status="technical_defaults_only",
            note="Synthetic regression-only coverage policy.",
        ),
        summary_metrics={
            "region_count": len(regions),
            "interval_bases": sum(region.end - region.start for region in regions),
        },
        regions=list(regions),
        target_bed_fingerprint=FileFingerprint(size_bytes=1, sha256="0" * 64),
        tool=ToolRecord(name="mosdepth", version="0.3.14"),
    )


def _span_event(event_type: EventType) -> GenomicEvent:
    return GenomicEvent(
        event_id=f"SYNTHETIC-{event_type.value}",
        event_type=event_type,
        primary=Locus(chromosome="chr1", start=100, end=900),
        length_bp=800,
    )


@pytest.mark.parametrize(
    "event_type",
    [EventType.DELETION, EventType.DUPLICATION, EventType.INVERSION],
)
def test_span_interior_overlap_does_not_count_as_breakpoint_observability(
    event_type: EventType,
) -> None:
    observed = apply_sv_observability(
        [_span_event(event_type)],
        assay_mode=AssayMode.ADAPTIVE_SAMPLING,
        coverage_report=_coverage(_region(start=400, end=600, mean_depth=30, region_id="INTERIOR")),
        minimum_mean_depth=10,
    )[0]

    assert observed.observability == SvObservability.OUTSIDE_TARGET
    assert observed.breakpoint_mean_depths == [None, None]
    assert observed.reportable is False


def test_span_breakpoints_use_start_and_final_covered_base() -> None:
    observed = apply_sv_observability(
        [_span_event(EventType.DELETION)],
        assay_mode=AssayMode.ADAPTIVE_SAMPLING,
        coverage_report=_coverage(
            _region(start=100, end=101, mean_depth=20, region_id="LEFT_ENDPOINT"),
            _region(start=899, end=900, mean_depth=30, region_id="RIGHT_ENDPOINT"),
        ),
        minimum_mean_depth=10,
    )[0]

    assert observed.observability == SvObservability.OBSERVED_ADEQUATELY
    assert observed.breakpoint_mean_depths == [20.0, 30.0]


def test_one_simple_sv_endpoint_inside_is_partial_not_adequate() -> None:
    observed = apply_sv_observability(
        [_span_event(EventType.INVERSION)],
        assay_mode=AssayMode.ADAPTIVE_SAMPLING,
        coverage_report=_coverage(
            _region(start=100, end=101, mean_depth=20, region_id="LEFT_ONLY")
        ),
        minimum_mean_depth=10,
    )[0]

    assert observed.observability == SvObservability.PARTIALLY_OBSERVED
    assert observed.breakpoint_mean_depths == [20.0, None]


def test_insertion_remains_a_single_anchor_locus() -> None:
    event = GenomicEvent(
        event_id="SYNTHETIC-insertion",
        event_type=EventType.INSERTION,
        primary=Locus(chromosome="chr1", start=500, end=501),
        length_bp=100,
    )
    observed = apply_sv_observability(
        [event],
        assay_mode=AssayMode.ADAPTIVE_SAMPLING,
        coverage_report=_coverage(
            _region(start=490, end=510, mean_depth=20, region_id="INS_ANCHOR")
        ),
        minimum_mean_depth=10,
    )[0]

    assert observed.observability == SvObservability.OBSERVED_ADEQUATELY
    assert observed.breakpoint_mean_depths == [20.0]

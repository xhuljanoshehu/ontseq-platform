from __future__ import annotations

from .models import AssayMode, GenomicEvent, Locus, SvObservability, TargetBedRole
from .target_coverage import TargetCoverageRegion, TargetCoverageReport


def _matching_regions(
    event_locus: Locus, regions: list[TargetCoverageRegion]
) -> list[TargetCoverageRegion]:
    canonical = event_locus.chromosome.removeprefix("chr")
    return [
        region
        for region in regions
        if region.chromosome.removeprefix("chr") == canonical
        and region.start < event_locus.end
        and event_locus.start < region.end
    ]


def apply_sv_observability(
    events: list[GenomicEvent],
    *,
    assay_mode: AssayMode,
    coverage_report: TargetCoverageReport | None,
    minimum_mean_depth: float,
) -> list[GenomicEvent]:
    """Describe breakpoint observability using an explicit unvalidated technical depth floor."""
    if minimum_mean_depth < 0:
        raise ValueError("minimum_mean_depth must be non-negative")
    if assay_mode != AssayMode.ADAPTIVE_SAMPLING:
        return [
            event.model_copy(update={"observability": SvObservability.NOT_APPLICABLE})
            for event in events
        ]
    if coverage_report is None:
        raise ValueError("Adaptive Sampling SV observability requires a coverage report")

    result: list[GenomicEvent] = []
    for event in events:
        loci = [event.primary, *([event.secondary] if event.secondary is not None else [])]
        hits = [_matching_regions(locus, coverage_report.regions) for locus in loci]
        inside = [bool(items) for items in hits]
        mean_depths: list[float | None] = [
            max(item.mean_depth for item in items) if items else None for items in hits
        ]
        adequate = [depth is not None and depth >= minimum_mean_depth for depth in mean_depths]
        if not any(inside):
            status = SvObservability.OUTSIDE_TARGET
        elif not all(inside):
            status = SvObservability.PARTIALLY_OBSERVED
        elif not all(adequate):
            status = SvObservability.INSUFFICIENT_COVERAGE
        elif coverage_report.target_bed_role == TargetBedRole.SELECTION_PANEL_BUFFERED:
            status = SvObservability.PARTIALLY_OBSERVED
        else:
            status = SvObservability.OBSERVED_ADEQUATELY
        notes = list(event.notes)
        notes.append(
            f"Adaptive Sampling observability={status.value} using an unvalidated technical "
            f"mean-depth floor of {minimum_mean_depth:g}x and target role "
            f"{coverage_report.target_bed_role.value}."
        )
        result.append(
            event.model_copy(
                update={
                    "observability": status,
                    "breakpoint_mean_depths": mean_depths,
                    "observability_target_role": coverage_report.target_bed_role,
                    "notes": notes,
                    "reportable": False,
                }
            )
        )
    return result

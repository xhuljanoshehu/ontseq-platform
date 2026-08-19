from __future__ import annotations

from .cutesv import CuteSVCallReport, cutesv_observations
from .models import ModuleRunStatus, SnifflesCallReport
from .sv_concordance import (
    SVCallerObservation,
    SVConcordancePolicy,
    SVConcordanceReport,
    compare_sv_caller_observations,
)


def sniffles_observations(report: SnifflesCallReport) -> list[SVCallerObservation]:
    observations: list[SVCallerObservation] = []
    for event in report.events:
        if len(event.evidence) != 1:
            raise ValueError("Sniffles2 normalized event must contain exactly one evidence record")
        observations.append(
            SVCallerObservation(
                observation_id=event.event_id,
                caller="Sniffles2",
                caller_version=report.tool.version,
                source_event_id=event.event_id,
                event_type=event.event_type,
                primary=event.primary,
                secondary=event.secondary,
                evidence=event.evidence[0],
            )
        )
    return observations


def compare_sniffles_and_cutesv(
    sniffles: SnifflesCallReport,
    cutesv: CuteSVCallReport,
    policy: SVConcordancePolicy,
) -> SVConcordanceReport:
    """Compare normalized caller evidence without promoting agreement to biological truth."""

    if sniffles.sample_id != cutesv.sample_id:
        raise ValueError("Sniffles2 and cuteSV reports must refer to the same sample")
    if sniffles.genome_build != cutesv.genome_build:
        raise ValueError("Sniffles2 and cuteSV reports must use the same genome build")

    report = compare_sv_caller_observations(
        sniffles_observations(sniffles),
        cutesv_observations(cutesv),
        policy,
    )
    warnings = list(report.warnings)
    if sniffles.status == ModuleRunStatus.NO_CALL or cutesv.status == ModuleRunStatus.NO_CALL:
        warnings.append(
            "A caller-level NO_CALL is not a validated negative result and must not be used as "
            "evidence of absence."
        )
    return report.model_copy(update={"warnings": list(dict.fromkeys(warnings))})

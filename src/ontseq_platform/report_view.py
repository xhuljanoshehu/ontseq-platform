from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .models import (
    Evidence,
    EventAnnotation,
    GenomicEvent,
    ModuleOutcome,
    ModuleRunStatus,
    PipelineResult,
    Verdict,
)


AlertLevel = Literal["critical", "warning", "info"]


@dataclass(frozen=True)
class ReviewAlert:
    level: AlertLevel
    title: str
    detail: str


@dataclass(frozen=True)
class ModuleView:
    name: str
    status: ModuleRunStatus
    reason: str
    meaning: str
    css_class: str


@dataclass(frozen=True)
class EvidenceView:
    caller: str
    caller_version: str
    support_reads: int | None
    local_coverage: float | None
    variant_allele_fraction: float | None
    quality: float | None
    filters: tuple[str, ...]
    supporting_read_strands: str | None
    precise: bool | None


@dataclass(frozen=True)
class AnnotationView:
    source_id: str
    source_release: str
    record_id: str
    assertion: str
    assertion_vocabulary: str
    record_origin: str
    scope_alignment: str
    scope_note: str
    caveats: tuple[str, ...]


@dataclass(frozen=True)
class EventView:
    event_id: str
    event_type: str
    length_bp: int | None
    copy_number: float | None
    primary_locus: str
    secondary_locus: str | None
    cytobands: str | None
    genes: tuple[str, ...]
    confidence: str
    reportable: bool
    reportability_text: str
    evidence: tuple[EvidenceView, ...]
    notes: tuple[str, ...]
    annotations: tuple[AnnotationView, ...]


@dataclass(frozen=True)
class ReportView:
    sample_id: str
    run_id: str
    assay_mode: str
    genome_build: str
    reference_id: str
    analysis_profile: str
    analysis_intent: str
    qc_verdict: str
    release_status: str
    pipeline_version: str
    git_commit: str
    created_at: str
    target_bed_version: str | None
    modules: tuple[ModuleView, ...]
    alerts: tuple[ReviewAlert, ...]
    events: tuple[EventView, ...]
    qc_metrics: tuple[tuple[str, float | int | str | None], ...]
    qc_failed_gates: tuple[str, ...]
    warnings: tuple[str, ...]
    reference_checksums: tuple[tuple[str, str], ...]


_STATUS_MEANING = {
    ModuleRunStatus.COMPLETED: (
        "Analysis completed. Any biological interpretation remains bounded by the module "
        "contract, observability and validation status."
    ),
    ModuleRunStatus.NO_CALL: (
        "Analysis ran but did not produce an interpretable call. This is not a biological "
        "negative result."
    ),
    ModuleRunStatus.FAILED: (
        "Execution was attempted and failed. Downstream absence of findings must not be "
        "interpreted biologically."
    ),
    ModuleRunStatus.NOT_RUN: (
        "The module did not run or was not applicable. This is not a negative result."
    ),
}

_STATUS_CLASS = {
    ModuleRunStatus.COMPLETED: "state-completed",
    ModuleRunStatus.NO_CALL: "state-no-call",
    ModuleRunStatus.FAILED: "state-failed",
    ModuleRunStatus.NOT_RUN: "state-not-run",
}


def _locus_text(event: GenomicEvent, secondary: bool = False) -> str | None:
    locus = event.secondary if secondary else event.primary
    if locus is None:
        return None
    return f"{locus.chromosome}:{locus.start}-{locus.end}"


def _cytoband_text(event: GenomicEvent) -> str | None:
    parts: list[str] = []
    for locus in (event.primary, event.secondary):
        if locus is None or locus.cytoband_start is None:
            continue
        band = locus.cytoband_start
        if locus.cytoband_end and locus.cytoband_end != locus.cytoband_start:
            band = f"{band}-{locus.cytoband_end}"
        parts.append(f"{locus.chromosome} {band}")
    return "; ".join(parts) if parts else None


def _evidence_view(item: Evidence) -> EvidenceView:
    return EvidenceView(
        caller=item.caller,
        caller_version=item.caller_version,
        support_reads=item.support_reads,
        local_coverage=item.local_coverage,
        variant_allele_fraction=item.variant_allele_fraction,
        quality=item.quality,
        filters=tuple(item.filters),
        supporting_read_strands=item.supporting_read_strands,
        precise=item.precise,
    )


def _annotation_view(item: EventAnnotation) -> AnnotationView:
    return AnnotationView(
        source_id=item.source_id,
        source_release=item.source_release,
        record_id=item.record_id,
        assertion=item.assertion,
        assertion_vocabulary=item.assertion_vocabulary,
        record_origin=item.record_origin,
        scope_alignment=item.scope_alignment,
        scope_note=item.scope_note,
        caveats=tuple(item.caveats),
    )


def _event_view(event: GenomicEvent) -> EventView:
    reportability_text = "false"
    if event.reportable:
        reportability_text = (
            "true — pipeline flag only; this RUO report is not clinically validated"
        )
    return EventView(
        event_id=event.event_id,
        event_type=event.event_type.value,
        length_bp=event.length_bp,
        copy_number=event.copy_number,
        primary_locus=_locus_text(event) or "not available",
        secondary_locus=_locus_text(event, secondary=True),
        cytobands=_cytoband_text(event),
        genes=tuple(event.genes),
        confidence=event.confidence,
        reportable=event.reportable,
        reportability_text=reportability_text,
        evidence=tuple(_evidence_view(item) for item in event.evidence),
        notes=tuple(event.notes),
        annotations=tuple(_annotation_view(item) for item in event.annotations),
    )


def _module_view(module: ModuleOutcome) -> ModuleView:
    return ModuleView(
        name=module.module.value,
        status=module.status,
        reason=module.reason,
        meaning=_STATUS_MEANING[module.status],
        css_class=_STATUS_CLASS[module.status],
    )


def _deduplicate(items: list[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(item for item in items if item))


def _alerts(result: PipelineResult) -> tuple[ReviewAlert, ...]:
    alerts: list[ReviewAlert] = []
    if result.qc.verdict == Verdict.FAIL:
        alerts.append(
            ReviewAlert(
                level="critical",
                title="QC failed",
                detail=(
                    "The result is technically blocked. Genomic absence must not be read as a "
                    "negative finding."
                ),
            )
        )
    elif result.qc.verdict == Verdict.WARN:
        alerts.append(
            ReviewAlert(
                level="warning",
                title="QC warning",
                detail="Review QC warnings and failed gates before interpreting genomic findings.",
            )
        )

    for module in result.modules:
        if module.status == ModuleRunStatus.FAILED:
            alerts.append(
                ReviewAlert(
                    level="critical",
                    title=f"{module.module.value}: FAILED",
                    detail=module.reason or _STATUS_MEANING[module.status],
                )
            )
        elif module.status == ModuleRunStatus.NO_CALL:
            alerts.append(
                ReviewAlert(
                    level="warning",
                    title=f"{module.module.value}: NO_CALL",
                    detail=(
                        f"{module.reason} {_STATUS_MEANING[module.status]}".strip()
                        if module.reason
                        else _STATUS_MEANING[module.status]
                    ),
                )
            )

    if not result.events:
        alerts.append(
            ReviewAlert(
                level="info",
                title="No normalized events in the result contract",
                detail=(
                    "This statement describes the result object only. Review module status and "
                    "observability before making any biological inference."
                ),
            )
        )
    return tuple(alerts)


def build_report_view(result: PipelineResult) -> ReportView:
    warnings = _deduplicate(result.warnings + result.qc.warnings + result.iscn.warnings)
    analysis_intent = (
        result.manifest.analysis.intent.value
        if result.manifest.analysis.intent is not None
        else "not declared"
    )
    return ReportView(
        sample_id=result.manifest.sample_id,
        run_id=result.manifest.run_id,
        assay_mode=result.manifest.assay.mode.value,
        genome_build=result.manifest.assay.genome_build.value,
        reference_id=result.manifest.assay.reference_id,
        analysis_profile=result.manifest.analysis.profile,
        analysis_intent=analysis_intent,
        qc_verdict=result.qc.verdict.value,
        release_status=result.release_status.value,
        pipeline_version=result.provenance.pipeline_version,
        git_commit=result.provenance.git_commit,
        created_at=result.provenance.created_at.isoformat(),
        target_bed_version=result.manifest.assay.target_bed_version,
        modules=tuple(_module_view(item) for item in result.modules),
        alerts=_alerts(result),
        events=tuple(_event_view(item) for item in result.events),
        qc_metrics=tuple(result.qc.metrics.items()),
        qc_failed_gates=tuple(result.qc.failed_gates),
        warnings=warnings,
        reference_checksums=tuple(sorted(result.provenance.reference_checksums.items())),
    )

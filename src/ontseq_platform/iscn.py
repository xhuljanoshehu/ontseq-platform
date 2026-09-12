from __future__ import annotations

from .iscn_syntax import render_supported_cnv_fragment
from .models import (
    AnalysisModule,
    EventType,
    GenomeBuild,
    GenomicEvent,
    ISCNAssessmentBlocker,
    ISCNDispositionOutcome,
    ISCNEventDisposition,
    ISCNEventFragment,
    ISCNProposal,
    ISCNProposalStatus,
    ISCNResourceProvenance,
    ISCNSelectionPolicy,
    ModuleOutcome,
    ModuleRunStatus,
    ResolvedResourceContext,
    SvObservability,
)

ISCN_RULE_PROFILE = "event-fragments-v0.2-unvalidated"

_TECHNICAL_CNV_TYPES = {
    EventType.CHROMOSOME_GAIN,
    EventType.CHROMOSOME_LOSS,
    EventType.DELETION,
    EventType.DUPLICATION,
}
_UNOBSERVABLE = {
    SvObservability.PARTIALLY_OBSERVED,
    SvObservability.INSUFFICIENT_COVERAGE,
    SvObservability.OUTSIDE_TARGET,
}


def _chromosome_number(chromosome: str) -> int:
    value = chromosome.removeprefix("chr")
    return {"X": 23, "Y": 24}.get(value, int(value) if value.isdigit() else 99)


def _fragment(event: GenomicEvent) -> tuple[str | None, str | None]:
    fragment = render_supported_cnv_fragment(
        event_type=event.event_type.value,
        chromosome=event.primary.chromosome,
        copy_number=event.copy_number,
        cytoband_start=event.primary.cytoband_start,
        cytoband_end=event.primary.cytoband_end,
        whole_chromosome_span_confirmed=event.whole_chromosome_span_confirmed,
    )
    if fragment is not None:
        return fragment, None
    if event.event_type in {EventType.DELETION, EventType.DUPLICATION}:
        return None, "no usable named cytoband range"
    if event.event_type in {EventType.CHROMOSOME_GAIN, EventType.CHROMOSOME_LOSS}:
        return None, "an exact complete chromosome span was not confirmed"
    return None, "event construct is outside the implemented ISCN subset"


def _selection_reason(
    event: GenomicEvent,
    selection_policy: ISCNSelectionPolicy,
    cnv_source_event_ids: set[str],
) -> tuple[str | None, str | None]:
    if event.observability in _UNOBSERVABLE:
        return None, f"observability is {event.observability.value}"
    chromosomes = {
        event.primary.chromosome.removeprefix("chr"),
        *([event.secondary.chromosome.removeprefix("chr")] if event.secondary is not None else []),
    }
    if chromosomes.intersection({"X", "Y"}):
        return None, "sex-chromosome copy number and complement were not assessed"
    if event.event_type not in _TECHNICAL_CNV_TYPES:
        return (
            None,
            "automatic SV-to-ISCN conversion is disabled until reciprocal/balanced and derivative "
            "structure semantics are validated",
        )
    if event.copy_number is None:
        return None, "typed CNV event lacks an absolute copy-number estimate"
    if event.event_id not in cnv_source_event_ids:
        return None, "event did not originate in the typed CNV call report"
    if selection_policy == ISCNSelectionPolicy.REPORTABLE_ONLY_V1:
        if event.reportable:
            return "REPORTABLE_CNV_EVENT", "CNV event explicitly marked reportable"
        return None, "CNV event is not marked reportable"
    return "TECHNICAL_CNV_CALL", "technical CNV call with an absolute copy-number estimate"


def resource_provenance_for_iscn(
    context: ResolvedResourceContext | None,
) -> ISCNResourceProvenance | None:
    """Extract only checksum-pinned resources required by the proposal contract."""

    if context is None:
        return None
    reference_lock_key = "reference.reference_lock"
    annotation_cache_key = "reference.annotation_cache"
    cytoband_key = "reference.cytobands"
    if reference_lock_key not in context.resource_checksums:
        return None
    if cytoband_key not in context.resource_checksums:
        return None
    if annotation_cache_key not in context.resource_checksums:
        return None
    return ISCNResourceProvenance(
        genome_build=context.genome_build,
        reference_dictionary_contract=context.reference_dictionary_contract,
        reference_bundle_id=context.reference_bundle_id,
        reference_bundle_version=context.reference_bundle_version,
        reference_lock_sha256=context.resource_checksums[reference_lock_key],
        annotation_cache_sha256=context.resource_checksums[annotation_cache_key],
        cytoband_release=context.resource_releases.get(cytoband_key, "unspecified"),
        cytoband_sha256=context.resource_checksums[cytoband_key],
    )


def _event_sort_key(event: GenomicEvent) -> tuple[int, int, int, int, str, str]:
    secondary_chromosome = (
        _chromosome_number(event.secondary.chromosome) if event.secondary is not None else 99
    )
    secondary_start = event.secondary.start if event.secondary is not None else 0
    return (
        _chromosome_number(event.primary.chromosome),
        event.primary.start,
        secondary_chromosome,
        secondary_start,
        event.event_type.value,
        event.event_id,
    )


def build_iscn_proposal(
    events: list[GenomicEvent],
    *,
    genome_build: GenomeBuild,
    selection_policy: ISCNSelectionPolicy = ISCNSelectionPolicy.REPORTABLE_ONLY_V1,
    requested: bool = True,
    upstream_evidence_assessed: bool = True,
    resource_provenance: ISCNResourceProvenance | None = None,
    policy_parameters: dict[str, str | int | float | bool] | None = None,
    technical_assumptions: list[str] | None = None,
    assessment_blockers: list[ISCNAssessmentBlocker] | None = None,
    cnv_source_event_ids: set[str] | None = None,
) -> ISCNProposal:
    """Render a build-bound, evidence-linked and deliberately partial proposal.

    The routine never infers chromosome count, sex-chromosome complement, clonality,
    phase, derivative structure or a normal result. It returns event fragments only.
    """

    ordered = sorted(events, key=_event_sort_key)
    considered_event_ids = [event.event_id for event in ordered]
    typed_cnv_ids = set(cnv_source_event_ids or set())
    unknown_cnv_ids = typed_cnv_ids.difference(considered_event_ids)
    if unknown_cnv_ids:
        raise ValueError(
            "CNV source IDs are absent from the proposal event set: "
            + ", ".join(sorted(unknown_cnv_ids))
        )
    common_warnings = [
        "Technical ISCN proposal only; not clinically validated; expert cytogenetic review "
        "is mandatory.",
        "The renderer implements only a small event-level subset and does not establish a "
        "complete karyotype.",
        "Chromosome count, sex-chromosome complement, clonality, phase, derivative structure "
        "and normality were not inferred.",
    ]
    common = {
        "conformance_profile": ISCN_RULE_PROFILE,
        "genome_build": genome_build,
        "selection_policy": selection_policy,
        "resource_provenance": resource_provenance,
        "policy_parameters": policy_parameters or {},
        "technical_assumptions": technical_assumptions or [],
        "considered_event_ids": considered_event_ids,
        "cnv_source_event_ids": sorted(typed_cnv_ids),
    }
    if not requested:
        return ISCNProposal(
            proposal_status=ISCNProposalStatus.NOT_REQUESTED,
            event_dispositions=[
                ISCNEventDisposition(
                    event_id=event.event_id,
                    outcome=ISCNDispositionOutcome.NOT_REQUESTED,
                    reason_code="MODULE_NOT_REQUESTED",
                    detail="the ISCN module was not requested in the analysis manifest",
                )
                for event in ordered
            ],
            warnings=[*common_warnings, "ISCN proposal generation was not requested."],
            **common,
        )
    blockers = list(assessment_blockers or [])
    if not upstream_evidence_assessed:
        blockers.append(
            ISCNAssessmentBlocker(
                reason_code="UPSTREAM_EVIDENCE_NOT_ASSESSED",
                detail="no completed eligible CNV evidence assessment was available",
            )
        )
    if resource_provenance is None and not any(
        item.reason_code.startswith("RESOURCE_") for item in blockers
    ):
        blockers.append(
            ISCNAssessmentBlocker(
                reason_code="RESOURCE_PROVENANCE_UNAVAILABLE",
                detail="checksum-pinned reference-lock/cytoband provenance was unavailable",
            )
        )
    if blockers:
        return ISCNProposal(
            proposal_status=ISCNProposalStatus.NOT_ASSESSED,
            event_dispositions=[
                ISCNEventDisposition(
                    event_id=event.event_id,
                    outcome=ISCNDispositionOutcome.ASSESSMENT_BLOCKED,
                    reason_code="PROPOSAL_NOT_ASSESSED",
                    detail="; ".join(item.detail for item in blockers),
                )
                for event in ordered
            ],
            assessment_blockers=blockers,
            warnings=[
                *common_warnings,
                "An unassessed proposal is not a biological negative and does not establish "
                "a normal karyotype.",
                *(
                    f"ISCN assessment blocker {item.reason_code}: {item.detail}."
                    for item in blockers
                ),
            ],
            **common,
        )
    fragments: list[ISCNEventFragment] = []
    dispositions: list[ISCNEventDisposition] = []
    warnings = list(common_warnings)
    seen_fragments: set[str] = set()

    for event in ordered:
        reason_code, selection_detail = _selection_reason(event, selection_policy, typed_cnv_ids)
        if reason_code is None:
            dispositions.append(
                ISCNEventDisposition(
                    event_id=event.event_id,
                    outcome=ISCNDispositionOutcome.EXCLUDED_BY_POLICY,
                    reason_code="EXCLUDED_BY_SELECTION_POLICY",
                    detail=selection_detail or "event did not satisfy the selection policy",
                )
            )
            continue
        fragment, omission_detail = _fragment(event)
        if fragment is None:
            dispositions.append(
                ISCNEventDisposition(
                    event_id=event.event_id,
                    outcome=ISCNDispositionOutcome.OMITTED_UNSUPPORTED,
                    reason_code="UNSUPPORTED_OR_INCOMPLETE_CONSTRUCT",
                    detail=omission_detail or "event could not be rendered",
                )
            )
            warnings.append(f"Event {event.event_id} was omitted: {omission_detail}.")
            continue
        if fragment in seen_fragments:
            dispositions.append(
                ISCNEventDisposition(
                    event_id=event.event_id,
                    outcome=ISCNDispositionOutcome.OMITTED_UNSUPPORTED,
                    reason_code="DUPLICATE_FRAGMENT",
                    detail=f"fragment {fragment!r} is already represented by another event",
                )
            )
            warnings.append(
                f"Event {event.event_id} duplicated proposed fragment {fragment!r} and was omitted."
            )
            continue
        seen_fragments.add(fragment)
        fragments.append(
            ISCNEventFragment(
                event_id=event.event_id,
                event_type=event.event_type,
                fragment=fragment,
                confidence=event.confidence,
                reportable=event.reportable,
                selection_reason=selection_detail or reason_code,
            )
        )
        dispositions.append(
            ISCNEventDisposition(
                event_id=event.event_id,
                outcome=ISCNDispositionOutcome.RENDERED,
                reason_code=reason_code,
                detail=selection_detail or reason_code,
            )
        )

    if fragments:
        proposal_status = ISCNProposalStatus.PARTIAL_EVENT_LEVEL
        notation = ",".join(item.fragment for item in fragments)
    else:
        proposal_status = ISCNProposalStatus.NO_RENDERABLE_CANDIDATE
        notation = None
        warnings.append(
            "No supported ISCN fragment was generated; this is not evidence of a normal karyotype."
        )

    if selection_policy == ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1:
        warnings.append(
            "Technical candidate selection does not change reportable=false, analytical "
            "validation, or clinical release status."
        )

    return ISCNProposal(
        notation=notation,
        conformance_profile=ISCN_RULE_PROFILE,
        proposal_status=proposal_status,
        genome_build=genome_build,
        selection_policy=selection_policy,
        resource_provenance=resource_provenance,
        policy_parameters=policy_parameters or {},
        event_fragments=fragments,
        event_dispositions=dispositions,
        considered_event_ids=considered_event_ids,
        cnv_source_event_ids=sorted(typed_cnv_ids),
        source_event_ids=[item.event_id for item in fragments],
        technical_assumptions=technical_assumptions or [],
        warnings=warnings,
    )


def iscn_module_outcome(proposal: ISCNProposal) -> ModuleOutcome:
    """Map proposal completeness to module execution without implying a biological verdict."""

    if proposal.proposal_status == ISCNProposalStatus.PARTIAL_EVENT_LEVEL:
        return ModuleOutcome(
            module=AnalysisModule.ISCN,
            status=ModuleRunStatus.COMPLETED,
            reason=(
                f"Generated {len(proposal.event_fragments)} traceable, unvalidated ISCN event "
                "fragment(s); this is not a complete or clinically releasable karyotype."
            ),
        )
    if proposal.proposal_status == ISCNProposalStatus.NO_RENDERABLE_CANDIDATE:
        return ModuleOutcome(
            module=AnalysisModule.ISCN,
            status=ModuleRunStatus.NO_CALL,
            reason=(
                "ISCN proposal evaluation completed, but no supported event fragment was "
                "generated; this is not a biological negative or a normal karyotype."
            ),
        )
    reason = {
        ISCNProposalStatus.NOT_REQUESTED: "ISCN proposal generation was not requested.",
        ISCNProposalStatus.NOT_ASSESSED: (
            "ISCN proposal generation was not assessed because required upstream evidence or "
            "checksum-pinned cytoband provenance was unavailable."
        ),
        ISCNProposalStatus.LEGACY_UNSPECIFIED: (
            "Legacy result does not record structured ISCN proposal execution state."
        ),
    }[proposal.proposal_status]
    return ModuleOutcome(
        module=AnalysisModule.ISCN,
        status=ModuleRunStatus.NOT_RUN,
        reason=reason,
    )

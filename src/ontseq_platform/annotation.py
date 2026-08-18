"""Wiring for ``ontseq annotate``: a locked ClinVar release applied to a result contract.

The matching and the scope logic live in :mod:`ontseq_platform.knowledge`, which carries no
dependencies and is unit tested on its own. This module knows only how to turn a result
contract's events into intervals, hand them over, and put the answers back — without
touching ``confidence`` or ``reportable``, which is the whole point.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .knowledge.annotate import (
    DEFAULT_EXACT_TOLERANCE_BP,
    DEFAULT_MINIMUM_OVERLAP,
    AnnotationSource,
    annotate_finding,
)
from .knowledge.clinvar import ClinVarRecord, load, release_sha256
from .knowledge.scope import Interval, Origin
from .models import (
    AnalysisIntent,
    EventAnnotation,
    GenomeBuild,
    KnowledgeResourceLock,
    PipelineResult,
)

CLINVAR_VOCABULARY = "acmg_germline"

#: Manifest intent mapped onto the origin vocabulary. ``BOTH`` resolves to ``UNKNOWN``
#: rather than to either: an analysis asking both questions cannot settle which one a
#: particular assertion answers, and claiming alignment would be the guess this avoids.
_INTENT: dict[AnalysisIntent, Origin] = {
    AnalysisIntent.SOMATIC: Origin.SOMATIC,
    AnalysisIntent.GERMLINE: Origin.GERMLINE,
    AnalysisIntent.BOTH: Origin.UNKNOWN,
}


@dataclass(frozen=True)
class AnnotationOutcome:
    """An annotated result and the resource it was annotated against."""

    result: PipelineResult
    lock: KnowledgeResourceLock
    events_annotated: int
    annotations_added: int


def assay_origin(result: PipelineResult) -> Origin:
    """What the manifest says this analysis is looking for, or ``UNKNOWN`` when it is silent."""
    intent = result.manifest.analysis.intent
    return Origin.UNKNOWN if intent is None else _INTENT[intent]


def _interval(chromosome: str, start: int, end: int) -> Interval | None:
    """Build an interval, or ``None`` for a locus that cannot carry one."""
    try:
        return Interval(contig=chromosome, start=start, end=end)
    except ValueError:
        return None


def annotate_result(
    result: PipelineResult,
    records: list[ClinVarRecord],
    *,
    lock: KnowledgeResourceLock,
    minimum_reciprocal_overlap: float = DEFAULT_MINIMUM_OVERLAP,
    exact_tolerance_bp: int = DEFAULT_EXACT_TOLERANCE_BP,
) -> AnnotationOutcome:
    """Attach matching records to every event, leaving every judgement field untouched.

    Both loci of a paired event are annotated: a translocation's partner region is where the
    knowledge base is most likely to have something to say, and annotating only the primary
    would silently halve the coverage for exactly the events that matter most.
    """
    source = AnnotationSource(
        source_id=lock.source_id,
        release=lock.release,
        sha256=lock.sha256,
        vocabulary=lock.assertion_vocabulary,
    )
    intent = assay_origin(result)

    annotated_events = []
    events_annotated = 0
    annotations_added = 0
    for event in result.events:
        found: list[EventAnnotation] = []
        seen: set[tuple[str, str]] = set()
        for locus in (event.primary, event.secondary):
            if locus is None:
                continue
            interval = _interval(locus.chromosome, locus.start, locus.end)
            if interval is None:
                continue
            for annotation in annotate_finding(
                interval,
                records,
                source=source,
                assay_intent=intent,
                minimum_reciprocal_overlap=minimum_reciprocal_overlap,
                exact_tolerance_bp=exact_tolerance_bp,
            ):
                key = (annotation.source_id, annotation.record_id)
                if key in seen:
                    # A record spanning both loci of a paired event is one piece of
                    # evidence, not two.
                    continue
                seen.add(key)
                found.append(
                    EventAnnotation(
                        source_id=annotation.source_id,
                        source_release=annotation.source_release,
                        source_sha256=annotation.source_sha256,
                        record_id=annotation.record_id,
                        record_type=annotation.record_type,
                        assertion=annotation.assertion,
                        assertion_vocabulary=annotation.assertion_vocabulary,
                        record_origin=annotation.record_origin.value,
                        scope_alignment=annotation.scope_alignment.value,
                        scope_note=annotation.scope_note,
                        match_type=annotation.match_type.value,
                        reciprocal_overlap=annotation.reciprocal_overlap,
                        review_status=annotation.review_status,
                        review_stars=annotation.review_stars,
                        genes=list(annotation.genes),
                        conditions=list(annotation.conditions),
                        caveats=list(annotation.caveats),
                    )
                )
        if found:
            events_annotated += 1
            annotations_added += len(found)
        # Only `annotations` changes. `confidence` and `reportable` are copied through
        # untouched, and this is the line that has to stay that way.
        annotated_events.append(event.model_copy(update={"annotations": found}))

    updated = result.model_copy(update={"events": annotated_events})
    return AnnotationOutcome(
        result=updated,
        lock=lock,
        events_annotated=events_annotated,
        annotations_added=annotations_added,
    )


def load_clinvar(
    path: Path, *, genome_build: GenomeBuild, release: str
) -> tuple[list[ClinVarRecord], KnowledgeResourceLock]:
    """Read and lock a ClinVar release for one genome build."""
    records, summary = load(path, assembly=genome_build.value)
    lock = KnowledgeResourceLock(
        source_id="clinvar",
        release=release,
        sha256=release_sha256(path),
        genome_build=genome_build,
        assertion_vocabulary=CLINVAR_VOCABULARY,
        records_loaded=len(records),
        load_summary=summary.describe(),
    )
    return records, lock


def describe(outcome: AnnotationOutcome) -> list[str]:
    """Lines a person can read without opening the JSON."""
    lines = [
        f"{outcome.lock.source_id} {outcome.lock.release} "
        f"({outcome.lock.genome_build.value}, {outcome.lock.sha256[:16]}…)",
        f"  {outcome.lock.load_summary}",
        f"  {outcome.annotations_added} annotation(s) on "
        f"{outcome.events_annotated}/{len(outcome.result.events)} event(s)",
    ]
    if outcome.result.manifest.analysis.intent is None:
        lines.append(
            "  the manifest does not declare somatic or germline intent, so no annotation's "
            "scope could be checked against the question this assay was asking"
        )
    lines.append(
        "  annotations are classifications of database records; none of them makes a "
        "finding reportable"
    )
    return lines

"""Attaching knowledge-base evidence to a finding without turning it into a verdict.

The whole point of this module is what it refuses to do. It matches findings to records,
records how well they matched and whether the record even concerns the kind of variation the
assay was asking about — and then stops. It does not set ``reportable``, it does not raise
``confidence``, and it does not translate ClinVar's germline vocabulary into a statement
about a somatic finding.

That restraint is the design, not a limitation of it. Deciding that a finding is reportable
in AML needs somatic criteria — ELN, ICC, a locally agreed gene list — and this repository
has none of them. Code that promoted a ClinVar *Pathogenic* into a reportable finding would
be inventing a clinical rule nobody agreed to, in the one place where that is least
recoverable: inside a report a physician is about to sign.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from .clinvar import ClinVarRecord
from .scope import (
    Interval,
    MatchType,
    Origin,
    ScopeAlignment,
    align,
    canonical_contig,
    classify_match,
)

#: Default reciprocal overlap for calling a region record a match. An engineering default,
#: pre-registerable like every other comparison threshold in this repository, and not a
#: validated concordance criterion.
DEFAULT_MINIMUM_OVERLAP = 0.5
#: Breakpoints from a read-depth caller and from a clinical array are never bit-identical.
DEFAULT_EXACT_TOLERANCE_BP = 10_000


@dataclass(frozen=True)
class Annotation:
    """One knowledge-base record attached to one finding, with everything needed to judge it.

    Deliberately carries no field that could be mistaken for a conclusion about the sample.
    ``assertion`` is the source's own words; ``assertion_vocabulary`` names the rule set
    those words belong to, so *Pathogenic* cannot be read as a somatic driver claim.
    """

    source_id: str
    source_release: str
    source_sha256: str
    record_id: str
    record_type: str
    assertion: str
    #: The classification system ``assertion`` is expressed in. ClinVar's is ACMG germline.
    assertion_vocabulary: str
    record_origin: Origin
    scope_alignment: ScopeAlignment
    scope_note: str
    match_type: MatchType
    reciprocal_overlap: float
    review_status: str
    review_stars: int | None
    genes: tuple[str, ...]
    conditions: tuple[str, ...]

    @property
    def caveats(self) -> tuple[str, ...]:
        """Everything a reader must know before using this annotation."""
        notes = [self.scope_note] if self.scope_alignment is not ScopeAlignment.ALIGNED else []
        if self.review_stars is not None and self.review_stars <= 1:
            notes.append(
                f"the submitting evidence is weak ({self.review_stars}★: {self.review_status}); "
                "a single-submitter or conflicting assertion is not a settled classification"
            )
        if self.review_stars is None:
            notes.append(
                f"review status {self.review_status!r} is not in this code's mapping, so the "
                "assertion's evidential weight could not be rated"
            )
        if self.match_type is not MatchType.EXACT:
            notes.append(
                f"matched by {self.match_type.value.replace('_', ' ')} at "
                f"{self.reciprocal_overlap:.2f} reciprocal overlap, not by identical "
                "coordinates; the record may describe a different event in the same region"
            )
        notes.append(
            "this is a classification of a database record, not a finding about this "
            "sample, and it does not make anything reportable"
        )
        return tuple(notes)


@dataclass(frozen=True)
class AnnotationSource:
    """The locked release an annotation came from."""

    source_id: str
    release: str
    sha256: str
    vocabulary: str = "acmg_germline"


def annotate_finding(
    finding: Interval,
    records: Sequence[ClinVarRecord],
    *,
    source: AnnotationSource,
    assay_intent: Origin,
    minimum_reciprocal_overlap: float = DEFAULT_MINIMUM_OVERLAP,
    exact_tolerance_bp: int = DEFAULT_EXACT_TOLERANCE_BP,
) -> list[Annotation]:
    """Attach every record that matches this finding, strongest match first.

    Records whose origin does not match the assay's question are **kept**, not filtered.
    A germline pathogenic deletion underlying a somatic finding is a secondary finding a
    reviewer needs to see; dropping it would be a clinical decision disguised as a filter.
    It is returned with its mismatch stated instead.
    """
    annotations: list[Annotation] = []
    wanted = canonical_contig(finding.contig)
    for record in records:
        if canonical_contig(record.interval.contig) != wanted:
            continue
        matched = classify_match(
            finding,
            record.interval,
            minimum_reciprocal_overlap=minimum_reciprocal_overlap,
            exact_tolerance_bp=exact_tolerance_bp,
        )
        if matched is None:
            continue
        match_type, overlap = matched
        alignment, note = align(record.origin, assay_intent)
        annotations.append(
            Annotation(
                source_id=source.source_id,
                source_release=source.release,
                source_sha256=source.sha256,
                record_id=record.variation_id,
                record_type=record.record_type,
                assertion=record.assertion,
                assertion_vocabulary=source.vocabulary,
                record_origin=record.origin,
                scope_alignment=alignment,
                scope_note=note,
                match_type=match_type,
                reciprocal_overlap=overlap,
                review_status=record.review_status,
                review_stars=record.stars,
                genes=record.genes,
                conditions=record.conditions,
            )
        )

    order = {
        MatchType.EXACT: 0,
        MatchType.FINDING_WITHIN_RECORD: 1,
        MatchType.RECORD_WITHIN_FINDING: 2,
        MatchType.OVERLAP: 3,
    }
    annotations.sort(
        key=lambda item: (
            order[item.match_type],
            -(item.review_stars if item.review_stars is not None else -1),
            -item.reciprocal_overlap,
            item.record_id,
        )
    )
    return annotations


def summarize(annotations: Sequence[Annotation]) -> str:
    """One line a reviewer can read without opening the JSON."""
    if not annotations:
        return "no knowledge-base record matched this finding"
    aligned = sum(1 for item in annotations if item.scope_alignment is ScopeAlignment.ALIGNED)
    mismatched = sum(1 for item in annotations if item.scope_alignment is ScopeAlignment.MISMATCHED)
    unknown = len(annotations) - aligned - mismatched
    parts = [f"{len(annotations)} record(s) matched"]
    if aligned:
        parts.append(f"{aligned} of matching origin")
    if mismatched:
        parts.append(f"{mismatched} of a different origin (secondary findings)")
    if unknown:
        parts.append(f"{unknown} whose origin could not be checked")
    return "; ".join(parts) + ". None of this makes a finding reportable."

from __future__ import annotations

from collections import defaultdict
from enum import StrEnum
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .models import (
    Evidence,
    GenomeBuild,
    GenomicEvent,
    ModuleRunStatus,
    SnifflesCallReport,
    StrictModel,
)


class ObservabilityStatus(StrEnum):
    OBSERVABLE = "observable"
    LIMITED = "limited"
    NOT_ASSESSABLE = "not_assessable"
    UNKNOWN = "unknown"


class FusionClassification(StrEnum):
    GENE_GENE = "gene_gene"
    GENE_INTERGENIC = "gene_intergenic"
    INTERGENIC_INTERGENIC = "intergenic_intergenic"
    UNRESOLVED = "unresolved"


class GeneFeature(StrictModel):
    chromosome: str = Field(pattern=r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    gene: str = Field(min_length=1)
    strand: Literal["+", "-"]
    transcript_ids: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def end_after_start(self) -> GeneFeature:
        if self.end <= self.start:
            raise ValueError("gene feature end must be greater than start")
        return self


class ObservabilityRegion(StrictModel):
    chromosome: str = Field(pattern=r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    status: ObservabilityStatus
    reason: str = Field(min_length=1)

    @model_validator(mode="after")
    def end_after_start(self) -> ObservabilityRegion:
        if self.end <= self.start:
            raise ValueError("observability region end must be greater than start")
        return self


class BreakpointGeneHit(StrictModel):
    gene: str
    strand: Literal["+", "-"]
    transcript_ids: list[str] = Field(default_factory=list)
    distance_bp: int = Field(ge=0)


class FusionBreakpoint(StrictModel):
    chromosome: str
    position: int = Field(ge=0)
    genes: list[BreakpointGeneHit] = Field(default_factory=list)
    observability: ObservabilityStatus = ObservabilityStatus.UNKNOWN
    observability_reason: str | None = None
    breakend_orientation: Literal["left", "right", "unknown"] = "unknown"


class FusionGenePair(StrictModel):
    gene_5prime: str | None = None
    gene_3prime: str | None = None
    orientation_resolved: bool = False
    known_pair: bool = False


class FusionCandidate(StrictModel):
    candidate_id: str
    source_event_id: str
    primary: FusionBreakpoint
    secondary: FusionBreakpoint
    classification: FusionClassification
    gene_pairs: list[FusionGenePair] = Field(default_factory=list)
    evidence: list[Evidence] = Field(default_factory=list)
    confidence: Literal["unclassified"] = "unclassified"
    reportable: Literal[False] = False
    research_only: Literal[True] = True
    limitations: list[str] = Field(default_factory=list)


class FusionInterpretationReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    status: ModuleRunStatus
    annotation_resource_id: str
    annotation_resource_version: str
    candidates: list[FusionCandidate]
    source_translocation_count: int = Field(ge=0)
    unresolved_source_event_ids: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def counts_are_consistent(self) -> FusionInterpretationReport:
        if self.source_translocation_count != len(self.candidates) + len(
            self.unresolved_source_event_ids
        ):
            raise ValueError("fusion source-event accounting is inconsistent")
        return self


class GeneAnnotationIndex:
    """Minimal build-locked gene interval index for breakpoint annotation.

    The parser intentionally consumes a simple BED6/BED7 contract instead of a presentation
    file. BED columns: chrom, start, end, gene, score, strand, optional transcript_id.
    Multiple transcript rows for the same gene interval are collapsed at annotation time.
    """

    def __init__(
        self,
        features: list[GeneFeature],
        *,
        resource_id: str,
        resource_version: str,
    ) -> None:
        self.resource_id = resource_id
        self.resource_version = resource_version
        by_chromosome: dict[str, list[GeneFeature]] = defaultdict(list)
        for feature in features:
            by_chromosome[_canonical_chromosome(feature.chromosome)].append(feature)
        self._features = {
            chromosome: sorted(items, key=lambda item: (item.start, item.end, item.gene))
            for chromosome, items in by_chromosome.items()
        }

    @classmethod
    def from_bed(
        cls,
        path: Path,
        *,
        resource_id: str,
        resource_version: str,
    ) -> GeneAnnotationIndex:
        features: list[GeneFeature] = []
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line or line.startswith("#") or line.startswith("track "):
                    continue
                fields = line.split("\t")
                if len(fields) < 6:
                    raise ValueError(
                        f"gene annotation BED line {line_number} requires at least 6 columns"
                    )
                chromosome, raw_start, raw_end, gene, _score, strand = fields[:6]
                transcript_id = fields[6].strip() if len(fields) >= 7 else ""
                try:
                    start = int(raw_start)
                    end = int(raw_end)
                except ValueError as exc:
                    raise ValueError(
                        f"gene annotation BED line {line_number} has invalid coordinates"
                    ) from exc
                features.append(
                    GeneFeature(
                        chromosome=chromosome,
                        start=start,
                        end=end,
                        gene=gene,
                        strand=strand,
                        transcript_ids=[transcript_id] if transcript_id else [],
                    )
                )
        if not features:
            raise ValueError("gene annotation BED contains no usable features")
        return cls(features, resource_id=resource_id, resource_version=resource_version)

    def annotate(self, chromosome: str, position: int, *, flank_bp: int = 0) -> list[BreakpointGeneHit]:
        if flank_bp < 0:
            raise ValueError("flank_bp must be non-negative")
        canonical = _canonical_chromosome(chromosome)
        grouped: dict[tuple[str, str], dict[str, object]] = {}
        for feature in self._features.get(canonical, []):
            if feature.start - flank_bp <= position < feature.end + flank_bp:
                distance = _distance_to_interval(position, feature.start, feature.end)
                key = (feature.gene, feature.strand)
                entry = grouped.setdefault(
                    key,
                    {"distance": distance, "transcripts": set()},
                )
                entry["distance"] = min(int(entry["distance"]), distance)
                transcripts = entry["transcripts"]
                assert isinstance(transcripts, set)
                transcripts.update(feature.transcript_ids)
        hits = [
            BreakpointGeneHit(
                gene=gene,
                strand=strand,
                transcript_ids=sorted(str(item) for item in entry["transcripts"] if item),
                distance_bp=int(entry["distance"]),
            )
            for (gene, strand), entry in grouped.items()
        ]
        return sorted(hits, key=lambda item: (item.distance_bp, item.gene, item.strand))


def _canonical_chromosome(chromosome: str) -> str:
    return chromosome[3:] if chromosome.startswith("chr") else chromosome


def _distance_to_interval(position: int, start: int, end: int) -> int:
    if start <= position < end:
        return 0
    if position < start:
        return start - position
    return position - end + 1


def _observability_at(
    chromosome: str,
    position: int,
    regions: list[ObservabilityRegion],
) -> tuple[ObservabilityStatus, str | None]:
    canonical = _canonical_chromosome(chromosome)
    matching = [
        region
        for region in regions
        if _canonical_chromosome(region.chromosome) == canonical
        and region.start <= position < region.end
    ]
    if not matching:
        return ObservabilityStatus.UNKNOWN, None
    priority = {
        ObservabilityStatus.NOT_ASSESSABLE: 3,
        ObservabilityStatus.LIMITED: 2,
        ObservabilityStatus.OBSERVABLE: 1,
        ObservabilityStatus.UNKNOWN: 0,
    }
    selected = max(matching, key=lambda item: priority[item.status])
    return selected.status, selected.reason


def _classification(primary: FusionBreakpoint, secondary: FusionBreakpoint) -> FusionClassification:
    if primary.genes and secondary.genes:
        return FusionClassification.GENE_GENE
    if primary.genes or secondary.genes:
        return FusionClassification.GENE_INTERGENIC
    if not primary.genes and not secondary.genes:
        return FusionClassification.INTERGENIC_INTERGENIC
    return FusionClassification.UNRESOLVED


def _canonical_pair(gene_a: str, gene_b: str) -> tuple[str, str]:
    return tuple(sorted((gene_a.upper(), gene_b.upper())))  # type: ignore[return-value]


def _gene_pairs(
    primary: FusionBreakpoint,
    secondary: FusionBreakpoint,
    known_pairs: set[tuple[str, str]],
) -> list[FusionGenePair]:
    pairs: list[FusionGenePair] = []
    for first in primary.genes:
        for second in secondary.genes:
            canonical = _canonical_pair(first.gene, second.gene)
            pairs.append(
                FusionGenePair(
                    gene_5prime=first.gene,
                    gene_3prime=second.gene,
                    orientation_resolved=False,
                    known_pair=canonical in known_pairs,
                )
            )
    return pairs


def _breakpoint(
    chromosome: str,
    position: int,
    annotation: GeneAnnotationIndex,
    observability: list[ObservabilityRegion],
    *,
    flank_bp: int,
) -> FusionBreakpoint:
    status, reason = _observability_at(chromosome, position, observability)
    return FusionBreakpoint(
        chromosome=chromosome,
        position=position,
        genes=annotation.annotate(chromosome, position, flank_bp=flank_bp),
        observability=status,
        observability_reason=reason,
        breakend_orientation="unknown",
    )


def candidate_from_event(
    event: GenomicEvent,
    annotation: GeneAnnotationIndex,
    *,
    observability: list[ObservabilityRegion] | None = None,
    known_pairs: set[tuple[str, str]] | None = None,
    flank_bp: int = 0,
) -> FusionCandidate:
    if event.event_type.value != "translocation" or event.secondary is None:
        raise ValueError("fusion interpretation requires a paired translocation/BND event")
    regions = observability or []
    normalized_known_pairs = {
        _canonical_pair(first, second) for first, second in (known_pairs or set())
    }
    primary = _breakpoint(
        event.primary.chromosome,
        event.primary.start,
        annotation,
        regions,
        flank_bp=flank_bp,
    )
    secondary = _breakpoint(
        event.secondary.chromosome,
        event.secondary.start,
        annotation,
        regions,
        flank_bp=flank_bp,
    )
    limitations = [
        "DNA breakend evidence is not equivalent to an expressed or functional fusion transcript.",
        "Breakend orientation is unresolved because the current normalized Sniffles2 event contract does not retain VCF ALT breakend orientation.",
        "Gene overlap alone does not establish transcript compatibility, reading frame or oncogenic relevance.",
        "Candidate remains research-only and non-reportable until assay-specific analytical validation and expert review.",
    ]
    if primary.observability != ObservabilityStatus.OBSERVABLE or secondary.observability != ObservabilityStatus.OBSERVABLE:
        limitations.append(
            "One or both breakpoints lack confirmed observable status; absence of additional evidence cannot be interpreted as a biological negative."
        )
    return FusionCandidate(
        candidate_id=f"FUSION-{event.event_id}",
        source_event_id=event.event_id,
        primary=primary,
        secondary=secondary,
        classification=_classification(primary, secondary),
        gene_pairs=_gene_pairs(primary, secondary, normalized_known_pairs),
        evidence=event.evidence,
        limitations=limitations,
    )


def interpret_sniffles_fusions(
    report: SnifflesCallReport,
    annotation: GeneAnnotationIndex,
    *,
    observability: list[ObservabilityRegion] | None = None,
    known_pairs: set[tuple[str, str]] | None = None,
    flank_bp: int = 0,
) -> FusionInterpretationReport:
    source_events = [event for event in report.events if event.event_type.value == "translocation"]
    candidates: list[FusionCandidate] = []
    unresolved: list[str] = []
    for event in source_events:
        try:
            candidates.append(
                candidate_from_event(
                    event,
                    annotation,
                    observability=observability,
                    known_pairs=known_pairs,
                    flank_bp=flank_bp,
                )
            )
        except ValueError:
            unresolved.append(event.event_id)
    warnings = [
        "Fusion candidates are derived from DNA structural-variant evidence and require independent biological interpretation.",
        "No candidate is clinically reportable under this research-only profile.",
    ]
    if not source_events:
        warnings.append(
            "NO_CALL means no normalized translocation/BND source event was available; it is not a validated negative fusion result."
        )
    return FusionInterpretationReport(
        sample_id=report.sample_id,
        genome_build=report.genome_build,
        status=ModuleRunStatus.COMPLETED if candidates else ModuleRunStatus.NO_CALL,
        annotation_resource_id=annotation.resource_id,
        annotation_resource_version=annotation.resource_version,
        candidates=candidates,
        source_translocation_count=len(source_events),
        unresolved_source_event_ids=unresolved,
        warnings=warnings,
    )

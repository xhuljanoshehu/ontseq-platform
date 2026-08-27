from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .models import GenomeBuild, GenomicEvent, IntervalResourceLock, Locus
from .reference import sha256_file


@dataclass(frozen=True)
class IntervalAnnotation:
    chromosome: str
    start: int
    end: int
    label: str


def _canonical_chromosome(value: str) -> str:
    return value if value.startswith("chr") else f"chr{value}"


def load_interval_resource(
    path: Path, lock: IntervalResourceLock
) -> dict[str, tuple[IntervalAnnotation, ...]]:
    """Load a locked BED-like four-column annotation without external interval dependencies."""
    if not path.is_file():
        raise ValueError(f"annotation resource is missing: {path}")
    observed = sha256_file(path)
    if observed != lock.sha256:
        raise ValueError(
            f"{lock.resource_id} checksum mismatch: expected {lock.sha256}, observed {observed}"
        )
    records: dict[str, list[IntervalAnnotation]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"{lock.resource_id} line {line_number}: expected chrom/start/end/label"
                )
            chromosome = _canonical_chromosome(fields[0])
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{lock.resource_id} line {line_number}: non-integer interval"
                ) from exc
            if start < 0 or end <= start or not fields[3].strip():
                raise ValueError(f"{lock.resource_id} line {line_number}: invalid interval")
            records.setdefault(chromosome, []).append(
                IntervalAnnotation(chromosome, start, end, fields[3].strip())
            )
    return {
        chromosome: tuple(sorted(items, key=lambda item: (item.start, item.end, item.label)))
        for chromosome, items in records.items()
    }


def _overlaps(
    locus: Locus, resource: dict[str, tuple[IntervalAnnotation, ...]]
) -> list[IntervalAnnotation]:
    chromosome = _canonical_chromosome(locus.chromosome)
    return [
        item
        for item in resource.get(chromosome, ())
        if item.start < locus.end and locus.start < item.end
    ]


def _nearest(
    locus: Locus, resource: dict[str, tuple[IntervalAnnotation, ...]]
) -> tuple[IntervalAnnotation, int] | None:
    chromosome = _canonical_chromosome(locus.chromosome)
    candidates = resource.get(chromosome, ())
    if not candidates:
        return None
    nearest = min(
        candidates,
        key=lambda item: min(abs(locus.start - item.end), abs(item.start - locus.end)),
    )
    distance = max(0, min(abs(locus.start - nearest.end), abs(nearest.start - locus.end)))
    return nearest, distance


def _annotate_locus(
    locus: Locus,
    *,
    genes: dict[str, tuple[IntervalAnnotation, ...]] | None,
    cytobands: dict[str, tuple[IntervalAnnotation, ...]] | None,
) -> tuple[Locus, list[str], str | None]:
    gene_hits = _overlaps(locus, genes) if genes is not None else []
    band_hits = _overlaps(locus, cytobands) if cytobands is not None else []
    gene_labels = sorted({item.label for item in gene_hits})
    nearest_note: str | None = None
    if genes is not None and not gene_labels:
        nearest = _nearest(locus, genes)
        if nearest is not None:
            item, distance = nearest
            nearest_note = f"nearest_gene={item.label}; distance_bp={distance}"
    band_labels = sorted({item.label for item in band_hits})
    annotated = locus.model_copy(
        update={
            "gene": gene_labels[0] if gene_labels else None,
            "cytoband_start": band_labels[0] if band_labels else None,
            "cytoband_end": band_labels[-1] if band_labels else None,
        }
    )
    return annotated, gene_labels, nearest_note


def annotate_sv_events(
    events: list[GenomicEvent],
    *,
    genome_build: GenomeBuild,
    gene_resource: tuple[Path, IntervalResourceLock] | None = None,
    cytoband_resource: tuple[Path, IntervalResourceLock] | None = None,
    context_resources: list[tuple[Path, IntervalResourceLock]] | None = None,
) -> list[GenomicEvent]:
    """Attach build-locked genes, cytobands, and artifact-context flags to both breakpoints."""
    resources = [item for item in [gene_resource, cytoband_resource] if item is not None]
    resources.extend(context_resources or [])
    for _path, lock in resources:
        if lock.genome_build != genome_build:
            raise ValueError(
                f"Refusing {lock.resource_id}: {lock.genome_build.value} does not match "
                f"event build {genome_build.value}"
            )
    genes = load_interval_resource(*gene_resource) if gene_resource else None
    cytobands = load_interval_resource(*cytoband_resource) if cytoband_resource else None
    contexts = [
        (lock.resource_type, load_interval_resource(path, lock))
        for path, lock in (context_resources or [])
    ]
    annotated_events: list[GenomicEvent] = []
    for event in events:
        primary, primary_genes, primary_nearest = _annotate_locus(
            event.primary, genes=genes, cytobands=cytobands
        )
        secondary: Locus | None = None
        secondary_genes: list[str] = []
        secondary_nearest: str | None = None
        if event.secondary is not None:
            secondary, secondary_genes, secondary_nearest = _annotate_locus(
                event.secondary, genes=genes, cytobands=cytobands
            )
        flags = list(event.technical_flags)
        for resource_type, resource in contexts:
            if _overlaps(event.primary, resource):
                flags.append(f"primary:{resource_type}")
            if event.secondary is not None and _overlaps(event.secondary, resource):
                flags.append(f"secondary:{resource_type}")
        notes = list(event.notes)
        if primary_nearest:
            notes.append(f"Primary breakpoint {primary_nearest}.")
        if secondary_nearest:
            notes.append(f"Secondary breakpoint {secondary_nearest}.")
        annotated_events.append(
            event.model_copy(
                update={
                    "primary": primary,
                    "secondary": secondary,
                    "genes": sorted(set([*primary_genes, *secondary_genes])),
                    "technical_flags": sorted(set(flags)),
                    "notes": notes,
                    "reportable": False,
                }
            )
        )
    return annotated_events

"""Compile build-native panel selection and analysis derivatives from locked caches."""

from __future__ import annotations

import csv
import math
import os
import sqlite3
import tempfile
from collections.abc import Iterator
from contextlib import closing, contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final, TextIO

import yaml

from .io import load_mapping
from .models import (
    GenomeBuild,
    PanelBundle,
    PanelCoordinateMappingLock,
    ReferenceDictionaryContract,
)
from .panel_bundle import sha256_file
from .reference import (
    grch37_ucsc_hg19_canonical_25_contigs,
    grch38_canonical_25_contigs,
)

TRANSCRIPT_CACHE_HEADER: Final = (
    "target_label",
    "preferred",
    "rank_tier",
    "gene_id",
    "transcript_id",
    "transcript_name",
    "chrom",
    "start",
    "end",
    "strand",
    "transcript_type",
    "mane_status",
    "mane_refseq_id",
    "appris",
    "is_canonical",
    "is_basic",
    "cds_length",
    "transcript_length",
)
TARGET_GENE_MAP_HEADER: Final = (
    "source_target_label",
    "target_label",
    "declared_chromosome",
    "gene_id",
    "native_gene_name",
    "resolution",
)
_TARGET_GENE_RESOLUTIONS: Final = {
    "direct",
    "curated_historical_symbol",
    "unresolved",
    "selection_unmapped",
}


class PanelCompilerError(ValueError):
    """The annotation cache cannot produce a deterministic panel derivative."""


@dataclass(frozen=True)
class PanelTarget:
    chromosome: str
    start: int
    end: int
    label: str
    gene_id: str | None = None


@dataclass(frozen=True)
class TargetGeneMapping:
    source_target_label: str
    target_label: str
    declared_chromosome: str
    gene_id: str | None
    native_gene_name: str | None
    resolution: str


@dataclass(frozen=True)
class CachedGene:
    gene_id: str
    gene_name: str
    chromosome: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class CachedTranscript:
    gene_id: str
    transcript_id: str
    transcript_name: str | None
    chromosome: str
    start: int
    end: int
    strand: str
    transcript_type: str | None
    mane_status: str | None
    mane_refseq_id: str | None
    appris: str | None
    is_canonical: bool
    is_basic: bool
    cds_length: int
    transcript_length: int

    @property
    def rank_tier(self) -> int:
        mane = (self.mane_status or "").lower().replace("_", " ").replace("-", " ")
        appris = (self.appris or "").lower()
        transcript_type = (self.transcript_type or "").lower().replace("_", "-")
        if mane == "select" or mane == "mane select":
            return 1
        if mane in {"plus clinical", "mane plus clinical"}:
            return 2
        if self.is_canonical or appris.startswith(("appris_principal", "principal")):
            return 3
        if transcript_type == "protein-coding" or self.is_basic:
            return 4
        return 5

    @property
    def sort_key(self) -> tuple[int, int, int, str]:
        return (self.rank_tier, -self.cds_length, -self.transcript_length, self.transcript_id)


@dataclass(frozen=True)
class PanelCompilationSummary:
    target_count: int
    resolved_target_count: int
    unresolved_targets: tuple[str, ...]
    roi_interval_count: int
    transcript_count: int
    preferred_transcript_count: int


@dataclass(frozen=True)
class MaterializedPanelSummary:
    compilation: PanelCompilationSummary
    analysis_roi_sha256: str
    transcript_cache_sha256: str
    manifest_sha256: str


def _read_panel_targets(path: Path) -> tuple[PanelTarget, ...]:
    targets: list[PanelTarget] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise PanelCompilerError(f"{path.name} line {line_number}: target label is missing")
        chromosome, raw_start, raw_end, label = fields[:4]
        try:
            start, end = int(raw_start), int(raw_end)
        except ValueError as exc:
            raise PanelCompilerError(
                f"{path.name} line {line_number}: non-integer coordinates"
            ) from exc
        if start < 0 or end <= start or not label.strip():
            raise PanelCompilerError(f"{path.name} line {line_number}: invalid BED interval")
        gene_id = fields[4].strip() if len(fields) >= 5 else ""
        targets.append(
            PanelTarget(
                chromosome=chromosome,
                start=start,
                end=end,
                label=label.strip(),
                gene_id=gene_id if gene_id and gene_id != "." else None,
            )
        )
    if not targets:
        raise PanelCompilerError(f"{path.name}: no target intervals")
    labels = [target.label for target in targets]
    if len(labels) != len(set(labels)):
        raise PanelCompilerError(f"{path.name}: target labels must be unique")
    return tuple(targets)


def read_target_gene_map(path: Path) -> tuple[TargetGeneMapping, ...]:
    """Read an exact, coordinate-free target-to-native-gene curation."""

    records: list[TargetGeneMapping] = []
    header_seen = False
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise PanelCompilerError(f"target gene map cannot be read: {path}") from exc
    for line_number, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = tuple(line.split("\t"))
        if not header_seen:
            if fields != TARGET_GENE_MAP_HEADER:
                raise PanelCompilerError(
                    f"{path.name} line {line_number}: unexpected target-gene-map header"
                )
            header_seen = True
            continue
        if len(fields) != len(TARGET_GENE_MAP_HEADER):
            raise PanelCompilerError(
                f"{path.name} line {line_number}: expected six tab-separated fields"
            )
        source_label, label, chromosome, raw_gene_id, raw_gene_name, resolution = (
            field.strip() for field in fields
        )
        if (
            not source_label
            or not label
            or not chromosome
            or resolution not in _TARGET_GENE_RESOLUTIONS
        ):
            raise PanelCompilerError(f"{path.name} line {line_number}: invalid target gene mapping")
        gene_id = None if raw_gene_id == "." else raw_gene_id
        native_gene_name = None if raw_gene_name == "." else raw_gene_name
        if resolution in {"unresolved", "selection_unmapped"}:
            if gene_id is not None or native_gene_name is not None:
                raise PanelCompilerError(
                    f"{path.name} line {line_number}: unresolved target carries a gene mapping"
                )
        elif gene_id is None or native_gene_name is None:
            raise PanelCompilerError(
                f"{path.name} line {line_number}: resolved target lacks an exact gene mapping"
            )
        elif resolution == "direct" and label != native_gene_name:
            raise PanelCompilerError(
                f"{path.name} line {line_number}: direct mapping changes the target symbol"
            )
        elif resolution == "curated_historical_symbol" and label == native_gene_name:
            raise PanelCompilerError(
                f"{path.name} line {line_number}: historical-symbol mapping does not change symbol"
            )
        records.append(
            TargetGeneMapping(
                source_target_label=source_label,
                target_label=label,
                declared_chromosome=chromosome,
                gene_id=gene_id,
                native_gene_name=native_gene_name,
                resolution=resolution,
            )
        )
    if not header_seen or not records:
        raise PanelCompilerError(f"{path.name}: no target gene mappings")
    labels = [record.target_label for record in records]
    if len(labels) != len(set(labels)):
        raise PanelCompilerError(f"{path.name}: target labels must be unique")
    source_labels = [record.source_target_label for record in records]
    if len(source_labels) != len(set(source_labels)):
        raise PanelCompilerError(f"{path.name}: source target labels must be unique")
    resolved_keys = [
        (record.gene_id, record.declared_chromosome)
        for record in records
        if record.gene_id is not None
    ]
    if len(resolved_keys) != len(set(resolved_keys)):
        raise PanelCompilerError(f"{path.name}: resolved gene mappings must be unique")
    return tuple(records)


def _metadata(connection: sqlite3.Connection) -> dict[str, str]:
    try:
        return {
            str(key): str(value)
            for key, value in connection.execute("SELECT key, value FROM metadata")
        }
    except sqlite3.Error as exc:
        raise PanelCompilerError("annotation cache has no readable metadata table") from exc


def _genes(connection: sqlite3.Connection, label: str, chromosome: str) -> tuple[CachedGene, ...]:
    try:
        rows = connection.execute(
            """
            SELECT gene_id, gene_name, chrom, start, end, strand
            FROM genes
            WHERE gene_name = ? AND chrom = ?
            ORDER BY gene_id
            """,
            (label, chromosome),
        ).fetchall()
    except sqlite3.Error as exc:
        raise PanelCompilerError(
            "annotation cache genes table does not satisfy its contract"
        ) from exc
    return tuple(
        CachedGene(
            gene_id=str(row[0]),
            gene_name=str(row[1]),
            chromosome=str(row[2]),
            start=int(row[3]),
            end=int(row[4]),
            strand=str(row[5]),
        )
        for row in rows
    )


def _gene_by_id(connection: sqlite3.Connection, gene_id: str) -> tuple[CachedGene, ...]:
    try:
        rows = connection.execute(
            """
            SELECT gene_id, gene_name, chrom, start, end, strand
            FROM genes
            WHERE gene_id = ?
            ORDER BY gene_id
            """,
            (gene_id,),
        ).fetchall()
    except sqlite3.Error as exc:
        raise PanelCompilerError(
            "annotation cache genes table does not satisfy its contract"
        ) from exc
    return tuple(
        CachedGene(
            gene_id=str(row[0]),
            gene_name=str(row[1]),
            chromosome=str(row[2]),
            start=int(row[3]),
            end=int(row[4]),
            strand=str(row[5]),
        )
        for row in rows
    )


def _transcripts(connection: sqlite3.Connection, gene_id: str) -> tuple[CachedTranscript, ...]:
    try:
        rows = connection.execute(
            """
            SELECT gene_id, transcript_id, transcript_name, chrom, start, end, strand,
                   transcript_type, mane_status, mane_refseq_id, appris, is_canonical,
                   is_basic, cds_length, transcript_length
            FROM transcripts
            WHERE gene_id = ?
            """,
            (gene_id,),
        ).fetchall()
    except sqlite3.Error as exc:
        raise PanelCompilerError(
            "annotation cache transcripts table does not satisfy its contract"
        ) from exc
    transcripts = tuple(
        CachedTranscript(
            gene_id=str(row[0]),
            transcript_id=str(row[1]),
            transcript_name=str(row[2]) if row[2] is not None else None,
            chromosome=str(row[3]),
            start=int(row[4]),
            end=int(row[5]),
            strand=str(row[6]),
            transcript_type=str(row[7]) if row[7] is not None else None,
            mane_status=str(row[8]) if row[8] is not None else None,
            mane_refseq_id=str(row[9]) if row[9] is not None else None,
            appris=str(row[10]) if row[10] is not None else None,
            is_canonical=bool(row[11]),
            is_basic=bool(row[12]),
            cds_length=int(row[13] or 0),
            transcript_length=int(row[14] or 0),
        )
        for row in rows
    )
    return tuple(sorted(transcripts, key=lambda transcript: transcript.sort_key))


def _write_tsv_row(handle: TextIO, values: tuple[object, ...]) -> None:
    csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(values)


def _atomic_output(path: Path) -> tuple[TextIO, Path]:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    handle = os.fdopen(descriptor, mode="w", encoding="utf-8", newline="")
    return handle, Path(temporary_name)


@contextmanager
def _adjacent_staging_paths(*destinations: Path) -> Iterator[tuple[Path, ...]]:
    """Create publish candidates beside their destinations and always clean leftovers.

    A same-volume Windows rename preserves the source ACL.  Creating a publish candidate in a
    private temporary directory can therefore move a non-inheriting ACL onto a WSL-visible NTFS
    resource.  Adjacent candidates inherit the intended destination-parent ACL before the final
    atomic rename and retain ordinary POSIX behavior on other filesystems.
    """

    staged: list[Path] = []
    try:
        for destination in destinations:
            destination.parent.mkdir(parents=True, exist_ok=True)
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f".{destination.stem}.panel-compile.",
                suffix=destination.suffix or ".tmp",
                dir=destination.parent,
            )
            os.close(descriptor)
            staged.append(Path(temporary_name))
        yield tuple(staged)
    finally:
        for temporary in staged:
            temporary.unlink(missing_ok=True)


def compile_panel_derivatives(
    selection_bed: Path,
    annotation_cache: Path,
    analysis_roi_output: Path,
    transcript_cache_output: Path,
    *,
    expected_genome_build: GenomeBuild = GenomeBuild.GRCH38,
    target_gene_map: Path | None = None,
) -> PanelCompilationSummary:
    """Compile exact gene-body ROI and ranked transcripts from a build-matched cache.

    A target is resolved by its explicitly curated gene ID when a target gene map is supplied;
    otherwise its label must map to exactly one GENCODE gene on the declared chromosome.
    Review-required, curated-unresolved and ambiguous/missing labels remain visible in the
    summary and create neither an ROI nor a negative-observability target.
    """
    targets = _read_panel_targets(selection_bed)
    gene_mappings = read_target_gene_map(target_gene_map) if target_gene_map is not None else ()
    selected_mappings = tuple(
        item for item in gene_mappings if item.resolution != "selection_unmapped"
    )
    if gene_mappings and tuple(item.target_label for item in selected_mappings) != tuple(
        target.label for target in targets
    ):
        raise PanelCompilerError(
            "target gene map labels must exactly match selection BED labels and order"
        )
    mappings_by_label = {item.target_label: item for item in selected_mappings}
    if analysis_roi_output.resolve() == transcript_cache_output.resolve():
        raise PanelCompilerError("ROI and transcript cache outputs must be different files")
    roi_handle, roi_temporary = _atomic_output(analysis_roi_output)
    transcript_handle, transcript_temporary = _atomic_output(transcript_cache_output)
    resolved = 0
    transcript_count = 0
    preferred_count = 0
    unresolved = [
        item.target_label for item in gene_mappings if item.resolution == "selection_unmapped"
    ]
    try:
        with closing(sqlite3.connect(annotation_cache)) as connection:
            connection.execute("PRAGMA query_only = ON")
            metadata = _metadata(connection)
            if metadata.get("genome_build") != expected_genome_build.value:
                raise PanelCompilerError(
                    f"panel derivatives require a {expected_genome_build.value} annotation "
                    "cache; observed "
                    f"{metadata.get('genome_build', 'unspecified')!r}"
                )
            _write_tsv_row(transcript_handle, TRANSCRIPT_CACHE_HEADER)
            for target in targets:
                mapping = mappings_by_label.get(target.label)
                if mapping is not None and mapping.resolution == "unresolved":
                    if target.gene_id is not None:
                        raise PanelCompilerError(
                            f"unresolved target {target.label!r} carries a selection gene ID"
                        )
                    unresolved.append(target.label)
                    continue
                if mapping is None and target.label.endswith("_REVIEW_REQUIRED"):
                    unresolved.append(target.label)
                    continue
                if mapping is not None:
                    assert mapping.gene_id is not None
                    assert mapping.native_gene_name is not None
                    if target.gene_id != mapping.gene_id:
                        raise PanelCompilerError(
                            f"target {target.label!r} selection gene ID disagrees with its "
                            "curated target gene map"
                        )
                    genes = _gene_by_id(connection, mapping.gene_id)
                    if len(genes) != 1:
                        raise PanelCompilerError(
                            f"target {target.label!r} does not resolve to exactly one curated "
                            f"gene ID {mapping.gene_id!r}"
                        )
                    gene = genes[0]
                    if (
                        mapping.declared_chromosome != target.chromosome
                        or gene.chromosome != target.chromosome
                        or mapping.native_gene_name != gene.gene_name
                    ):
                        raise PanelCompilerError(
                            f"target {target.label!r} curated gene mapping disagrees with the "
                            "selection or annotation cache"
                        )
                else:
                    genes = (
                        _gene_by_id(connection, target.gene_id)
                        if target.gene_id is not None
                        else _genes(connection, target.label, target.chromosome)
                    )
                    if len(genes) != 1:
                        unresolved.append(target.label)
                        continue
                    gene = genes[0]
                    if gene.chromosome != target.chromosome:
                        unresolved.append(target.label)
                        continue
                if gene.start >= target.end or gene.end <= target.start:
                    raise PanelCompilerError(
                        f"target {target.label!r} resolves to gene {gene.gene_id!r}, but the "
                        "native gene body does not overlap the mapped selection interval"
                    )
                resolved += 1
                _write_tsv_row(
                    roi_handle,
                    (gene.chromosome, gene.start, gene.end, target.label, gene.gene_id),
                )
                ranked = _transcripts(connection, gene.gene_id)
                for rank_index, transcript in enumerate(ranked):
                    preferred = rank_index == 0
                    transcript_count += 1
                    preferred_count += int(preferred)
                    _write_tsv_row(
                        transcript_handle,
                        (
                            target.label,
                            str(preferred).lower(),
                            transcript.rank_tier,
                            transcript.gene_id,
                            transcript.transcript_id,
                            transcript.transcript_name or "",
                            transcript.chromosome,
                            transcript.start,
                            transcript.end,
                            transcript.strand,
                            transcript.transcript_type or "",
                            transcript.mane_status or "",
                            transcript.mane_refseq_id or "",
                            transcript.appris or "",
                            str(transcript.is_canonical).lower(),
                            str(transcript.is_basic).lower(),
                            transcript.cds_length,
                            transcript.transcript_length,
                        ),
                    )
        roi_handle.close()
        transcript_handle.close()
        os.replace(roi_temporary, analysis_roi_output)
        os.replace(transcript_temporary, transcript_cache_output)
    except Exception:
        roi_handle.close()
        transcript_handle.close()
        roi_temporary.unlink(missing_ok=True)
        transcript_temporary.unlink(missing_ok=True)
        raise
    return PanelCompilationSummary(
        target_count=len(gene_mappings) if gene_mappings else len(targets),
        resolved_target_count=resolved,
        unresolved_targets=tuple(unresolved),
        roi_interval_count=resolved,
        transcript_count=transcript_count,
        preferred_transcript_count=preferred_count,
    )


def transcript_rank_key(transcript: CachedTranscript) -> tuple[int, int, int, str]:
    """Public deterministic key used by compiler and breakpoint annotation tests."""
    return transcript.sort_key


def _bundle_path(bundle_directory: Path, relative_path: str) -> Path:
    root = bundle_directory.resolve()
    candidate = (root / relative_path).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise PanelCompilerError(f"panel resource escapes bundle: {relative_path}") from exc
    return candidate


def _pin_resource(
    document: dict[str, Any],
    resource_id: str,
    artifact: Path,
) -> None:
    resources = document.get("resources")
    if not isinstance(resources, list):
        raise PanelCompilerError("panel bundle resources must be a list")
    for resource in resources:
        if isinstance(resource, dict) and resource.get("resource_id") == resource_id:
            resource["sha256"] = sha256_file(artifact)
            resource["size_bytes"] = artifact.stat().st_size
            return
    raise PanelCompilerError(f"panel bundle does not declare resource {resource_id!r}")


def _canonical_contract_contigs(bundle: PanelBundle) -> set[str] | None:
    compatible: list[set[str]] = []
    if ReferenceDictionaryContract.GRCH38_CANONICAL_25 in bundle.reference_dictionary_contracts:
        compatible.append({name for name, _length in grch38_canonical_25_contigs()})
    if (
        ReferenceDictionaryContract.GRCH37_UCSC_HG19_CANONICAL_25
        in bundle.reference_dictionary_contracts
    ):
        compatible.append({name for name, _length in grch37_ucsc_hg19_canonical_25_contigs()})
    if not compatible:
        return None
    return set.intersection(*compatible)


def _validate_selection_contigs(selection: Path, bundle: PanelBundle) -> None:
    allowed = _canonical_contract_contigs(bundle)
    if allowed is None:
        return
    incompatible = sorted(
        {target.chromosome for target in _read_panel_targets(selection)}.difference(allowed)
    )
    if incompatible:
        raise PanelCompilerError(
            f"panel {bundle.bundle_id!r} declares Canonical-25 compatibility but its "
            f"selection uses unavailable contigs: {', '.join(incompatible)}"
        )


def validate_coordinate_mapping_provenance(
    bundle_directory: Path,
    bundle: PanelBundle,
    selection_path: Path,
    target_gene_map_path: Path,
) -> PanelCoordinateMappingLock:
    """Verify a precomputed build-time mapping without executing or fetching anything."""

    lock_id = bundle.coordinate_mapping_lock_resource_id
    source_id = bundle.coordinate_mapping_source_resource_id
    output_id = bundle.coordinate_mapping_output_resource_id
    unmapped_id = bundle.coordinate_mapping_unmapped_resource_id
    roundtrip_id = bundle.coordinate_mapping_roundtrip_output_resource_id
    roundtrip_unmapped_id = bundle.coordinate_mapping_roundtrip_unmapped_resource_id
    if None in {
        lock_id,
        source_id,
        output_id,
        unmapped_id,
        roundtrip_id,
        roundtrip_unmapped_id,
    }:
        raise PanelCompilerError("panel coordinate mapping resources are incomplete")
    assert lock_id is not None
    assert source_id is not None
    assert output_id is not None
    assert unmapped_id is not None
    assert roundtrip_id is not None
    assert roundtrip_unmapped_id is not None
    lock_resource = bundle.resource(lock_id)
    source_resource = bundle.resource(source_id)
    output_resource = bundle.resource(output_id)
    unmapped_resource = bundle.resource(unmapped_id)
    roundtrip_resource = bundle.resource(roundtrip_id)
    roundtrip_unmapped_resource = bundle.resource(roundtrip_unmapped_id)
    lock_path = _bundle_path(bundle_directory, lock_resource.path)
    source_path = _bundle_path(bundle_directory, source_resource.path)
    output_path = _bundle_path(bundle_directory, output_resource.path)
    unmapped_path = _bundle_path(bundle_directory, unmapped_resource.path)
    roundtrip_path = _bundle_path(bundle_directory, roundtrip_resource.path)
    roundtrip_unmapped_path = _bundle_path(bundle_directory, roundtrip_unmapped_resource.path)
    for label, resource, path in (
        ("mapping lock", lock_resource, lock_path),
        ("source selection", source_resource, source_path),
        ("mapped output", output_resource, output_path),
        ("unmapped output", unmapped_resource, unmapped_path),
        ("roundtrip mapped output", roundtrip_resource, roundtrip_path),
        (
            "roundtrip unmapped output",
            roundtrip_unmapped_resource,
            roundtrip_unmapped_path,
        ),
    ):
        if not path.is_file() or resource.sha256 is None:
            raise PanelCompilerError(f"panel {label} is missing or not checksum-pinned")
        if sha256_file(path) != resource.sha256:
            raise PanelCompilerError(f"panel {label} checksum mismatch")
    try:
        lock = PanelCoordinateMappingLock.model_validate(load_mapping(lock_path))
    except (OSError, ValueError) as exc:
        raise PanelCompilerError("panel coordinate mapping lock is invalid") from exc
    if lock.target_genome_build != bundle.genome_build:
        raise PanelCompilerError("panel mapping target build disagrees with the panel bundle")
    if (
        lock.source_panel_bundle_id != bundle.coordinate_mapping_origin_bundle_id
        or lock.source_panel_resource_id != bundle.coordinate_mapping_origin_resource_id
    ):
        raise PanelCompilerError("panel mapping source identity disagrees with its bundle")
    if lock.source_selection_sha256 != source_resource.sha256:
        raise PanelCompilerError("panel mapping source checksum disagrees with its lock")
    if lock.mapped_output_sha256 != output_resource.sha256:
        raise PanelCompilerError("panel mapped-output checksum disagrees with its lock")
    if lock.unmapped_output_sha256 != unmapped_resource.sha256:
        raise PanelCompilerError("panel unmapped-output checksum disagrees with its lock")
    if (
        lock.roundtrip_mapped_output_sha256 != roundtrip_resource.sha256
        or lock.roundtrip_mapped_output_size_bytes != roundtrip_path.stat().st_size
    ):
        raise PanelCompilerError("panel roundtrip mapped-output identity disagrees with its lock")
    if (
        lock.roundtrip_unmapped_output_sha256 != roundtrip_unmapped_resource.sha256
        or lock.roundtrip_unmapped_output_size_bytes != roundtrip_unmapped_path.stat().st_size
    ):
        raise PanelCompilerError("panel roundtrip unmapped-output identity disagrees with its lock")
    selection_resource = bundle.resource(bundle.selection_panel_resource_id)
    if (
        selection_resource.sha256 is None
        or lock.final_selection_sha256 != selection_resource.sha256
        or sha256_file(selection_path) != selection_resource.sha256
    ):
        raise PanelCompilerError("panel final-selection checksum disagrees with its mapping lock")
    gene_map_id = bundle.target_gene_map_resource_id
    if gene_map_id is None:
        raise PanelCompilerError("coordinate-mapped panel has no target gene map")
    gene_map_resource = bundle.resource(gene_map_id)
    if (
        gene_map_resource.sha256 is None
        or _bundle_path(bundle_directory, gene_map_resource.path).resolve()
        != target_gene_map_path.resolve()
        or sha256_file(target_gene_map_path) != gene_map_resource.sha256
    ):
        raise PanelCompilerError("panel target gene map checksum mismatch")

    source_targets = _read_panel_targets(source_path)
    mapped_targets = _read_panel_targets(output_path)
    unmapped_targets = _read_panel_targets(unmapped_path)
    roundtrip_targets = _read_panel_targets(roundtrip_path)
    roundtrip_unmapped_targets = _read_panel_targets(roundtrip_unmapped_path)
    selection_targets = _read_panel_targets(selection_path)
    gene_mappings = read_target_gene_map(target_gene_map_path)
    if len(source_targets) != lock.source_interval_count:
        raise PanelCompilerError("panel mapping source interval count disagrees with its lock")
    if len(mapped_targets) != lock.mapped_interval_count:
        raise PanelCompilerError("panel mapped interval count disagrees with its lock")
    if tuple(target.label for target in unmapped_targets) != tuple(lock.unmapped_target_labels):
        raise PanelCompilerError("panel unmapped targets disagree with its mapping lock")
    if len(roundtrip_targets) != lock.roundtrip_mapped_interval_count:
        raise PanelCompilerError("panel roundtrip mapped count disagrees with its mapping lock")
    if tuple(target.label for target in roundtrip_unmapped_targets) != tuple(
        lock.roundtrip_unmapped_target_labels
    ):
        raise PanelCompilerError("panel roundtrip unmapped targets disagree with its mapping lock")
    if tuple(target.label for target in source_targets) != tuple(
        mapping.source_target_label for mapping in gene_mappings
    ):
        raise PanelCompilerError(
            "target gene map source labels must exactly match the source selection and order"
        )
    unmapped_labels = set(lock.unmapped_target_labels)
    expected_mapped_sources = tuple(
        target for target in source_targets if target.label not in unmapped_labels
    )
    selected_mappings = tuple(
        mapping for mapping in gene_mappings if mapping.resolution != "selection_unmapped"
    )
    if tuple(target.label for target in mapped_targets) != tuple(
        target.label for target in expected_mapped_sources
    ) or tuple(target.label for target in mapped_targets) != tuple(
        mapping.source_target_label for mapping in selected_mappings
    ):
        raise PanelCompilerError("mapped target labels or order disagree with source provenance")
    if len(selection_targets) != len(mapped_targets):
        raise PanelCompilerError("final selection count differs from mapped coordinate output")
    missing_roundtrip_review = set(lock.roundtrip_review_required_target_labels).difference(
        target.label for target in mapped_targets
    )
    if missing_roundtrip_review:
        raise PanelCompilerError(
            "roundtrip-review targets are absent from the mapped coordinate output"
        )

    same_chromosome = 0
    span_identical = 0
    maximum_span_delta = 0
    maximum_span_delta_fraction = 0.0
    for source, mapped, selection, mapping in zip(
        expected_mapped_sources,
        mapped_targets,
        selection_targets,
        selected_mappings,
        strict=True,
    ):
        same_chromosome += int(source.chromosome == mapped.chromosome)
        source_span = source.end - source.start
        mapped_span = mapped.end - mapped.start
        span_identical += int(source_span == mapped_span)
        span_delta = abs(mapped_span - source_span)
        maximum_span_delta = max(maximum_span_delta, span_delta)
        maximum_span_delta_fraction = max(maximum_span_delta_fraction, span_delta / source_span)
        expected_gene_id = mapping.gene_id
        if (
            (selection.chromosome, selection.start, selection.end)
            != (mapped.chromosome, mapped.start, mapped.end)
            or selection.label != mapping.target_label
            or selection.gene_id != expected_gene_id
        ):
            raise PanelCompilerError(
                f"final selection row {selection.label!r} changes mapped coordinates or curation"
            )
    if same_chromosome != lock.same_chromosome_mapped_count:
        raise PanelCompilerError("same-chromosome mapping count disagrees with its lock")
    if span_identical != lock.span_identical_count:
        raise PanelCompilerError("span-identical mapping count disagrees with its lock")
    if maximum_span_delta != lock.maximum_span_delta_bases:
        raise PanelCompilerError("maximum mapped span delta disagrees with its lock")
    if not math.isclose(
        maximum_span_delta_fraction,
        lock.maximum_span_delta_fraction,
        rel_tol=1e-12,
        abs_tol=1e-15,
    ):
        raise PanelCompilerError("maximum mapped span-delta fraction disagrees with its lock")
    return lock


def materialize_and_pin_panel_derivatives(
    bundle_directory: Path,
    annotation_cache: Path,
) -> MaterializedPanelSummary:
    """Compile pending panel artifacts and activate their manifest last.

    Compilation happens through temporary files adjacent to each destination so Windows/WSL
    ACL inheritance is retained. ROI and transcript artifacts are atomically replaced first;
    the checksum-pinned, model-validated manifest is replaced last. A failure before activation
    leaves the original manifest and existing derivatives untouched, so the registry cannot
    resolve a partial generation.
    """
    bundle_directory = bundle_directory.resolve()
    manifest_path = bundle_directory / "bundle.yaml"
    try:
        document = load_mapping(manifest_path)
        bundle = PanelBundle.model_validate(document)
    except (OSError, ValueError) as exc:
        raise PanelCompilerError(f"panel bundle manifest is invalid: {manifest_path}") from exc
    selection_resource = bundle.resource(bundle.selection_panel_resource_id)
    roi_resource = bundle.resource(bundle.analysis_roi_resource_id)
    transcript_resource = bundle.resource(bundle.transcript_cache_resource_id)
    selection_path = _bundle_path(bundle_directory, selection_resource.path)
    roi_path = _bundle_path(bundle_directory, roi_resource.path)
    transcript_path = _bundle_path(bundle_directory, transcript_resource.path)
    target_gene_map_path: Path | None = None
    if bundle.target_gene_map_resource_id is not None:
        target_gene_map_resource = bundle.resource(bundle.target_gene_map_resource_id)
        target_gene_map_path = _bundle_path(bundle_directory, target_gene_map_resource.path)
        if not target_gene_map_path.is_file() or target_gene_map_resource.sha256 is None:
            raise PanelCompilerError("target gene map is missing or not checksum-pinned")
        if sha256_file(target_gene_map_path) != target_gene_map_resource.sha256:
            raise PanelCompilerError("target gene map checksum mismatch")
    if not annotation_cache.is_file():
        raise PanelCompilerError(f"annotation cache is missing: {annotation_cache}")
    if not selection_path.is_file() or selection_resource.sha256 is None:
        raise PanelCompilerError("normalized selection panel is missing or not checksum-pinned")
    if sha256_file(selection_path) != selection_resource.sha256:
        raise PanelCompilerError("normalized selection panel checksum mismatch")
    if bundle.coordinate_mapping_lock_resource_id is not None:
        if target_gene_map_path is None:
            raise PanelCompilerError(
                "coordinate-mapped panel requires a checksum-pinned target gene map"
            )
        validate_coordinate_mapping_provenance(
            bundle_directory,
            bundle,
            selection_path,
            target_gene_map_path,
        )

    with _adjacent_staging_paths(
        roi_path,
        transcript_path,
        manifest_path,
    ) as (staged_roi, staged_transcripts, staged_manifest):
        _validate_selection_contigs(selection_path, bundle)
        compilation = compile_panel_derivatives(
            selection_path,
            annotation_cache,
            staged_roi,
            staged_transcripts,
            expected_genome_build=bundle.genome_build,
            target_gene_map=target_gene_map_path,
        )
        document["unresolved_targets"] = list(
            dict.fromkeys((*bundle.unresolved_targets, *compilation.unresolved_targets))
        )
        _pin_resource(document, bundle.analysis_roi_resource_id, staged_roi)
        _pin_resource(document, bundle.transcript_cache_resource_id, staged_transcripts)
        staged_manifest.write_text(
            yaml.safe_dump(document, sort_keys=False),
            encoding="utf-8",
            newline="\n",
        )
        try:
            PanelBundle.model_validate(load_mapping(staged_manifest))
        except (OSError, ValueError) as exc:
            raise PanelCompilerError("materialized panel manifest is invalid") from exc
        roi_sha256 = sha256_file(staged_roi)
        transcript_sha256 = sha256_file(staged_transcripts)
        os.replace(staged_roi, roi_path)
        os.replace(staged_transcripts, transcript_path)
        os.replace(staged_manifest, manifest_path)

    return MaterializedPanelSummary(
        compilation=compilation,
        analysis_roi_sha256=roi_sha256,
        transcript_cache_sha256=transcript_sha256,
        manifest_sha256=sha256_file(manifest_path),
    )

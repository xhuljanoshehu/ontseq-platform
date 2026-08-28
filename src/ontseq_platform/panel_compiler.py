"""Compile panel analysis derivatives from the locked GRCh38 annotation cache."""

from __future__ import annotations

import csv
import os
import sqlite3
import tempfile
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final, TextIO

import yaml

from .io import load_mapping
from .models import PanelBundle
from .panel_bundle import sha256_file

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


class PanelCompilerError(ValueError):
    """The annotation cache cannot produce a deterministic panel derivative."""


@dataclass(frozen=True)
class PanelTarget:
    chromosome: str
    start: int
    end: int
    label: str


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
        targets.append(PanelTarget(chromosome, start, end, label.strip()))
    if not targets:
        raise PanelCompilerError(f"{path.name}: no target intervals")
    labels = [target.label for target in targets]
    if len(labels) != len(set(labels)):
        raise PanelCompilerError(f"{path.name}: target labels must be unique")
    return tuple(targets)


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


def compile_panel_derivatives(
    selection_bed: Path,
    annotation_cache: Path,
    analysis_roi_output: Path,
    transcript_cache_output: Path,
) -> PanelCompilationSummary:
    """Compile exact gene-body ROI and ranked transcripts from a GRCh38 cache.

    A label is resolved only if it maps to exactly one GENCODE gene on the target's declared
    chromosome. Review-required and ambiguous/missing labels remain visible in the summary and
    create neither an ROI nor a negative-observability target.
    """
    targets = _read_panel_targets(selection_bed)
    if analysis_roi_output.resolve() == transcript_cache_output.resolve():
        raise PanelCompilerError("ROI and transcript cache outputs must be different files")
    roi_handle, roi_temporary = _atomic_output(analysis_roi_output)
    transcript_handle, transcript_temporary = _atomic_output(transcript_cache_output)
    resolved = 0
    transcript_count = 0
    preferred_count = 0
    unresolved: list[str] = []
    try:
        with closing(sqlite3.connect(annotation_cache)) as connection:
            connection.execute("PRAGMA query_only = ON")
            metadata = _metadata(connection)
            if metadata.get("genome_build") != "GRCh38":
                raise PanelCompilerError(
                    "panel derivatives require a GRCh38 annotation cache; observed "
                    f"{metadata.get('genome_build', 'unspecified')!r}"
                )
            _write_tsv_row(transcript_handle, TRANSCRIPT_CACHE_HEADER)
            for target in targets:
                if target.label.endswith("_REVIEW_REQUIRED"):
                    unresolved.append(target.label)
                    continue
                genes = _genes(connection, target.label, target.chromosome)
                if len(genes) != 1:
                    unresolved.append(target.label)
                    continue
                gene = genes[0]
                if gene.chromosome != target.chromosome:
                    unresolved.append(target.label)
                    continue
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
        target_count=len(targets),
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


def materialize_and_pin_panel_derivatives(
    bundle_directory: Path,
    annotation_cache: Path,
) -> MaterializedPanelSummary:
    """Compile pending panel artifacts and activate their manifest last.

    Compilation happens entirely in a bundle-local staging directory. ROI and transcript
    artifacts are atomically replaced first; the checksum-pinned, model-validated manifest
    is replaced last. A failure before activation leaves the original manifest and existing
    derivatives untouched, so the registry cannot resolve a partial generation.
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
    if not selection_path.is_file() or selection_resource.sha256 is None:
        raise PanelCompilerError("normalized selection panel is missing or not checksum-pinned")
    if sha256_file(selection_path) != selection_resource.sha256:
        raise PanelCompilerError("normalized selection panel checksum mismatch")
    if not annotation_cache.is_file():
        raise PanelCompilerError(f"annotation cache is missing: {annotation_cache}")

    with tempfile.TemporaryDirectory(prefix=".panel-compile-", dir=bundle_directory) as raw:
        staging = Path(raw)
        staged_roi = staging / "analysis_roi.bed"
        staged_transcripts = staging / "transcripts.tsv"
        staged_manifest = staging / "bundle.yaml"
        compilation = compile_panel_derivatives(
            selection_path,
            annotation_cache,
            staged_roi,
            staged_transcripts,
        )
        document["unresolved_targets"] = list(compilation.unresolved_targets)
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
        roi_path.parent.mkdir(parents=True, exist_ok=True)
        transcript_path.parent.mkdir(parents=True, exist_ok=True)
        os.replace(staged_roi, roi_path)
        os.replace(staged_transcripts, transcript_path)
        os.replace(staged_manifest, manifest_path)

    return MaterializedPanelSummary(
        compilation=compilation,
        analysis_roi_sha256=roi_sha256,
        transcript_cache_sha256=transcript_sha256,
        manifest_sha256=sha256_file(manifest_path),
    )

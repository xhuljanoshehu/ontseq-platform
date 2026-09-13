from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from .marlin_contracts import MarlinArtifactLock
from .models import GenomeBuild
from .reference import sha256_file

_EXPECTED_FEATURE_COUNT = 357340


@dataclass(frozen=True)
class MarlinArtifactPaths:
    model: Path
    feature_rdata: Path
    canonical_feature_list: Path
    class_annotations: Path
    probe_bed: Path


def _require_file(path: Path, *, label: str) -> Path:
    candidate = Path(path)
    try:
        resolved = candidate.resolve(strict=True)
    except FileNotFoundError as exc:
        raise ValueError(f"MARLIN {label} artifact is missing: {candidate}") from exc
    if not resolved.is_file():
        raise ValueError(f"MARLIN {label} artifact is not a regular file: {candidate}")
    return resolved


def load_exported_feature_ids(path: Path) -> tuple[str, ...]:
    """Load the deterministic text export of MARLIN v1 feature identifiers."""
    resolved = _require_file(path, label="canonical feature-list")
    feature_ids: list[str] = []
    seen: set[str] = set()
    with resolved.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            feature_id = raw_line.rstrip("\r\n")
            if not feature_id:
                raise ValueError(
                    f"MARLIN feature list contains a blank identifier at line {line_number}"
                )
            if feature_id != feature_id.strip():
                raise ValueError(
                    f"MARLIN feature list contains surrounding whitespace at line {line_number}"
                )
            if "\x00" in feature_id:
                raise ValueError(f"MARLIN feature list contains NUL at line {line_number}")
            if feature_id in seen:
                raise ValueError(
                    f"MARLIN feature list contains duplicate identifier {feature_id!r}"
                )
            seen.add(feature_id)
            feature_ids.append(feature_id)
    if len(feature_ids) != _EXPECTED_FEATURE_COUNT:
        raise ValueError(
            f"MARLIN v1 canonical feature list requires exactly {_EXPECTED_FEATURE_COUNT} IDs; "
            f"found {len(feature_ids)}"
        )
    return tuple(feature_ids)


def _artifact_hash(path: Path, *, label: str) -> str:
    return sha256_file(_require_file(path, label=label))


def create_marlin_artifact_lock(
    paths: MarlinArtifactPaths,
    *,
    lock_id: str,
    code_version_or_commit: str,
    code_manifest_sha256: str,
    genome_build: GenomeBuild,
    R_version: str,
    keras_version: str,
    tensorflow_version: str,
    python_version_if_used: str | None,
    execution_backend: str,
    created_at: datetime,
    model_source_uri: str = "https://doi.org/10.5281/zenodo.15565404",
    code_source_uri: str = "https://github.com/hovestadt/MARLIN",
    feature_source_uri: str = "MARLIN_realtime/files/marlin_v1.features.RData",
    class_annotation_source_uri: str = "MARLIN_realtime/files/marlin_v1.class_annotations.xlsx",
    probe_resource_source_uri: str = "MARLIN_realtime/files/marlin_v1.probes_hg19.bed.gz",
) -> MarlinArtifactLock:
    """Create a content lock without downloading or modifying any MARLIN artifact."""
    if genome_build is not GenomeBuild.GRCH37:
        raise ValueError("MARLIN v1 ONTSeq artifact lock currently requires GRCh37/hg19")
    load_exported_feature_ids(paths.canonical_feature_list)
    return MarlinArtifactLock(
        lock_id=lock_id,
        model_version="1.0.0",
        model_source_uri=model_source_uri,
        model_sha256=_artifact_hash(paths.model, label="model"),
        code_source_uri=code_source_uri,
        code_version_or_commit=code_version_or_commit,
        code_manifest_sha256=code_manifest_sha256,
        feature_source_uri=feature_source_uri,
        feature_sha256=_artifact_hash(paths.feature_rdata, label="feature RData"),
        canonical_feature_list_sha256=_artifact_hash(
            paths.canonical_feature_list, label="canonical feature-list"
        ),
        class_annotation_source_uri=class_annotation_source_uri,
        class_annotation_sha256=_artifact_hash(paths.class_annotations, label="class annotation"),
        probe_resource_source_uri=probe_resource_source_uri,
        probe_resource_sha256=_artifact_hash(paths.probe_bed, label="hg19 probe BED"),
        genome_build=genome_build,
        R_version=R_version,
        keras_version=keras_version,
        tensorflow_version=tensorflow_version,
        python_version_if_used=python_version_if_used,
        execution_backend=execution_backend,
        created_at=created_at,
    )


def _verify_hash(path: Path, *, expected: str, label: str) -> None:
    observed = _artifact_hash(path, label=label)
    if observed != expected:
        raise ValueError(
            f"MARLIN {label} SHA-256 mismatch: expected {expected}, observed {observed}"
        )


def verify_marlin_artifact_lock(lock: MarlinArtifactLock, paths: MarlinArtifactPaths) -> None:
    """Rehash every artifact and revalidate the canonical feature-list semantics."""
    if lock.genome_build is not GenomeBuild.GRCH37:
        raise ValueError("MARLIN v1 ONTSeq artifact lock currently requires GRCh37/hg19")
    if lock.expected_feature_count != _EXPECTED_FEATURE_COUNT:
        raise ValueError("MARLIN v1 lock must declare 357340 model features")
    if lock.expected_model_unit_count != 42:
        raise ValueError("MARLIN v1 lock must declare 42 model units")

    load_exported_feature_ids(paths.canonical_feature_list)
    _verify_hash(paths.model, expected=lock.model_sha256, label="model")
    _verify_hash(paths.feature_rdata, expected=lock.feature_sha256, label="feature RData")
    _verify_hash(
        paths.canonical_feature_list,
        expected=lock.canonical_feature_list_sha256,
        label="canonical feature-list",
    )
    _verify_hash(
        paths.class_annotations,
        expected=lock.class_annotation_sha256,
        label="class annotation",
    )
    _verify_hash(paths.probe_bed, expected=lock.probe_resource_sha256, label="hg19 probe BED")

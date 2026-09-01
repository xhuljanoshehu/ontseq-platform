"""Metadata contracts for large tabular artifacts kept outside PipelineResult."""

from __future__ import annotations

import gzip
from pathlib import Path, PurePosixPath
from typing import TextIO

from .models import SidecarArtifact
from .reference import sha256_file


def _open_text(path: Path) -> TextIO:
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def tabular_sidecar(
    *,
    artifact_id: str,
    envelope_root: Path,
    relative_path: str,
    schema_version: str = "1.0.0",
) -> SidecarArtifact:
    """Fingerprint a TSV/TSV.GZ and count rows without loading it into memory."""

    relative = PurePosixPath(relative_path)
    path = envelope_root.joinpath(*relative.parts).resolve()
    root = envelope_root.resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ValueError("sidecar path escapes the run envelope") from exc
    if not path.is_file():
        raise FileNotFoundError(f"sidecar not found: {path}")

    columns: list[str] = []
    row_count = 0
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                continue
            if not columns:
                columns = stripped.split("\t")
            else:
                row_count += 1
    if not columns:
        raise ValueError(f"sidecar has no tabular header: {path}")
    return SidecarArtifact(
        artifact_id=artifact_id,
        relative_path=relative.as_posix(),
        schema_version=schema_version,
        sha256=sha256_file(path),
        row_count=row_count,
        size_bytes=path.stat().st_size,
        media_type=(
            "application/gzip" if path.suffix.casefold() == ".gz" else "text/tab-separated-values"
        ),
        columns=columns,
    )


__all__ = ["tabular_sidecar"]

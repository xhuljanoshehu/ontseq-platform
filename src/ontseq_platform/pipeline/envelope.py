"""The run envelope: directory layout, atomic writes, fingerprints and resume signatures.

The layout follows ``docs/ARCHITECTURE.md`` so that a run on disk is self-describing and
archivable. Three properties matter more than the layout itself:

**Every recorded path is relative to the envelope root.** Absolute paths leak the source
BAM location and the local directory structure into reviewer artifacts, which the data
boundary in ``docs/DATA_SECURITY.md`` forbids. The API makes the relative form the only
one that can be recorded.

**Every write is atomic.** A stage writes to a temporary name in the same directory and
renames it into place. A run interrupted mid-write therefore leaves either the previous
artifact or none, never a truncated one that a later resume would happily accept.

**Resume is content-addressed, not timestamp-based.** A stage may be skipped only when its
declared inputs, parameters and tool versions hash to the value recorded last time *and*
every output it claims to have produced is still present with a matching checksum. A
mtime comparison would silently accept an artifact produced by different parameters.

This module is dependency-free so the resume logic can be tested without pydantic.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

CHUNK_SIZE = 1024 * 1024

#: Subdirectories created for every run, mirroring docs/ARCHITECTURE.md.
SUBDIRECTORIES: tuple[str, ...] = (
    "manifest",
    "qc",
    "evidence/cnv",
    "evidence/sv",
    "evidence/fusion",
    "alignment",
    "normalized",
    "reports",
    "provenance",
    "release",
    "work",
)

#: Directories holding intermediate data that must never leave the execution system.
NON_EXPORTABLE_DIRECTORIES: frozenset[str] = frozenset({"work", "alignment"})

#: Raw genomic formats are never exportable, wherever they sit. Mirrors the banned
#: suffixes in ``scripts/check_repository_safety.py``: the rule that keeps them out of Git
#: is the same rule that must keep them out of a release bundle.
NON_EXPORTABLE_SUFFIXES: frozenset[str] = frozenset(
    {
        ".bam",
        ".bai",
        ".cram",
        ".crai",
        ".sam",
        ".pod5",
        ".fast5",
        ".fastq",
        ".fq",
        ".vcf",
        ".bcf",
        ".tbi",
        ".csi",
        ".gz",
        ".bedmethyl",
        ".bigwig",
        ".bw",
    }
)


def is_exportable(relative_path: str) -> bool:
    """Return whether an artifact may appear in a release bundle.

    Fails closed twice over: by directory for intermediates, and by suffix for raw
    genomic formats. A file needs to clear both to be exportable.
    """
    if relative_path.split("/", 1)[0] in NON_EXPORTABLE_DIRECTORIES:
        return False
    return not any(relative_path.lower().endswith(suffix) for suffix in NON_EXPORTABLE_SUFFIXES)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path, *, chunk_size: int = CHUNK_SIZE) -> str:
    """Stream a file through SHA-256 without loading it into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class Artifact:
    """A file produced by a stage, addressed relative to the envelope root."""

    relative_path: str
    size_bytes: int
    sha256: str
    exportable: bool

    @property
    def directory(self) -> str:
        return self.relative_path.split("/", 1)[0]


class EnvelopeError(RuntimeError):
    """Raised when the run envelope cannot be created or used safely."""


@dataclass(frozen=True)
class RunEnvelope:
    """A single run's directory, with safe accessors."""

    root: Path
    run_id: str
    sample_id: str

    @classmethod
    def create(cls, base_dir: Path, *, run_id: str, sample_id: str) -> RunEnvelope:
        """Create (or reuse) the envelope for a run.

        Reuse is intentional and is what makes resume possible. The caller decides whether
        an existing envelope is resumed or rejected; this function only guarantees the
        directory skeleton exists.
        """
        for component in (run_id, sample_id):
            if not component or "/" in component or component in {".", ".."}:
                raise EnvelopeError(f"unsafe path component: {component!r}")
        root = base_dir / run_id / sample_id
        for name in SUBDIRECTORIES:
            (root / name).mkdir(parents=True, exist_ok=True)
        return cls(root=root, run_id=run_id, sample_id=sample_id)

    def path(self, relative_path: str) -> Path:
        """Resolve a path inside the envelope, refusing to escape it."""
        if relative_path.startswith("/"):
            raise EnvelopeError("envelope paths must be relative")
        candidate = (self.root / relative_path).resolve()
        root = self.root.resolve()
        if candidate != root and root not in candidate.parents:
            raise EnvelopeError(f"path escapes the run envelope: {relative_path}")
        return candidate

    def relative(self, path: Path) -> str:
        """Return a path's envelope-relative form, which is the only recordable form."""
        try:
            return path.resolve().relative_to(self.root.resolve()).as_posix()
        except ValueError as error:
            raise EnvelopeError(
                "refusing to record a path outside the run envelope; absolute source "
                "paths must not appear in reviewer artifacts"
            ) from error

    def exists(self, relative_path: str) -> bool:
        return self.path(relative_path).is_file()

    def atomic_write_text(self, relative_path: str, text: str) -> Artifact:
        return self.atomic_write_bytes(relative_path, text.encode("utf-8"))

    def atomic_write_bytes(self, relative_path: str, payload: bytes) -> Artifact:
        """Write a file atomically and return its fingerprint."""
        target = self.path(relative_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        descriptor, staged_name = tempfile.mkstemp(
            dir=target.parent, prefix=f".{target.name}.", suffix=".tmp"
        )
        try:
            with os.fdopen(descriptor, "wb") as handle:
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(staged_name, target)
        except BaseException:
            Path(staged_name).unlink(missing_ok=True)
            raise
        return self.fingerprint(relative_path)

    def adopt(self, source: Path, relative_path: str) -> Artifact:
        """Move a file produced elsewhere into the envelope atomically.

        Tools write where they are told to write; when that is a scratch location, this
        brings the result into the envelope without a window in which a partially copied
        file is visible under its final name.
        """
        target = self.path(relative_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        staged = target.parent / f".{target.name}.adopt.tmp"
        staged.unlink(missing_ok=True)
        try:
            shutil.move(str(source), str(staged))
            os.replace(staged, target)
        except BaseException:
            staged.unlink(missing_ok=True)
            raise
        return self.fingerprint(relative_path)

    def fingerprint(self, relative_path: str) -> Artifact:
        """Fingerprint an artifact that already exists in the envelope."""
        target = self.path(relative_path)
        if not target.is_file():
            raise EnvelopeError(f"artifact is missing: {relative_path}")
        normalized = self.relative(target)
        return Artifact(
            relative_path=normalized,
            size_bytes=target.stat().st_size,
            sha256=sha256_file(target),
            exportable=is_exportable(normalized),
        )

    def scratch(self, name: str) -> Path:
        """Return a path under ``work/`` for intermediate files."""
        return self.path(f"work/{name}")

    def verify(self, artifacts: Sequence[Artifact]) -> list[str]:
        """Return a list of artifacts that no longer match their recorded checksum."""
        problems: list[str] = []
        for artifact in artifacts:
            target = self.path(artifact.relative_path)
            if not target.is_file():
                problems.append(f"{artifact.relative_path}: missing")
                continue
            if target.stat().st_size != artifact.size_bytes:
                problems.append(f"{artifact.relative_path}: size changed")
                continue
            if sha256_file(target) != artifact.sha256:
                problems.append(f"{artifact.relative_path}: checksum changed")
        return problems


def canonical_signature(payload: Mapping[str, object]) -> str:
    """Hash a mapping deterministically.

    Used for stage resume signatures. Keys are sorted and separators fixed so that the
    same logical inputs always hash identically regardless of dictionary ordering, and
    ``default=str`` keeps the function total rather than raising on an unexpected value
    type at the worst possible moment, mid-run.
    """
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return sha256_bytes(encoded.encode("utf-8"))


def stage_signature(
    *,
    stage: str,
    upstream: Sequence[Artifact],
    parameters: Mapping[str, object],
    tool_versions: Mapping[str, str],
    external_inputs: Sequence[tuple[str, str]] = (),
) -> str:
    """Compute the resume signature of a stage.

    Includes the checksums of every upstream artifact, the stage's own parameters, the
    resolved tool versions and the fingerprints of any inputs from outside the envelope.
    Change any one of them and the stage re-runs, which is the only safe default: a
    resume that ignores a parameter change would silently mix results from two different
    configurations inside one envelope.
    """
    return canonical_signature(
        {
            "stage": stage,
            "upstream": sorted((item.relative_path, item.sha256) for item in upstream),
            "parameters": dict(sorted(parameters.items())),
            "tool_versions": dict(sorted(tool_versions.items())),
            "external_inputs": sorted(external_inputs),
        }
    )

"""Catalog and transactional installation services for ONTSeq reference bundles.

Only an activated ``references/<bundle-id>/bundle.yaml`` is discoverable by the registry.  An
installation is assembled under ``references/.staging`` and moves into place only after every
source checksum, derived artifact and bundle invariant has passed.  Normal analyses only consume
local activated files; this module performs network I/O exclusively in explicit install/repair
operations.
"""

from __future__ import annotations

import gzip
import json
import os
import shutil
import tempfile
import urllib.request
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager, suppress
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path, PurePosixPath
from typing import BinaryIO, cast
from urllib.parse import urlparse

import yaml

from .annotation_cache import compile_annotation_cache, validate_annotation_cache
from .io import load_mapping
from .models import GenomeBuild, ReferenceBundle, ReferenceLock, ResourceFile
from .path_safety import (
    assert_plain_tree,
    assert_safe_descendant,
    is_link_like,
    path_lexists,
)
from .reference import (
    reference_lock_from_fai,
    sha256_file,
    validate_canonical_reference,
)
from .resource_registry import resource_root_from_environment

BUNDLE_MANIFEST_NAME = "bundle.yaml"


class ResourceValidationState(StrEnum):
    OK = "ok"
    MISSING = "missing"
    UNRESOLVED_CHECKSUM = "unresolved_checksum"
    SIZE_MISMATCH = "size_mismatch"
    CHECKSUM_MISMATCH = "checksum_mismatch"
    UNSAFE_PATH = "unsafe_path"
    INVALID_DERIVED_ARTIFACT = "invalid_derived_artifact"


@dataclass(frozen=True)
class ResourceValidation:
    resource_id: str
    role: str
    path: Path
    required: bool
    state: ResourceValidationState
    message: str
    observed_size_bytes: int | None = None
    observed_sha256: str | None = None

    @property
    def valid(self) -> bool:
        return self.state == ResourceValidationState.OK or (
            not self.required and self.state == ResourceValidationState.MISSING
        )


@dataclass(frozen=True)
class BundleValidationReport:
    bundle_id: str
    bundle_path: Path
    manifest_valid: bool
    resources: tuple[ResourceValidation, ...]
    errors: tuple[str, ...] = ()

    @property
    def valid(self) -> bool:
        return (
            self.manifest_valid
            and not self.errors
            and all(resource.valid for resource in self.resources)
        )

    @property
    def missing_resource_ids(self) -> tuple[str, ...]:
        return tuple(
            item.resource_id
            for item in self.resources
            if item.state == ResourceValidationState.MISSING
        )


@dataclass(frozen=True)
class InstalledReferenceBundle:
    bundle: ReferenceBundle
    path: Path
    validation: BundleValidationReport


DownloadOpener = Callable[[str], BinaryIO]


def _default_opener(url: str) -> BinaryIO:
    request = urllib.request.Request(url, headers={"User-Agent": "ONTSeq-reference-installer/1"})
    return cast(BinaryIO, urllib.request.urlopen(request, timeout=120))


def _load_reference_bundle(path: Path) -> ReferenceBundle:
    bundle = ReferenceBundle.model_validate(load_mapping(path))
    if bundle.genome_build != GenomeBuild.GRCH38:
        raise ValueError(
            f"active implementation accepts only GRCh38 reference bundles, got "
            f"{bundle.genome_build.value}"
        )
    return bundle


def _resource_path(bundle_path: Path, resource: ResourceFile) -> Path:
    candidate = bundle_path.joinpath(*PurePosixPath(resource.path).parts)
    try:
        assert_safe_descendant(
            bundle_path, candidate, label=f"resource {resource.resource_id!r} path"
        )
    except ValueError as exc:
        raise ValueError(f"resource {resource.resource_id!r} has an unsafe path: {exc}") from exc
    return candidate


def _write_bundle_manifest(bundle: ReferenceBundle, path: Path) -> None:
    payload = bundle.model_dump(mode="json", exclude_none=True)
    rendered = yaml.safe_dump(
        payload,
        allow_unicode=True,
        default_flow_style=False,
        sort_keys=False,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    ) as handle:
        handle.write(rendered)
        temporary = Path(handle.name)
    os.replace(temporary, path)


class ReferenceCatalog:
    """A deterministic collection of non-activated reference bundle recipes."""

    def __init__(self, bundles: Mapping[str, ReferenceBundle]) -> None:
        self._bundles = dict(sorted(bundles.items()))

    @classmethod
    def discover(cls, catalog_root: Path) -> ReferenceCatalog:
        if not catalog_root.is_dir():
            raise FileNotFoundError(f"reference catalog directory not found: {catalog_root}")
        bundles: dict[str, ReferenceBundle] = {}
        candidates = sorted(
            {
                *catalog_root.rglob("bundle.yaml"),
                *catalog_root.rglob("bundle.recipe.yaml"),
            }
        )
        for manifest_path in candidates:
            bundle = _load_reference_bundle(manifest_path)
            if bundle.bundle_id in bundles:
                raise ValueError(f"duplicate reference catalog bundle ID: {bundle.bundle_id}")
            bundles[bundle.bundle_id] = bundle
        return cls(bundles)

    @classmethod
    def from_manifests(cls, manifests: Sequence[Path]) -> ReferenceCatalog:
        bundles: dict[str, ReferenceBundle] = {}
        for path in manifests:
            bundle = _load_reference_bundle(path)
            if bundle.bundle_id in bundles:
                raise ValueError(f"duplicate reference catalog bundle ID: {bundle.bundle_id}")
            bundles[bundle.bundle_id] = bundle
        return cls(bundles)

    def list(self) -> tuple[ReferenceBundle, ...]:
        return tuple(self._bundles.values())

    def get(self, bundle_id: str) -> ReferenceBundle:
        try:
            return self._bundles[bundle_id]
        except KeyError as exc:
            available = ", ".join(self._bundles) or "none"
            raise KeyError(
                f"unknown reference bundle {bundle_id!r}; available: {available}"
            ) from exc

    def find(self, bundle_id: str) -> ReferenceBundle | None:
        """Return a catalog authority when one exists, without claiming custom bundle IDs."""

        return self._bundles.get(bundle_id)


_REFERENCE_CONTRACT_FIELDS = (
    "schema_version",
    "bundle_type",
    "bundle_id",
    "version",
    "genome_build",
    "reference_lock_resource_id",
    "fasta_resource_id",
    "fai_resource_id",
    "annotation_cache_resource_id",
)
_RESOURCE_CONTRACT_FIELDS = (
    "resource_id",
    "role",
    "path",
    "source_url",
    "release",
    "source_date",
    "coordinate_system",
    "generated",
    "derived_from",
    "required",
)


def validate_reference_recipe_contract(
    candidate: ReferenceBundle,
    authority: ReferenceBundle,
    *,
    allow_materialized_generated_pins: bool = True,
) -> None:
    """Require one immutable source/generator contract for a bundle ID and version.

    Activated manifests legitimately add the observed SHA256 and size of generated artifacts.
    Those two materialized values are accepted only when the catalog recipe left them unresolved;
    all source pins, generator dependencies, paths, roles and bundle pointers remain exact.
    """

    differences: list[str] = []
    for field in _REFERENCE_CONTRACT_FIELDS:
        if getattr(candidate, field) != getattr(authority, field):
            differences.append(field)
    candidate_ids = [resource.resource_id for resource in candidate.resources]
    authority_ids = [resource.resource_id for resource in authority.resources]
    if candidate_ids != authority_ids:
        differences.append("resources.order_or_membership")
    else:
        for candidate_resource, authority_resource in zip(
            candidate.resources, authority.resources, strict=True
        ):
            for field in _RESOURCE_CONTRACT_FIELDS:
                if getattr(candidate_resource, field) != getattr(authority_resource, field):
                    differences.append(f"{candidate_resource.resource_id}.{field}")
            if (
                authority_resource.sha256 is not None
                and candidate_resource.sha256 != authority_resource.sha256
            ) or (
                authority_resource.generated
                and not allow_materialized_generated_pins
                and candidate_resource.sha256 is not None
            ):
                differences.append(f"{candidate_resource.resource_id}.sha256")
            if (
                authority_resource.size_bytes is not None
                and candidate_resource.size_bytes != authority_resource.size_bytes
            ) or (
                authority_resource.generated
                and not allow_materialized_generated_pins
                and candidate_resource.size_bytes is not None
            ):
                differences.append(f"{candidate_resource.resource_id}.size_bytes")
    if differences:
        examples = ", ".join(differences[:5])
        suffix = " ..." if len(differences) > 5 else ""
        raise ValueError(
            f"reference bundle {candidate.bundle_id!r} does not match the immutable catalog "
            f"source/generator contract ({examples}{suffix}); publish a new bundle ID/version "
            "for changed sources or derivation rules"
        )


def _validate_resource(
    bundle_path: Path,
    resource: ResourceFile,
    *,
    verify_checksum: bool,
) -> ResourceValidation:
    try:
        path = _resource_path(bundle_path, resource)
    except ValueError as exc:
        return ResourceValidation(
            resource.resource_id,
            resource.role,
            bundle_path / resource.path,
            resource.required,
            ResourceValidationState.UNSAFE_PATH,
            str(exc),
        )
    if not path.is_file():
        return ResourceValidation(
            resource.resource_id,
            resource.role,
            path,
            resource.required,
            ResourceValidationState.MISSING,
            "required file is missing" if resource.required else "optional file is missing",
        )
    if path.is_symlink():
        return ResourceValidation(
            resource.resource_id,
            resource.role,
            path,
            resource.required,
            ResourceValidationState.UNSAFE_PATH,
            "bundle resources may not be symbolic links",
        )
    size = path.stat().st_size
    if resource.size_bytes is not None and size != resource.size_bytes:
        return ResourceValidation(
            resource.resource_id,
            resource.role,
            path,
            resource.required,
            ResourceValidationState.SIZE_MISMATCH,
            f"expected {resource.size_bytes} bytes, found {size}",
            observed_size_bytes=size,
        )
    if resource.sha256 is None:
        return ResourceValidation(
            resource.resource_id,
            resource.role,
            path,
            resource.required,
            ResourceValidationState.UNRESOLVED_CHECKSUM,
            "activated resources require a resolved SHA256 checksum",
            observed_size_bytes=size,
        )
    observed_sha256 = sha256_file(path) if verify_checksum else None
    if observed_sha256 is not None and observed_sha256 != resource.sha256:
        return ResourceValidation(
            resource.resource_id,
            resource.role,
            path,
            resource.required,
            ResourceValidationState.CHECKSUM_MISMATCH,
            f"expected SHA256 {resource.sha256}, found {observed_sha256}",
            observed_size_bytes=size,
            observed_sha256=observed_sha256,
        )
    if resource.role == "annotation_cache" and verify_checksum:
        try:
            validate_annotation_cache(path)
        except (OSError, ValueError) as exc:
            return ResourceValidation(
                resource.resource_id,
                resource.role,
                path,
                resource.required,
                ResourceValidationState.INVALID_DERIVED_ARTIFACT,
                str(exc),
                observed_size_bytes=size,
                observed_sha256=observed_sha256,
            )
    return ResourceValidation(
        resource.resource_id,
        resource.role,
        path,
        resource.required,
        ResourceValidationState.OK,
        "checksum and size verified" if verify_checksum else "present with expected size",
        observed_size_bytes=size,
        observed_sha256=observed_sha256,
    )


def validate_reference_bundle_directory(
    bundle_path: Path,
    *,
    verify_checksums: bool = True,
) -> BundleValidationReport:
    """Validate an activated directory without reaching the network or mutating it."""

    if is_link_like(bundle_path):
        return BundleValidationReport(
            bundle_id=bundle_path.name,
            bundle_path=bundle_path,
            manifest_valid=False,
            resources=(),
            errors=(f"reference bundle directory is a symbolic link or junction: {bundle_path}",),
        )
    manifest_path = bundle_path / BUNDLE_MANIFEST_NAME
    if not manifest_path.is_file():
        return BundleValidationReport(
            bundle_id=bundle_path.name,
            bundle_path=bundle_path,
            manifest_valid=False,
            resources=(),
            errors=(f"bundle manifest is missing: {manifest_path}",),
        )
    try:
        bundle = _load_reference_bundle(manifest_path)
    except (OSError, ValueError) as exc:
        return BundleValidationReport(
            bundle_id=bundle_path.name,
            bundle_path=bundle_path,
            manifest_valid=False,
            resources=(),
            errors=(f"invalid bundle manifest: {exc}",),
        )
    resources = tuple(
        _validate_resource(bundle_path, item, verify_checksum=verify_checksums)
        for item in bundle.resources
    )
    errors: list[str] = []
    try:
        fasta_path = _resource_path(bundle_path, bundle.resource(bundle.fasta_resource_id))
        fai_path = _resource_path(bundle_path, bundle.resource(bundle.fai_resource_id))
        if fasta_path.is_file() and fai_path.is_file():
            if verify_checksums:
                validate_fasta_fai_consistency(fasta_path, fai_path)
            expected_lock = reference_lock_from_fai(
                fai_path,
                reference_id=bundle.bundle_id,
                genome_build=bundle.genome_build,
                allow_extra_contigs=True,
                require_canonical_assembly=True,
            )
            validate_canonical_reference(
                ((item.name, item.length) for item in expected_lock.contigs),
                bundle.genome_build,
            )
            lock_resource = bundle.resource(bundle.reference_lock_resource_id)
            lock_path = _resource_path(bundle_path, lock_resource)
            if lock_path.is_file():
                declared_lock = ReferenceLock.model_validate(load_mapping(lock_path))
                if declared_lock != expected_lock:
                    raise ValueError("reference lock content does not match the pinned FASTA index")
    except (OSError, UnicodeError, ValueError) as exc:
        errors.append(f"FASTA/index/build validation failed: {exc}")
    return BundleValidationReport(
        bundle_id=bundle.bundle_id,
        bundle_path=bundle_path,
        manifest_valid=True,
        resources=resources,
        errors=tuple(errors),
    )


def _fasta_index_records(fasta_path: Path) -> tuple[tuple[str, int, int, int, int], ...]:
    """Scan FASTA bytes and return the exact five-column samtools-faidx records.

    This intentionally reads the complete FASTA.  Comparing only sequence names and declared
    lengths would allow a synthetic or stale FAI to make an unrelated FASTA look like GRCh38.
    """

    records: list[tuple[str, int, int, int, int]] = []
    with fasta_path.open("rb") as handle:
        name: str | None = None
        sequence_length = 0
        sequence_offset = 0
        line_bases: int | None = None
        line_width: int | None = None
        saw_short_line = False
        while raw_line := handle.readline():
            if raw_line.startswith(b">"):
                if name is not None:
                    if line_bases is None or line_width is None:
                        raise ValueError(f"FASTA sequence {name!r} has no bases")
                    records.append((name, sequence_length, sequence_offset, line_bases, line_width))
                header = raw_line[1:].strip().split(maxsplit=1)
                if not header:
                    raise ValueError("FASTA contains an empty sequence header")
                name = header[0].decode("utf-8")
                sequence_length = 0
                sequence_offset = handle.tell()
                line_bases = None
                line_width = None
                saw_short_line = False
                continue
            if name is None:
                if raw_line.strip():
                    raise ValueError("FASTA sequence data occurs before the first header")
                continue
            stripped = raw_line.rstrip(b"\r\n")
            if not stripped:
                raise ValueError(f"FASTA sequence {name!r} contains a blank line")
            if saw_short_line:
                raise ValueError(f"FASTA sequence {name!r} has data after a short final line")
            bases, width = len(stripped), len(raw_line)
            if line_bases is None:
                line_bases, line_width = bases, width
            elif bases > line_bases or (bases == line_bases and width != line_width):
                raise ValueError(f"FASTA sequence {name!r} has inconsistent wrapping")
            elif bases < line_bases:
                saw_short_line = True
            sequence_length += bases
        if name is not None:
            if line_bases is None or line_width is None:
                raise ValueError(f"FASTA sequence {name!r} has no bases")
            records.append((name, sequence_length, sequence_offset, line_bases, line_width))
    if not records:
        raise ValueError(f"FASTA contains no sequences: {fasta_path}")
    return tuple(records)


def _fasta_index(fasta_path: Path, output: Path) -> None:
    records = _fasta_index_records(fasta_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "".join("\t".join(map(str, record)) + "\n" for record in records),
        encoding="utf-8",
        newline="\n",
    )


def _read_fasta_index_records(fai_path: Path) -> tuple[tuple[str, int, int, int, int], ...]:
    records: list[tuple[str, int, int, int, int]] = []
    with fai_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                raise ValueError(f"FAI line {line_number}: empty lines are not permitted")
            fields = line.split("\t")
            if len(fields) != 5:
                raise ValueError(
                    f"FAI line {line_number}: expected exactly five FASTA index columns"
                )
            name = fields[0]
            if not name:
                raise ValueError(f"FAI line {line_number}: contig name is empty")
            try:
                length, offset, line_bases, line_width = map(int, fields[1:])
            except ValueError as exc:
                raise ValueError(f"FAI line {line_number}: index fields must be integers") from exc
            if length <= 0 or offset < 0 or line_bases <= 0 or line_width < line_bases:
                raise ValueError(f"FAI line {line_number}: invalid FASTA index geometry")
            records.append((name, length, offset, line_bases, line_width))
    if not records:
        raise ValueError(f"FASTA index contains no contigs: {fai_path}")
    return tuple(records)


def validate_fasta_fai_consistency(fasta_path: Path, fai_path: Path) -> None:
    """Prove that ``fai_path`` indexes the exact bytes in ``fasta_path``.

    The full record includes sequence order, length, byte offset and wrapping geometry.  It is
    therefore stronger than the FAI-to-lock comparison and closes the gap where a canonical
    GRCh38-shaped FAI could previously be paired with a tiny unrelated FASTA.
    """

    observed = _read_fasta_index_records(fai_path)
    expected = _fasta_index_records(fasta_path)
    if observed == expected:
        return
    mismatch = next(
        (
            index
            for index, (observed_record, expected_record) in enumerate(
                zip(observed, expected, strict=False), start=1
            )
            if observed_record != expected_record
        ),
        min(len(observed), len(expected)) + 1,
    )
    observed_record = observed[mismatch - 1] if mismatch <= len(observed) else None
    expected_record = expected[mismatch - 1] if mismatch <= len(expected) else None
    raise ValueError(
        "FASTA index is not derived from the installed FASTA bytes: "
        f"first mismatch at record {mismatch}; FAI={observed_record!r}, "
        f"FASTA-derived={expected_record!r}"
    )


def _chromosome_sizes(fai_path: Path, output: Path) -> None:
    rows: list[str] = []
    with fai_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"FAI line {line_number}: expected at least two columns")
            int(fields[1])
            rows.append(f"{fields[0]}\t{fields[1]}\n")
    if not rows:
        raise ValueError("cannot derive chromosome sizes from an empty FASTA index")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(rows), encoding="utf-8", newline="\n")


def _gunzip(source: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rb") as source_handle, output.open("wb") as output_handle:
        shutil.copyfileobj(source_handle, output_handle, length=1024 * 1024)


def _normalize_ucsc_context_table(source: Path, output: Path, *, role: str) -> None:
    """Compile selected UCSC database dumps into canonical 0-based BED.

    UCSC database tables already use zero-based, half-open coordinates, but the coordinate
    columns differ by table. Keeping this adapter in the installer makes the active context
    resources uniform and lets analyses consume only BED rather than repeatedly interpreting
    raw database schemas.
    """

    layouts = {
        "repeatmasker": (5, 6, 7, (10, 11, 12)),
        "simple_repeats": (1, 2, 3, (4,)),
        "segmental_duplication": (1, 2, 3, (4,)),
    }
    try:
        chrom_index, start_index, end_index, label_indices = layouts[role]
    except KeyError as exc:
        raise ValueError(f"no UCSC context-table layout is registered for {role!r}") from exc
    output.parent.mkdir(parents=True, exist_ok=True)
    with (
        gzip.open(source, "rt", encoding="utf-8") as source_handle,
        output.open("w", encoding="utf-8", newline="\n") as output_handle,
    ):
        for line_number, raw in enumerate(source_handle, start=1):
            fields = raw.rstrip("\r\n").split("\t")
            required_index = max(chrom_index, start_index, end_index, *label_indices)
            if len(fields) <= required_index:
                raise ValueError(
                    f"{role} source line {line_number} has {len(fields)} columns; "
                    f"expected at least {required_index + 1}"
                )
            try:
                start, end = int(fields[start_index]), int(fields[end_index])
            except ValueError as exc:
                raise ValueError(
                    f"{role} source line {line_number} has non-integer coordinates"
                ) from exc
            chromosome = fields[chrom_index]
            if start < 0 or end <= start or not chromosome:
                raise ValueError(f"{role} source line {line_number} has an invalid interval")
            label = "|".join(fields[index] for index in label_indices if fields[index]) or role
            output_handle.write(f"{chromosome}\t{start}\t{end}\t{label}\n")


class ReferenceBundleInstaller:
    """Install, validate, import and repair manifest-pinned reference bundles."""

    def __init__(
        self,
        resource_root: Path | str | None = None,
        *,
        opener: DownloadOpener = _default_opener,
    ) -> None:
        self.resource_root = resource_root_from_environment(resource_root).resolve()
        self.references_root = self.resource_root / "references"
        self._opener = opener

    def _assert_safe_reference_path(self, path: Path, *, label: str) -> None:
        try:
            assert_safe_descendant(self.resource_root, path, label=label)
        except ValueError as exc:
            raise ValueError(str(exc)) from exc

    def _staging_directory(self, bundle_id: str) -> Path:
        staging_parent = self.references_root / ".staging"
        self._assert_safe_reference_path(staging_parent, label="reference staging directory")
        staging_parent.mkdir(parents=True, exist_ok=True)
        self._assert_safe_reference_path(staging_parent, label="reference staging directory")
        staging = Path(tempfile.mkdtemp(prefix=f"{bundle_id}.", dir=staging_parent))
        self._assert_safe_reference_path(staging, label="reference staging transaction")
        return staging

    @contextmanager
    def _lock(self, bundle_id: str) -> Iterator[None]:
        self._assert_safe_reference_path(self.references_root, label="reference bundle namespace")
        self.references_root.mkdir(parents=True, exist_ok=True)
        self._assert_safe_reference_path(self.references_root, label="reference bundle namespace")
        lock_path = self.references_root / f".{bundle_id}.install.lock"
        self._assert_safe_reference_path(lock_path, label="reference installer lock")
        try:
            descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        except FileExistsError as exc:
            raise RuntimeError(
                f"another install or repair is already active for {bundle_id}"
            ) from exc
        try:
            os.write(descriptor, str(os.getpid()).encode("ascii"))
            os.close(descriptor)
            yield
        finally:
            with suppress(OSError):
                os.close(descriptor)
            lock_path.unlink(missing_ok=True)

    def _download(self, resource: ResourceFile, target: Path, *, offline: bool) -> None:
        if offline:
            raise RuntimeError(
                f"offline mode cannot retrieve missing source {resource.resource_id!r}"
            )
        if resource.source_url is None:
            raise ValueError(
                f"source resource {resource.resource_id!r} has no source_url; use import"
            )
        scheme = urlparse(resource.source_url).scheme.casefold()
        if scheme not in {"https", "file"}:
            raise ValueError(
                f"resource {resource.resource_id!r} uses unsupported URL scheme {scheme!r}"
            )
        target.parent.mkdir(parents=True, exist_ok=True)
        with self._opener(resource.source_url) as response, target.open("wb") as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)
        observed_size = target.stat().st_size
        if resource.size_bytes is not None and observed_size != resource.size_bytes:
            raise ValueError(
                f"downloaded {resource.resource_id!r}: expected {resource.size_bytes} bytes, "
                f"found {observed_size}"
            )
        if resource.sha256 is None:
            raise ValueError(f"source resource {resource.resource_id!r} has no pinned SHA256")
        observed_hash = sha256_file(target)
        if observed_hash != resource.sha256:
            raise ValueError(
                f"downloaded {resource.resource_id!r}: expected SHA256 {resource.sha256}, "
                f"found {observed_hash}"
            )

    @staticmethod
    def _role_path(
        bundle: ReferenceBundle, bundle_path: Path, role: str
    ) -> tuple[ResourceFile, Path]:
        matches = [item for item in bundle.resources if item.role == role]
        if len(matches) != 1:
            raise ValueError(
                f"bundle {bundle.bundle_id!r} requires exactly one resource with role {role!r}"
            )
        return matches[0], _resource_path(bundle_path, matches[0])

    def _generate(
        self,
        bundle: ReferenceBundle,
        resource: ResourceFile,
        bundle_path: Path,
    ) -> None:
        output = _resource_path(bundle_path, resource)
        output.parent.mkdir(parents=True, exist_ok=True)
        if resource.role == "genome_fasta":
            if len(resource.derived_from) != 1:
                raise ValueError("derived genome FASTA requires exactly one gzip source")
            source = _resource_path(bundle_path, bundle.resource(resource.derived_from[0]))
            _gunzip(source, output)
        elif resource.role == "fasta_index":
            _, fasta_path = self._role_path(bundle, bundle_path, "genome_fasta")
            _fasta_index(fasta_path, output)
        elif resource.role == "chromosome_sizes":
            _, fai_path = self._role_path(bundle, bundle_path, "fasta_index")
            _chromosome_sizes(fai_path, output)
        elif resource.role == "reference_lock":
            _, fai_path = self._role_path(bundle, bundle_path, "fasta_index")
            lock = reference_lock_from_fai(
                fai_path,
                reference_id=bundle.bundle_id,
                genome_build=bundle.genome_build,
                allow_extra_contigs=True,
                require_canonical_assembly=True,
            )
            output.write_text(
                json.dumps(
                    lock.model_dump(mode="json"),
                    indent=2,
                    ensure_ascii=False,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
                newline="\n",
            )
        elif resource.role == "annotation_cache":
            gencode_resource, gencode = self._role_path(bundle, bundle_path, "gencode_gtf")
            mane_resource, mane = self._role_path(bundle, bundle_path, "mane_gff3")
            cytoband_resource, cytobands = self._role_path(bundle, bundle_path, "cytobands")
            compile_annotation_cache(
                gencode,
                mane,
                cytobands,
                output,
                metadata={
                    "bundle_id": bundle.bundle_id,
                    "bundle_version": bundle.version,
                    "genome_build": bundle.genome_build.value,
                    "gencode_release": gencode_resource.release or "unspecified",
                    "gencode_sha256": sha256_file(gencode),
                    "mane_release": mane_resource.release or "unspecified",
                    "mane_sha256": sha256_file(mane),
                    "cytoband_release": cytoband_resource.release or "unspecified",
                    "cytoband_sha256": sha256_file(cytobands),
                },
            )
        elif resource.role in {"repeatmasker", "simple_repeats", "segmental_duplication"}:
            if len(resource.derived_from) != 1:
                raise ValueError(f"derived {resource.role} requires exactly one UCSC source table")
            source = _resource_path(bundle_path, bundle.resource(resource.derived_from[0]))
            _normalize_ucsc_context_table(source, output, role=resource.role)
        else:
            raise ValueError(f"no deterministic generator is registered for role {resource.role!r}")

    def _materialize(
        self,
        recipe: ReferenceBundle,
        bundle_path: Path,
        *,
        offline: bool,
        only_resource_ids: set[str] | None = None,
        local_source_root: Path | None = None,
    ) -> ReferenceBundle:
        resolved: dict[str, ResourceFile] = {}
        for resource in recipe.resources:
            if resource.generated:
                continue
            target = _resource_path(bundle_path, resource)
            if only_resource_ids is None or resource.resource_id in only_resource_ids:
                download_target = target
                if only_resource_ids is not None:
                    download_target = target.with_name(f".{target.name}.repair.tmp")
                    download_target.unlink(missing_ok=True)
                try:
                    if local_source_root is None:
                        self._download(resource, download_target, offline=offline)
                    else:
                        source = _resource_path(local_source_root, resource)
                        if not source.is_file() or source.is_symlink():
                            raise FileNotFoundError(
                                f"pinned import source is absent or unsafe: {source}"
                            )
                        if (
                            resource.size_bytes is not None
                            and (observed_size := source.stat().st_size) != resource.size_bytes
                        ):
                            raise ValueError(
                                f"import source {resource.resource_id!r}: expected "
                                f"{resource.size_bytes} bytes, found {observed_size}"
                            )
                        if resource.sha256 is None:
                            raise ValueError(
                                f"import source {resource.resource_id!r} has no pinned SHA256"
                            )
                        observed_hash = sha256_file(source)
                        if observed_hash != resource.sha256:
                            raise ValueError(
                                f"import source {resource.resource_id!r}: expected SHA256 "
                                f"{resource.sha256}, found {observed_hash}"
                            )
                        download_target.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copyfile(source, download_target)
                except BaseException:
                    if download_target != target:
                        download_target.unlink(missing_ok=True)
                    raise
                if download_target != target:
                    os.replace(download_target, target)
            elif not target.is_file():
                raise FileNotFoundError(f"retained source is missing during repair: {target}")
            resolved[resource.resource_id] = resource.model_copy(
                update={"size_bytes": target.stat().st_size, "sha256": sha256_file(target)}
            )

        pending = {item.resource_id: item for item in recipe.resources if item.generated}
        while pending:
            progressed = False
            for resource_id in list(pending):
                resource = pending[resource_id]
                if not set(resource.derived_from).issubset(resolved):
                    continue
                target = _resource_path(bundle_path, resource)
                should_generate = only_resource_ids is None or resource_id in only_resource_ids
                if only_resource_ids is not None and set(resource.derived_from).intersection(
                    only_resource_ids
                ):
                    should_generate = True
                if should_generate:
                    temporary = target.with_name(f".{target.name}.repair.tmp")
                    temporary.unlink(missing_ok=True)
                    temporary_resource = resource.model_copy(
                        update={"path": temporary.relative_to(bundle_path).as_posix()}
                    )
                    self._generate(recipe, temporary_resource, bundle_path)
                    expected_hash = resource.sha256
                    observed_hash = sha256_file(temporary)
                    if expected_hash is not None and observed_hash != expected_hash:
                        temporary.unlink(missing_ok=True)
                        raise ValueError(
                            f"derived {resource_id!r}: expected SHA256 {expected_hash}, "
                            f"found {observed_hash}"
                        )
                    target.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(temporary, target)
                elif not target.is_file():
                    raise FileNotFoundError(
                        f"retained derived resource is missing during repair: {target}"
                    )
                observed_hash = sha256_file(target)
                resolved[resource_id] = resource.model_copy(
                    update={"sha256": observed_hash, "size_bytes": target.stat().st_size}
                )
                pending.pop(resource_id)
                progressed = True
            if not progressed:
                raise ValueError(
                    "generated resource dependencies contain a cycle or unresolved source"
                )

        ordered = [resolved[item.resource_id] for item in recipe.resources]
        return recipe.model_copy(update={"resources": ordered})

    def install(
        self, recipe: ReferenceBundle, *, offline: bool = False
    ) -> InstalledReferenceBundle:
        destination = self.references_root / recipe.bundle_id
        with self._lock(recipe.bundle_id):
            self._assert_safe_reference_path(destination, label="reference bundle destination")
            if path_lexists(destination):
                raise FileExistsError(f"reference bundle already exists: {destination}; use repair")
            staging = self._staging_directory(recipe.bundle_id)
            try:
                activated_bundle = self._materialize(
                    recipe, staging, offline=offline, only_resource_ids=None
                )
                _write_bundle_manifest(activated_bundle, staging / BUNDLE_MANIFEST_NAME)
                report = validate_reference_bundle_directory(staging)
                if not report.valid:
                    invalid_messages = (item.message for item in report.resources if not item.valid)
                    messages = "; ".join([*report.errors, *invalid_messages])
                    raise ValueError(f"staged reference bundle is invalid: {messages}")
                self._assert_safe_reference_path(destination, label="reference bundle destination")
                os.replace(staging, destination)
            except BaseException:
                shutil.rmtree(staging, ignore_errors=True)
                raise
        validation = validate_reference_bundle_directory(destination)
        return InstalledReferenceBundle(activated_bundle, destination, validation)

    def validate(self, bundle_id: str) -> BundleValidationReport:
        return validate_reference_bundle_directory(self.references_root / bundle_id)

    def status(self) -> tuple[BundleValidationReport, ...]:
        if not self.references_root.is_dir():
            return ()
        return tuple(
            validate_reference_bundle_directory(path, verify_checksums=False)
            for path in sorted(self.references_root.iterdir())
            if path.is_dir() and not path.name.startswith(".")
        )

    def repair(
        self,
        recipe: ReferenceBundle,
        *,
        offline: bool = False,
    ) -> InstalledReferenceBundle:
        destination = self.references_root / recipe.bundle_id
        with self._lock(recipe.bundle_id):
            self._assert_safe_reference_path(destination, label="reference bundle destination")
            if not path_lexists(destination) or not destination.is_dir():
                raise FileNotFoundError(
                    f"reference bundle is not installed: {destination}; use install"
                )
            active_bundle = _load_reference_bundle(destination / BUNDLE_MANIFEST_NAME)
            if (
                active_bundle.bundle_id,
                active_bundle.version,
                active_bundle.genome_build,
            ) != (recipe.bundle_id, recipe.version, recipe.genome_build):
                raise ValueError(
                    "repair requires the exact installed bundle ID, version and build; "
                    "install a new bundle ID for an annotation release change"
                )
            validate_reference_recipe_contract(active_bundle, recipe)
            report = validate_reference_bundle_directory(destination)
            broken = {
                item.resource_id
                for item in report.resources
                if item.state != ResourceValidationState.OK
            }
            if not report.manifest_valid:
                raise ValueError("cannot repair a bundle whose activated manifest is invalid")
            if report.errors:
                broken.update(
                    resource.resource_id
                    for resource in recipe.resources
                    if resource.role
                    in {
                        "genome_fasta_archive",
                        "genome_fasta",
                        "fasta_index",
                        "reference_lock",
                    }
                )
            if broken:
                activated_bundle = self._materialize(
                    recipe,
                    destination,
                    offline=offline,
                    only_resource_ids=broken,
                )
                _write_bundle_manifest(activated_bundle, destination / BUNDLE_MANIFEST_NAME)
            else:
                activated_bundle = _load_reference_bundle(destination / BUNDLE_MANIFEST_NAME)
            repaired_report = validate_reference_bundle_directory(destination)
            if not repaired_report.valid:
                raise ValueError("reference repair completed but validation still fails")
        return InstalledReferenceBundle(activated_bundle, destination, repaired_report)

    def import_bundle(
        self,
        source: Path,
        *,
        authority_catalog: ReferenceCatalog | None = None,
    ) -> InstalledReferenceBundle:
        """Import an activated bundle or a pinned recipe/source tree without network access."""

        if not source.is_dir():
            raise FileNotFoundError(f"reference import directory not found: {source}")
        assert_plain_tree(source, label="reference import tree")
        recipe_path = source / "bundle.recipe.yaml"
        if not (source / BUNDLE_MANIFEST_NAME).is_file() and recipe_path.is_file():
            recipe = _load_reference_bundle(recipe_path)
            if (
                authority_catalog is not None
                and (authority := authority_catalog.find(recipe.bundle_id)) is not None
            ):
                validate_reference_recipe_contract(
                    recipe,
                    authority,
                    allow_materialized_generated_pins=False,
                )
            destination = self.references_root / recipe.bundle_id
            with self._lock(recipe.bundle_id):
                self._assert_safe_reference_path(destination, label="reference bundle destination")
                if path_lexists(destination):
                    raise FileExistsError(f"reference bundle already exists: {destination}")
                staging = self._staging_directory(recipe.bundle_id)
                try:
                    bundle = self._materialize(
                        recipe,
                        staging,
                        offline=True,
                        local_source_root=source,
                    )
                    _write_bundle_manifest(bundle, staging / BUNDLE_MANIFEST_NAME)
                    imported_report = validate_reference_bundle_directory(staging)
                    if not imported_report.valid:
                        raise ValueError("materialized recipe import failed validation")
                    self._assert_safe_reference_path(
                        destination, label="reference bundle destination"
                    )
                    os.replace(staging, destination)
                except BaseException:
                    shutil.rmtree(staging, ignore_errors=True)
                    raise
            final_report = validate_reference_bundle_directory(destination)
            return InstalledReferenceBundle(bundle, destination, final_report)

        source_report = validate_reference_bundle_directory(source)
        if not source_report.valid:
            raise ValueError("import source is not a valid activated reference bundle")
        bundle = _load_reference_bundle(source / BUNDLE_MANIFEST_NAME)
        authority = (
            authority_catalog.find(bundle.bundle_id) if authority_catalog is not None else None
        )
        if authority is not None:
            validate_reference_recipe_contract(bundle, authority)
        destination = self.references_root / bundle.bundle_id
        with self._lock(bundle.bundle_id):
            self._assert_safe_reference_path(destination, label="reference bundle destination")
            if path_lexists(destination):
                raise FileExistsError(f"reference bundle already exists: {destination}")
            staging = self._staging_directory(bundle.bundle_id)
            try:
                if authority is None:
                    shutil.rmtree(staging)
                    shutil.copytree(source, staging, symlinks=False)
                    imported_bundle = bundle
                else:
                    imported_bundle = self._materialize(
                        authority,
                        staging,
                        offline=True,
                        local_source_root=source,
                    )
                    mismatched_generated = [
                        resource.resource_id
                        for resource in imported_bundle.resources
                        if resource.generated
                        and (
                            resource.sha256 != bundle.resource(resource.resource_id).sha256
                            or resource.size_bytes
                            != bundle.resource(resource.resource_id).size_bytes
                        )
                    ]
                    if mismatched_generated:
                        raise ValueError(
                            "official reference import contains generated artifacts that do not "
                            "match deterministic catalog derivation: "
                            + ", ".join(mismatched_generated)
                        )
                    _write_bundle_manifest(imported_bundle, staging / BUNDLE_MANIFEST_NAME)
                imported_report = validate_reference_bundle_directory(staging)
                if not imported_report.valid:
                    raise ValueError("copied reference bundle failed validation")
                self._assert_safe_reference_path(destination, label="reference bundle destination")
                os.replace(staging, destination)
            except BaseException:
                shutil.rmtree(staging, ignore_errors=True)
                raise
        final_report = validate_reference_bundle_directory(destination)
        return InstalledReferenceBundle(imported_bundle, destination, final_report)


__all__ = [
    "BUNDLE_MANIFEST_NAME",
    "BundleValidationReport",
    "InstalledReferenceBundle",
    "ReferenceBundleInstaller",
    "ReferenceCatalog",
    "ResourceValidation",
    "ResourceValidationState",
    "validate_fasta_fai_consistency",
    "validate_reference_recipe_contract",
    "validate_reference_bundle_directory",
]

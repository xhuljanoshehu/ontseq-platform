"""Activate the curated GRCh38 panel, knowledge and analysis-profile resources.

The reference installer deliberately owns only ``references/``.  This module is the small
transaction boundary that follows a successful reference installation: it copies the curated
repository resources byte-for-byte into staging, compiles the panel derivatives against the
installed annotation SQLite cache, validates every checksum, and exposes profiles last.

There is intentionally no build inference, liftover or cross-build fallback here.  The current
bootstrap contract is the single, pinned GRCh38 resource family.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager, suppress
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path, PurePosixPath
from typing import cast

from pydantic import ValidationError

from .io import load_mapping
from .models import (
    AnalysisProfile,
    AssayMode,
    GenomeBuild,
    KnowledgeBundle,
    PanelBundle,
    ReferenceBundle,
    ReferenceDictionaryContract,
    ResourceBundle,
)
from .panel_compiler import MaterializedPanelSummary, materialize_and_pin_panel_derivatives
from .path_safety import (
    assert_plain_tree,
    assert_safe_descendant,
    is_link_like,
    path_lexists,
)
from .reference import sha256_file
from .reference_catalog import validate_reference_bundle_directory
from .resource_registry import resource_root_from_environment

REFERENCE_BUNDLE_ID = "GRCh38_GENCODE50_MANE1.5_v1"
PANEL_BUNDLE_ID = "AML_AS_111_GRCh38_v1"
KNOWLEDGE_BUNDLE_ID = "HEMATOLOGY_v3"
PROFILE_CONTRACTS: dict[
    str,
    tuple[AssayMode, ReferenceDictionaryContract, str | None],
] = {
    "AML_LCWGS_GRCh38": (
        AssayMode.LOW_COVERAGE_WGS,
        ReferenceDictionaryContract.EXACT_FULL,
        None,
    ),
    "AML_AS_111_GRCh38": (
        AssayMode.ADAPTIVE_SAMPLING,
        ReferenceDictionaryContract.EXACT_FULL,
        PANEL_BUNDLE_ID,
    ),
    "AML_LCWGS_GRCh38_CANONICAL25": (
        AssayMode.LOW_COVERAGE_WGS,
        ReferenceDictionaryContract.GRCH38_CANONICAL_25,
        None,
    ),
    "AML_AS_111_GRCh38_CANONICAL25": (
        AssayMode.ADAPTIVE_SAMPLING,
        ReferenceDictionaryContract.GRCH38_CANONICAL_25,
        PANEL_BUNDLE_ID,
    ),
}
PROFILE_IDS = tuple(PROFILE_CONTRACTS)


class ResourceBootstrapError(ValueError):
    """Curated resources cannot be activated without weakening a pinned invariant."""


@dataclass(frozen=True)
class BootstrapResult:
    reference_bundle_id: str
    panel_bundle_id: str
    knowledge_bundle_id: str
    profile_ids: tuple[str, ...]
    panel_summary: MaterializedPanelSummary
    activated_paths: tuple[Path, ...]
    repaired_paths: tuple[Path, ...]
    already_active_paths: tuple[Path, ...]


class _DestinationState(StrEnum):
    ABSENT = "absent"
    IDENTICAL = "identical"
    REPLACE = "replace"


class _RollbackIncompleteError(Exception):
    """Signal that transaction-local recovery data must not be cleaned up."""

    def __init__(self, errors: tuple[str, ...]) -> None:
        super().__init__("; ".join(errors))
        self.errors = errors


def default_packaged_config_root() -> Path:
    """Return the repository's curated config root or fail with an actionable message.

    Callers that install ONTSeq without its source tree should pass an explicit config root.
    Keeping this lookup explicit prevents an arbitrary working-directory ``configs`` folder from
    being trusted as packaged provenance.
    """

    candidate = Path(__file__).resolve().parents[2] / "configs"
    if not candidate.is_dir():
        raise FileNotFoundError(
            "curated ONTSeq configs are not present beside this installation; pass "
            "packaged_config_root explicitly"
        )
    return candidate


def _bundle_file(bundle_directory: Path, relative_path: str) -> Path:
    root = bundle_directory.resolve()
    candidate = bundle_directory.joinpath(*PurePosixPath(relative_path).parts).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ResourceBootstrapError(
            f"bundle resource escapes its directory: {relative_path!r}"
        ) from exc
    return candidate


def _copy_tree_byte_exact(source: Path, destination: Path) -> None:
    try:
        assert_plain_tree(source, label="curated resource tree")
    except (OSError, ValueError) as exc:
        raise ResourceBootstrapError(
            f"curated resource directory is absent or unsafe: {source}: {exc}"
        ) from exc
    shutil.copytree(source, destination, copy_function=shutil.copyfile)


def _validate_bundle_files(bundle: ResourceBundle, directory: Path) -> None:
    for resource in bundle.resources:
        path = _bundle_file(directory, resource.path)
        if not path.is_file():
            if resource.required:
                raise ResourceBootstrapError(
                    f"required resource {bundle.bundle_id}:{resource.resource_id} is missing"
                )
            continue
        if resource.sha256 is None:
            raise ResourceBootstrapError(
                f"resource {bundle.bundle_id}:{resource.resource_id} is not checksum-pinned"
            )
        observed_sha256 = sha256_file(path)
        if observed_sha256 != resource.sha256:
            raise ResourceBootstrapError(
                f"resource {bundle.bundle_id}:{resource.resource_id} checksum mismatch: "
                f"expected {resource.sha256}, observed {observed_sha256}"
            )
        if resource.size_bytes is not None and path.stat().st_size != resource.size_bytes:
            raise ResourceBootstrapError(
                f"resource {bundle.bundle_id}:{resource.resource_id} size mismatch: "
                f"expected {resource.size_bytes}, observed {path.stat().st_size}"
            )


def _tree_digest(directory: Path) -> str:
    """Hash every plain entry and byte; reject links, special entries and type changes."""

    assert_plain_tree(directory, label="resource activation tree")
    digest = hashlib.sha256()
    for path in sorted(directory.rglob("*")):
        relative = path.relative_to(directory).as_posix().encode("utf-8")
        if is_link_like(path):
            raise ResourceBootstrapError(f"resource activation tree contains a link: {path}")
        if path.is_dir():
            kind = b"D"
        elif path.is_file():
            kind = b"F"
        else:
            raise ResourceBootstrapError(
                f"resource activation tree contains a special entry: {path}"
            )
        digest.update(kind)
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        if kind == b"D":
            continue
        digest.update(path.stat().st_size.to_bytes(8, "big"))
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


class GRCh38ResourceBootstrapper:
    """Stage and activate the curated GRCh38 non-reference resource family."""

    def __init__(
        self,
        resource_root: Path | str | None = None,
        *,
        packaged_config_root: Path | None = None,
    ) -> None:
        self.resource_root = resource_root_from_environment(resource_root).resolve()
        self.packaged_config_root = (
            packaged_config_root.resolve()
            if packaged_config_root is not None
            else default_packaged_config_root().resolve()
        )

    @contextmanager
    def _lock(self) -> Iterator[None]:
        self.resource_root.mkdir(parents=True, exist_ok=True)
        lock_path = self.resource_root / ".grch38-resource-bootstrap.lock"
        try:
            assert_safe_descendant(self.resource_root, lock_path, label="resource bootstrap lock")
        except ValueError as exc:
            raise ResourceBootstrapError(str(exc)) from exc
        descriptor: int | None = None
        try:
            descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
            os.write(descriptor, str(os.getpid()).encode("ascii"))
            os.close(descriptor)
            descriptor = None
        except FileExistsError as exc:
            raise RuntimeError("another GRCh38 resource activation is already running") from exc
        try:
            yield
        finally:
            if descriptor is not None:
                with suppress(OSError):
                    os.close(descriptor)
            lock_path.unlink(missing_ok=True)

    def _reference_annotation_cache(self) -> Path:
        reference_directory = self.resource_root / "references" / REFERENCE_BUNDLE_ID
        report = validate_reference_bundle_directory(reference_directory)
        if not report.valid:
            messages = [*report.errors]
            messages.extend(item.message for item in report.resources if not item.valid)
            rendered = "; ".join(messages) or "bundle manifest is absent"
            raise ResourceBootstrapError(
                f"reference bundle {REFERENCE_BUNDLE_ID!r} must be installed and valid before "
                f"resource activation: {rendered}"
            )
        try:
            reference = ReferenceBundle.model_validate(
                load_mapping(reference_directory / "bundle.yaml")
            )
        except (OSError, ValueError, ValidationError) as exc:
            raise ResourceBootstrapError("installed reference manifest is invalid") from exc
        if reference.bundle_id != REFERENCE_BUNDLE_ID:
            raise ResourceBootstrapError("installed reference directory and bundle ID disagree")
        if reference.genome_build != GenomeBuild.GRCH38:
            raise ResourceBootstrapError("resource activation supports GRCh38 references only")
        if reference.annotation_cache_resource_id is None:
            raise ResourceBootstrapError("installed reference does not declare an annotation cache")
        cache_resource = reference.resource(reference.annotation_cache_resource_id)
        cache_path = _bundle_file(reference_directory, cache_resource.path)
        if not cache_path.is_file():
            raise ResourceBootstrapError(f"installed annotation cache is missing: {cache_path}")
        return cache_path

    def _load_curated_profiles(self) -> tuple[AnalysisProfile, ...]:
        profiles: list[AnalysisProfile] = []
        for profile_id, expected in PROFILE_CONTRACTS.items():
            path = self.packaged_config_root / "profiles" / f"{profile_id}.yaml"
            try:
                profile = AnalysisProfile.model_validate(load_mapping(path))
            except (OSError, ValueError, ValidationError) as exc:
                raise ResourceBootstrapError(f"curated profile is invalid: {path}") from exc
            if profile.profile_id != profile_id:
                raise ResourceBootstrapError(
                    f"curated profile filename and ID disagree: {path.name}"
                )
            if profile.genome_build != GenomeBuild.GRCH38:
                raise ResourceBootstrapError(f"curated profile {profile_id!r} is not GRCh38")
            if profile.reference_bundle != REFERENCE_BUNDLE_ID:
                raise ResourceBootstrapError(
                    f"curated profile {profile_id!r} does not pin {REFERENCE_BUNDLE_ID!r}"
                )
            if profile.knowledge_bundle != KNOWLEDGE_BUNDLE_ID:
                raise ResourceBootstrapError(
                    f"curated profile {profile_id!r} does not pin {KNOWLEDGE_BUNDLE_ID!r}"
                )
            expected_assay, expected_dictionary, expected_panel = expected
            if (
                profile.assay_mode != expected_assay
                or profile.reference_dictionary_contract != expected_dictionary
                or profile.panel_bundle != expected_panel
            ):
                raise ResourceBootstrapError(
                    f"curated profile {profile_id!r} does not match its published assay, "
                    "dictionary and panel contract"
                )
            profiles.append(profile)
        return tuple(profiles)

    @staticmethod
    def _require_curated_bundle(
        directory: Path,
        model: type[PanelBundle] | type[KnowledgeBundle],
        expected_id: str,
    ) -> PanelBundle | KnowledgeBundle:
        try:
            bundle = model.model_validate(load_mapping(directory / "bundle.yaml"))
        except (OSError, ValueError, ValidationError) as exc:
            raise ResourceBootstrapError(f"curated bundle is invalid: {directory}") from exc
        if bundle.bundle_id != expected_id:
            raise ResourceBootstrapError(
                f"curated directory {directory.name!r} declares bundle {bundle.bundle_id!r}"
            )
        if bundle.genome_build != GenomeBuild.GRCH38:
            raise ResourceBootstrapError(f"curated bundle {expected_id!r} is not GRCh38")
        return bundle

    @staticmethod
    def _destination_state(
        staged: Path,
        destination: Path,
        *,
        resource_root: Path,
        directory: bool,
        repair_existing: bool,
    ) -> _DestinationState:
        try:
            assert_safe_descendant(
                resource_root, destination, label="resource activation destination"
            )
        except ValueError as exc:
            raise ResourceBootstrapError(str(exc)) from exc
        if not path_lexists(destination):
            return _DestinationState.ABSENT
        safe_type = destination.is_dir() if directory else destination.is_file()
        label = "resource activation" if directory else "profile activation"
        if not safe_type or is_link_like(destination):
            raise FileExistsError(f"{label} destination is unsafe: {destination}")
        identical = (
            _tree_digest(staged) == _tree_digest(destination)
            if directory
            else staged.read_bytes() == destination.read_bytes()
        )
        if identical:
            return _DestinationState.IDENTICAL
        if repair_existing:
            return _DestinationState.REPLACE
        kind = "bundle" if directory else "profile"
        raise FileExistsError(f"{label} would overwrite a different {kind}: {destination}")

    @staticmethod
    def _apply_destinations(
        plans: tuple[tuple[Path, Path, bool], ...],
        *,
        resource_root: Path,
        transaction: Path,
        repair_existing: bool,
    ) -> tuple[list[Path], list[Path], list[Path]]:
        """Activate or repair a fully validated family, rolling back a partial replacement.

        Every destination is classified before the first mutation. Divergent trees remain a hard
        error for normal activation. Explicit repair moves each old path into transaction-local
        backup before atomically publishing its staged replacement; profiles are ordered after
        bundles by the caller. If a later move fails, completed moves are reversed.
        """

        states = tuple(
            GRCh38ResourceBootstrapper._destination_state(
                staged,
                destination,
                resource_root=resource_root,
                directory=directory,
                repair_existing=repair_existing,
            )
            for staged, destination, directory in plans
        )
        activated: list[Path] = []
        repaired: list[Path] = []
        already_active: list[Path] = []
        completed: list[tuple[_DestinationState, Path, Path | None]] = []
        try:
            for index, ((staged, destination, _directory), state) in enumerate(
                zip(plans, states, strict=True)
            ):
                if state == _DestinationState.IDENTICAL:
                    already_active.append(destination)
                    continue
                try:
                    assert_safe_descendant(
                        resource_root, destination, label="resource activation destination"
                    )
                except ValueError as exc:
                    raise ResourceBootstrapError(str(exc)) from exc
                destination.parent.mkdir(parents=True, exist_ok=True)
                if state == _DestinationState.REPLACE:
                    replacement_backup = transaction / "backups" / f"{index:02d}"
                    replacement_backup.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(destination, replacement_backup)
                    completed.append((state, destination, replacement_backup))
                    os.replace(staged, destination)
                    repaired.append(destination)
                    continue
                os.replace(staged, destination)
                completed.append((state, destination, None))
                activated.append(destination)
        except BaseException as original_error:
            rollback_errors: list[str] = []
            rollback_root = transaction / "rollback-new"
            for index, (state, destination, backup_path) in enumerate(reversed(completed)):
                try:
                    if destination.exists():
                        displaced = rollback_root / f"{index:02d}"
                        displaced.parent.mkdir(parents=True, exist_ok=True)
                        os.replace(destination, displaced)
                    if state == _DestinationState.REPLACE and backup_path is not None:
                        os.replace(backup_path, destination)
                except OSError as rollback_error:
                    rollback_errors.append(f"{destination}: {rollback_error}")
            if rollback_errors:
                raise _RollbackIncompleteError(tuple(rollback_errors)) from original_error
            raise
        return activated, repaired, already_active

    def _apply(self, *, repair_existing: bool) -> BootstrapResult:
        """Compile, validate and atomically expose the curated GRCh38 resources.

        Bundle directories are moved into place atomically. Profiles are activated last, so a
        failure can leave at most an unreferenced valid bundle and never a profile that selects a
        partial panel. Explicit repair may replace divergent pinned resources, with transaction-
        local backups and rollback; normal activation still refuses every divergent collision.
        """

        with self._lock():
            annotation_cache = self._reference_annotation_cache()
            profiles = self._load_curated_profiles()
            staging_root = self.resource_root / ".bootstrap-staging"
            try:
                assert_safe_descendant(
                    self.resource_root, staging_root, label="resource bootstrap staging"
                )
            except ValueError as exc:
                raise ResourceBootstrapError(str(exc)) from exc
            staging_root.mkdir(parents=True, exist_ok=True)
            raw_transaction = tempfile.mkdtemp(prefix="grch38.", dir=staging_root)
            transaction = Path(raw_transaction)
            retain_transaction = False
            activated: list[Path] = []
            repaired: list[Path] = []
            already_active: list[Path] = []
            try:
                staged_knowledge = transaction / "knowledge" / KNOWLEDGE_BUNDLE_ID
                staged_panel = transaction / "panels" / PANEL_BUNDLE_ID
                staged_profiles = transaction / "profiles"
                _copy_tree_byte_exact(
                    self.packaged_config_root / "knowledge_bundles" / KNOWLEDGE_BUNDLE_ID,
                    staged_knowledge,
                )
                _copy_tree_byte_exact(
                    self.packaged_config_root / "panels" / PANEL_BUNDLE_ID,
                    staged_panel,
                )
                staged_profiles.mkdir(parents=True)
                for profile in profiles:
                    source = self.packaged_config_root / "profiles" / f"{profile.profile_id}.yaml"
                    shutil.copyfile(source, staged_profiles / source.name)

                knowledge = cast(
                    KnowledgeBundle,
                    self._require_curated_bundle(
                        staged_knowledge, KnowledgeBundle, KNOWLEDGE_BUNDLE_ID
                    ),
                )
                panel = cast(
                    PanelBundle,
                    self._require_curated_bundle(staged_panel, PanelBundle, PANEL_BUNDLE_ID),
                )
                if (
                    panel.resource(panel.analysis_roi_resource_id).derived_from.count(
                        f"{REFERENCE_BUNDLE_ID}:annotation_cache"
                    )
                    != 1
                ):
                    raise ResourceBootstrapError(
                        "curated panel ROI does not pin the active reference annotation cache"
                    )
                if (
                    panel.resource(panel.transcript_cache_resource_id).derived_from.count(
                        f"{REFERENCE_BUNDLE_ID}:annotation_cache"
                    )
                    != 1
                ):
                    raise ResourceBootstrapError(
                        "curated panel transcript cache does not pin the active reference cache"
                    )
                panel_summary = materialize_and_pin_panel_derivatives(
                    staged_panel, annotation_cache
                )
                panel = PanelBundle.model_validate(load_mapping(staged_panel / "bundle.yaml"))
                _validate_bundle_files(knowledge, staged_knowledge)
                _validate_bundle_files(panel, staged_panel)

                destinations = (
                    (
                        staged_knowledge,
                        self.resource_root / "knowledge" / KNOWLEDGE_BUNDLE_ID,
                        True,
                    ),
                    (
                        staged_panel,
                        self.resource_root / "panels" / PANEL_BUNDLE_ID,
                        True,
                    ),
                )
                profile_destinations = tuple(
                    (
                        staged_profiles / f"{profile.profile_id}.yaml",
                        self.resource_root / "profiles" / f"{profile.profile_id}.yaml",
                        False,
                    )
                    for profile in profiles
                )
                activated, repaired, already_active = self._apply_destinations(
                    destinations + profile_destinations,
                    resource_root=self.resource_root,
                    transaction=transaction,
                    repair_existing=repair_existing,
                )
            except _RollbackIncompleteError as exc:
                retain_transaction = True
                raise ResourceBootstrapError(
                    "resource-family repair failed and rollback was incomplete; recovery data "
                    f"was retained at {transaction}: {'; '.join(exc.errors)}"
                ) from exc.__cause__
            finally:
                if not retain_transaction:
                    shutil.rmtree(transaction, ignore_errors=True)

        return BootstrapResult(
            reference_bundle_id=REFERENCE_BUNDLE_ID,
            panel_bundle_id=PANEL_BUNDLE_ID,
            knowledge_bundle_id=KNOWLEDGE_BUNDLE_ID,
            profile_ids=PROFILE_IDS,
            panel_summary=panel_summary,
            activated_paths=tuple(activated),
            repaired_paths=tuple(repaired),
            already_active_paths=tuple(already_active),
        )

    def activate(self) -> BootstrapResult:
        """Activate absent resources and verify existing resources byte-for-byte."""

        return self._apply(repair_existing=False)

    def repair(self) -> BootstrapResult:
        """Restore the pinned panel, knowledge and profiles without manual deletion."""

        return self._apply(repair_existing=True)


__all__ = [
    "BootstrapResult",
    "GRCh38ResourceBootstrapper",
    "KNOWLEDGE_BUNDLE_ID",
    "PANEL_BUNDLE_ID",
    "PROFILE_CONTRACTS",
    "PROFILE_IDS",
    "REFERENCE_BUNDLE_ID",
    "ResourceBootstrapError",
    "default_packaged_config_root",
]

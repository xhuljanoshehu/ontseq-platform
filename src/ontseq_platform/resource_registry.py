"""Manifest-driven resource discovery and fail-closed profile resolution."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Generic, TypeVar

from pydantic import ValidationError

from .io import load_mapping
from .models import (
    AnalysisProfile,
    GenomeBuild,
    KnowledgeBundle,
    PanelBundle,
    PanelCoordinateMappingLock,
    PanelResolutionSummary,
    ReferenceBundle,
    ReferenceDictionaryContract,
    ReferenceLock,
    ResolvedResourceContext,
    ResourceBundle,
)
from .reference import (
    GRCH37_UCSC_HG19_CANONICAL_25_FAI_ROLE,
    GRCH37_UCSC_HG19_CANONICAL_25_FASTA_ROLE,
    reference_lock_for_dictionary_contract,
    reference_lock_from_fai,
    sha256_file,
    validate_grch37_ucsc_hg19_canonical_25,
)

DEFAULT_RESOURCE_ROOT = Path("/opt/ontseq")
RESOURCE_ROOT_ENV = "ONTSEQ_RESOURCE_ROOT"

BundleT = TypeVar("BundleT", bound=ResourceBundle, covariant=True)


def _bed_interval_count(path: Path) -> int:
    """Count non-comment BED rows after basic coordinate validation."""

    count = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path.name} line {line_number} has fewer than three BED fields")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path.name} line {line_number} has non-integer BED coordinates"
                ) from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path.name} line {line_number} has an invalid BED interval")
            count += 1
    return count


@dataclass(frozen=True)
class RegistryDiagnostic:
    path: Path
    code: str
    message: str


@dataclass(frozen=True)
class _RegisteredBundle(Generic[BundleT]):
    bundle: BundleT
    directory: Path


def resource_root_from_environment(explicit: Path | str | None = None) -> Path:
    """Resolve ``--resource-root`` semantics without coupling the registry to a CLI."""

    if explicit is not None:
        return Path(explicit).expanduser()
    configured = os.environ.get(RESOURCE_ROOT_ENV)
    return Path(configured).expanduser() if configured else DEFAULT_RESOURCE_ROOT


class ResourceRegistry:
    """Index one build's activated bundle manifests below a resource root.

    Separate registry instances can target separate builds. A directory without a valid
    ``bundle.yaml`` is never interpreted from its filename or loose contents.
    """

    def __init__(
        self,
        resource_root: Path | str | None = None,
        *,
        active_build: GenomeBuild = GenomeBuild.GRCH38,
    ) -> None:
        self.resource_root = resource_root_from_environment(resource_root).resolve()
        self.active_build = active_build
        self._references: dict[str, _RegisteredBundle[ReferenceBundle]] = {}
        self._panels: dict[str, _RegisteredBundle[PanelBundle]] = {}
        self._knowledge: dict[str, _RegisteredBundle[KnowledgeBundle]] = {}
        self._profiles: dict[str, AnalysisProfile] = {}
        self._diagnostics: list[RegistryDiagnostic] = []
        self.refresh()

    @property
    def references(self) -> MappingProxyType[str, ReferenceBundle]:
        return MappingProxyType(
            {bundle_id: registered.bundle for bundle_id, registered in self._references.items()}
        )

    @property
    def panels(self) -> MappingProxyType[str, PanelBundle]:
        return MappingProxyType(
            {bundle_id: registered.bundle for bundle_id, registered in self._panels.items()}
        )

    @property
    def knowledge(self) -> MappingProxyType[str, KnowledgeBundle]:
        return MappingProxyType(
            {bundle_id: registered.bundle for bundle_id, registered in self._knowledge.items()}
        )

    @property
    def profiles(self) -> MappingProxyType[str, AnalysisProfile]:
        return MappingProxyType(dict(self._profiles))

    @property
    def diagnostics(self) -> tuple[RegistryDiagnostic, ...]:
        return tuple(self._diagnostics)

    def refresh(self) -> None:
        self._references.clear()
        self._panels.clear()
        self._knowledge.clear()
        self._profiles.clear()
        self._diagnostics.clear()
        self._scan_bundle_category("references", ReferenceBundle, self._references)
        self._scan_bundle_category("panels", PanelBundle, self._panels)
        self._scan_bundle_category("knowledge", KnowledgeBundle, self._knowledge)
        self._scan_profiles()

    def _diagnose(self, path: Path, code: str, message: str) -> None:
        self._diagnostics.append(RegistryDiagnostic(path=path, code=code, message=message))

    def _scan_bundle_category(
        self,
        category: str,
        model: type[BundleT],
        destination: dict[str, _RegisteredBundle[BundleT]],
    ) -> None:
        category_root = self.resource_root / category
        if not category_root.is_dir():
            return
        for directory in sorted(category_root.iterdir(), key=lambda item: item.name):
            if not directory.is_dir():
                continue
            manifest_path = directory / "bundle.yaml"
            if not manifest_path.is_file():
                continue
            try:
                bundle = model.model_validate(load_mapping(manifest_path))
            except (OSError, ValueError, ValidationError) as exc:
                self._diagnose(manifest_path, "invalid_manifest", str(exc))
                continue
            if bundle.bundle_id != directory.name:
                self._diagnose(
                    manifest_path,
                    "bundle_id_directory_mismatch",
                    f"bundle_id {bundle.bundle_id!r} does not match directory {directory.name!r}",
                )
                continue
            if bundle.genome_build not in {None, self.active_build}:
                observed_build = bundle.genome_build
                assert observed_build is not None
                self._diagnose(
                    manifest_path,
                    "inactive_build",
                    f"{observed_build.value} bundle is not active in the "
                    f"{self.active_build.value} registry",
                )
                continue
            destination[bundle.bundle_id] = _RegisteredBundle(
                bundle=bundle,
                directory=directory.resolve(),
            )

    def _scan_profiles(self) -> None:
        profiles_root = self.resource_root / "profiles"
        if not profiles_root.is_dir():
            return
        for path in sorted(profiles_root.glob("*.yaml"), key=lambda item: item.name):
            try:
                profile = AnalysisProfile.model_validate(load_mapping(path))
            except (OSError, ValueError, ValidationError) as exc:
                self._diagnose(path, "invalid_profile", str(exc))
                continue
            if path.stem != profile.profile_id:
                self._diagnose(
                    path,
                    "profile_id_filename_mismatch",
                    f"profile_id {profile.profile_id!r} does not match filename {path.stem!r}",
                )
                continue
            if profile.genome_build != self.active_build:
                self._diagnose(
                    path,
                    "inactive_build",
                    f"{profile.genome_build.value} profile is not active in the "
                    f"{self.active_build.value} registry",
                )
                continue
            self._profiles[profile.profile_id] = profile

    def resolve_profile(
        self,
        profile_id: str,
        *,
        verify_files: bool = True,
    ) -> ResolvedResourceContext:
        """Resolve one profile to absolute paths after enforcing every build invariant.

        Required files, declared sizes and checksum pins are always checked.  ``verify_files``
        controls only the expensive act of reading every resource byte and recomputing SHA256.
        This distinction keeps an interactive run preflight fast even when the reference FASTA
        is several gigabytes, while ``ontseq references validate`` can still perform the full
        integrity audit explicitly.
        """

        try:
            profile = self._profiles[profile_id]
        except KeyError as exc:
            raise KeyError(
                f"analysis profile {profile_id!r} is not active for {self.active_build.value}"
            ) from exc
        reference = self._require_bundle(self._references, profile.reference_bundle, "reference")
        knowledge = self._require_bundle(self._knowledge, profile.knowledge_bundle, "knowledge")
        panel = (
            self._require_bundle(self._panels, profile.panel_bundle, "panel")
            if profile.panel_bundle is not None
            else None
        )

        self._require_matching_build(profile, reference.bundle)
        self._require_matching_build(profile, knowledge.bundle)
        if panel is not None:
            self._require_matching_build(profile, panel.bundle)
            if panel.bundle.assay_mode != profile.assay_mode:
                raise ValueError(
                    f"profile {profile.profile_id!r} uses assay {profile.assay_mode.value}, but "
                    f"panel {panel.bundle.bundle_id!r} declares {panel.bundle.assay_mode.value}"
                )
            if (
                profile.reference_dictionary_contract
                not in panel.bundle.reference_dictionary_contracts
            ):
                raise ValueError(
                    f"profile {profile.profile_id!r} uses dictionary contract "
                    f"{profile.reference_dictionary_contract.value!r}, but panel "
                    f"{panel.bundle.bundle_id!r} is not declared compatible with it"
                )

        paths: dict[str, str] = {}
        checksums: dict[str, str] = {}
        releases: dict[str, str] = {}
        registered: list[tuple[str, _RegisteredBundle[ResourceBundle]]] = [
            ("reference", reference),
            ("knowledge", knowledge),
        ]
        if panel is not None:
            registered.append(("panel", panel))
        self._validate_cross_bundle_derivations([item.bundle for _, item in registered])
        for namespace, item in registered:
            self._resolve_bundle_files(
                namespace,
                item,
                paths=paths,
                checksums=checksums,
                releases=releases,
                verify_files=verify_files,
            )
        self._validate_reference_lock(profile, reference, paths)
        panel_resolution = (
            self._resolve_panel_resolution(panel, paths, checksums) if panel is not None else None
        )

        return ResolvedResourceContext(
            profile_id=profile.profile_id,
            profile_version=profile.version,
            genome_build=profile.genome_build,
            reference_dictionary_contract=profile.reference_dictionary_contract,
            reference_bundle_id=reference.bundle.bundle_id,
            reference_bundle_version=reference.bundle.version,
            panel_bundle_id=panel.bundle.bundle_id if panel else None,
            panel_bundle_version=panel.bundle.version if panel else None,
            panel_resolution=panel_resolution,
            knowledge_bundle_id=knowledge.bundle.bundle_id,
            knowledge_bundle_version=knowledge.bundle.version,
            resource_root=str(self.resource_root),
            resource_paths=paths,
            resource_checksums=checksums,
            resource_releases=releases,
        )

    @staticmethod
    def _resolve_panel_resolution(
        panel: _RegisteredBundle[PanelBundle],
        paths: dict[str, str],
        checksums: dict[str, str],
    ) -> PanelResolutionSummary:
        """Bind reportable panel-resolution state to the already validated resource set."""

        bundle = panel.bundle
        selection = bundle.resource(bundle.selection_panel_resource_id)
        analysis_roi = bundle.resource(bundle.analysis_roi_resource_id)
        assert selection.sha256 is not None
        assert analysis_roi.sha256 is not None
        selection_count = _bed_interval_count(Path(paths["panel.selection_panel_buffered"]))
        analysis_roi_count = _bed_interval_count(Path(paths["panel.analysis_roi_unbuffered"]))
        if selection_count < 1:
            raise ValueError("resolved panel selection contains no intervals")
        lock_id = bundle.coordinate_mapping_lock_resource_id
        if lock_id is None:
            return PanelResolutionSummary(
                mapping_status="native_build_not_required",
                source_genome_build=bundle.genome_build,
                target_genome_build=bundle.genome_build,
                selection_interval_count=selection_count,
                analysis_roi_interval_count=analysis_roi_count,
                unresolved_target_labels=bundle.unresolved_targets,
                selection_panel_sha256=selection.sha256,
                analysis_roi_sha256=analysis_roi.sha256,
            )

        source_id = bundle.coordinate_mapping_source_resource_id
        output_id = bundle.coordinate_mapping_output_resource_id
        unmapped_id = bundle.coordinate_mapping_unmapped_resource_id
        roundtrip_id = bundle.coordinate_mapping_roundtrip_output_resource_id
        roundtrip_unmapped_id = bundle.coordinate_mapping_roundtrip_unmapped_resource_id
        assert source_id is not None
        assert output_id is not None
        assert unmapped_id is not None
        assert roundtrip_id is not None
        assert roundtrip_unmapped_id is not None
        lock = PanelCoordinateMappingLock.model_validate(
            load_mapping(Path(paths["panel.coordinate_mapping_lock"]))
        )
        source = bundle.resource(source_id)
        output = bundle.resource(output_id)
        unmapped = bundle.resource(unmapped_id)
        roundtrip = bundle.resource(roundtrip_id)
        roundtrip_unmapped = bundle.resource(roundtrip_unmapped_id)
        expected_values = {
            "target genome build": (lock.target_genome_build, bundle.genome_build),
            "source selection SHA256": (lock.source_selection_sha256, source.sha256),
            "mapped output SHA256": (lock.mapped_output_sha256, output.sha256),
            "unmapped output SHA256": (lock.unmapped_output_sha256, unmapped.sha256),
            "final selection SHA256": (lock.final_selection_sha256, selection.sha256),
            "source panel bundle ID": (
                lock.source_panel_bundle_id,
                bundle.coordinate_mapping_origin_bundle_id,
            ),
            "source panel resource ID": (
                lock.source_panel_resource_id,
                bundle.coordinate_mapping_origin_resource_id,
            ),
            "roundtrip mapped output SHA256": (
                lock.roundtrip_mapped_output_sha256,
                roundtrip.sha256,
            ),
            "roundtrip unmapped output SHA256": (
                lock.roundtrip_unmapped_output_sha256,
                roundtrip_unmapped.sha256,
            ),
        }
        mismatches = [
            label
            for label, (lock_value, resource_value) in expected_values.items()
            if lock_value != resource_value
        ]
        if mismatches:
            raise ValueError(
                "panel coordinate mapping lock does not match resolved resources: "
                + ", ".join(mismatches)
            )
        if _bed_interval_count(Path(paths["panel.coordinate_mapping_roundtrip_output"])) != (
            lock.roundtrip_mapped_interval_count
        ):
            raise ValueError("panel roundtrip mapped interval count disagrees with its lock")
        if selection_count != lock.mapped_interval_count:
            raise ValueError("panel final selection interval count disagrees with its lock")
        lock_checksum = checksums["panel.coordinate_mapping_lock"]
        return PanelResolutionSummary(
            mapping_status=lock.status,
            mapping_id=lock.mapping_id,
            mapping_method=lock.mapping_method,
            mapping_tool=lock.mapping_tool,
            source_genome_build=lock.source_genome_build,
            target_genome_build=lock.target_genome_build,
            source_panel_bundle_id=lock.source_panel_bundle_id,
            source_panel_resource_id=lock.source_panel_resource_id,
            selection_interval_count=selection_count,
            analysis_roi_interval_count=analysis_roi_count,
            source_interval_count=lock.source_interval_count,
            mapped_interval_count=lock.mapped_interval_count,
            unmapped_target_labels=lock.unmapped_target_labels,
            roundtrip_mapped_interval_count=lock.roundtrip_mapped_interval_count,
            reciprocal_exact_interval_count=lock.reciprocal_exact_interval_count,
            roundtrip_review_required_target_labels=(lock.roundtrip_review_required_target_labels),
            unresolved_target_labels=bundle.unresolved_targets,
            selection_panel_sha256=selection.sha256,
            analysis_roi_sha256=analysis_roi.sha256,
            coordinate_mapping_lock_sha256=lock_checksum,
            roundtrip_mapping_output_sha256=roundtrip.sha256,
            roundtrip_mapping_unmapped_sha256=roundtrip_unmapped.sha256,
            runtime_mapping_allowed=lock.runtime_mapping_allowed,
        )

    @staticmethod
    def _validate_reference_lock(
        profile: AnalysisProfile,
        reference: _RegisteredBundle[ReferenceBundle],
        paths: dict[str, str],
    ) -> None:
        lock_path = Path(paths["reference.reference_lock"])
        reference_lock = ReferenceLock.model_validate(load_mapping(lock_path))
        if reference_lock.reference_id != reference.bundle.bundle_id:
            raise ValueError(
                f"reference lock ID {reference_lock.reference_id!r} does not match bundle "
                f"{reference.bundle.bundle_id!r}"
            )
        if reference_lock.genome_build != profile.genome_build:
            raise ValueError(
                f"reference lock is {reference_lock.genome_build.value}, but profile "
                f"{profile.profile_id!r} requires {profile.genome_build.value}"
            )
        fai_resource = reference.bundle.resource(reference.bundle.fai_resource_id)
        if reference_lock.source_fai_sha256 != fai_resource.sha256:
            raise ValueError(
                "reference lock source_fai_sha256 does not match the pinned FASTA index"
            )
        if (
            profile.reference_dictionary_contract
            == ReferenceDictionaryContract.GRCH37_UCSC_HG19_CANONICAL_25
        ):
            fasta_key = f"reference.{GRCH37_UCSC_HG19_CANONICAL_25_FASTA_ROLE}"
            fai_key = f"reference.{GRCH37_UCSC_HG19_CANONICAL_25_FAI_ROLE}"
            if fasta_key not in paths or fai_key not in paths:
                raise ValueError(
                    "UCSC hg19 Canonical-25 profile requires its exact analysis FASTA and index"
                )
            contract_fai = Path(paths[fai_key])
            contract_lock = reference_lock_from_fai(
                contract_fai,
                reference_id=reference.bundle.bundle_id,
                genome_build=GenomeBuild.GRCH37,
                allow_extra_contigs=False,
            )
            validate_grch37_ucsc_hg19_canonical_25(
                (contig.name, contig.length) for contig in contract_lock.contigs
            )
            return
        reference_lock_for_dictionary_contract(
            reference_lock,
            profile.reference_dictionary_contract,
        )

    @staticmethod
    def _validate_cross_bundle_derivations(bundles: list[ResourceBundle]) -> None:
        selected: dict[str, ResourceBundle] = {}
        for bundle in bundles:
            if bundle.bundle_id in selected:
                raise ValueError(
                    f"profile resolves duplicate bundle ID {bundle.bundle_id!r}; bundle IDs "
                    "must be globally unambiguous"
                )
            selected[bundle.bundle_id] = bundle
        for bundle in bundles:
            for resource in bundle.resources:
                for reference in resource.derived_from:
                    if ":" not in reference:
                        continue
                    source_bundle_id, source_resource_id = reference.split(":", maxsplit=1)
                    try:
                        source_bundle = selected[source_bundle_id]
                    except KeyError as exc:
                        raise ValueError(
                            f"resource {bundle.bundle_id}:{resource.resource_id} derives from "
                            f"unselected bundle {source_bundle_id!r}"
                        ) from exc
                    source_bundle.resource(source_resource_id)

    @staticmethod
    def _require_bundle(
        bundles: dict[str, _RegisteredBundle[BundleT]],
        bundle_id: str,
        kind: str,
    ) -> _RegisteredBundle[BundleT]:
        try:
            return bundles[bundle_id]
        except KeyError as exc:
            raise ValueError(f"profile pins unavailable {kind} bundle {bundle_id!r}") from exc

    @staticmethod
    def _require_matching_build(
        profile: AnalysisProfile,
        bundle: ResourceBundle,
    ) -> None:
        if bundle.genome_build is not None and bundle.genome_build != profile.genome_build:
            raise ValueError(
                f"profile {profile.profile_id!r} is {profile.genome_build.value}, but bundle "
                f"{bundle.bundle_id!r} is {bundle.genome_build.value}; cross-build resolution "
                "is prohibited"
            )

    def _resolve_bundle_files(
        self,
        namespace: str,
        registered: _RegisteredBundle[ResourceBundle],
        *,
        paths: dict[str, str],
        checksums: dict[str, str],
        releases: dict[str, str],
        verify_files: bool,
    ) -> None:
        bundle = registered.bundle
        for resource in bundle.resources:
            key = f"{namespace}.{resource.role}"
            candidate = (registered.directory / resource.path).resolve()
            try:
                candidate.relative_to(registered.directory)
                candidate.relative_to(self.resource_root)
            except ValueError as exc:
                raise ValueError(
                    f"resource {bundle.bundle_id}:{resource.resource_id} escapes its bundle"
                ) from exc

            if not candidate.is_file():
                if not resource.required:
                    continue
                raise FileNotFoundError(
                    f"required resource {bundle.bundle_id}:{resource.resource_id} is missing: "
                    f"{candidate}"
                )
            if resource.sha256 is None:
                raise ValueError(
                    f"materialized resource {bundle.bundle_id}:{resource.resource_id} has no "
                    "pinned SHA256; the bundle is not activated completely"
                )
            observed_size = candidate.stat().st_size
            if verify_files:
                observed_checksum = sha256_file(candidate)
                if observed_checksum != resource.sha256:
                    raise ValueError(
                        f"resource {bundle.bundle_id}:{resource.resource_id} checksum mismatch: "
                        f"expected {resource.sha256}, observed {observed_checksum}"
                    )
            if resource.size_bytes is not None and observed_size != resource.size_bytes:
                raise ValueError(
                    f"resource {bundle.bundle_id}:{resource.resource_id} size mismatch: "
                    f"expected {resource.size_bytes}, observed {observed_size}"
                )
            paths[key] = str(candidate)
            checksums[key] = resource.sha256
            if resource.release is not None:
                releases[key] = resource.release


__all__ = [
    "DEFAULT_RESOURCE_ROOT",
    "RESOURCE_ROOT_ENV",
    "RegistryDiagnostic",
    "ResourceRegistry",
    "resource_root_from_environment",
]

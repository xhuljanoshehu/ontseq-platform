"""CLI surface for explicit, network-bounded reference lifecycle operations."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .io import load_mapping
from .profile_analysis import configuration_root
from .reference_catalog import (
    BundleValidationReport,
    ReferenceBundleInstaller,
    ReferenceCatalog,
)
from .resource_bootstrap import REFERENCE_BUNDLE_ID, GRCh38ResourceBootstrapper
from .resource_registry import ResourceRegistry


def _resource_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--resource-root", type=Path)
    parser.add_argument("--config-root", type=Path)


def add_references_parser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    references = subparsers.add_parser(
        "references", help="Install and validate manifest-locked GRCh38 resources"
    )
    commands = references.add_subparsers(dest="references_command", required=True)

    listing = commands.add_parser("list", help="List installable and release-lock inventories")
    _resource_options(listing)
    listing.add_argument("--json", action="store_true", dest="as_json")

    status = commands.add_parser("status", help="Show fast local bundle and profile status")
    _resource_options(status)
    status.add_argument("--json", action="store_true", dest="as_json")

    validate = commands.add_parser("validate", help="Checksum-validate installed references")
    _resource_options(validate)
    validate.add_argument("bundle_id", nargs="?")
    validate.add_argument("--json", action="store_true", dest="as_json")

    install = commands.add_parser("install", help="Install a pinned GRCh38 reference recipe")
    install.add_argument("bundle_id")
    _resource_options(install)
    install.add_argument("--offline", action="store_true")

    imported = commands.add_parser("import", help="Import a pinned local reference tree")
    imported.add_argument("path", type=Path)
    _resource_options(imported)

    repair = commands.add_parser("repair", help="Repair missing or checksum-failed resources")
    repair.add_argument("bundle_id")
    _resource_options(repair)
    repair.add_argument("--offline", action="store_true")


def _config_root(args: argparse.Namespace) -> Path:
    return configuration_root(getattr(args, "config_root", None))


def _catalog(args: argparse.Namespace) -> ReferenceCatalog:
    return ReferenceCatalog.discover(_config_root(args) / "reference_bundles")


def _packaged_config_root() -> Path:
    """Locate shipped configs without trusting cwd or a caller-supplied ``--config-root``."""

    package = Path(__file__).resolve()
    candidates = [
        package.parents[2] / "configs",
        Path(sys.prefix) / "share" / "ontseq" / "configs",
        *(parent / "share" / "ontseq" / "configs" for parent in package.parents),
    ]
    for candidate in candidates:
        if (candidate / "reference_bundles" / REFERENCE_BUNDLE_ID / "bundle.recipe.yaml").is_file():
            return candidate.resolve()
    rendered = ", ".join(str(candidate) for candidate in candidates)
    raise FileNotFoundError(
        "the shipped ONTSeq reference authority catalog was not found; checked: " + rendered
    )


def _packaged_authority_catalog() -> ReferenceCatalog:
    return ReferenceCatalog.discover(_packaged_config_root() / "reference_bundles")


def _declared_import_bundle_id(source: Path) -> str | None:
    manifest = source / "bundle.yaml"
    if not manifest.is_file():
        manifest = source / "bundle.recipe.yaml"
    if not manifest.is_file():
        return None
    value = load_mapping(manifest).get("bundle_id")
    return value if isinstance(value, str) else None


def _authority_catalog_for_import(source: Path) -> ReferenceCatalog | None:
    """Claim the official ID from shipped provenance; leave custom IDs catalog-independent."""

    if _declared_import_bundle_id(source) != REFERENCE_BUNDLE_ID:
        return None
    catalog = _packaged_authority_catalog()
    catalog.get(REFERENCE_BUNDLE_ID)
    return catalog


def _report_payload(report: BundleValidationReport) -> dict[str, object]:
    return {
        "bundle_id": report.bundle_id,
        "path": str(report.bundle_path),
        "valid": report.valid,
        "manifest_valid": report.manifest_valid,
        "errors": list(report.errors),
        "resources": [
            {
                "resource_id": item.resource_id,
                "role": item.role,
                "state": item.state.value,
                "required": item.required,
                "message": item.message,
            }
            for item in report.resources
        ],
    }


def _print_reports(reports: tuple[BundleValidationReport, ...], *, as_json: bool) -> None:
    payload = [_report_payload(report) for report in reports]
    if as_json:
        print(json.dumps(payload, indent=2))
        return
    if not reports:
        print("No activated reference bundles.")
        return
    for report in reports:
        marker = "OK" if report.valid else "INVALID"
        print(f"{report.bundle_id}\t{marker}\t{report.bundle_path}")
        for error in report.errors:
            print(f"  ERROR\t{error}")
        for item in report.resources:
            if not item.valid:
                print(f"  {item.state.value.upper()}\t{item.resource_id}\t{item.message}")


def _bootstrap_if_official(
    args: argparse.Namespace, bundle_id: str, *, repair_existing: bool = False
) -> None:
    if bundle_id != REFERENCE_BUNDLE_ID:
        return
    bootstrapper = GRCh38ResourceBootstrapper(
        args.resource_root,
        packaged_config_root=_packaged_config_root(),
    )
    result = bootstrapper.repair() if repair_existing else bootstrapper.activate()
    action = "repaired" if repair_existing else "activated"
    print(
        f"{action} profile resources: "
        + ", ".join(result.profile_ids)
        + f"; panel intervals: {result.panel_summary.compilation.roi_interval_count} ROI"
        + f"; replaced paths: {len(result.repaired_paths)}"
    )


def _release_lock_inventories(root: Path) -> list[dict[str, object]]:
    inventories: list[dict[str, object]] = []
    for path in sorted(root.glob("*/release-lock-required.yaml")):
        payload = load_mapping(path)
        inventories.append(
            {
                "bundle_id": payload.get("target_bundle_id", path.parent.name),
                "version": payload.get("target_bundle_version", "unspecified"),
                "genome_build": payload.get("genome_build", "unspecified"),
                "status": payload.get("status", "RELEASE_LOCK_REQUIRED"),
                "installable": False,
                "path": str(path),
            }
        )
    return inventories


def _profile_status(
    registry: ResourceRegistry,
    *,
    verify_checksums: bool,
) -> list[dict[str, object]]:
    """Validate each profile's pinned reference, knowledge and optional panel context."""

    status: list[dict[str, object]] = []
    for profile_id in sorted(registry.profiles):
        try:
            context = registry.resolve_profile(profile_id, verify_files=verify_checksums)
        except (KeyError, OSError, ValueError) as error:
            status.append(
                {
                    "profile_id": profile_id,
                    "valid": False,
                    "error": str(error),
                }
            )
            continue
        status.append(
            {
                "profile_id": profile_id,
                "valid": True,
                "reference_bundle": context.reference_bundle_id,
                "panel_bundle": context.panel_bundle_id,
                "knowledge_bundle": context.knowledge_bundle_id,
            }
        )
    return status


def handle_references_command(args: argparse.Namespace) -> bool:
    """Execute ``ontseq references`` and return whether the command was handled."""

    if args.command != "references":
        return False
    installer = ReferenceBundleInstaller(args.resource_root)
    command = args.references_command
    if command == "list":
        catalog_root = _config_root(args) / "reference_bundles"
        catalog = ReferenceCatalog.discover(catalog_root)
        payload: list[dict[str, object]] = [
            {
                "bundle_id": bundle.bundle_id,
                "version": bundle.version,
                "genome_build": bundle.genome_build.value,
                "status": "installable",
                "installable": True,
            }
            for bundle in catalog.list()
        ]
        payload.extend(_release_lock_inventories(catalog_root))
        if args.as_json:
            print(json.dumps(payload, indent=2))
        elif payload:
            for item in payload:
                print(f"{item['bundle_id']}\t{item['genome_build']}\t{item['status']}")
        else:
            print("No reference recipes or release-lock inventories found.")
        return True
    if command == "status":
        reports = installer.status()
        registry = ResourceRegistry(args.resource_root)
        profile_status = _profile_status(registry, verify_checksums=False)
        ready_profiles = [
            str(item["profile_id"]) for item in profile_status if item["valid"] is True
        ]
        if args.as_json:
            print(
                json.dumps(
                    {
                        "references": [_report_payload(report) for report in reports],
                        "profiles": ready_profiles,
                        "profile_status": profile_status,
                        "diagnostics": [
                            {
                                "path": str(item.path),
                                "code": item.code,
                                "message": item.message,
                            }
                            for item in registry.diagnostics
                        ],
                    },
                    indent=2,
                )
            )
        else:
            _print_reports(reports, as_json=False)
            print("profiles\t" + (", ".join(ready_profiles) or "none"))
            for item in profile_status:
                if item["valid"] is not True:
                    print(f"  INVALID\t{item['profile_id']}\t{item['error']}")
        return True
    if command == "validate":
        reports = (
            (installer.validate(args.bundle_id),)
            if args.bundle_id is not None
            else tuple(installer.validate(report.bundle_id) for report in installer.status())
        )
        _print_reports(reports, as_json=args.as_json)
        validation_profile_status: list[dict[str, object]] = []
        if args.bundle_id is None:
            registry = ResourceRegistry(args.resource_root)
            validation_profile_status = _profile_status(registry, verify_checksums=True)
            if not args.as_json:
                for item in validation_profile_status:
                    marker = "OK" if item["valid"] is True else "INVALID"
                    detail = "" if item["valid"] is True else f"\t{item['error']}"
                    print(f"profile\t{marker}\t{item['profile_id']}{detail}")
        if (
            not reports
            or any(not report.valid for report in reports)
            or any(item["valid"] is not True for item in validation_profile_status)
        ):
            raise SystemExit(2)
        return True
    if command in {"install", "repair"}:
        if args.bundle_id == REFERENCE_BUNDLE_ID:
            recipe = _packaged_authority_catalog().get(REFERENCE_BUNDLE_ID)
        else:
            try:
                recipe = _catalog(args).get(args.bundle_id)
            except KeyError as exc:
                inventory = (
                    _config_root(args)
                    / "reference_bundles"
                    / args.bundle_id
                    / "release-lock-required.yaml"
                )
                if inventory.is_file():
                    raise ValueError(
                        f"{args.bundle_id} is RELEASE_LOCK_REQUIRED and cannot be installed until "
                        "its exact source sizes and SHA256 values have been independently pinned; "
                        f"see {inventory}"
                    ) from exc
                raise
        installed = (
            installer.install(recipe, offline=args.offline)
            if command == "install"
            else installer.repair(recipe, offline=args.offline)
        )
        print(f"{command}ed {installed.bundle.bundle_id}: {installed.path}")
        _bootstrap_if_official(
            args,
            installed.bundle.bundle_id,
            repair_existing=command == "repair",
        )
        return True
    if command == "import":
        imported = installer.import_bundle(
            args.path,
            authority_catalog=_authority_catalog_for_import(args.path),
        )
        print(f"imported {imported.bundle.bundle_id}: {imported.path}")
        _bootstrap_if_official(args, imported.bundle.bundle_id)
        return True
    raise AssertionError(f"unhandled references command: {command}")


__all__ = ["add_references_parser", "handle_references_command"]

from __future__ import annotations

import hashlib
import json
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import yaml

from ontseq_platform.models import GenomeBuild, ReferenceContig, ReferenceLock
from ontseq_platform.reference import grch38_canonical_25_contigs
from ontseq_platform.resource_registry import (
    DEFAULT_RESOURCE_ROOT,
    RESOURCE_ROOT_ENV,
    ResourceRegistry,
    resource_root_from_environment,
)


def _write_file(bundle_dir: Path, relative_path: str, content: str) -> dict[str, object]:
    path = bundle_dir / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content.encode())
    return {
        "path": relative_path,
        "sha256": hashlib.sha256(content.encode()).hexdigest(),
        "size_bytes": len(content.encode()),
    }


def _write_manifest(directory: Path, payload: dict[str, object]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "bundle.yaml").write_text(
        yaml.safe_dump(payload, sort_keys=False),
        encoding="utf-8",
    )


def _resource(
    directory: Path,
    resource_id: str,
    role: str,
    filename: str,
) -> dict[str, object]:
    return {
        "resource_id": resource_id,
        "role": role,
        **_write_file(directory, filename, f"{resource_id}\n"),
    }


def _reference_bundle(
    root: Path,
    bundle_id: str,
    build: str = "GRCh38",
    *,
    contigs: tuple[tuple[str, int], ...] = (("chr1", 1),),
) -> None:
    directory = root / "references" / bundle_id
    fai_content = "".join(f"{name}\t{length}\t0\t1\t2\n" for name, length in contigs)
    fai = {
        "resource_id": "fai",
        "role": "fasta_index",
        **_write_file(directory, "genome.fa.fai", fai_content),
    }
    reference_lock = ReferenceLock(
        reference_id=bundle_id,
        genome_build=GenomeBuild(build),
        contigs=[ReferenceContig(name=name, length=length) for name, length in contigs],
        source_fai_sha256=str(fai["sha256"]),
    )
    resources = [
        _resource(directory, "fasta", "genome_fasta", "genome.fa"),
        fai,
        {
            "resource_id": "lock",
            "role": "reference_lock",
            **_write_file(
                directory,
                "reference.lock.json",
                reference_lock.model_dump_json() + "\n",
            ),
        },
    ]
    _write_manifest(
        directory,
        {
            "schema_version": "1.0.0",
            "bundle_type": "reference",
            "bundle_id": bundle_id,
            "version": "1",
            "genome_build": build,
            "resources": resources,
            "reference_lock_resource_id": "lock",
            "fasta_resource_id": "fasta",
            "fai_resource_id": "fai",
        },
    )


def _knowledge_bundle(root: Path, bundle_id: str = "HEMATOLOGY_v1") -> None:
    directory = root / "knowledge" / bundle_id
    _write_manifest(
        directory,
        {
            "schema_version": "1.0.0",
            "bundle_type": "knowledge",
            "bundle_id": bundle_id,
            "version": "1",
            "genome_build": "GRCh38",
            "coordinate_bearing": False,
            "resources": [
                _resource(
                    directory,
                    "rearrangements",
                    "recurrent_rearrangements",
                    "rearrangements.tsv",
                )
            ],
        },
    )


def _panel_bundle(root: Path, bundle_id: str = "AML_AS_111_GRCh38_v1") -> None:
    directory = root / "panels" / bundle_id
    selection = _resource(directory, "selection", "selection_panel_buffered", "selection.bed")
    selection["coordinate_system"] = "zero_based_half_open"
    roi = _resource(directory, "roi", "analysis_roi_unbuffered", "roi.bed")
    roi["coordinate_system"] = "zero_based_half_open"
    _write_manifest(
        directory,
        {
            "schema_version": "1.0.0",
            "bundle_type": "panel",
            "bundle_id": bundle_id,
            "version": "1",
            "genome_build": "GRCh38",
            "assay_mode": "adaptive_sampling",
            "resources": [
                selection,
                roi,
                _resource(
                    directory,
                    "transcripts",
                    "transcript_cache",
                    "transcripts.tsv",
                ),
            ],
            "selection_panel_resource_id": "selection",
            "analysis_roi_resource_id": "roi",
            "transcript_cache_resource_id": "transcripts",
        },
    )


def _profile(
    root: Path,
    profile_id: str,
    *,
    adaptive: bool,
    dictionary_contract: str = "exact_full",
) -> None:
    directory = root / "profiles"
    directory.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": "1.0.0",
        "profile_id": profile_id,
        "version": "1",
        "genome_build": "GRCh38",
        "assay_mode": "adaptive_sampling" if adaptive else "lcwgs",
        "reference_bundle": "GRCh38_TEST_v1",
        "reference_dictionary_contract": dictionary_contract,
        "knowledge_bundle": "HEMATOLOGY_v1",
        "panel_bundle": "AML_AS_111_GRCh38_v1" if adaptive else None,
        "adaptive_sampling": "enabled" if adaptive else "disabled",
    }
    (directory / f"{profile_id}.yaml").write_text(
        yaml.safe_dump(payload, sort_keys=False),
        encoding="utf-8",
    )


class ResourceRegistryTests(unittest.TestCase):
    def _complete_root(self, root: Path) -> None:
        _reference_bundle(root, "GRCh38_TEST_v1")
        _knowledge_bundle(root)
        _panel_bundle(root)
        _profile(root, "AML_LCWGS_GRCh38", adaptive=False)
        _profile(root, "AML_AS_111_GRCh38", adaptive=True)

    def test_manifested_grch38_bundles_resolve_to_absolute_pinned_paths(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            registry = ResourceRegistry(root, active_build=GenomeBuild.GRCH38)

            context = registry.resolve_profile("AML_AS_111_GRCh38")

            self.assertEqual(context.reference_bundle_id, "GRCh38_TEST_v1")
            self.assertEqual(context.panel_bundle_id, "AML_AS_111_GRCh38_v1")
            self.assertEqual(context.knowledge_bundle_id, "HEMATOLOGY_v1")
            self.assertIn("reference.genome_fasta", context.resource_paths)
            self.assertIn("panel.selection_panel_buffered", context.resource_paths)
            self.assertTrue(
                all(Path(path).is_absolute() for path in context.resource_paths.values())
            )
            self.assertEqual(
                set(context.resource_paths),
                set(context.resource_checksums),
            )

    def test_loose_files_and_unmanifested_directories_are_not_resources(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            loose = root / "references" / "LOOSE"
            loose.mkdir(parents=True)
            (loose / "genome.fa").write_text("not active", encoding="utf-8")
            (root / "references" / "also-loose.bed").write_text("chr1\t0\t1\n")

            registry = ResourceRegistry(root)

            self.assertEqual(set(registry.references), {"GRCh38_TEST_v1"})

    def test_other_build_bundle_is_not_active_in_grch38_registry(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            _reference_bundle(root, "GRCh37_TEST_v1", build="GRCh37")

            registry = ResourceRegistry(root, active_build=GenomeBuild.GRCH38)

            self.assertNotIn("GRCh37_TEST_v1", registry.references)
            self.assertTrue(
                any(
                    item.code == "inactive_build" and "GRCh37" in item.message
                    for item in registry.diagnostics
                )
            )

    def test_checksum_mismatch_fails_before_context_is_returned(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            fasta = root / "references" / "GRCh38_TEST_v1" / "genome.fa"
            fasta.write_text("modified after activation\n", encoding="utf-8")
            registry = ResourceRegistry(root)

            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                registry.resolve_profile("AML_LCWGS_GRCh38")

    def test_fast_resolution_checks_pins_and_sizes_without_hashing_file_contents(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            registry = ResourceRegistry(root)

            with patch(
                "ontseq_platform.resource_registry.sha256_file",
                side_effect=AssertionError("fast resolution must not hash resource bytes"),
            ):
                context = registry.resolve_profile(
                    "AML_AS_111_GRCh38",
                    verify_files=False,
                )

            self.assertEqual(context.panel_bundle_id, "AML_AS_111_GRCh38_v1")
            self.assertEqual(context.knowledge_bundle_id, "HEMATOLOGY_v1")

    def test_fast_resolution_still_rejects_a_wrong_declared_size(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            panel = root / "panels" / "AML_AS_111_GRCh38_v1" / "selection.bed"
            panel.write_text("changed-size\n", encoding="utf-8")
            registry = ResourceRegistry(root)

            with self.assertRaisesRegex(ValueError, "size mismatch"):
                registry.resolve_profile("AML_AS_111_GRCh38", verify_files=False)

    def test_profile_rejects_reference_lock_for_another_build_in_fast_mode(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            lock_path = root / "references" / "GRCh38_TEST_v1" / "reference.lock.json"
            payload = json.loads(lock_path.read_text(encoding="utf-8"))
            payload["genome_build"] = "GRCh37"
            lock_path.write_text(
                ReferenceLock.model_validate(payload).model_dump_json() + "\n",
                encoding="utf-8",
                newline="\n",
            )

            with self.assertRaisesRegex(ValueError, "reference lock is GRCh37"):
                ResourceRegistry(root).resolve_profile(
                    "AML_LCWGS_GRCh38",
                    verify_files=False,
                )

    def test_profile_rejects_reference_lock_for_a_different_fai_pin(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            lock_path = root / "references" / "GRCh38_TEST_v1" / "reference.lock.json"
            payload = json.loads(lock_path.read_text(encoding="utf-8"))
            payload["source_fai_sha256"] = "0" * 64
            lock_path.write_text(
                ReferenceLock.model_validate(payload).model_dump_json() + "\n",
                encoding="utf-8",
                newline="\n",
            )

            with self.assertRaisesRegex(ValueError, "pinned FASTA index"):
                ResourceRegistry(root).resolve_profile(
                    "AML_LCWGS_GRCh38",
                    verify_files=False,
                )

    def test_canonical25_profile_requires_a_compatible_full_reference_lock(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            _profile(
                root,
                "AML_LCWGS_GRCh38_CANONICAL25",
                adaptive=False,
                dictionary_contract="grch38_canonical_25",
            )

            with self.assertRaisesRegex(ValueError, "complete Canonical-25 dictionary"):
                ResourceRegistry(root).resolve_profile(
                    "AML_LCWGS_GRCh38_CANONICAL25",
                    verify_files=False,
                )

    def test_canonical25_profile_resolves_against_a_superset_reference_lock(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            contigs = (*grch38_canonical_25_contigs(), ("GL000008.2", 209709))
            _reference_bundle(
                root,
                "GRCh38_TEST_v1",
                contigs=contigs,
            )
            _knowledge_bundle(root)
            _profile(
                root,
                "AML_LCWGS_GRCh38_CANONICAL25",
                adaptive=False,
                dictionary_contract="grch38_canonical_25",
            )

            context = ResourceRegistry(root).resolve_profile(
                "AML_LCWGS_GRCh38_CANONICAL25",
                verify_files=False,
            )

            self.assertEqual(
                context.reference_dictionary_contract.value,
                "grch38_canonical_25",
            )

    def test_profile_cannot_fall_back_to_an_unavailable_bundle(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            self._complete_root(root)
            profile_path = root / "profiles" / "AML_LCWGS_GRCh38.yaml"
            payload = yaml.safe_load(profile_path.read_text(encoding="utf-8"))
            payload["reference_bundle"] = "GRCh37_TEST_v1"
            profile_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
            _reference_bundle(root, "GRCh37_TEST_v1", build="GRCh37")
            registry = ResourceRegistry(root, active_build=GenomeBuild.GRCH38)

            with self.assertRaisesRegex(ValueError, "unavailable reference bundle"):
                registry.resolve_profile("AML_LCWGS_GRCh38")

    def test_invalid_manifest_is_diagnosed_not_activated(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            bad = root / "references" / "BROKEN"
            _write_manifest(bad, {"bundle_type": "reference", "bundle_id": "BROKEN"})

            registry = ResourceRegistry(root)

            self.assertFalse(registry.references)
            self.assertEqual(registry.diagnostics[0].code, "invalid_manifest")

    def test_explicit_root_precedes_environment_and_environment_precedes_default(self) -> None:
        with (
            TemporaryDirectory() as raw,
            TemporaryDirectory() as explicit_raw,
            patch.dict(os.environ, {RESOURCE_ROOT_ENV: raw}),
        ):
            self.assertEqual(resource_root_from_environment(), Path(raw))
            self.assertEqual(
                resource_root_from_environment(Path(explicit_raw)),
                Path(explicit_raw),
            )

    def test_core_default_resource_root_remains_opt_ontseq(self) -> None:
        with patch.dict(os.environ, {RESOURCE_ROOT_ENV: ""}):
            self.assertEqual(resource_root_from_environment(), DEFAULT_RESOURCE_ROOT)
            self.assertEqual(DEFAULT_RESOURCE_ROOT, Path("/opt/ontseq"))


if __name__ == "__main__":
    unittest.main()

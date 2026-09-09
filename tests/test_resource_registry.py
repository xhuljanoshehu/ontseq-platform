from __future__ import annotations

import hashlib
import json
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import yaml

from ontseq_platform.models import (
    GenomeBuild,
    ReferenceContig,
    ReferenceDictionaryContract,
    ReferenceLock,
)
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
    selection.update(**_write_file(directory, "selection.bed", "chr1\t0\t10\tA\n"))
    selection["coordinate_system"] = "zero_based_half_open"
    roi = _resource(directory, "roi", "analysis_roi_unbuffered", "roi.bed")
    roi.update(**_write_file(directory, "roi.bed", "chr1\t1\t9\tA\n"))
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


def _mapped_panel_bundle(root: Path, bundle_id: str = "AML_AS_MAPPED_GRCh38_v1") -> None:
    directory = root / "panels" / bundle_id
    source = _resource(directory, "source", "coordinate_mapping_source", "source.bed")
    source.update(
        **_write_file(
            directory,
            "source.bed",
            "chr1\t0\t10\tACACA\nchrX\t0\t10\tCT45A2\n",
        )
    )
    source["coordinate_system"] = "zero_based_half_open"
    output = _resource(directory, "mapped", "coordinate_mapping_output", "mapped.bed")
    output.update(**_write_file(directory, "mapped.bed", "chr1\t0\t10\tACACA\n"))
    output.update(
        coordinate_system="zero_based_half_open",
        generated=True,
        derived_from=["source", "mapping_lock"],
    )
    unmapped = _resource(
        directory,
        "unmapped",
        "coordinate_mapping_unmapped",
        "unmapped.txt",
    )
    unmapped.update(
        **_write_file(directory, "unmapped.txt", "#Split in new\nchrX\t0\t10\tCT45A2\n")
    )
    unmapped.update(
        coordinate_system="zero_based_half_open",
        generated=True,
        derived_from=["source", "mapping_lock"],
    )
    gene_map = _resource(directory, "gene_map", "target_gene_map", "gene-map.tsv")
    selection = _resource(
        directory,
        "selection",
        "selection_panel_buffered",
        "selection.bed",
    )
    selection.update(**_write_file(directory, "selection.bed", "chr1\t0\t10\tACACA\n"))
    selection.update(
        coordinate_system="zero_based_half_open",
        generated=True,
        derived_from=["mapped", "gene_map"],
    )
    roi = _resource(directory, "roi", "analysis_roi_unbuffered", "roi.bed")
    roi.update(**_write_file(directory, "roi.bed", "chr1\t1\t9\tACACA\n"))
    roi.update(
        coordinate_system="zero_based_half_open",
        generated=True,
        derived_from=["selection", "gene_map"],
    )
    transcripts = _resource(directory, "transcripts", "transcript_cache", "transcripts.tsv")
    transcripts.update(generated=True, derived_from=["selection", "gene_map"])
    roundtrip = _resource(
        directory,
        "roundtrip",
        "coordinate_mapping_roundtrip_output",
        "roundtrip.bed",
    )
    roundtrip.update(
        **_write_file(directory, "roundtrip.bed", ""),
        coordinate_system="zero_based_half_open",
        generated=True,
        derived_from=["mapped", "mapping_lock"],
    )
    roundtrip_unmapped = _resource(
        directory,
        "roundtrip_unmapped",
        "coordinate_mapping_roundtrip_unmapped",
        "roundtrip.unmapped.txt",
    )
    roundtrip_unmapped.update(
        **_write_file(
            directory,
            "roundtrip.unmapped.txt",
            "#Split in new\nchr1\t0\t10\tACACA\n",
        ),
        coordinate_system="zero_based_half_open",
        generated=True,
        derived_from=["mapped", "mapping_lock"],
    )
    lock_payload = {
        "schema_version": "1.0.0",
        "mapping_id": "TEST_GRCh37_to_GRCh38_v1",
        "status": "controlled_build_time_unvalidated",
        "source_genome_build": "GRCh37",
        "target_genome_build": "GRCh38",
        "source_panel_bundle_id": "SOURCE_PANEL_GRCh37_v1",
        "source_panel_resource_id": "selection",
        "source_selection_sha256": source["sha256"],
        "source_interval_count": 2,
        "mapping_method": "ucsc_liftover_chain",
        "mapping_tool": "UCSC liftOver",
        "mapping_tool_url": "https://example.invalid/liftOver",
        "mapping_tool_sha256": "1" * 64,
        "mapping_tool_size_bytes": 1,
        "mapping_tool_build_id": "test-build",
        "mapping_tool_last_modified": "2026-01-01T00:00:00Z",
        "chain_name": "test.chain.gz",
        "chain_url": "https://example.invalid/test.chain.gz",
        "chain_md5": "2" * 32,
        "chain_sha256": "3" * 64,
        "min_match": 0.99,
        "min_blocks": 1.0,
        "multiple": False,
        "mapped_output_sha256": output["sha256"],
        "unmapped_output_sha256": unmapped["sha256"],
        "mapped_interval_count": 1,
        "same_chromosome_mapped_count": 1,
        "span_identical_count": 1,
        "maximum_span_delta_bases": 0,
        "maximum_span_delta_fraction": 0,
        "unmapped_target_labels": ["CT45A2"],
        "split_target_labels": ["CT45A2"],
        "roundtrip_chain_name": "reverse.chain.gz",
        "roundtrip_chain_url": "https://example.invalid/reverse.chain.gz",
        "roundtrip_chain_md5": "4" * 32,
        "roundtrip_chain_sha256": "5" * 64,
        "roundtrip_min_match": 0.99,
        "roundtrip_min_blocks": 1.0,
        "roundtrip_multiple": False,
        "roundtrip_mapped_output_sha256": roundtrip["sha256"],
        "roundtrip_mapped_output_size_bytes": roundtrip["size_bytes"],
        "roundtrip_unmapped_output_sha256": roundtrip_unmapped["sha256"],
        "roundtrip_unmapped_output_size_bytes": roundtrip_unmapped["size_bytes"],
        "roundtrip_mapped_interval_count": 0,
        "reciprocal_exact_interval_count": 0,
        "roundtrip_unmapped_target_labels": ["ACACA"],
        "roundtrip_split_target_labels": ["ACACA"],
        "roundtrip_review_required_target_labels": ["ACACA"],
        "final_selection_sha256": selection["sha256"],
        "runtime_mapping_allowed": False,
        "note": "Synthetic RUO mapping fixture.",
    }
    lock = {
        "resource_id": "mapping_lock",
        "role": "coordinate_mapping_lock",
        **_write_file(
            directory,
            "mapping.lock.yaml",
            yaml.safe_dump(lock_payload, sort_keys=False),
        ),
    }
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
                source,
                lock,
                output,
                unmapped,
                roundtrip,
                roundtrip_unmapped,
                gene_map,
                selection,
                roi,
                transcripts,
            ],
            "selection_panel_resource_id": "selection",
            "analysis_roi_resource_id": "roi",
            "transcript_cache_resource_id": "transcripts",
            "target_gene_map_resource_id": "gene_map",
            "coordinate_mapping_lock_resource_id": "mapping_lock",
            "coordinate_mapping_source_resource_id": "source",
            "coordinate_mapping_output_resource_id": "mapped",
            "coordinate_mapping_unmapped_resource_id": "unmapped",
            "coordinate_mapping_roundtrip_output_resource_id": "roundtrip",
            "coordinate_mapping_roundtrip_unmapped_resource_id": "roundtrip_unmapped",
            "coordinate_mapping_origin_bundle_id": "SOURCE_PANEL_GRCh37_v1",
            "coordinate_mapping_origin_resource_id": "selection",
            "unresolved_targets": ["CT45A2_REVIEW_REQUIRED", "IGH_REVIEW_REQUIRED"],
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
            self.assertIsNotNone(context.panel_resolution)
            assert context.panel_resolution is not None
            self.assertEqual(
                context.panel_resolution.mapping_status,
                "native_build_not_required",
            )
            self.assertIn("reference.genome_fasta", context.resource_paths)
            self.assertIn("panel.selection_panel_buffered", context.resource_paths)
            self.assertTrue(
                all(Path(path).is_absolute() for path in context.resource_paths.values())
            )
            self.assertEqual(
                set(context.resource_paths),
                set(context.resource_checksums),
            )

    def test_mapped_panel_resolution_is_bound_to_validated_lock_and_resources(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            _reference_bundle(root, "GRCh38_TEST_v1")
            _knowledge_bundle(root)
            _mapped_panel_bundle(root)
            _profile(root, "AML_AS_MAPPED_GRCh38", adaptive=True)
            profile_path = root / "profiles" / "AML_AS_MAPPED_GRCh38.yaml"
            profile_payload = yaml.safe_load(profile_path.read_text(encoding="utf-8"))
            profile_payload["panel_bundle"] = "AML_AS_MAPPED_GRCh38_v1"
            profile_path.write_text(
                yaml.safe_dump(profile_payload, sort_keys=False),
                encoding="utf-8",
            )

            context = ResourceRegistry(root).resolve_profile("AML_AS_MAPPED_GRCh38")

            summary = context.panel_resolution
            self.assertIsNotNone(summary)
            assert summary is not None
            self.assertEqual(summary.mapping_status, "controlled_build_time_unvalidated")
            self.assertEqual(summary.source_panel_bundle_id, "SOURCE_PANEL_GRCh37_v1")
            self.assertEqual(summary.source_panel_resource_id, "selection")
            self.assertEqual(summary.selection_interval_count, 1)
            self.assertEqual(summary.analysis_roi_interval_count, 1)
            self.assertEqual(summary.mapped_interval_count, 1)
            self.assertEqual(summary.unmapped_target_labels, ["CT45A2"])
            self.assertEqual(summary.roundtrip_review_required_target_labels, ["ACACA"])
            self.assertEqual(
                summary.unresolved_target_labels,
                ["CT45A2_REVIEW_REQUIRED", "IGH_REVIEW_REQUIRED"],
            )
            self.assertEqual(
                summary.coordinate_mapping_lock_sha256,
                context.resource_checksums["panel.coordinate_mapping_lock"],
            )
            self.assertEqual(
                summary.roundtrip_mapping_output_sha256,
                context.resource_checksums["panel.coordinate_mapping_roundtrip_output"],
            )

    def test_mapped_panel_resolution_rejects_lock_that_disagrees_with_selection(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            _reference_bundle(root, "GRCh38_TEST_v1")
            _knowledge_bundle(root)
            _mapped_panel_bundle(root)
            _profile(root, "AML_AS_MAPPED_GRCh38", adaptive=True)
            profile_path = root / "profiles" / "AML_AS_MAPPED_GRCh38.yaml"
            profile_payload = yaml.safe_load(profile_path.read_text(encoding="utf-8"))
            profile_payload["panel_bundle"] = "AML_AS_MAPPED_GRCh38_v1"
            profile_path.write_text(
                yaml.safe_dump(profile_payload, sort_keys=False),
                encoding="utf-8",
            )
            panel_directory = root / "panels" / "AML_AS_MAPPED_GRCh38_v1"
            lock_path = panel_directory / "mapping.lock.yaml"
            lock_payload = yaml.safe_load(lock_path.read_text(encoding="utf-8"))
            lock_payload["final_selection_sha256"] = "0" * 64
            lock_content = yaml.safe_dump(lock_payload, sort_keys=False)
            lock_path.write_bytes(lock_content.encode())
            manifest_path = panel_directory / "bundle.yaml"
            manifest_payload = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
            lock_resource = next(
                item
                for item in manifest_payload["resources"]
                if item["resource_id"] == "mapping_lock"
            )
            lock_resource["sha256"] = hashlib.sha256(lock_content.encode()).hexdigest()
            lock_resource["size_bytes"] = len(lock_content.encode())
            manifest_path.write_text(
                yaml.safe_dump(manifest_payload, sort_keys=False),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "final selection SHA256"):
                ResourceRegistry(root).resolve_profile(
                    "AML_AS_MAPPED_GRCh38",
                    verify_files=False,
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

    def test_adaptive_profile_requires_panel_dictionary_compatibility(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            contigs = (*grch38_canonical_25_contigs(), ("GL000008.2", 209709))
            _reference_bundle(root, "GRCh38_TEST_v1", contigs=contigs)
            _knowledge_bundle(root)
            _panel_bundle(root)
            _profile(
                root,
                "AML_AS_111_GRCh38_CANONICAL25",
                adaptive=True,
                dictionary_contract="grch38_canonical_25",
            )

            registry = ResourceRegistry(root)
            with self.assertRaisesRegex(ValueError, "not declared compatible"):
                registry.resolve_profile(
                    "AML_AS_111_GRCh38_CANONICAL25",
                    verify_files=False,
                )

            panel_manifest = root / "panels" / "AML_AS_111_GRCh38_v1" / "bundle.yaml"
            payload = yaml.safe_load(panel_manifest.read_text(encoding="utf-8"))
            payload["reference_dictionary_contracts"] = [
                "exact_full",
                "grch38_canonical_25",
            ]
            panel_manifest.write_text(
                yaml.safe_dump(payload, sort_keys=False),
                encoding="utf-8",
            )
            registry.refresh()

            context = registry.resolve_profile(
                "AML_AS_111_GRCh38_CANONICAL25",
                verify_files=False,
            )
            self.assertEqual(
                context.reference_dictionary_contract,
                ReferenceDictionaryContract.GRCH38_CANONICAL_25,
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

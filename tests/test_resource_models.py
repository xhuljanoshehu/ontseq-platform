from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from pydantic import ValidationError

from ontseq_platform.models import (
    AnalysisProfile,
    AssayMode,
    CoordinateSystem,
    GenomeBuild,
    PanelBundle,
    ReferenceBundle,
    ReferenceDictionaryContract,
    ResolvedResourceContext,
    ResourceFile,
    SidecarArtifact,
)

SHA = "a" * 64


def _resource(resource_id: str, role: str, path: str) -> ResourceFile:
    return ResourceFile(
        resource_id=resource_id,
        role=role,
        path=path,
        sha256=SHA,
        coordinate_system=CoordinateSystem.ZERO_BASED_HALF_OPEN,
    )


class ResourceContractTests(unittest.TestCase):
    def test_reference_bundle_required_resource_ids_are_distinct(self) -> None:
        resources = [
            _resource("fasta", "genome_fasta", "genome.fa"),
            _resource("fai", "fasta_index", "genome.fa.fai"),
            _resource("lock", "reference_lock", "reference.lock.json"),
            _resource("cache", "annotation_cache", "annotation.sqlite"),
        ]
        with self.assertRaisesRegex(ValidationError, "must be distinct"):
            ReferenceBundle(
                bundle_id="GRCh38_TEST_v1",
                version="1",
                genome_build=GenomeBuild.GRCH38,
                resources=resources,
                reference_lock_resource_id="fasta",
                fasta_resource_id="fasta",
                fai_resource_id="fai",
                annotation_cache_resource_id="cache",
            )

    def test_source_resource_requires_checksum_and_relative_path(self) -> None:
        with self.assertRaises(ValidationError):
            ResourceFile(resource_id="genes", role="genes", path="genes.bed")
        with self.assertRaises(ValidationError):
            ResourceFile(
                resource_id="genes",
                role="genes",
                path="../genes.bed",
                sha256=SHA,
            )
        with self.assertRaises(ValidationError):
            ResourceFile(
                resource_id="genes",
                role="genes",
                path="C:/references/genes.bed",
                sha256=SHA,
            )

    def test_pending_generated_resource_retains_derivation(self) -> None:
        resource = ResourceFile(
            resource_id="annotation_cache",
            role="annotation_cache",
            path="derived/annotation.sqlite",
            generated=True,
            derived_from=["gencode", "mane"],
        )
        self.assertIsNone(resource.sha256)

    def test_generated_resource_can_name_cross_bundle_provenance(self) -> None:
        resource = ResourceFile(
            resource_id="analysis_roi",
            role="analysis_roi_unbuffered",
            path="derived/analysis_roi.bed",
            generated=True,
            derived_from=["GRCh38_GENCODE50_MANE1.5_v1:annotation_cache"],
        )
        self.assertEqual(
            resource.derived_from,
            ["GRCh38_GENCODE50_MANE1.5_v1:annotation_cache"],
        )
        with self.assertRaises(ValidationError):
            ResourceFile(
                resource_id="bad",
                role="analysis_roi_unbuffered",
                path="derived/bad.bed",
                generated=True,
                derived_from=["bundle:resource:extra"],
            )

    def test_panel_bundle_references_typed_roles(self) -> None:
        bundle = PanelBundle(
            bundle_id="AML_AS_111_GRCh38_v1",
            version="1",
            genome_build=GenomeBuild.GRCH38,
            assay_mode=AssayMode.ADAPTIVE_SAMPLING,
            resources=[
                _resource("selection", "selection_panel_buffered", "selection.bed"),
                _resource("roi", "analysis_roi_unbuffered", "roi.bed"),
                _resource("transcripts", "transcript_cache", "transcripts.tsv"),
            ],
            selection_panel_resource_id="selection",
            analysis_roi_resource_id="roi",
            transcript_cache_resource_id="transcripts",
            unresolved_targets=["IGH_REVIEW_REQUIRED"],
        )
        self.assertEqual(bundle.genome_build, GenomeBuild.GRCH38)

    def test_panel_bundle_does_not_accept_swapped_selection_and_roi(self) -> None:
        with self.assertRaisesRegex(ValidationError, "must have role"):
            PanelBundle(
                bundle_id="AML_AS_111_GRCh38_v1",
                version="1",
                genome_build=GenomeBuild.GRCH38,
                assay_mode=AssayMode.ADAPTIVE_SAMPLING,
                resources=[
                    _resource("selection", "analysis_roi_unbuffered", "selection.bed"),
                    _resource("roi", "selection_panel_buffered", "roi.bed"),
                    _resource("transcripts", "transcript_cache", "transcripts.tsv"),
                ],
                selection_panel_resource_id="selection",
                analysis_roi_resource_id="roi",
                transcript_cache_resource_id="transcripts",
            )

    def test_panel_beds_require_explicit_zero_based_half_open_coordinates(self) -> None:
        for resource_id, role in (
            ("selection", "selection_panel_buffered"),
            ("roi", "analysis_roi_unbuffered"),
        ):
            for invalid_coordinate in (None, CoordinateSystem.ONE_BASED_INCLUSIVE):
                with (
                    self.subTest(role=role, coordinate_system=invalid_coordinate),
                    self.assertRaisesRegex(ValidationError, "zero_based_half_open"),
                ):
                    resources = [
                        _resource("selection", "selection_panel_buffered", "selection.bed"),
                        _resource("roi", "analysis_roi_unbuffered", "roi.bed"),
                        _resource("transcripts", "transcript_cache", "transcripts.tsv"),
                    ]
                    resource_index = {"selection": 0, "roi": 1}[resource_id]
                    resources[resource_index] = resources[resource_index].model_copy(
                        update={"coordinate_system": invalid_coordinate}
                    )
                    PanelBundle(
                        bundle_id="AML_AS_111_GRCh38_v1",
                        version="1",
                        genome_build=GenomeBuild.GRCH38,
                        assay_mode=AssayMode.ADAPTIVE_SAMPLING,
                        resources=resources,
                        selection_panel_resource_id="selection",
                        analysis_roi_resource_id="roi",
                        transcript_cache_resource_id="transcripts",
                    )

    def test_profiles_keep_assay_modes_separate(self) -> None:
        with self.assertRaises(ValidationError):
            AnalysisProfile(
                profile_id="AML_LCWGS_GRCh38",
                version="1",
                genome_build=GenomeBuild.GRCH38,
                assay_mode=AssayMode.LOW_COVERAGE_WGS,
                reference_bundle="REF",
                knowledge_bundle="HEMATOLOGY",
                panel_bundle="PANEL",
                adaptive_sampling="disabled",
            )

    def test_profiles_default_to_full_dictionary_and_canonical_25_is_grch38_only(self) -> None:
        profile = AnalysisProfile(
            profile_id="AML_LCWGS_GRCh38",
            version="1",
            genome_build=GenomeBuild.GRCH38,
            assay_mode=AssayMode.LOW_COVERAGE_WGS,
            reference_bundle="REF",
            knowledge_bundle="HEMATOLOGY",
            adaptive_sampling="disabled",
        )
        self.assertEqual(
            profile.reference_dictionary_contract,
            ReferenceDictionaryContract.EXACT_FULL,
        )
        with self.assertRaisesRegex(ValidationError, "valid only for GRCh38"):
            AnalysisProfile(
                profile_id="AML_LCWGS_GRCh37_CANONICAL25",
                version="1",
                genome_build=GenomeBuild.GRCH37,
                assay_mode=AssayMode.LOW_COVERAGE_WGS,
                reference_bundle="REF",
                reference_dictionary_contract=(ReferenceDictionaryContract.GRCH38_CANONICAL_25),
                knowledge_bundle="HEMATOLOGY",
                adaptive_sampling="disabled",
            )

    def test_pre_052_resolved_context_defaults_to_exact_full(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            context = ResolvedResourceContext.model_validate(
                {
                    "profile_id": "AML_LCWGS_GRCh38",
                    "profile_version": "v1",
                    "genome_build": "GRCh38",
                    "reference_bundle_id": "GRCh38_TEST_v1",
                    "reference_bundle_version": "v1",
                    "knowledge_bundle_id": "HEMATOLOGY_v1",
                    "knowledge_bundle_version": "v1",
                    "resource_root": str(root),
                    "resource_paths": {"reference.genome_fasta": str(root / "genome.fa")},
                    "resource_checksums": {"reference.genome_fasta": SHA},
                    "resource_releases": {},
                }
            )

            self.assertEqual(
                context.reference_dictionary_contract,
                ReferenceDictionaryContract.EXACT_FULL,
            )

    def test_resolved_context_rejects_canonical_25_outside_grch38(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            with self.assertRaisesRegex(ValidationError, "valid only for GRCh38"):
                ResolvedResourceContext.model_validate(
                    {
                        "profile_id": "AML_LCWGS_GRCh37_CANONICAL25",
                        "profile_version": "v1",
                        "genome_build": "GRCh37",
                        "reference_dictionary_contract": "grch38_canonical_25",
                        "reference_bundle_id": "GRCh37_TEST_v1",
                        "reference_bundle_version": "v1",
                        "knowledge_bundle_id": "HEMATOLOGY_v1",
                        "knowledge_bundle_version": "v1",
                        "resource_root": str(root),
                        "resource_paths": {"reference.genome_fasta": str(root / "genome.fa")},
                        "resource_checksums": {"reference.genome_fasta": SHA},
                        "resource_releases": {},
                    }
                )

    def test_sidecar_path_is_relative_and_checksum_pinned(self) -> None:
        artifact = SidecarArtifact(
            artifact_id="cnv_bins",
            relative_path="sidecars/cnv_bins.tsv.gz",
            schema_version="1.0.0",
            sha256=SHA,
            row_count=10,
        )
        self.assertEqual(artifact.row_count, 10)
        with self.assertRaises(ValidationError):
            SidecarArtifact(
                artifact_id="cnv_bins",
                relative_path="../outside.tsv",
                schema_version="1.0.0",
                sha256=SHA,
                row_count=10,
            )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import yaml

from ontseq_platform.models import (
    GenomeBuild,
    PanelBundle,
    ReferenceDictionaryContract,
)
from ontseq_platform.panel_bundle import sha256_file
from ontseq_platform.panel_compiler import (
    PanelCompilerError,
    read_target_gene_map,
    validate_coordinate_mapping_provenance,
)

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "configs" / "panels" / "AML_AS_111_GRCh37_v1"


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def _synthetic_mapped_bundle(root: Path) -> PanelBundle:
    source = root / "source" / "selection.grch38.bed"
    mapped = root / "derived" / "mapped.grch37.bed"
    unmapped = root / "derived" / "unmapped.txt"
    roundtrip = root / "derived" / "roundtrip.grch38.bed"
    roundtrip_unmapped = root / "derived" / "roundtrip.unmapped.txt"
    gene_map = root / "source" / "target_gene_map.tsv"
    selection = root / "derived" / "selection.bed"
    _write(source, "chr1\t10\t20\tA\nchr2\t30\t50\tB\nchrX\t70\t90\tSPLIT\n")
    _write(mapped, "chr1\t11\t21\tA\nchr2\t31\t53\tB\n")
    _write(unmapped, "#Split in new\nchrX\t70\t90\tSPLIT\n")
    _write(roundtrip, "chr1\t10\t20\tA\n")
    _write(roundtrip_unmapped, "#Split in new\nchr2\t31\t53\tB\n")
    _write(
        gene_map,
        "source_target_label\ttarget_label\tdeclared_chromosome\tgene_id\t"
        "native_gene_name\tresolution\n"
        "A\tA\tchr1\tENSG_A\tA\tdirect\n"
        "B\tHIST\tchr2\tENSG_B\tNATIVE_B\tcurated_historical_symbol\n"
        "SPLIT\tSPLIT_REVIEW_REQUIRED\tchrX\t.\t.\tselection_unmapped\n",
    )
    _write(selection, "chr1\t11\t21\tA\tENSG_A\nchr2\t31\t53\tHIST\tENSG_B\n")
    lock_path = root / "mapping.lock.yaml"
    lock = {
        "schema_version": "1.0.0",
        "mapping_id": "SYNTHETIC_GRCh38_to_GRCh37_v1",
        "status": "controlled_build_time_unvalidated",
        "source_genome_build": "GRCh38",
        "target_genome_build": "GRCh37",
        "source_panel_bundle_id": "SYNTHETIC_GRCh38_v1",
        "source_panel_resource_id": "source_selection",
        "source_selection_sha256": sha256_file(source),
        "source_interval_count": 3,
        "mapping_method": "ucsc_liftover_chain",
        "mapping_tool": "UCSC liftOver",
        "mapping_tool_url": "https://example.invalid/liftOver",
        "mapping_tool_sha256": "1" * 64,
        "mapping_tool_size_bytes": 1,
        "mapping_tool_build_id": "synthetic",
        "mapping_tool_last_modified": "2026-01-01T00:00:00Z",
        "chain_name": "synthetic.forward.chain.gz",
        "chain_url": "https://example.invalid/forward.chain.gz",
        "chain_md5": "2" * 32,
        "chain_sha256": "3" * 64,
        "min_match": 0.99,
        "min_blocks": 1.0,
        "multiple": False,
        "mapped_output_sha256": sha256_file(mapped),
        "unmapped_output_sha256": sha256_file(unmapped),
        "mapped_interval_count": 2,
        "same_chromosome_mapped_count": 2,
        "span_identical_count": 1,
        "maximum_span_delta_bases": 2,
        "maximum_span_delta_fraction": 0.1,
        "unmapped_target_labels": ["SPLIT"],
        "split_target_labels": ["SPLIT"],
        "roundtrip_chain_name": "synthetic.reverse.chain.gz",
        "roundtrip_chain_url": "https://example.invalid/reverse.chain.gz",
        "roundtrip_chain_md5": "4" * 32,
        "roundtrip_chain_sha256": "5" * 64,
        "roundtrip_min_match": 0.99,
        "roundtrip_min_blocks": 1.0,
        "roundtrip_multiple": False,
        "roundtrip_mapped_output_sha256": sha256_file(roundtrip),
        "roundtrip_mapped_output_size_bytes": roundtrip.stat().st_size,
        "roundtrip_unmapped_output_sha256": sha256_file(roundtrip_unmapped),
        "roundtrip_unmapped_output_size_bytes": roundtrip_unmapped.stat().st_size,
        "roundtrip_mapped_interval_count": 1,
        "reciprocal_exact_interval_count": 1,
        "roundtrip_unmapped_target_labels": ["B"],
        "roundtrip_split_target_labels": ["B"],
        "roundtrip_review_required_target_labels": ["B"],
        "final_selection_sha256": sha256_file(selection),
        "runtime_mapping_allowed": False,
        "note": "synthetic test fixture",
    }
    _write(lock_path, yaml.safe_dump(lock, sort_keys=False))
    resources = [
        {
            "resource_id": "source_selection",
            "role": "coordinate_mapping_source",
            "path": "source/selection.grch38.bed",
            "sha256": sha256_file(source),
            "coordinate_system": "zero_based_half_open",
        },
        {
            "resource_id": "mapping_lock",
            "role": "coordinate_mapping_lock",
            "path": "mapping.lock.yaml",
            "sha256": sha256_file(lock_path),
        },
        {
            "resource_id": "mapped",
            "role": "coordinate_mapping_output",
            "path": "derived/mapped.grch37.bed",
            "sha256": sha256_file(mapped),
            "coordinate_system": "zero_based_half_open",
            "generated": True,
            "derived_from": ["source_selection", "mapping_lock"],
        },
        {
            "resource_id": "unmapped",
            "role": "coordinate_mapping_unmapped",
            "path": "derived/unmapped.txt",
            "sha256": sha256_file(unmapped),
            "coordinate_system": "zero_based_half_open",
            "generated": True,
            "derived_from": ["source_selection", "mapping_lock"],
        },
        {
            "resource_id": "roundtrip",
            "role": "coordinate_mapping_roundtrip_output",
            "path": "derived/roundtrip.grch38.bed",
            "sha256": sha256_file(roundtrip),
            "coordinate_system": "zero_based_half_open",
            "generated": True,
            "derived_from": ["mapped", "mapping_lock"],
        },
        {
            "resource_id": "roundtrip_unmapped",
            "role": "coordinate_mapping_roundtrip_unmapped",
            "path": "derived/roundtrip.unmapped.txt",
            "sha256": sha256_file(roundtrip_unmapped),
            "coordinate_system": "zero_based_half_open",
            "generated": True,
            "derived_from": ["mapped", "mapping_lock"],
        },
        {
            "resource_id": "gene_map",
            "role": "target_gene_map",
            "path": "source/target_gene_map.tsv",
            "sha256": sha256_file(gene_map),
        },
        {
            "resource_id": "selection",
            "role": "selection_panel_buffered",
            "path": "derived/selection.bed",
            "sha256": sha256_file(selection),
            "coordinate_system": "zero_based_half_open",
            "generated": True,
            "derived_from": ["mapped", "gene_map"],
        },
        {
            "resource_id": "roi",
            "role": "analysis_roi_unbuffered",
            "path": "derived/roi.bed",
            "coordinate_system": "zero_based_half_open",
            "generated": True,
            "derived_from": ["selection", "gene_map", "REF:annotation_cache"],
        },
        {
            "resource_id": "transcripts",
            "role": "transcript_cache",
            "path": "derived/transcripts.tsv",
            "generated": True,
            "derived_from": ["selection", "gene_map", "REF:annotation_cache"],
        },
    ]
    return PanelBundle.model_validate(
        {
            "bundle_type": "panel",
            "bundle_id": "SYNTHETIC_GRCh37_v1",
            "version": "v1",
            "genome_build": "GRCh37",
            "assay_mode": "adaptive_sampling",
            "reference_dictionary_contracts": ["exact_full"],
            "resources": resources,
            "selection_panel_resource_id": "selection",
            "analysis_roi_resource_id": "roi",
            "transcript_cache_resource_id": "transcripts",
            "target_gene_map_resource_id": "gene_map",
            "coordinate_mapping_lock_resource_id": "mapping_lock",
            "coordinate_mapping_source_resource_id": "source_selection",
            "coordinate_mapping_output_resource_id": "mapped",
            "coordinate_mapping_unmapped_resource_id": "unmapped",
            "coordinate_mapping_roundtrip_output_resource_id": "roundtrip",
            "coordinate_mapping_roundtrip_unmapped_resource_id": "roundtrip_unmapped",
            "coordinate_mapping_origin_bundle_id": "SYNTHETIC_GRCh38_v1",
            "coordinate_mapping_origin_resource_id": "source_selection",
            "unresolved_targets": ["SPLIT_REVIEW_REQUIRED"],
        }
    )


class Grch37PanelBundleTests(unittest.TestCase):
    def test_official_bundle_is_fully_pinned_and_reports_partial_mapping(self) -> None:
        bundle = PanelBundle.model_validate(
            yaml.safe_load((PANEL / "bundle.yaml").read_text(encoding="utf-8"))
        )
        self.assertEqual(bundle.genome_build, GenomeBuild.GRCH37)
        self.assertEqual(
            bundle.reference_dictionary_contracts,
            [
                ReferenceDictionaryContract.EXACT_FULL,
                ReferenceDictionaryContract.GRCH37_UCSC_HG19_CANONICAL_25,
            ],
        )
        for resource in bundle.resources:
            path = PANEL / resource.path
            self.assertEqual(resource.sha256, sha256_file(path))
            self.assertEqual(resource.size_bytes, path.stat().st_size)

        selection_resource = bundle.resource(bundle.selection_panel_resource_id)
        gene_map_resource = bundle.resource(bundle.target_gene_map_resource_id or "")
        lock = validate_coordinate_mapping_provenance(
            PANEL,
            bundle,
            PANEL / selection_resource.path,
            PANEL / gene_map_resource.path,
        )
        self.assertEqual(lock.source_interval_count, 111)
        self.assertEqual(lock.mapped_interval_count, 110)
        self.assertEqual(lock.same_chromosome_mapped_count, 110)
        self.assertEqual(lock.span_identical_count, 102)
        self.assertEqual(lock.unmapped_target_labels, ["CT45A2"])
        self.assertEqual(lock.reciprocal_exact_interval_count, 109)
        self.assertEqual(lock.roundtrip_unmapped_target_labels, ["ACACA"])
        self.assertEqual(lock.roundtrip_split_target_labels, ["ACACA"])
        self.assertEqual(lock.roundtrip_review_required_target_labels, ["ACACA"])
        self.assertFalse(lock.runtime_mapping_allowed)

        rows = [
            line.split("\t")
            for line in (PANEL / selection_resource.path).read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(len(rows), 110)
        self.assertEqual(sum(row[4] != "." for row in rows), 108)
        labels = {row[3] for row in rows}
        self.assertNotIn("CT45A2", labels)
        self.assertIn("IGH_REVIEW_REQUIRED", labels)
        self.assertIn("GPR128_REVIEW_REQUIRED", labels)

        mappings = read_target_gene_map(PANEL / gene_map_resource.path)
        self.assertEqual(len(mappings), 111)
        aliases = {
            item.target_label: item.native_gene_name
            for item in mappings
            if item.resolution == "curated_historical_symbol"
        }
        self.assertEqual(
            aliases,
            {
                "SEPTIN2": "SEPT2",
                "CIP2A": "KIAA1524",
                "SEPTIN11": "SEPT11",
                "AFDN": "MLLT4",
                "KNL1": "CASC5",
                "SEPTIN9": "SEPT9",
                "SEPTIN5": "SEPT5",
                "SEPTIN3": "SEPT3",
                "SEPTIN6": "SEPT6",
            },
        )
        roi_rows = (
            (PANEL / bundle.resource(bundle.analysis_roi_resource_id).path)
            .read_text(encoding="utf-8")
            .splitlines()
        )
        self.assertEqual(len(roi_rows), 108)
        transcript_rows = (
            (PANEL / bundle.resource(bundle.transcript_cache_resource_id).path)
            .read_text(encoding="utf-8")
            .splitlines()
        )
        self.assertEqual(sum(row.split("\t")[1] == "true" for row in transcript_rows[1:]), 108)

    def test_synthetic_mapping_provenance_fails_closed_after_selection_tamper(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            bundle = _synthetic_mapped_bundle(root)
            selection = root / bundle.resource(bundle.selection_panel_resource_id).path
            gene_map = root / bundle.resource(bundle.target_gene_map_resource_id or "").path
            lock = validate_coordinate_mapping_provenance(root, bundle, selection, gene_map)
            self.assertEqual(lock.reciprocal_exact_interval_count, 1)

            _write(selection, "chr1\t12\t21\tA\tENSG_A\nchr2\t31\t53\tHIST\tENSG_B\n")
            with self.assertRaisesRegex(PanelCompilerError, "final-selection checksum"):
                validate_coordinate_mapping_provenance(root, bundle, selection, gene_map)

    def test_synthetic_mapping_provenance_binds_origin_and_roundtrip_artifacts(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw)
            bundle = _synthetic_mapped_bundle(root)
            selection = root / bundle.resource(bundle.selection_panel_resource_id).path
            gene_map = root / bundle.resource(bundle.target_gene_map_resource_id or "").path

            wrong_origin = bundle.model_copy(
                update={"coordinate_mapping_origin_bundle_id": "OTHER_SOURCE_PANEL"}
            )
            with self.assertRaisesRegex(PanelCompilerError, "source identity"):
                validate_coordinate_mapping_provenance(
                    root,
                    wrong_origin,
                    selection,
                    gene_map,
                )

            roundtrip = (
                root
                / bundle.resource(bundle.coordinate_mapping_roundtrip_output_resource_id or "").path
            )
            _write(roundtrip, roundtrip.read_text(encoding="utf-8") + "chr2\t30\t50\tB\n")
            with self.assertRaisesRegex(PanelCompilerError, "roundtrip mapped output checksum"):
                validate_coordinate_mapping_provenance(root, bundle, selection, gene_map)


if __name__ == "__main__":
    unittest.main()

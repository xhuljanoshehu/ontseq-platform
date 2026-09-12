from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from openpyxl import load_workbook

from ontseq_platform.demo import build_demo_result
from ontseq_platform.models import (
    GenomeBuild,
    PanelResolutionSummary,
    ReferenceDictionaryContract,
    ResolvedResourceContext,
)
from ontseq_platform.report import render_html
from ontseq_platform.workbook import render_workbook


class ReportResourceProvenanceTests(unittest.TestCase):
    def test_resource_releases_and_checksums_render_in_methods_provenance(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            resource = root / "reference.fa"
            resource.write_text(">chr1\nA\n", encoding="utf-8")
            context = ResolvedResourceContext(
                profile_id="AML_LCWGS_GRCh38",
                profile_version="v1",
                genome_build=GenomeBuild.GRCH38,
                reference_dictionary_contract=(ReferenceDictionaryContract.GRCH38_CANONICAL_25),
                reference_bundle_id="GRCh38_GENCODE50_MANE1.5_v1",
                reference_bundle_version="v1",
                knowledge_bundle_id="HEMATOLOGY_v1",
                knowledge_bundle_version="v1",
                resource_root=str(root),
                resource_paths={"reference.genome_fasta": str(resource)},
                resource_checksums={"reference.genome_fasta": "a" * 64},
                resource_releases={
                    "reference.genome_fasta": "GRCh38.p14",
                    "reference.gencode_gtf": "GENCODE 50",
                    "reference.mane_gff3": "MANE 1.5",
                    "reference.cytobands": "UCSC hg38",
                },
            )
            result = build_demo_result().model_copy(update={"reference_context": context})
            html_path = render_html(result, root / "report.html")
            html = html_path.read_text(encoding="utf-8")
            for expected in (
                "GRCh38.p14",
                "GRCh38_GENCODE50_MANE1.5_v1",
                "grch38_canonical_25",
                "GENCODE 50",
                "MANE 1.5",
                "UCSC hg38",
                "HEMATOLOGY_v1",
                "Resource SHA256 provenance",
            ):
                self.assertIn(expected, html)

            workbook_path = render_workbook(result, root / "report.xlsx")
            workbook = load_workbook(workbook_path, read_only=True)
            try:
                values = {
                    row[0].value: row[1].value
                    for row in workbook["00_Summary"].iter_rows(min_row=2)
                    if row[0].value is not None
                }
                self.assertEqual(values["Genome assembly"], "GRCh38.p14")
                self.assertIn("GRCh38_GENCODE50", values["ReferenceBundle"])
                self.assertEqual(
                    values["BAM dictionary contract"],
                    "grch38_canonical_25",
                )
            finally:
                workbook.close()

    def test_panel_mapping_and_unresolved_states_render_with_ruo_interpretation(self) -> None:
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            paths = {
                "panel.selection_panel_buffered": root / "selection.bed",
                "panel.analysis_roi_unbuffered": root / "roi.bed",
                "panel.coordinate_mapping_lock": root / "mapping.lock.yaml",
                "panel.coordinate_mapping_roundtrip_output": root / "roundtrip.bed",
                "panel.coordinate_mapping_roundtrip_unmapped": root / "roundtrip.unmapped.txt",
            }
            for path in paths.values():
                path.write_text("synthetic\n", encoding="utf-8")
            checksums = {
                "panel.selection_panel_buffered": "a" * 64,
                "panel.analysis_roi_unbuffered": "b" * 64,
                "panel.coordinate_mapping_lock": "c" * 64,
                "panel.coordinate_mapping_roundtrip_output": "d" * 64,
                "panel.coordinate_mapping_roundtrip_unmapped": "e" * 64,
            }
            panel_resolution = PanelResolutionSummary(
                mapping_status="controlled_build_time_unvalidated",
                mapping_id="AML_AS_TEST_GRCh37_to_GRCh38_v1",
                mapping_method="ucsc_liftover_chain",
                mapping_tool="UCSC liftOver",
                source_genome_build=GenomeBuild.GRCH37,
                target_genome_build=GenomeBuild.GRCH38,
                source_panel_bundle_id="AML_AS_TEST_GRCh37_v1",
                source_panel_resource_id="selection_panel_normalized",
                selection_interval_count=1,
                analysis_roi_interval_count=1,
                source_interval_count=2,
                mapped_interval_count=1,
                unmapped_target_labels=["CT45A2"],
                roundtrip_mapped_interval_count=1,
                reciprocal_exact_interval_count=1,
                roundtrip_review_required_target_labels=["ACACA"],
                unresolved_target_labels=["CT45A2_REVIEW_REQUIRED", "IGH_REVIEW_REQUIRED"],
                selection_panel_sha256=checksums["panel.selection_panel_buffered"],
                analysis_roi_sha256=checksums["panel.analysis_roi_unbuffered"],
                coordinate_mapping_lock_sha256=checksums["panel.coordinate_mapping_lock"],
                roundtrip_mapping_output_sha256=checksums[
                    "panel.coordinate_mapping_roundtrip_output"
                ],
                roundtrip_mapping_unmapped_sha256=checksums[
                    "panel.coordinate_mapping_roundtrip_unmapped"
                ],
            )
            context = ResolvedResourceContext(
                profile_id="AML_AS_111_GRCh38",
                profile_version="v1",
                genome_build=GenomeBuild.GRCH38,
                reference_bundle_id="GRCh38_GENCODE50_MANE1.5_v1",
                reference_bundle_version="v1",
                panel_bundle_id="AML_AS_111_GRCh38_v1",
                panel_bundle_version="v1",
                panel_resolution=panel_resolution,
                knowledge_bundle_id="HEMATOLOGY_v1",
                knowledge_bundle_version="v1",
                resource_root=str(root),
                resource_paths={name: str(path) for name, path in paths.items()},
                resource_checksums=checksums,
            )
            result = build_demo_result().model_copy(update={"reference_context": context})

            html = render_html(result, root / "report.html").read_text(encoding="utf-8")
            for expected in (
                "Panel coordinate resolution",
                "controlled_build_time_unvalidated",
                "1 of 2",
                "Native analysis ROI intervals",
                "CT45A2",
                "ACACA",
                "mapped does not mean analytically validated",
                "not negative findings for this sample",
            ):
                self.assertIn(expected, html)

            workbook_path = render_workbook(result, root / "report.xlsx")
            workbook = load_workbook(workbook_path, read_only=True)
            try:
                values = {
                    row[0].value: row[1].value
                    for row in workbook["00_Summary"].iter_rows(min_row=2)
                    if row[0].value is not None
                }
                self.assertEqual(
                    values["Panel mapping status"],
                    "controlled_build_time_unvalidated",
                )
                self.assertEqual(values["Panel forward-mapped intervals"], "1 of 2")
                self.assertIn("not negative findings", values["Panel resolution interpretation"])
            finally:
                workbook.close()


if __name__ == "__main__":
    unittest.main()

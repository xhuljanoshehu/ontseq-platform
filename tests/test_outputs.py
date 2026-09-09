from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from openpyxl import load_workbook

from ontseq_platform import __version__
from ontseq_platform.demo import build_demo_result
from ontseq_platform.io import write_json
from ontseq_platform.iscn import build_iscn_proposal, iscn_module_outcome
from ontseq_platform.models import (
    AnalysisModule,
    EventType,
    FusionAnnotation,
    FusionPartnerAnnotation,
    ISCNProposalStatus,
    PathologyAssociation,
    PipelineResult,
)
from ontseq_platform.report import render_html
from ontseq_platform.workbook import render_workbook


class OutputTests(unittest.TestCase):
    def test_demo_round_trip_and_reports(self) -> None:
        result = build_demo_result()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            json_path = write_json(result, root / "result.json")
            html_path = render_html(result, root / "report.html")
            xlsx_path = render_workbook(result, root / "report.xlsx")

            restored = PipelineResult.model_validate(json.loads(json_path.read_text()))
            self.assertEqual(restored.manifest.sample_id, "SYNTHETIC_AML_001")
            self.assertEqual(restored.provenance.pipeline_version, __version__)
            html = html_path.read_text(encoding="utf-8")
            self.assertIn("RESEARCH USE ONLY", html)
            self.assertIn("Technical ISCN proposal", html)
            self.assertIn("PARTIAL_EVENT_LEVEL", html)
            self.assertIn("fachzytogenetisch zu prüfen", html)
            self.assertIn("Source event", html)

            workbook = load_workbook(xlsx_path, read_only=True)
            self.assertEqual(
                workbook.sheetnames,
                [
                    "00_Summary",
                    "01_Key_Findings",
                    "12_AS_Coverage",
                    "01_QC",
                    "02_CNV_Chromosomes",
                    "03_CNV_Segments",
                    "04_SV",
                    "05_Fusions",
                    "06_ISCN",
                    "07_Warnings",
                    "08_Methods_Versions",
                    "09_Run_Log",
                    "10_Module_Status",
                ],
            )
            iscn_rows = list(workbook["06_ISCN"].iter_rows(values_only=True))
            self.assertTrue(
                any(row[:2] == ("Proposal status", "PARTIAL_EVENT_LEVEL") for row in iscn_rows)
            )
            self.assertTrue(any(row[:2] == ("Genome build", "GRCh38") for row in iscn_rows))
            self.assertTrue(
                any(row[:2] == ("Clinical release allowed", False) for row in iscn_rows)
            )
            self.assertTrue(any(row[0] == "CNV-001" and row[1] == "+8" for row in iscn_rows))
            workbook.close()

    def test_fragmentless_workbook_rows_preserve_the_exact_proposal_status(self) -> None:
        baseline = build_demo_result()
        build = baseline.manifest.assay.genome_build
        provenance = baseline.iscn.resource_provenance
        self.assertIsNotNone(provenance)
        cnv_ids = {"CNV-001", "CNV-002"}
        proposals = [
            build_iscn_proposal(
                baseline.events,
                genome_build=build,
                requested=False,
                resource_provenance=provenance,
                cnv_source_event_ids=cnv_ids,
            ),
            build_iscn_proposal(
                baseline.events,
                genome_build=build,
                upstream_evidence_assessed=False,
                resource_provenance=provenance,
                cnv_source_event_ids=cnv_ids,
            ),
            build_iscn_proposal(
                baseline.events,
                genome_build=build,
                resource_provenance=provenance,
                cnv_source_event_ids=set(),
            ),
        ]

        expected_statuses = {
            ISCNProposalStatus.NOT_REQUESTED,
            ISCNProposalStatus.NOT_ASSESSED,
            ISCNProposalStatus.NO_RENDERABLE_CANDIDATE,
        }
        self.assertEqual({proposal.proposal_status for proposal in proposals}, expected_statuses)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for proposal in proposals:
                with self.subTest(proposal_status=proposal.proposal_status):
                    modules = [
                        iscn_module_outcome(proposal)
                        if item.module == AnalysisModule.ISCN
                        else item
                        for item in baseline.modules
                    ]
                    result = PipelineResult.model_validate(
                        baseline.model_copy(
                            update={"iscn": proposal, "modules": modules}
                        ).model_dump(mode="python")
                    )
                    path = render_workbook(
                        result,
                        root / f"{proposal.proposal_status.value}.xlsx",
                    )
                    html_path = render_html(
                        result,
                        root / f"{proposal.proposal_status.value}.html",
                    )
                    workbook = load_workbook(path, read_only=True)
                    rows = list(workbook["06_ISCN"].iter_rows(values_only=True))
                    summary_rows = list(workbook["00_Summary"].iter_rows(values_only=True))
                    workbook.close()

                    self.assertTrue(any(row[0] == proposal.proposal_status.value for row in rows))
                    self.assertFalse(any(row[0] == "NO_CALL" for row in rows))
                    self.assertFalse(
                        any(
                            "NO RENDERABLE EVENT FRAGMENT" in str(cell)
                            for row in rows
                            for cell in row
                        )
                    )
                    self.assertFalse(
                        any(
                            "NO RENDERABLE EVENT FRAGMENT" in str(cell)
                            for row in summary_rows
                            for cell in row
                        )
                    )
                    html = html_path.read_text(encoding="utf-8")
                    if proposal.proposal_status in {
                        ISCNProposalStatus.NOT_REQUESTED,
                        ISCNProposalStatus.NOT_ASSESSED,
                    }:
                        self.assertNotIn("No event-level fragment was renderable", html)

    def test_annotated_translocation_is_included_in_fusion_sheet(self) -> None:
        result = build_demo_result()
        candidate = next(event for event in result.events if event.event_type == EventType.FUSION)
        annotated_bnd = candidate.model_copy(
            update={
                "event_type": EventType.TRANSLOCATION,
                "fusion_evidence": FusionAnnotation(
                    gene_a=FusionPartnerAnnotation(
                        gene="RUNX1T1", preferred_transcript="ENST_RUNX1T1"
                    ),
                    gene_b=FusionPartnerAnnotation(gene="RUNX1", preferred_transcript="ENST_RUNX1"),
                    orientation="+-",
                    frame_status="unknown",
                ),
                "known_pathologies": [
                    PathologyAssociation(
                        disease_id="DOID:0081093",
                        name="Acute Myeloid Leukemia With T(8;21); (q22; Q22.1)",
                        source_id="CIVIC-2026-08-29",
                        source_record_id="CIVIC-FUSION-TEST/DISEASE-TEST",
                        source_url="https://civicdb.org/features/test",
                    )
                ],
            }
        )
        result = result.model_copy(
            update={
                "events": [
                    annotated_bnd if event.event_id == candidate.event_id else event
                    for event in result.events
                ]
            }
        )

        with tempfile.TemporaryDirectory() as temporary:
            path = render_workbook(result, Path(temporary) / "report.xlsx")
            workbook = load_workbook(path, read_only=True)
            rows = list(workbook["05_Fusions"].iter_rows(values_only=True))
            workbook.close()

        self.assertEqual(rows[1][2], "FUS-001")
        self.assertEqual(rows[1][9], "+-")
        self.assertIn("DOID:0081093", rows[1][15])


if __name__ == "__main__":
    unittest.main()

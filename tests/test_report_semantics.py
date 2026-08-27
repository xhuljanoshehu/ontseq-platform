from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.demo import build_demo_result
from ontseq_platform.models import ModuleRunStatus, Verdict
from ontseq_platform.report import render_html
from ontseq_platform.report_view import build_report_view


class ReportSemanticTests(unittest.TestCase):
    def _render(self) -> str:
        result = build_demo_result()
        with tempfile.TemporaryDirectory() as temporary:
            path = render_html(result, Path(temporary) / "report.html")
            return path.read_text(encoding="utf-8")

    def test_portable_report_has_no_remote_runtime_dependency(self) -> None:
        document = self._render()
        self.assertNotIn("https://", document)
        self.assertNotIn("http://", document)
        self.assertNotIn("<script src=", document)
        self.assertIn("offline/self-contained presentation", document)

    def test_no_call_is_explicitly_not_a_negative_result(self) -> None:
        result = build_demo_result()
        result.events = []
        result.modules[2].status = ModuleRunStatus.NO_CALL
        result.modules[2].reason = "Insufficient interpretable evidence in the assessed scope."
        with tempfile.TemporaryDirectory() as temporary:
            document = render_html(result, Path(temporary) / "report.html").read_text(
                encoding="utf-8"
            )
        self.assertIn("NO_CALL", document)
        self.assertIn("not a biological negative result", document)
        self.assertIn("No normalized events were produced", document)

    def test_failed_module_and_qc_failure_are_front_loaded(self) -> None:
        result = build_demo_result()
        result.qc.verdict = Verdict.FAIL
        result.modules[1].status = ModuleRunStatus.FAILED
        result.modules[1].reason = "Synthetic execution failure."
        view = build_report_view(result)
        self.assertEqual(view.alerts[0].title, "QC failed")
        self.assertTrue(any(item.title == "cnv: FAILED" for item in view.alerts))

    def test_zero_evidence_values_are_not_rendered_as_missing(self) -> None:
        result = build_demo_result()
        evidence = result.events[-1].evidence[0]
        evidence.support_reads = 0
        evidence.local_coverage = 0
        evidence.variant_allele_fraction = 0
        evidence.quality = 0
        evidence.precise = False
        view = build_report_view(result)
        normalized = view.events[-1].evidence[0]
        self.assertEqual(normalized.support_reads, 0)
        self.assertEqual(normalized.local_coverage, 0)
        self.assertEqual(normalized.variant_allele_fraction, 0)
        self.assertEqual(normalized.quality, 0)
        self.assertFalse(normalized.precise)
        with tempfile.TemporaryDirectory() as temporary:
            document = render_html(result, Path(temporary) / "report.html").read_text(
                encoding="utf-8"
            )
        self.assertIn("0.0%", document)
        self.assertIn(">False</td>", document)

    def test_reportable_true_keeps_ruo_boundary_visible(self) -> None:
        result = build_demo_result()
        view = build_report_view(result)
        self.assertTrue(view.events[0].reportable)
        self.assertIn("pipeline flag only", view.events[0].reportability_text)
        document = self._render()
        self.assertIn("RESEARCH USE ONLY", document)
        self.assertIn("pipeline flag only; this RUO report is not clinically validated", document)

    def test_fusion_label_does_not_claim_expressed_transcript(self) -> None:
        document = self._render()
        self.assertIn(
            "does not establish an expressed, in-frame or functional fusion transcript",
            document,
        )

    def test_warning_content_is_html_escaped_and_deduplicated(self) -> None:
        result = build_demo_result()
        result.warnings.extend(
            ["<script>alert(1)</script>", "duplicate warning", "duplicate warning"]
        )
        view = build_report_view(result)
        self.assertEqual(view.warnings.count("duplicate warning"), 1)
        with tempfile.TemporaryDirectory() as temporary:
            document = render_html(result, Path(temporary) / "report.html").read_text(
                encoding="utf-8"
            )
        self.assertNotIn("<script>alert(1)</script>", document)
        self.assertIn("&lt;script&gt;alert(1)&lt;/script&gt;", document)

    def test_reference_and_analysis_context_are_visible(self) -> None:
        document = self._render()
        self.assertIn("GRCh38-demo-not-for-analysis", document)
        self.assertIn("synthetic-v1", document)
        self.assertIn("not declared", document)
        self.assertIn("Reference checksums", document)


if __name__ == "__main__":
    unittest.main()

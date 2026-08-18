from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from openpyxl import load_workbook

from ontseq_platform.demo import build_demo_result
from ontseq_platform.io import write_json
from ontseq_platform.models import EventAnnotation, PipelineResult
from ontseq_platform.report import render_html
from ontseq_platform.workbook import ANNOTATION_BANNER, render_workbook

SOURCE_SHA = "0" * 64


def _annotation(**overrides: object) -> EventAnnotation:
    """A minimal annotation. The overrides name exactly what a given test is about."""
    base = EventAnnotation(
        source_id="clinvar",
        source_release="TEST-2026-01-01",
        source_sha256=SOURCE_SHA,
        record_id="VCV000000001",
        record_type="copy number loss",
        assertion="Pathogenic",
        assertion_vocabulary="acmg_germline",
        record_origin="germline",
        scope_alignment="aligned",
        scope_note="record and assay agree on origin",
        match_type="overlap",
        reciprocal_overlap=0.75,
        review_status="criteria provided, single submitter",
        review_stars=1,
        genes=["SYNTHETIC_GENE"],
        conditions=["Synthetic condition"],
        caveats=["This is a classification of a database record, not a finding."],
    )
    return base.model_copy(update=overrides)


def _with_annotations(result: PipelineResult, *annotations: EventAnnotation) -> PipelineResult:
    """Attach annotations to the first event, leaving the rest of the result alone."""
    first = result.events[0].model_copy(update={"annotations": list(annotations)})
    return result.model_copy(update={"events": [first, *result.events[1:]]})


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
            self.assertIn("RESEARCH USE ONLY", html_path.read_text(encoding="utf-8"))

            workbook = load_workbook(xlsx_path, read_only=True)
            self.assertEqual(
                workbook.sheetnames,
                [
                    "00_Summary",
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
                    "11_Annotations",
                ],
            )


class AnnotationSheetTests(unittest.TestCase):
    """The workbook is where a knowledge-base assertion is most easily misread.

    In the HTML report the reading rule sits in prose above the table. A spreadsheet has no
    prose: a reviewer sorts by the assertion column and sees "Pathogenic" beside an event
    identifier. Everything asserted here exists to keep that from reading as a verdict.
    """

    def _render(self, result: PipelineResult, root: Path) -> Path:
        return render_workbook(result, root / "report.xlsx")

    def test_the_reading_rule_is_in_the_sheet_not_only_in_the_html(self) -> None:
        result = _with_annotations(build_demo_result(), _annotation())
        with tempfile.TemporaryDirectory() as temporary:
            sheet = load_workbook(self._render(result, Path(temporary)))["11_Annotations"]
            self.assertEqual(sheet["A1"].value, ANNOTATION_BANNER)
            self.assertIn("not findings about this sample", str(sheet["A1"].value))

    def test_the_header_is_row_two_because_row_one_carries_the_banner(self) -> None:
        result = _with_annotations(build_demo_result(), _annotation())
        with tempfile.TemporaryDirectory() as temporary:
            sheet = load_workbook(self._render(result, Path(temporary)))["11_Annotations"]
            headers = [cell.value for cell in sheet[2]]
            self.assertEqual(headers[0], "event_id")
            self.assertIn("assertion", headers)
            self.assertIn("scope_alignment", headers)
            self.assertIn("caveats", headers)
            self.assertEqual(sheet.freeze_panes, "A3")

    def test_the_annotation_reaches_the_sheet_with_its_caveats(self) -> None:
        result = _with_annotations(build_demo_result(), _annotation())
        with tempfile.TemporaryDirectory() as temporary:
            sheet = load_workbook(self._render(result, Path(temporary)))["11_Annotations"]
            headers = [cell.value for cell in sheet[2]]
            row = {header: cell.value for header, cell in zip(headers, sheet[3], strict=True)}
            self.assertEqual(row["event_id"], result.events[0].event_id)
            self.assertEqual(row["assertion"], "Pathogenic")
            self.assertEqual(row["assertion_vocabulary"], "acmg_germline")
            self.assertTrue(str(row["caveats"]).strip())

    def test_an_unknown_star_rating_is_not_written_as_zero(self) -> None:
        """Zero stars and "this build does not know the status" are different facts."""
        result = _with_annotations(build_demo_result(), _annotation(review_stars=None))
        with tempfile.TemporaryDirectory() as temporary:
            sheet = load_workbook(self._render(result, Path(temporary)))["11_Annotations"]
            headers = [cell.value for cell in sheet[2]]
            row = {header: cell.value for header, cell in zip(headers, sheet[3], strict=True)}
            self.assertEqual(row["review_stars"], "unknown")

    def test_a_scope_mismatch_is_marked_in_the_row_not_only_named_in_a_column(self) -> None:
        """A reviewer scanning a long sheet reads colour before they read column ten."""
        result = _with_annotations(
            build_demo_result(),
            _annotation(record_id="VCV000000002", scope_alignment="mismatched"),
        )
        with tempfile.TemporaryDirectory() as temporary:
            sheet = load_workbook(self._render(result, Path(temporary)))["11_Annotations"]
            self.assertEqual(sheet["A3"].fill.fgColor.rgb, "00FEE2E2")

    def test_an_aligned_row_is_not_marked(self) -> None:
        """Marking everything would be the same as marking nothing."""
        result = _with_annotations(build_demo_result(), _annotation())
        with tempfile.TemporaryDirectory() as temporary:
            sheet = load_workbook(self._render(result, Path(temporary)))["11_Annotations"]
            self.assertNotEqual(sheet["A3"].fill.fgColor.rgb, "00FEE2E2")

    def test_a_result_with_no_annotations_still_gets_the_sheet_and_the_banner(self) -> None:
        """An absent sheet would be ambiguous: nothing matched, or nothing was run?"""
        with tempfile.TemporaryDirectory() as temporary:
            book = load_workbook(self._render(build_demo_result(), Path(temporary)))
            sheet = book["11_Annotations"]
            self.assertEqual(sheet["A1"].value, ANNOTATION_BANNER)
            self.assertEqual(sheet.max_row, 2)

    def test_the_event_sheets_say_how_many_records_matched(self) -> None:
        """Otherwise a reviewer has no reason to open sheet eleven at all."""
        result = _with_annotations(
            build_demo_result(),
            _annotation(),
            _annotation(record_id="VCV000000003"),
        )
        target = result.events[0].event_id
        with tempfile.TemporaryDirectory() as temporary:
            book = load_workbook(self._render(result, Path(temporary)))
            found = False
            for name in ("02_CNV_Chromosomes", "03_CNV_Segments", "04_SV", "05_Fusions"):
                sheet = book[name]
                headers = [cell.value for cell in sheet[1]]
                self.assertIn("db_records_matched", headers)
                column = headers.index("db_records_matched")
                for row in sheet.iter_rows(min_row=2, values_only=True):
                    if row[0] == target:
                        self.assertEqual(row[column], 2)
                        found = True
                    else:
                        self.assertEqual(row[column], 0)
            self.assertTrue(found, "the annotated event appeared on no event sheet")


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.guideline_criteria import COMPUTABLE, DRAFT, Criterion, load_for_review
from ontseq_platform.panel_reachability import (
    DEFAULT_CRITERIA,
    DEFAULT_PANEL,
    GENERATED_REPORT,
    NO_GENES_NAMED,
    NOT_PANEL_GATED,
    PARTIAL,
    REACHABLE,
    UNREACHABLE,
    Panel,
    PanelReachabilityError,
    PanelTarget,
    assess_criterion,
    assess_panel,
    build_default_report,
    format_report,
    load_panel,
)

BED = "\n".join(
    (
        "# a comment",
        "chr1\t100\t200\tAAA",
        "",
        "chr2\t300\t900\tBBB",
    )
)


def _criterion(
    record_id: str = "C1",
    *,
    genes: tuple[str, ...] = ("AAA",),
    detectable_by: tuple[str, ...] = ("small_variant",),
) -> Criterion:
    return Criterion(
        record_id=record_id,
        category="eln2022_adverse",
        display_name=f"criterion {record_id}",
        pattern_type="gene_mutation",
        detectable_by=detectable_by,
        assay_status=COMPUTABLE,
        verification=DRAFT,
        guideline_reference=None,
        reviewer_note="",
        caveat="",
        genes=genes,
    )


def _panel(*labels: str) -> Panel:
    targets = tuple(
        PanelTarget(chrom="chr1", start=index * 1000, end=index * 1000 + 500, label=label)
        for index, label in enumerate(labels)
    )
    return Panel(path=Path("synthetic.bed"), targets=targets)


class LoadPanelTests(unittest.TestCase):
    def _write(self, text: str) -> Path:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / "panel.bed"
        path.write_text(text, encoding="utf-8")
        return path

    def test_skips_comments_and_blank_lines(self) -> None:
        panel = load_panel(self._write(BED))
        self.assertEqual(panel.labels, ("AAA", "BBB"))
        self.assertEqual(panel.total_span, 100 + 600)

    def test_rejects_row_without_label(self) -> None:
        with self.assertRaises(PanelReachabilityError):
            load_panel(self._write("chr1\t100\t200"))

    def test_rejects_non_numeric_coordinates(self) -> None:
        with self.assertRaises(PanelReachabilityError):
            load_panel(self._write("chr1\tstart\t200\tAAA"))

    def test_rejects_reversed_interval(self) -> None:
        with self.assertRaises(PanelReachabilityError):
            load_panel(self._write("chr1\t400\t200\tAAA"))

    def test_rejects_file_with_no_intervals(self) -> None:
        with self.assertRaises(PanelReachabilityError):
            load_panel(self._write("# only a comment\n"))

    def test_rejects_unreadable_path(self) -> None:
        with self.assertRaises(PanelReachabilityError):
            load_panel(Path("/nonexistent/panel.bed"))


class AssessCriterionTests(unittest.TestCase):
    def test_all_genes_present_is_reachable(self) -> None:
        reach = assess_criterion(_criterion(genes=("AAA", "BBB")), _panel("AAA", "BBB"))
        self.assertEqual(reach.status, REACHABLE)
        self.assertTrue(reach.reachable)
        self.assertEqual(reach.genes_missing, ())

    def test_no_gene_present_is_unreachable(self) -> None:
        reach = assess_criterion(_criterion(genes=("XXX",)), _panel("AAA"))
        self.assertEqual(reach.status, UNREACHABLE)
        self.assertFalse(reach.reachable)

    def test_partial_coverage_is_not_reachable(self) -> None:
        """The central pin: a criterion answered from a subset of its genes is wrong, not partial.

        Evaluating a nine-gene criterion against the one gene that happens to be targeted
        reports absence for eight genes that were never sequenced deeply enough to see, and
        that error runs one way -- towards favourable.
        """
        reach = assess_criterion(_criterion(genes=("AAA", "XXX", "YYY")), _panel("AAA"))
        self.assertEqual(reach.status, PARTIAL)
        self.assertFalse(reach.reachable)
        self.assertEqual(reach.genes_in_panel, ("AAA",))
        self.assertEqual(reach.genes_missing, ("XXX", "YYY"))
        self.assertIn("XXX", reach.reason())
        self.assertIn("YYY", reach.reason())

    def test_criterion_not_needing_small_variants_gets_no_verdict(self) -> None:
        reach = assess_criterion(
            _criterion(genes=("XXX",), detectable_by=("cnv", "sv_breakpoint")), _panel("AAA")
        )
        self.assertEqual(reach.status, NOT_PANEL_GATED)
        self.assertTrue(reach.reachable)
        self.assertFalse(reach.needs_small_variants)
        self.assertIn("no verdict", reach.reason().lower())

    def test_small_variant_criterion_naming_no_genes(self) -> None:
        reach = assess_criterion(_criterion(genes=()), _panel("AAA"))
        self.assertEqual(reach.status, NO_GENES_NAMED)
        self.assertFalse(reach.reachable)

    def test_genes_to_add_excludes_targeted_genes(self) -> None:
        report = assess_panel(
            [_criterion("C1", genes=("AAA", "XXX")), _criterion("C2", genes=("YYY",))],
            _panel("AAA"),
            bundle_id="B",
        )
        self.assertEqual(report.genes_to_add, ("XXX", "YYY"))

    def test_genes_to_add_ignores_criteria_not_gated_by_the_panel(self) -> None:
        report = assess_panel(
            [_criterion("C1", genes=("ZZZ",), detectable_by=("cnv",))], _panel("AAA"), bundle_id="B"
        )
        self.assertEqual(report.genes_to_add, ())


class ParserAgreementTests(unittest.TestCase):
    def test_labels_agree_with_panel_lock(self) -> None:
        """Two parsers read the same BED, so pin them to the same answer.

        ``panel_lock`` is the authority but imports yaml; this module is dependency-free so it
        can run wherever the criteria table can. The duplication is deliberate, so it needs a
        test that fails if the two ever disagree.
        """
        try:
            from ontseq_platform.pipeline.panel_lock import target_labels
        except ImportError as error:  # pragma: no cover - only when yaml is absent
            self.skipTest(f"panel_lock unavailable: {error}")
        self.assertEqual(load_panel(DEFAULT_PANEL).labels, target_labels(DEFAULT_PANEL))


class ShippedDesignTests(unittest.TestCase):
    """Pin what the shipped panel does and does not reach.

    This is the drift guard. If either the panel BED or the criteria bundle changes, these
    fail and somebody has to look at whether the change moved a criterion across the line.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls.report, cls.panel = build_default_report()

    def test_every_criterion_is_assessed(self) -> None:
        bundle = load_for_review(DEFAULT_CRITERIA)
        self.assertEqual(len(self.report.reaches), len(bundle.criteria))

    def test_seven_criteria_need_small_variants(self) -> None:
        self.assertEqual(len(self.report.small_variant_criteria), 7)

    def test_npm1_and_flt3_are_already_targeted(self) -> None:
        reachable = {
            item.record_id for item in self.report.small_variant_criteria if item.reachable
        }
        self.assertEqual(len(reachable), 3)
        for record_id in reachable:
            reach = next(i for i in self.report.reaches if i.record_id == record_id)
            self.assertEqual(reach.status, REACHABLE)
            self.assertIn("NPM1", reach.genes_in_panel)

    def test_mds_gene_set_is_only_partially_targeted(self) -> None:
        partial = self.report.by_status(PARTIAL)
        self.assertEqual(len(partial), 1)
        self.assertEqual(partial[0].genes_in_panel, ("RUNX1",))
        self.assertEqual(len(partial[0].genes_missing), 8)
        self.assertFalse(partial[0].reachable)

    def test_tp53_and_cebpa_are_not_targeted_at_all(self) -> None:
        unreachable = {
            gene for item in self.report.by_status(UNREACHABLE) for gene in item.genes_missing
        }
        self.assertEqual(unreachable, {"TP53", "CEBPA"})

    def test_missing_gene_list_is_the_one_the_laboratory_needs(self) -> None:
        self.assertEqual(
            self.report.genes_to_add,
            ("ASXL1", "BCOR", "CEBPA", "EZH2", "SF3B1", "SRSF2", "STAG2", "TP53", "U2AF1", "ZRSR2"),
        )
        self.assertNotIn("RUNX1", self.report.genes_to_add)
        self.assertNotIn("NPM1", self.report.genes_to_add)

    def test_report_names_the_blocked_criteria_and_the_missing_genes(self) -> None:
        text = format_report(self.report, self.panel)
        self.assertIn("TP53", text)
        self.assertIn("CEBPA", text)
        self.assertIn("AS_FUSION_PANEL_V1_UNCONFIRMED", text)
        self.assertIn("GENERATED FILE", text)

    def test_generated_report_on_disk_is_current(self) -> None:
        self.assertTrue(
            GENERATED_REPORT.exists(), "run python -m ontseq_platform.panel_reachability"
        )
        expected = format_report(self.report, self.panel) + "\n"
        self.assertEqual(
            GENERATED_REPORT.read_text(encoding="utf-8"),
            expected,
            "docs/PANEL_REACHABILITY.md is stale; regenerate it",
        )


if __name__ == "__main__":
    unittest.main()

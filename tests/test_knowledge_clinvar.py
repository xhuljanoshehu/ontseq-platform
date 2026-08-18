from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.knowledge.annotate import AnnotationSource, annotate_finding, summarize
from ontseq_platform.knowledge.clinvar import ClinVarError, load, release_sha256
from ontseq_platform.knowledge.scope import Interval, MatchType, Origin, ScopeAlignment

COLUMNS = [
    "#AlleleID", "Type", "Name", "GeneSymbol", "ClinicalSignificance", "OriginSimple",
    "Assembly", "Chromosome", "Start", "Stop", "ReviewStatus", "PhenotypeList", "VariationID",
]


def row(**overrides: str) -> list[str]:
    values = {
        "#AlleleID": "1001",
        "Type": "copy number loss",
        "Name": "GRCh38/hg38 7q22.1(chr7:98000000-102000000)x1",
        "GeneSymbol": "CUX1;EZH2",
        "ClinicalSignificance": "Pathogenic",
        "OriginSimple": "germline",
        "Assembly": "GRCh38",
        "Chromosome": "7",
        "Start": "98000001",
        "Stop": "102000000",
        "ReviewStatus": "reviewed by expert panel",
        "PhenotypeList": "Myelodysplastic syndrome;not provided",
        "VariationID": "50001",
    }
    values.update(overrides)
    return [values[name] for name in COLUMNS]


class ClinVarCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.path = Path(self._temporary.name) / "variant_summary.txt"

    def write(self, *rows: list[str], columns: list[str] | None = None) -> Path:
        header = columns if columns is not None else COLUMNS
        lines = ["\t".join(header)] + ["\t".join(item) for item in rows]
        self.path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return self.path


class LoadTests(ClinVarCase):
    def test_a_usable_record_is_kept(self) -> None:
        records, summary = load(self.write(row()), assembly="GRCh38")
        self.assertEqual(len(records), 1)
        self.assertEqual(summary.kept, 1)
        self.assertEqual(records[0].variation_id, "50001")

    def test_start_is_converted_to_half_open_coordinates(self) -> None:
        """ClinVar is 1-based inclusive; every interval in this codebase is half-open."""
        records, _ = load(self.write(row()), assembly="GRCh38")
        self.assertEqual(records[0].interval.start, 98_000_000)
        self.assertEqual(records[0].interval.end, 102_000_000)

    def test_rows_for_another_assembly_are_dropped_and_counted(self) -> None:
        """A silent assembly mismatch looks exactly like a sample nobody knows anything about."""
        records, summary = load(self.write(row(Assembly="GRCh37")), assembly="GRCh38")
        self.assertEqual(records, [])
        self.assertEqual(summary.wrong_assembly, 1)

    def test_variants_without_a_matchable_extent_are_counted_not_matched(self) -> None:
        records, summary = load(
            self.write(row(Type="single nucleotide variant")), assembly="GRCh38"
        )
        self.assertEqual(records, [])
        self.assertEqual(summary.not_a_region_type, 1)

    def test_unplaced_records_are_counted_not_crashed_on(self) -> None:
        """ClinVar writes -1 for records whose placement is unknown."""
        records, summary = load(self.write(row(Start="-1", Stop="-1")), assembly="GRCh38")
        self.assertEqual(records, [])
        self.assertEqual(summary.unusable_coordinates, 1)

    def test_the_summary_says_what_was_skipped(self) -> None:
        _, summary = load(
            self.write(row(), row(Assembly="GRCh37"), row(Type="deletion", VariationID="2")),
            assembly="GRCh38",
        )
        self.assertIn("3 row(s) read", summary.describe())
        self.assertIn("1 for another assembly", summary.describe())

    def test_multi_valued_fields_are_split_and_placeholders_dropped(self) -> None:
        records, _ = load(self.write(row()), assembly="GRCh38")
        self.assertEqual(records[0].genes, ("CUX1", "EZH2"))
        self.assertEqual(records[0].conditions, ("Myelodysplastic syndrome",))

    def test_origin_is_read_from_the_record(self) -> None:
        records, _ = load(self.write(row(OriginSimple="somatic")), assembly="GRCh38")
        self.assertIs(records[0].origin, Origin.SOMATIC)

    def test_an_ambiguous_origin_is_unknown(self) -> None:
        records, _ = load(self.write(row(OriginSimple="germline/somatic")), assembly="GRCh38")
        self.assertIs(records[0].origin, Origin.UNKNOWN)

    def test_an_unrecognised_review_status_is_recorded_for_inspection(self) -> None:
        _, summary = load(self.write(row(ReviewStatus="future wording")), assembly="GRCh38")
        self.assertIn("future wording", summary.unknown_review_status)

    def test_a_file_missing_required_columns_is_refused_by_name(self) -> None:
        columns = [name for name in COLUMNS if name != "OriginSimple"]
        values = [item for name, item in zip(COLUMNS, row(), strict=True) if name != "OriginSimple"]
        with self.assertRaises(ClinVarError) as caught:
            load(self.write(values, columns=columns), assembly="GRCh38")
        self.assertIn("OriginSimple", str(caught.exception))

    def test_an_unknown_assembly_is_refused(self) -> None:
        with self.assertRaises(ClinVarError):
            load(self.write(row()), assembly="T2T-CHM13")

    def test_a_missing_release_is_refused(self) -> None:
        with self.assertRaises(ClinVarError):
            load(self.path.parent / "absent.txt", assembly="GRCh38")

    def test_the_release_is_checksummed_so_an_annotation_can_name_it(self) -> None:
        digest = release_sha256(self.write(row()))
        self.assertEqual(len(digest), 64)
        self.assertEqual(digest, release_sha256(self.path))


class AnnotationTests(ClinVarCase):
    SOURCE = AnnotationSource(source_id="clinvar", release="2026-08-01", sha256="a" * 64)

    def annotate(self, finding: Interval, *rows: list[str], intent: Origin = Origin.SOMATIC):
        records, _ = load(self.write(*(rows or (row(),))), assembly="GRCh38")
        return annotate_finding(finding, records, source=self.SOURCE, assay_intent=intent)

    def test_a_matching_record_is_attached_with_its_provenance(self) -> None:
        found = self.annotate(Interval("chr7", 98_000_000, 102_000_000))
        self.assertEqual(len(found), 1)
        self.assertEqual(found[0].source_release, "2026-08-01")
        self.assertEqual(found[0].source_sha256, "a" * 64)
        self.assertIs(found[0].match_type, MatchType.EXACT)

    def test_the_assertion_stays_in_clinvars_own_vocabulary(self) -> None:
        """`Pathogenic` means germline causation; it must not read as a somatic driver."""
        found = self.annotate(Interval("chr7", 98_000_000, 102_000_000))
        self.assertEqual(found[0].assertion, "Pathogenic")
        self.assertEqual(found[0].assertion_vocabulary, "acmg_germline")

    def test_a_germline_record_under_a_somatic_assay_is_kept_and_marked(self) -> None:
        """Dropping it would be a clinical decision disguised as a filter."""
        found = self.annotate(Interval("chr7", 98_000_000, 102_000_000))
        self.assertIs(found[0].scope_alignment, ScopeAlignment.MISMATCHED)
        self.assertTrue(any("secondary finding" in note for note in found[0].caveats))

    def test_every_annotation_says_it_does_not_make_anything_reportable(self) -> None:
        found = self.annotate(Interval("chr7", 98_000_000, 102_000_000))
        self.assertTrue(any("not make anything reportable" in note for note in found[0].caveats))

    def test_a_weak_assertion_carries_that_caveat(self) -> None:
        found = self.annotate(
            Interval("chr7", 98_000_000, 102_000_000),
            row(ReviewStatus="criteria provided, single submitter"),
        )
        self.assertTrue(any("evidence is weak" in note for note in found[0].caveats))

    def test_a_strong_exact_match_carries_no_weakness_caveat(self) -> None:
        found = self.annotate(Interval("chr7", 98_000_000, 102_000_000))
        self.assertFalse(any("evidence is weak" in note for note in found[0].caveats))
        self.assertFalse(any("matched by" in note for note in found[0].caveats))

    def test_a_containment_match_says_it_may_be_a_different_event(self) -> None:
        found = self.annotate(Interval("chr7", 0, 159_000_000))
        self.assertIs(found[0].match_type, MatchType.RECORD_WITHIN_FINDING)
        self.assertTrue(any("different event in the same region" in n for n in found[0].caveats))

    def test_a_distant_finding_matches_nothing(self) -> None:
        self.assertEqual(self.annotate(Interval("chr7", 0, 1_000_000)), [])

    def test_exact_matches_are_ordered_before_weaker_ones(self) -> None:
        found = self.annotate(
            Interval("chr7", 98_000_000, 102_000_000),
            row(VariationID="weak", Start="1", Stop="159000000"),
            row(VariationID="exact"),
        )
        self.assertEqual([item.record_id for item in found], ["exact", "weak"])

    def test_the_summary_names_the_secondary_findings_and_the_limit(self) -> None:
        text = summarize(self.annotate(Interval("chr7", 98_000_000, 102_000_000)))
        self.assertIn("different origin (secondary findings)", text)
        self.assertIn("None of this makes a finding reportable", text)

    def test_no_match_says_so_plainly(self) -> None:
        self.assertIn("no knowledge-base record matched", summarize([]))

    def test_an_undeclared_assay_intent_leaves_the_scope_unchecked(self) -> None:
        found = self.annotate(Interval("chr7", 98_000_000, 102_000_000), intent=Origin.UNKNOWN)
        self.assertIs(found[0].scope_alignment, ScopeAlignment.UNKNOWN)


if __name__ == "__main__":
    unittest.main()

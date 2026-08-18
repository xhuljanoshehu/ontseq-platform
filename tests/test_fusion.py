from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.fusion import (
    FusionClassification,
    FusionGenePair,
    GeneAnnotationIndex,
    ObservabilityRegion,
    ObservabilityStatus,
    candidate_from_event,
    interpret_sniffles_fusions,
)
from ontseq_platform.models import (
    EventType,
    Evidence,
    FileFingerprint,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    SnifflesCallReport,
    SnifflesPolicy,
    ToolRecord,
)


def _event() -> GenomicEvent:
    return GenomicEvent(
        event_id="SNIFFLES2-000001",
        event_type=EventType.TRANSLOCATION,
        primary=Locus(chromosome="chr9", start=133_729_450, end=133_729_451),
        secondary=Locus(chromosome="chr22", start=23_632_500, end=23_632_501),
        evidence=[
            Evidence(
                caller="Sniffles2",
                caller_version="2.8.0",
                support_reads=12,
                variant_allele_fraction=0.35,
                precise=True,
            )
        ],
        confidence="unclassified",
        reportable=False,
    )


def _report(event: GenomicEvent | None = None) -> SnifflesCallReport:
    events = [event] if event is not None else []
    return SnifflesCallReport(
        sample_id="SYNTHETIC_FUSION_001",
        genome_build=GenomeBuild.GRCH38,
        status=ModuleRunStatus.COMPLETED if events else ModuleRunStatus.NO_CALL,
        policy=SnifflesPolicy(
            profile_id="synthetic-fusion",
            status="technical_defaults_only",
            note="Synthetic technical policy only.",
        ),
        events=events,
        raw_record_count=len(events),
        accepted_record_count=len(events),
        rejected_record_count=0,
        rejection_counts={},
        tool=ToolRecord(name="Sniffles2", version="2.8.0"),
        vcf_fingerprint=FileFingerprint(size_bytes=1, sha256="0" * 64),
    )


class FusionInterpretationTests(unittest.TestCase):
    def _annotation(self, directory: str) -> GeneAnnotationIndex:
        path = Path(directory) / "genes.bed"
        path.write_text(
            "chr9\t133729000\t133730000\tABL1\t0\t+\tENST_SYNTH_ABL1\n"
            "chr22\t23632000\t23633000\tBCR\t0\t+\tENST_SYNTH_BCR\n",
            encoding="utf-8",
        )
        return GeneAnnotationIndex.from_bed(
            path,
            resource_id="synthetic-genes",
            resource_version="v1",
            genome_build=GenomeBuild.GRCH38,
        )

    def test_gene_gene_candidate_is_research_only_and_orientation_unresolved(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            annotation = self._annotation(directory)
            candidate = candidate_from_event(
                _event(),
                annotation,
                observability=[
                    ObservabilityRegion(
                        chromosome="chr9",
                        start=133_729_000,
                        end=133_730_000,
                        status=ObservabilityStatus.OBSERVABLE,
                        reason="synthetic target observed",
                    ),
                    ObservabilityRegion(
                        chromosome="chr22",
                        start=23_632_000,
                        end=23_633_000,
                        status=ObservabilityStatus.OBSERVABLE,
                        reason="synthetic target observed",
                    ),
                ],
                known_pairs={("BCR", "ABL1")},
            )

        self.assertEqual(candidate.classification, FusionClassification.GENE_GENE)
        self.assertFalse(candidate.reportable)
        self.assertTrue(candidate.research_only)
        self.assertEqual(candidate.primary.observability, ObservabilityStatus.OBSERVABLE)
        self.assertEqual(candidate.secondary.observability, ObservabilityStatus.OBSERVABLE)
        self.assertEqual(len(candidate.gene_pairs), 1)
        pair = candidate.gene_pairs[0]
        self.assertEqual({pair.gene_a, pair.gene_b}, {"ABL1", "BCR"})
        self.assertTrue(pair.known_pair)
        self.assertFalse(pair.orientation_resolved)
        self.assertIsNone(pair.gene_5prime)
        self.assertIsNone(pair.gene_3prime)
        self.assertEqual(candidate.evidence[0].support_reads, 12)

    def test_unobserved_partner_is_not_interpreted_as_negative(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            annotation = self._annotation(directory)
            candidate = candidate_from_event(
                _event(),
                annotation,
                observability=[
                    ObservabilityRegion(
                        chromosome="chr9",
                        start=133_729_000,
                        end=133_730_000,
                        status=ObservabilityStatus.OBSERVABLE,
                        reason="synthetic target observed",
                    )
                ],
            )

        self.assertEqual(candidate.secondary.observability, ObservabilityStatus.UNKNOWN)
        negative_limitation = any(
            "cannot be interpreted as a biological negative" in item
            for item in candidate.limitations
        )
        self.assertTrue(negative_limitation)

    def test_report_preserves_annotation_provenance_and_counts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            annotation = self._annotation(directory)
            report = interpret_sniffles_fusions(_report(_event()), annotation)

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.source_translocation_count, 1)
        self.assertEqual(len(report.candidates), 1)
        self.assertIsNotNone(report.annotation_source_sha256)
        assert report.annotation_source_sha256 is not None
        self.assertEqual(len(report.annotation_source_sha256), 64)
        self.assertTrue(report.research_only)

    def test_no_bnd_source_is_no_call_not_validated_negative(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            annotation = self._annotation(directory)
            report = interpret_sniffles_fusions(_report(), annotation)

        self.assertEqual(report.status, ModuleRunStatus.NO_CALL)
        self.assertTrue(any("not a validated negative" in item for item in report.warnings))

    def test_annotation_build_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "genes.bed"
            path.write_text("chr9\t1\t10\tGENE1\t0\t+\n", encoding="utf-8")
            annotation = GeneAnnotationIndex.from_bed(
                path,
                resource_id="wrong-build",
                resource_version="v1",
                genome_build=GenomeBuild.GRCH37,
            )
            with self.assertRaisesRegex(ValueError, "different genome builds"):
                interpret_sniffles_fusions(_report(_event()), annotation)

    def test_unresolved_pair_cannot_claim_5prime_3prime_direction(self) -> None:
        with self.assertRaisesRegex(ValueError, "must not assign"):
            FusionGenePair(
                gene_a="BCR",
                gene_b="ABL1",
                gene_5prime="BCR",
                gene_3prime="ABL1",
                orientation_resolved=False,
            )


if __name__ == "__main__":
    unittest.main()

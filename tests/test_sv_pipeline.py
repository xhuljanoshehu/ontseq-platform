from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.aml_rearrangements import prioritize_aml_rearrangements
from ontseq_platform.cutesv import normalize_cutesv_vcf
from ontseq_platform.demo import build_demo_result
from ontseq_platform.models import (
    AmlKnowledgeLock,
    AssayMode,
    CuteSvPolicy,
    EventType,
    Evidence,
    FileFingerprint,
    GenomeBuild,
    GenomicEvent,
    IntervalResourceLock,
    Locus,
    ModuleRunStatus,
    SvConsensusPolicy,
    SvObservability,
    TargetBedRole,
    ToolRecord,
)
from ontseq_platform.reference import sha256_file
from ontseq_platform.report import render_html
from ontseq_platform.sv_annotation import annotate_sv_events
from ontseq_platform.sv_consensus import consolidate_sv_events
from ontseq_platform.sv_observability import apply_sv_observability
from ontseq_platform.target_coverage import (
    TargetCoveragePolicy,
    TargetCoverageRegion,
    TargetCoverageReport,
)


def _event(
    event_id: str,
    primary: tuple[str, int],
    secondary: tuple[str, int] | None,
    caller: str,
    orientation: str | None = None,
) -> GenomicEvent:
    return GenomicEvent(
        event_id=event_id,
        event_type=EventType.TRANSLOCATION,
        primary=Locus(chromosome=primary[0], start=primary[1], end=primary[1] + 1),
        secondary=(
            Locus(chromosome=secondary[0], start=secondary[1], end=secondary[1] + 1)
            if secondary
            else None
        ),
        evidence=[
            Evidence(
                caller=caller,
                caller_version="test",
                support_reads=20,
                supporting_read_strands=orientation,
            )
        ],
    )


class CuteSvAndConsensusTests(unittest.TestCase):
    def test_cutesv_vcf_normalizes_support_vaf_and_breakend(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "calls.vcf"
            path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr21\t9827581\tc1\tN\tN]chr2:133011913]\t60\tPASS\t"
                "SVTYPE=BND;RE=30;PRECISE\tGT:DR:DV\t0/1:20:30\n",
                encoding="utf-8",
            )
            policy = CuteSvPolicy(
                profile_id="test",
                status="technical_defaults_only",
                note="test-only technical defaults",
            )
            report = normalize_cutesv_vcf(
                path,
                sample_id="SYNTHETIC_B418",
                genome_build=GenomeBuild.GRCH37,
                policy=policy,
                tool=ToolRecord(name="cuteSV", version="2.1.3"),
            )
            self.assertEqual(report.accepted_record_count, 1)
            event = report.events[0]
            self.assertEqual(event.secondary.chromosome if event.secondary else None, "chr2")
            self.assertEqual(event.evidence[0].support_reads, 30)
            self.assertAlmostEqual(event.evidence[0].variant_allele_fraction or 0, 0.6)
            self.assertFalse(event.reportable)

    def test_reversed_bnd_representations_merge_across_callers(self) -> None:
        policy = SvConsensusPolicy(
            profile_id="test",
            status="technical_defaults_only",
            maximum_breakpoint_distance_bp=100,
            note="test-only technical defaults",
        )
        merged = consolidate_sv_events(
            [
                _event("sniffles", ("chr2", 133_011_912), ("chr21", 9_827_580), "Sniffles2"),
                _event("cutesv", ("chr21", 9_827_600), ("chr2", 133_011_950), "cuteSV"),
            ],
            policy,
        )
        self.assertEqual(len(merged), 1)
        self.assertEqual({item.caller for item in merged[0].evidence}, {"Sniffles2", "cuteSV"})
        self.assertEqual(merged[0].source_event_ids, ["cutesv", "sniffles"])
        self.assertFalse(merged[0].reportable)

    def test_reversed_bnd_orientation_is_compared_in_canonical_order(self) -> None:
        policy = SvConsensusPolicy(
            profile_id="orientation-test",
            status="technical_defaults_only",
            maximum_breakpoint_distance_bp=100,
            require_orientation_when_available=True,
            note="test-only technical defaults",
        )
        merged = consolidate_sv_events(
            [
                _event("sniffles", ("chr2", 1000), ("chr21", 2000), "Sniffles2", "+-"),
                _event("cutesv", ("21", 2010), ("2", 1010), "cuteSV", "-+"),
            ],
            policy,
        )
        self.assertEqual(len(merged), 1)

    def test_nearby_same_caller_representations_cluster_without_fake_consensus(self) -> None:
        policy = SvConsensusPolicy(
            profile_id="within-caller-test",
            status="technical_defaults_only",
            maximum_breakpoint_distance_bp=100,
            merge_within_caller=True,
            note="test-only technical defaults",
        )
        merged = consolidate_sv_events(
            [
                _event("b418-a", ("chr2", 1000), ("chr21", 2000), "Sniffles2"),
                _event("b418-b", ("chr2", 1040), ("chr21", 2040), "Sniffles2"),
                _event("b418-c", ("chr2", 1080), ("chr21", 2080), "Sniffles2"),
            ],
            policy,
        )
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0].source_event_ids, ["b418-a", "b418-b", "b418-c"])
        self.assertEqual(merged[0].validation_status.value, "detected")

    def test_two_nearby_but_distinct_events_remain_separate(self) -> None:
        policy = SvConsensusPolicy(
            profile_id="test",
            status="technical_defaults_only",
            maximum_breakpoint_distance_bp=100,
            note="test-only technical defaults",
        )
        merged = consolidate_sv_events(
            [
                _event("one", ("chr2", 1000), ("chr21", 2000), "Sniffles2"),
                _event("two", ("chr2", 1200), ("chr21", 2200), "cuteSV"),
            ],
            policy,
        )
        self.assertEqual(len(merged), 2)


class AnnotationAndPrioritizationTests(unittest.TestCase):
    def test_picalm_mllt10_match_is_order_independent_and_uses_standard_display(self) -> None:
        root = Path(__file__).parents[1] / "configs" / "knowledge_bundles" / "HEMATOLOGY_v3"
        resource = root / "hematology_rearrangements.v0.3.json"
        lock = AmlKnowledgeLock.model_validate_json(
            (root / "hematology_rearrangements.v0.3.lock.json").read_text(encoding="utf-8")
        )
        genomic_order = _event(
            "picalm-mllt10",
            ("chr10", 21_634_899),
            ("chr11", 85_975_045),
            "Sniffles2",
        ).model_copy(update={"genes": ["MLLT10", "PICALM"]})

        prioritized = prioritize_aml_rearrangements(
            [genomic_order], resource_path=resource, lock=lock
        )[0]

        self.assertEqual(prioritized.known_rearrangement, "PICALM::MLLT10")
        self.assertEqual(prioritized.fusion_status.value, "fusion_candidate")
        self.assertEqual(prioritized.validation_status.value, "biologically_prioritized")
        self.assertEqual(
            [(item.disease_id, item.name) for item in prioritized.known_pathologies],
            [("DOID:9119", "Acute Myeloid Leukemia")],
        )
        self.assertFalse(prioritized.reportable)

    def test_build_locked_annotation_context_and_aml_pattern(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genes = root / "genes.tsv"
            bands = root / "bands.tsv"
            repeats = root / "repeats.tsv"
            genes.write_text(
                "chr8\t92000000\t92100000\tRUNX1T1\nchr21\t36000000\t36100000\tRUNX1\n",
                encoding="utf-8",
            )
            bands.write_text(
                "chr8\t90000000\t95000000\tq21.3\nchr21\t35000000\t38000000\tq22.12\n",
                encoding="utf-8",
            )
            repeats.write_text("chr8\t92000000\t92000100\tLINE\n", encoding="utf-8")

            def lock(path: Path, kind: str) -> IntervalResourceLock:
                return IntervalResourceLock(
                    resource_id=path.stem,
                    resource_type=kind,
                    source="synthetic",
                    release="test",
                    genome_build=GenomeBuild.GRCH37,
                    sha256=sha256_file(path),
                )

            event = _event("candidate", ("chr8", 92_000_010), ("chr21", 36_000_010), "Sniffles2")
            annotated = annotate_sv_events(
                [event],
                genome_build=GenomeBuild.GRCH37,
                gene_resource=(genes, lock(genes, "genes")),
                cytoband_resource=(bands, lock(bands, "cytobands")),
                context_resources=[(repeats, lock(repeats, "repeatmasker"))],
            )[0]
            self.assertEqual(annotated.genes, ["RUNX1", "RUNX1T1"])
            self.assertEqual(annotated.primary.cytoband_start, "q21.3")
            self.assertIn("primary:repeatmasker", annotated.technical_flags)

            knowledge_path = (
                Path(__file__).parents[1] / "configs" / "knowledge" / "aml_rearrangements.v0.1.json"
            )
            lock_path = knowledge_path.with_name("aml_rearrangements.v0.1.lock.json")
            knowledge_lock = AmlKnowledgeLock.model_validate_json(
                lock_path.read_text(encoding="utf-8")
            )
            prioritized = prioritize_aml_rearrangements(
                [annotated], resource_path=knowledge_path, lock=knowledge_lock
            )[0]
            self.assertEqual(prioritized.known_rearrangement, "RUNX1::RUNX1T1")
            self.assertEqual(prioritized.fusion_status.value, "fusion_candidate")
            self.assertEqual(prioritized.validation_status.value, "biologically_prioritized")
            self.assertFalse(prioritized.reportable)

    def test_wrong_build_annotation_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "genes.tsv"
            path.write_text("chr1\t0\t100\tGENE1\n", encoding="utf-8")
            lock = IntervalResourceLock(
                resource_id="genes",
                resource_type="genes",
                source="synthetic",
                release="test",
                genome_build=GenomeBuild.GRCH38,
                sha256=sha256_file(path),
            )
            with self.assertRaisesRegex(ValueError, "does not match"):
                annotate_sv_events(
                    [_event("event", ("chr1", 10), ("chr2", 20), "Sniffles2")],
                    genome_build=GenomeBuild.GRCH37,
                    gene_resource=(path, lock),
                )

    def test_adaptive_sampling_observability_distinguishes_partial_coverage(self) -> None:
        policy = TargetCoveragePolicy(
            profile_id="test", status="technical_defaults_only", note="test"
        )
        coverage = TargetCoverageReport(
            sample_id="SYNTHETIC_B418",
            genome_build=GenomeBuild.GRCH37,
            target_bed_version="test",
            target_bed_role=TargetBedRole.ANALYSIS_ROI_UNBUFFERED,
            status=ModuleRunStatus.COMPLETED,
            policy=policy,
            summary_metrics={"region_count": 1, "interval_bases": 100},
            regions=[
                TargetCoverageRegion(
                    chromosome="chr2",
                    start=900,
                    end=1000,
                    region_id="one",
                    mean_depth=30,
                    bases_at_threshold={"1x": 100, "10x": 100, "20x": 100, "30x": 100},
                    fraction_at_threshold={"1x": 1, "10x": 1, "20x": 1, "30x": 1},
                )
            ],
            target_bed_fingerprint=FileFingerprint(size_bytes=1, sha256="0" * 64),
            tool=ToolRecord(name="mosdepth", version="0.3.14"),
        )
        event = _event("partial", ("chr2", 950), ("chr21", 2000), "Sniffles2")
        observed = apply_sv_observability(
            [event],
            assay_mode=AssayMode.ADAPTIVE_SAMPLING,
            coverage_report=coverage,
            minimum_mean_depth=10,
        )[0]
        self.assertEqual(observed.observability, SvObservability.PARTIALLY_OBSERVED)
        self.assertEqual(observed.breakpoint_mean_depths, [30.0, None])
        self.assertEqual(observed.observability_target_role, TargetBedRole.ANALYSIS_ROI_UNBUFFERED)

    def test_both_targeted_breakpoints_below_depth_are_insufficient(self) -> None:
        policy = TargetCoveragePolicy(
            profile_id="test", status="technical_defaults_only", note="test"
        )
        coverage = TargetCoverageReport(
            sample_id="SYNTHETIC_B418",
            genome_build=GenomeBuild.GRCH37,
            target_bed_version="test",
            target_bed_role=TargetBedRole.ANALYSIS_ROI_UNBUFFERED,
            status=ModuleRunStatus.COMPLETED,
            policy=policy,
            summary_metrics={"region_count": 2, "interval_bases": 400},
            regions=[
                TargetCoverageRegion(
                    chromosome=chromosome,
                    start=900 if chromosome == "chr2" else 1900,
                    end=1100 if chromosome == "chr2" else 2100,
                    region_id=chromosome,
                    mean_depth=4,
                    bases_at_threshold={"1x": 200, "10x": 0, "20x": 0, "30x": 0},
                    fraction_at_threshold={"1x": 1, "10x": 0, "20x": 0, "30x": 0},
                )
                for chromosome in ("chr2", "chr21")
            ],
            target_bed_fingerprint=FileFingerprint(size_bytes=1, sha256="0" * 64),
            tool=ToolRecord(name="mosdepth", version="0.3.14"),
        )
        observed = apply_sv_observability(
            [_event("low-depth", ("chr2", 1000), ("chr21", 2000), "Sniffles2")],
            assay_mode=AssayMode.ADAPTIVE_SAMPLING,
            coverage_report=coverage,
            minimum_mean_depth=10,
        )[0]
        self.assertEqual(observed.observability, SvObservability.INSUFFICIENT_COVERAGE)
        self.assertEqual(observed.breakpoint_mean_depths, [4.0, 4.0])

    def test_html_retains_full_table_and_exposes_review_fields(self) -> None:
        result = build_demo_result()
        with tempfile.TemporaryDirectory() as temporary:
            path = render_html(result, Path(temporary) / "report.html")
            rendered = path.read_text(encoding="utf-8")
        self.assertIn("SV review queue", rendered)
        self.assertIn("Technical appendix", rendered)
        self.assertIn("Observability", rendered)
        self.assertIn("Validation status", rendered)
        self.assertIn("Filter review queue", rendered)
        self.assertIn("BENCHMARK_REQUIRED", rendered)


if __name__ == "__main__":
    unittest.main()

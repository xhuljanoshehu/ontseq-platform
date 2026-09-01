from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.demo import build_demo_result
from ontseq_platform.models import (
    FileFingerprint,
    GenomeBuild,
    ModuleRunStatus,
    PipelineResult,
    TargetBedRole,
    ToolRecord,
)
from ontseq_platform.report import render_html
from ontseq_platform.target_coverage import (
    TargetCoveragePolicy,
    TargetCoverageRegion,
    TargetCoverageReport,
)


class ReportHardeningTests(unittest.TestCase):
    def _render(
        self,
        result: PipelineResult,
        *,
        target_coverage: TargetCoverageReport | None = None,
    ) -> str:
        with tempfile.TemporaryDirectory() as temporary:
            path = render_html(
                result,
                Path(temporary) / "report.html",
                target_coverage=target_coverage,
            )
            return path.read_text(encoding="utf-8")

    def test_execution_states_and_empty_events_are_not_negative_findings(self) -> None:
        result = build_demo_result()
        result.events = []
        modules = {item.module.value: item for item in result.modules}
        modules["sv"].status = ModuleRunStatus.NO_CALL
        modules["sv"].reason = r"No call from C:\patients\CASE-7.bam"
        modules["cnv"].status = ModuleRunStatus.FAILED
        modules["cnv"].reason = "Failure at /mnt/private/CASE-7.bam"
        modules["report"].status = ModuleRunStatus.NOT_RUN
        modules["report"].reason = "Unavailable config: configs/private/report.yaml."
        result.warnings.extend(
            [
                "Input file:///private/CASE-7.bam was withheld",
                r"UNC \\server\share\UNC_SECRET.bam was withheld",
                "input_path=SECRET_EQUALS",
                "input_path:SECRET_COLON",
                '{"input_path":"SECRET_JSON"}',
                "Technical terms CNV/SV, tumor/normal and BCR/ABL1 remain meaningful.",
                "<img src=x onerror=alert(1)>",
            ]
        )

        document = self._render(result)

        for expected in (
            "NO_CALL",
            "FAILED",
            "NOT_RUN",
            "Analysis ran without an interpretable call",
            "Module did not run; this is not a biological negative result",
            "No structural-variant event was produced",
            "[redacted path-like token]",
            "CNV/SV",
            "tumor/normal",
            "BCR/ABL1",
        ):
            self.assertIn(expected, document)
        self.assertIn("<tr class='critical'><td>cnv</td>", document)
        for secret in (
            "CASE-7",
            "SECRET_EQUALS",
            "SECRET_COLON",
            "SECRET_JSON",
            "UNC_SECRET",
            "server\\share",
            "C:\\patients",
            "/mnt/private",
            "configs/private",
            "file://",
        ):
            self.assertNotIn(secret, document)
        self.assertNotIn("<img src=x onerror=alert(1)>", document)
        self.assertIn("&lt;img src=x onerror=alert(1)&gt;", document)

    def test_missing_qc_and_coverage_values_do_not_become_zero(self) -> None:
        result = build_demo_result()
        result.qc.metrics.update(mean_coverage_x=None, n50_bp=0)
        policy = TargetCoveragePolicy(
            profile_id="synthetic-coverage",
            status="technical_defaults_only",
            note="Synthetic test policy.",
        )
        coverage = TargetCoverageReport(
            sample_id=result.manifest.sample_id,
            genome_build=GenomeBuild.GRCH38,
            target_bed_version="synthetic-v1",
            target_bed_role=TargetBedRole.ANALYSIS_ROI_UNBUFFERED,
            status=ModuleRunStatus.COMPLETED,
            policy=policy,
            summary_metrics={"region_count": 1, "interval_bases": 100},
            regions=[
                TargetCoverageRegion(
                    chromosome="chr1",
                    start=0,
                    end=100,
                    region_id="synthetic-zero-depth",
                    mean_depth=0,
                    bases_at_threshold={"1x": 0, "10x": 0, "20x": 0, "30x": 0},
                    fraction_at_threshold={"1x": 0, "10x": 0, "20x": 0, "30x": 0},
                )
            ],
            target_bed_fingerprint=FileFingerprint(size_bytes=1, sha256="0" * 64),
            tool=ToolRecord(name="mosdepth", version="0.3.14"),
        )

        document = self._render(result, target_coverage=coverage)

        self.assertIn("Mean Coverage X</span><strong>not available", document)
        self.assertIn("N50 Bp</span><strong>0", document)
        for label in ("Target-weighted mean", "Median target mean", "Least-covered target"):
            self.assertIn(f"{label}</span><strong>not available", document)
        self.assertIn("Target bases ≥20×</span><strong>not available", document)
        self.assertIn("Buffered selection mean</span><strong>not available", document)
        self.assertIn("synthetic-zero-depth", document)
        self.assertIn("0.00×", document)
        self.assertIn("0.0%", document)

    def test_event_release_and_fusion_language_preserve_the_ruo_boundary(self) -> None:
        result = build_demo_result()
        default_document = self._render(result)
        self.assertIn(
            "BENCHMARK_REQUIRED — release gate not satisfied; not a biological negative result",
            default_document,
        )
        self.assertIn(
            "Neither state asserts an expressed, in-frame or functional fusion transcript",
            default_document,
        )
        result.events = [
            event.model_copy(update={"reportable": True}) if event.event_id == "FUS-001" else event
            for event in result.events
        ]
        document = self._render(result)
        self.assertIn("RESEARCH USE ONLY", document)
        self.assertIn(
            "REPORTABLE — pipeline flag only; this RUO report is not clinically validated",
            document,
        )

    def test_tool_parameters_are_allowlisted_and_report_remains_offline(self) -> None:
        result = build_demo_result()
        qdnaseq_parameters = {
            "profile": "AML_LCWGS_GRCh38",
            "bin_sizes_kbp": [100, 500, 1000],
            "primary_bin_size_kbp": 500,
            "ace_penalty": 0.6,
            "ploidy_min": 1.0,
            "ploidy_max": 4.0,
            "ploidy_step": 0.1,
        }
        result.provenance.tools = [
            ToolRecord(
                name="Sniffles2",
                version="2.8.0",
                parameters={
                    "threads": 8,
                    "minsupport": 5,
                    "minsvlen": 50,
                    "mapq": 20,
                    "caller_pass_only": False,
                    "normalizer_pass_only": True,
                    "symbolic": True,
                    "mosaic": False,
                    "output_read_names": False,
                    "expected_version": "2.8.0",
                    "input_path": r"C:\patients\CASE-9.bam",
                },
            ),
            ToolRecord(
                name="cuteSV",
                version="2.1.3",
                parameters={
                    "threads": 8,
                    "min_support": 5,
                    "min_size": 50,
                    "max_cluster_bias_INS": 100,
                    "diff_ratio_merging_INS": 0.3,
                    "max_cluster_bias_DEL": 100,
                    "diff_ratio_merging_DEL": 0.5,
                    "expected_version": "2.1.3",
                    "normalizer_pass_only": True,
                },
            ),
            ToolRecord(
                name="mosdepth",
                version="0.3.14",
                parameters={
                    "threads": 4,
                    "no_per_base": True,
                    "thresholds": [1, 10, 20, 30],
                    "mapq": 0,
                    "exclude_flags": 1796,
                    "target_bed_role": "analysis_roi_unbuffered",
                    "expected_version": "0.3.14",
                },
            ),
            ToolRecord(name="QDNAseq", version="1.42.0", parameters=qdnaseq_parameters),
            ToolRecord(name="ACE", version="1.24.0", parameters=qdnaseq_parameters),
            ToolRecord(
                name="cramino",
                version="1.3.0",
                parameters={
                    "threads": 4,
                    "format": "json",
                    "read_length_histogram_requested": True,
                },
            ),
            ToolRecord(
                name="minimap2",
                version="2.30",
                parameters={
                    "preset": "map-ont",
                    "threads": 8,
                    "md_tag": True,
                    "soft_clip_supplementary": True,
                    "carry_fastq_tags": True,
                },
            ),
            ToolRecord(
                name="samtools",
                version="1.24",
                parameters={"checks": ["quickcheck", "view -H", "idxstats -X"]},
            ),
            ToolRecord(
                name="dorado",
                version="1.1.1",
                parameters={
                    "model": "dna_r10.4.1_e8.2_400bps_sup@v5.0.0",
                    "model_sha256": "a" * 64,
                    "device": "cuda:all",
                    "modified_bases": ["5mCG_5hmCG"],
                    "minimum_qscore": 10,
                    "emit_moves": True,
                    "adapter_verification": "unverified-runtime-adapter-v1",
                },
            ),
            ToolRecord(
                name="dorado",
                version="path-redaction-fixture",
                parameters={"model": r"C:\models\SECRET_DORADO"},
            ),
            ToolRecord(name="unknown-tool", version="1.0", parameters={"threads": 99}),
        ]

        document = self._render(result)

        for expected in (
            "Sniffles2",
            "&quot;minsupport&quot;: 5",
            "&quot;minsvlen&quot;: 50",
            "&quot;output_read_names&quot;: false",
            "cuteSV",
            "&quot;max_cluster_bias_INS&quot;: 100",
            "&quot;diff_ratio_merging_DEL&quot;: 0.5",
            "mosdepth",
            "&quot;thresholds&quot;: [1, 10, 20, 30]",
            "&quot;target_bed_role&quot;: &quot;analysis_roi_unbuffered&quot;",
            "QDNAseq",
            "ACE",
            "&quot;profile&quot;: &quot;AML_LCWGS_GRCh38&quot;",
            "&quot;bin_sizes_kbp&quot;: [100, 500, 1000]",
            "cramino",
            "&quot;format&quot;: &quot;json&quot;",
            "minimap2",
            "&quot;preset&quot;: &quot;map-ont&quot;",
            "samtools",
            "&quot;checks&quot;: [&quot;quickcheck&quot;, &quot;view -H&quot;, "
            "&quot;idxstats -X&quot;]",
            "dorado",
            "&quot;device&quot;: &quot;cuda:all&quot;",
            "unsafe or path-like parameters redacted",
            "SV review queue",
            "Technical appendix",
            "Adaptive-sampling target coverage",
        ):
            self.assertIn(expected, document)
        self.assertIn(
            "<tr><td>unknown-tool</td><td>1.0</td><td><code>{} · unsafe or path-like "
            "parameters redacted</code></td></tr>",
            document,
        )
        for secret in (
            "/secure-input/",
            "configs/panels/",
            "C:\\patients",
            "SECRET_DORADO",
            "input_path",
        ):
            self.assertNotIn(secret, document)
        self.assertNotIn("https://", document)
        self.assertNotIn("http://", document)
        self.assertNotIn("<script src=", document)
        self.assertNotIn("<link rel=", document)
        self.assertIn("self-contained HTML has no CDN or remote runtime dependency", document)


if __name__ == "__main__":
    unittest.main()

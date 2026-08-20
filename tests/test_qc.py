from __future__ import annotations

import json
import unittest
from collections.abc import Sequence

from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    QCPolicy,
    SampleManifest,
    Verdict,
)
from ontseq_platform.qc import evaluate_qc_metrics, parse_cramino_json, run_cramino_qc

CRAMINO_JSON = json.dumps(
    {
        "file_info": {
            "name": "must-not-be-copied.bam",
            "path": "/sensitive/must-not-be-copied.bam",
            "creation_time": "01/01/2026 00:00:00",
        },
        "alignment_stats": {
            "num_alignments": 120,
            "percent_from_total": 75.5,
            "num_reads": 100,
        },
        "read_stats": {
            "yield_gb": 1.2,
            "mean_coverage": 3.1,
            "yield_gb_long": 1.0,
            "n50": 20000,
            "n75": 10000,
            "median_length": 8000.0,
            "mean_length": 12000.0,
        },
        "identity_stats": {
            "median_identity": 96.0,
            "mean_identity": 95.0,
            "modal_identity": 97.0,
            "is_estimated": False,
        },
    }
)


class FakeCraminoRunner:
    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(argv)
        if normalized[1:] == ("--version",):
            return CommandResult(normalized, 0, "cramino 1.3.0\n", "")
        return CommandResult(normalized, 0, CRAMINO_JSON, "")


class FailingCraminoRunner:
    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(argv)
        if normalized[1:] == ("--version",):
            return CommandResult(normalized, 0, "cramino 1.3.0\n", "")
        return CommandResult(
            normalized,
            101,
            "",
            "thread 'main' panicked at missing NM/de alignment tag\nsynthetic diagnostic",
        )


def _policy(mean_coverage: float | None = None) -> QCPolicy:
    return QCPolicy(
        status="technical_defaults_only",
        hard_failures=[],
        numeric_gates={"mean_coverage_x": mean_coverage},
        note="Synthetic unvalidated policy",
    )


def _manifest() -> SampleManifest:
    return SampleManifest(
        sample_id="SYNTHETIC_001",
        run_id="SYNTHETIC_RUN_001",
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path="/secure/SYNTHETIC_001.bam",
            index_path="/secure/SYNTHETIC_001.bam.bai",
        ),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id="SYNTHETIC_REF",
        ),
        analysis=AnalysisSpec(profile="lcwgs", modules=[]),
    )


class CraminoQCTests(unittest.TestCase):
    def test_parser_normalizes_metrics_without_copying_source_path(self) -> None:
        metrics = parse_cramino_json(CRAMINO_JSON)
        self.assertEqual(metrics["number_of_reads"], 100)
        self.assertEqual(metrics["mean_coverage_x"], 3.1)
        self.assertEqual(metrics["median_identity_percent"], 96.0)
        self.assertNotIn("path", metrics)
        self.assertNotIn("name", metrics)

    def test_unvalidated_empty_policy_is_warning_not_pass(self) -> None:
        qc = evaluate_qc_metrics(parse_cramino_json(CRAMINO_JSON), _policy())
        self.assertEqual(qc.verdict, Verdict.WARN)

    def test_configured_minimum_can_fail(self) -> None:
        qc = evaluate_qc_metrics(parse_cramino_json(CRAMINO_JSON), _policy(4.0))
        self.assertEqual(qc.verdict, Verdict.FAIL)
        self.assertEqual(qc.failed_gates, ["mean_coverage_x"])

    def test_runner_records_version_and_normalized_output(self) -> None:
        report = run_cramino_qc(
            _manifest(),
            _policy(3.0),
            runner=FakeCraminoRunner(),
            threads=2,
        )
        self.assertEqual(report.tool.version, "1.3.0")
        self.assertEqual(report.qc.verdict, Verdict.WARN)
        self.assertNotIn("/sensitive/", report.model_dump_json())

    def test_runner_includes_bounded_stderr_on_tool_failure(self) -> None:
        with self.assertRaisesRegex(
            ValueError,
            r"Cramino failed with exit code 101: .*missing NM/de alignment tag",
        ):
            run_cramino_qc(
                _manifest(),
                _policy(),
                runner=FailingCraminoRunner(),
                threads=2,
            )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import json
import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

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
from ontseq_platform.qc import (
    evaluate_qc_metrics,
    parse_cramino_json,
    read_length_histogram,
    run_cramino_qc,
    write_read_length_histogram,
)

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
    def __init__(self) -> None:
        self.calls: list[tuple[str, ...]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        normalized = tuple(argv)
        self.calls.append(normalized)
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
    def test_histogram_count_output_is_explicit_so_stdout_stays_json(self) -> None:
        runner = FakeCraminoRunner()
        with tempfile.TemporaryDirectory() as raw:
            output = Path(raw) / "read_length_histogram.tsv"
            run_cramino_qc(
                _manifest(),
                _policy(),
                runner=runner,
                histogram_output=output,
            )

            argv = runner.calls[-1]
            option = argv.index("--hist-count")
            raw_histogram = Path(argv[option + 1])
            self.assertEqual(argv[option + 2 : option + 4], ("--format", "json"))
            self.assertEqual(raw_histogram.parent, output.parent)
            self.assertFalse(raw_histogram.exists())

    def test_optional_numeric_read_length_histogram_is_written_as_tsv(self) -> None:
        payload = json.loads(CRAMINO_JSON)
        payload["histograms"] = {
            "read_length": [
                {"start": 0, "end": 2000, "count": 3, "bases": 2500},
                {"start": 2000, "end": None, "count": 1, "bases": 5000},
            ]
        }
        bins = read_length_histogram(json.dumps(payload))
        self.assertEqual(bins[-1], (2000, None, 1, 5000))
        with tempfile.TemporaryDirectory() as raw:
            output = Path(raw) / "hist.tsv"
            self.assertEqual(write_read_length_histogram(bins, output), output)
            self.assertEqual(
                output.read_text(encoding="utf-8").splitlines()[0],
                "start_bp\tend_bp\tread_count\tbase_count",
            )

    def test_empty_histogram_removes_stale_sidecar(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            output = Path(raw) / "hist.tsv"
            output.write_text("stale histogram\n", encoding="utf-8")

            self.assertIsNone(write_read_length_histogram([], output))
            self.assertFalse(output.exists())

    def test_cramino_130_nested_histogram_object_is_supported(self) -> None:
        payload = json.loads(CRAMINO_JSON)
        payload["histograms"] = {
            "read_length": {
                "step": 2000,
                "max_value": 4000,
                "bins": [
                    {"start": 0, "end": 2000, "count": 3, "bases": 2500},
                    {"start": 4000, "count": 1, "bases": 5000},
                ],
            }
        }

        self.assertEqual(
            read_length_histogram(json.dumps(payload)),
            [(0, 2000, 3, 2500), (4000, None, 1, 5000)],
        )

    def test_parser_normalizes_metrics_without_copying_source_path(self) -> None:
        metrics = parse_cramino_json(CRAMINO_JSON)
        self.assertEqual(metrics["number_of_reads"], 100)
        self.assertEqual(metrics["mean_coverage_x"], 3.1)
        self.assertEqual(metrics["median_identity_percent"], 96.0)
        self.assertNotIn("path", metrics)
        self.assertNotIn("name", metrics)

    def test_missing_numeric_gates_warn_independent_of_research_status(self) -> None:
        qc = evaluate_qc_metrics(parse_cramino_json(CRAMINO_JSON), _policy())
        self.assertEqual(qc.verdict, Verdict.WARN)
        self.assertIn(
            "No validated numeric QC gates are configured; metrics are descriptive only.",
            qc.warnings,
        )
        self.assertNotIn("Synthetic unvalidated policy", qc.warnings)

    def test_configured_minimum_can_fail(self) -> None:
        qc = evaluate_qc_metrics(parse_cramino_json(CRAMINO_JSON), _policy(4.0))
        self.assertEqual(qc.verdict, Verdict.FAIL)
        self.assertEqual(qc.failed_gates, ["mean_coverage_x"])

    def test_runner_research_status_does_not_change_qc_verdict(self) -> None:
        report = run_cramino_qc(
            _manifest(),
            _policy(3.0),
            runner=FakeCraminoRunner(),
            threads=2,
        )
        self.assertEqual(report.tool.version, "1.3.0")
        self.assertEqual(report.qc.verdict, Verdict.PASS)
        self.assertNotIn("Synthetic unvalidated policy", report.qc.warnings)
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

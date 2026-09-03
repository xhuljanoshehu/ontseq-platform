from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.dilution import (
    DilutionDetection,
    DilutionPolicy,
    LodPolicy,
    classify_detection,
    evaluate_lod,
    execute_dilution_series,
    plan_dilution_series,
)
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import (
    BenchmarkKind,
    BenchmarkMetrics,
    BenchmarkReport,
    BenchmarkThresholds,
    GenomeBuild,
)


def _policy(**overrides: object) -> DilutionPolicy:
    defaults: dict[str, object] = {
        "profile_id": "synthetic",
        "status": "technical_defaults_only",
        "tumor_fractions": [1.0, 0.5, 0.1],
        "note": "Synthetic technical dilution design",
    }
    defaults.update(overrides)
    return DilutionPolicy.model_validate(defaults)


def _plan(**overrides: object):
    return plan_dilution_series(
        _policy(**overrides.pop("policy", {})),  # type: ignore[arg-type]
        series_id=str(overrides.pop("series_id", "SERIES_001")),
        tumor_sample_id="TUMOR_001",
        normal_sample_id="NORMAL_001",
        genome_build=GenomeBuild.GRCH38,
        tumor_read_count=int(overrides.pop("tumor_read_count", 1000)),
        normal_read_count=int(overrides.pop("normal_read_count", 1000)),
    )


def _report(
    fraction: float,
    replicate: int,
    *,
    true_positive: int,
    false_negative: int,
    false_positive: int = 0,
    case_suffix: str = "",
) -> BenchmarkReport:
    truth = true_positive + false_negative
    recall = true_positive / truth if truth else None
    called = true_positive + false_positive
    precision = true_positive / called if called else None
    return BenchmarkReport(
        case_id=f"CASE_{round(fraction * 1000):04d}_R{replicate}{case_suffix}",
        kind=BenchmarkKind.CNV,
        genome_build=GenomeBuild.GRCH38,
        thresholds=BenchmarkThresholds(),
        strata={
            "dilution_series_id": "SERIES_001",
            "tumor_fraction": fraction,
            "replicate": replicate,
        },
        metrics=BenchmarkMetrics(
            true_positive=true_positive,
            false_positive=false_positive,
            false_negative=false_negative,
            precision=precision,
            recall=recall,
        ),
        matches=[],
        unmatched_truth_event_ids=[],
        unmatched_query_event_ids=[],
    )


def _lod_policy(**overrides: object) -> LodPolicy:
    defaults: dict[str, object] = {
        "profile_id": "synthetic",
        "status": "technical_defaults_only",
        "minimum_detection_rate": 0.5,
        "note": "Synthetic technical detection criterion",
    }
    defaults.update(overrides)
    return LodPolicy.model_validate(defaults)


class PlanTests(unittest.TestCase):
    def test_every_level_gets_the_same_depth(self) -> None:
        """A titration varies tumour fraction; varying depth as well confounds the result."""
        plan = _plan()
        totals = {level.total_reads_target for level in plan.levels}
        self.assertEqual(len(totals), 1)
        self.assertEqual(plan.tumor_read_count, 1000)

    def test_the_budget_is_the_largest_every_level_can_fund(self) -> None:
        # The 1.0 level caps the budget at the tumour's 1000 reads; the 0.1 level would
        # need 1000/0.9 normal reads for a larger budget than the normal BAM holds.
        plan = _plan()
        self.assertEqual(plan.levels[0].total_reads_target, 1000)
        self.assertEqual(plan.levels[0].tumor_reads_target, 1000)
        self.assertEqual(plan.levels[0].normal_reads_target, 0)

    def test_a_negative_control_is_planned_and_marked(self) -> None:
        plan = _plan()
        control = [level for level in plan.levels if level.is_negative_control]
        self.assertEqual(len(control), 1)
        self.assertEqual(control[0].nominal_tumor_fraction, 0.0)
        self.assertEqual(control[0].tumor_reads_target, 0)
        self.assertIsNone(control[0].tumor_subsample_argument)

    def test_a_source_taken_whole_carries_no_subsample_argument(self) -> None:
        """samtools cannot express "all reads" as a fraction, so it is not asked to."""
        plan = _plan()
        undiluted = plan.levels[0]
        self.assertEqual(undiluted.tumor_subsample_fraction, 1.0)
        self.assertIsNone(undiluted.tumor_subsample_argument)

    def test_subsample_arguments_carry_seed_and_fraction(self) -> None:
        plan = _plan()
        half = next(level for level in plan.levels if level.nominal_tumor_fraction == 0.5)
        self.assertEqual(half.tumor_subsample_fraction, 0.5)
        self.assertIsNotNone(half.tumor_subsample_argument)
        seed, _, digits = (half.tumor_subsample_argument or "").partition(".")
        self.assertTrue(seed.isdigit())
        self.assertEqual(digits, "500000")

    def test_no_two_subsamples_of_one_source_share_a_seed(self) -> None:
        plan = _plan(policy={"replicates": 3, "tumor_fractions": [0.5, 0.25]})
        arguments = [
            level.tumor_subsample_argument
            for level in plan.levels
            if level.tumor_subsample_argument is not None
        ]
        self.assertEqual(len(arguments), len(set(arguments)))

    def test_planning_is_deterministic(self) -> None:
        self.assertEqual(_plan().model_dump(), _plan().model_dump())

    def test_a_declared_budget_beyond_the_sources_is_refused(self) -> None:
        with self.assertRaises(ValueError) as raised:
            _plan(policy={"total_read_target": 5000})
        self.assertIn("exceeds what the sources can fund", str(raised.exception))

    def test_a_derived_budget_says_so(self) -> None:
        plan = _plan()
        self.assertTrue(
            any("derived from the sources" in item for item in plan.warnings), plan.warnings
        )

    def test_too_few_replicates_are_flagged_rather_than_refused(self) -> None:
        plan = _plan()
        self.assertTrue(
            any("cannot characterise a detection rate" in item for item in plan.warnings),
            plan.warnings,
        )

    def test_fractions_must_be_unique_and_ordered(self) -> None:
        with self.assertRaises(ValueError):
            _policy(tumor_fractions=[0.1, 0.5])
        with self.assertRaises(ValueError):
            _policy(tumor_fractions=[0.5, 0.5])

    def test_a_tumour_and_normal_of_the_same_sample_are_refused(self) -> None:
        with self.assertRaises(ValueError):
            plan_dilution_series(
                _policy(),
                series_id="SERIES_001",
                tumor_sample_id="SAME",
                normal_sample_id="SAME",
                genome_build=GenomeBuild.GRCH38,
                tumor_read_count=1000,
                normal_read_count=1000,
            )


class DetectionTests(unittest.TestCase):
    def test_an_undefined_recall_is_a_no_call_not_a_failure(self) -> None:
        report = _report(0.1, 1, true_positive=0, false_negative=0)
        self.assertEqual(classify_detection(report, _lod_policy()), DilutionDetection.NO_CALL)

    def test_partial_recall_fails_the_strict_default(self) -> None:
        report = _report(0.1, 1, true_positive=1, false_negative=1)
        self.assertEqual(classify_detection(report, _lod_policy()), DilutionDetection.NOT_DETECTED)
        self.assertEqual(
            classify_detection(report, _lod_policy(minimum_recall=0.5)),
            DilutionDetection.DETECTED,
        )


class LodTests(unittest.TestCase):
    def _evaluate(self, reports, **policy_overrides: object):
        return evaluate_lod(reports, _lod_policy(**policy_overrides), series_id="SERIES_001")

    def test_the_limit_is_the_lowest_level_that_still_detects(self) -> None:
        reports = [
            _report(0.5, 1, true_positive=2, false_negative=0),
            _report(0.2, 1, true_positive=2, false_negative=0),
            _report(0.05, 1, true_positive=0, false_negative=2),
        ]
        result = self._evaluate(reports)
        self.assertEqual(result.detection_limit_fraction, 0.2)
        self.assertTrue(result.bracketed)

    def test_a_level_that_passes_below_a_failing_one_is_not_the_limit(self) -> None:
        """Monotonic reading: one lucky low level does not become the reported limit."""
        reports = [
            _report(0.5, 1, true_positive=2, false_negative=0),
            _report(0.2, 1, true_positive=0, false_negative=2),
            _report(0.05, 1, true_positive=2, false_negative=0),
        ]
        result = self._evaluate(reports)
        self.assertEqual(result.detection_limit_fraction, 0.5)
        self.assertEqual(
            self._evaluate(reports, require_monotonic=False).detection_limit_fraction, 0.05
        )

    def test_an_unbracketed_limit_says_so(self) -> None:
        reports = [
            _report(0.5, 1, true_positive=2, false_negative=0),
            _report(0.05, 1, true_positive=2, false_negative=0),
        ]
        result = self._evaluate(reports)
        self.assertEqual(result.detection_limit_fraction, 0.05)
        self.assertFalse(result.bracketed)
        self.assertTrue(
            any("bounded from above and not located" in item for item in result.warnings),
            result.warnings,
        )

    def test_no_passing_level_yields_no_limit_and_no_impossibility_claim(self) -> None:
        reports = [_report(0.5, 1, true_positive=0, false_negative=2)]
        result = self._evaluate(reports)
        self.assertIsNone(result.detection_limit_fraction)
        self.assertFalse(result.bracketed)
        self.assertTrue(
            any(
                "did not establish that detection is impossible" in item for item in result.warnings
            ),
            result.warnings,
        )

    def test_no_call_replicates_leave_the_detection_rate_rather_than_depress_it(self) -> None:
        reports = [
            _report(0.5, 1, true_positive=2, false_negative=0),
            _report(0.5, 2, true_positive=0, false_negative=0, case_suffix="_B"),
        ]
        level = self._evaluate(reports).levels[0]
        self.assertEqual(level.replicates_total, 2)
        self.assertEqual(level.replicates_evaluated, 1)
        self.assertEqual(level.replicates_no_call, 1)
        self.assertEqual(level.detection_rate, 1.0)

    def test_a_level_with_nothing_evaluable_cannot_meet_the_criterion(self) -> None:
        reports = [_report(0.5, 1, true_positive=0, false_negative=0)]
        result = self._evaluate(reports)
        self.assertIsNone(result.levels[0].detection_rate)
        self.assertFalse(result.levels[0].meets_criterion)
        self.assertIsNone(result.detection_limit_fraction)

    def test_a_repeated_level_is_not_an_extra_replicate(self) -> None:
        reports = [
            _report(0.5, 1, true_positive=2, false_negative=0),
            _report(0.5, 1, true_positive=2, false_negative=0, case_suffix="_DUP"),
        ]
        with self.assertRaises(ValueError) as raised:
            self._evaluate(reports)
        self.assertIn("not an extra replicate", str(raised.exception))

    def test_a_report_from_another_series_is_refused(self) -> None:
        report = _report(0.5, 1, true_positive=2, false_negative=0)
        report = report.model_copy(
            update={"strata": {**report.strata, "dilution_series_id": "SERIES_OTHER"}}
        )
        with self.assertRaises(ValueError):
            self._evaluate([report])

    def test_a_report_without_a_fraction_cannot_be_placed(self) -> None:
        report = _report(0.5, 1, true_positive=2, false_negative=0)
        report = report.model_copy(update={"strata": {}})
        with self.assertRaises(ValueError) as raised:
            self._evaluate([report])
        self.assertIn("cannot be placed in a dilution series", str(raised.exception))

    def test_a_minimum_replicate_count_gates_the_limit(self) -> None:
        reports = [_report(0.5, 1, true_positive=2, false_negative=0)]
        self.assertIsNone(self._evaluate(reports, minimum_replicates=3).detection_limit_fraction)


class _FakeSamtools:
    """A samtools stand-in over line-counted text files.

    Subsampling keeps a deterministic prefix of the "reads" rather than a random draw,
    which is what makes the drift assertions below meaningful: any deviation the executor
    reports comes from the arithmetic under test, not from a random number generator.
    """

    def __init__(self, *, version: str = "1.24", tumor_skew: float = 1.0) -> None:
        self.version = version
        self.tumor_skew = tumor_skew

    def run(self, argv, *, timeout_seconds: int = 300) -> CommandResult:
        argv = [str(item) for item in argv]
        if argv[1:2] == ["--version"]:
            return CommandResult(tuple(argv), 0, f"samtools {self.version}", "")
        if "-c" in argv:
            path = Path(argv[-1])
            return CommandResult(tuple(argv), 0, str(_reads(path)), "")
        if "merge" in argv:
            destination = Path(argv[argv.index("merge") + 3])
            total = sum(_reads(Path(item)) for item in argv[argv.index("merge") + 4 :])
            _write(destination, total)
            return CommandResult(tuple(argv), 0, "", "")
        if "index" in argv:
            Path(str(argv[-1]) + ".bai").write_text("", encoding="utf-8")
            return CommandResult(tuple(argv), 0, "", "")
        raise AssertionError(f"unexpected command: {argv}")

    def run_to_file(self, argv, output_path: Path, *, timeout_seconds: int = 300):
        argv = [str(item) for item in argv]
        source = Path(argv[-1])
        available = _reads(source)
        if "-s" in argv:
            fraction = float("0." + argv[argv.index("-s") + 1].split(".")[1])
            if "tumor" in source.name:
                fraction *= self.tumor_skew
            kept = round(available * fraction)
        else:
            kept = available
        _write(output_path, kept)
        return CommandResult(tuple(argv), 0, "", "")


def _write(path: Path, reads: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("read\n" * reads, encoding="utf-8")


def _reads(path: Path) -> int:
    return len(path.read_text(encoding="utf-8").splitlines()) if path.is_file() else 0


class ExecutionTests(unittest.TestCase):
    def _sources(self, root: Path) -> tuple[Path, Path]:
        tumor = root / "tumor.bam"
        normal = root / "normal.bam"
        _write(tumor, 1000)
        _write(normal, 1000)
        return tumor, normal

    def test_each_level_records_what_it_actually_contains(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor, normal = self._sources(root)
            report = execute_dilution_series(
                _plan(policy={"tumor_fractions": [1.0, 0.5, 0.1]}),
                tumor_bam=tumor,
                normal_bam=normal,
                output_dir=root / "levels",
                runner=_FakeSamtools(),
            )
        by_fraction = {item.nominal_tumor_fraction: item for item in report.levels}
        self.assertEqual(by_fraction[0.5].tumor_reads_observed, 500)
        self.assertEqual(by_fraction[0.5].normal_reads_observed, 500)
        self.assertEqual(by_fraction[0.5].observed_tumor_fraction, 0.5)
        self.assertEqual(by_fraction[0.0].tumor_reads_observed, 0)
        self.assertEqual(by_fraction[0.0].observed_tumor_fraction, 0.0)
        self.assertEqual(by_fraction[1.0].normal_reads_observed, 0)
        self.assertTrue(all(item.mixed_bam_fingerprint.sha256 for item in report.levels))

    def test_a_level_that_drifts_past_the_tolerance_fails_the_series(self) -> None:
        """A level labelled with a fraction it does not contain poisons every later number."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor, normal = self._sources(root)
            with self.assertRaises(ValueError) as raised:
                execute_dilution_series(
                    _plan(policy={"tumor_fractions": [0.5]}),
                    tumor_bam=tumor,
                    normal_bam=normal,
                    output_dir=root / "levels",
                    runner=_FakeSamtools(tumor_skew=0.5),
                )
        self.assertIn("exceeds the policy tolerance", str(raised.exception))

    def test_a_samtools_outside_the_lock_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor, normal = self._sources(root)
            with self.assertRaises(ValueError):
                execute_dilution_series(
                    _plan(),
                    tumor_bam=tumor,
                    normal_bam=normal,
                    output_dir=root / "levels",
                    runner=_FakeSamtools(version="1.21"),
                )

    def test_an_existing_level_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor, normal = self._sources(root)
            plan = _plan(policy={"tumor_fractions": [0.5]})
            _write(root / "levels" / f"{plan.levels[0].level_id}.bam", 1)
            with self.assertRaises(ValueError) as raised:
                execute_dilution_series(
                    plan,
                    tumor_bam=tumor,
                    normal_bam=normal,
                    output_dir=root / "levels",
                    runner=_FakeSamtools(),
                )
        self.assertIn("Refusing to overwrite", str(raised.exception))


if __name__ == "__main__":
    unittest.main()

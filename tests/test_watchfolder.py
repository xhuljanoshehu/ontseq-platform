from __future__ import annotations

import json
import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

from ontseq_platform.models import InputKind
from ontseq_platform.pipeline.watch import Outcome, sample_id_from_directory
from ontseq_platform.watchfolder import (
    DropRejected,
    WatchConfigurationError,
    WatchSettings,
    find_input,
    resolve,
    sweep,
)

TEMPLATE = {
    "schema_version": "0.1.0",
    "run_id": "WATCH_RUN_001",
    "assay": {"mode": "lcwgs", "genome_build": "GRCh38", "reference_id": "REF_V1"},
    "analysis": {"profile": "lcwgs", "modules": ["qc", "sv", "report"]},
}

REFERENCE_LOCK = {
    "reference_id": "REF_V1",
    "genome_build": "GRCh38",
    "contigs": [{"name": "chr1", "length": 1000}],
    "source_fai_sha256": "a" * 64,
}

QC_POLICY = {"status": "technical_defaults_only", "note": "test"}


class WatchfolderCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.root = Path(self._temporary.name)
        self.watch = self.root / "drop"
        self.watch.mkdir()
        self.config = self.root / "config"
        self.config.mkdir()
        self._write("manifest.json", TEMPLATE)
        self._write("reference.lock.json", REFERENCE_LOCK)
        self._write("qc.json", QC_POLICY)

    def _write(self, name: str, payload: dict) -> Path:
        path = self.config / name
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def _settings(self, **overrides: object) -> WatchSettings:
        values: dict[str, object] = {
            "watch_dir": self.watch,
            "output_dir": self.root / "out",
            "manifest_template": self.config / "manifest.json",
            "reference_lock": self.config / "reference.lock.json",
            "qc_policy": self.config / "qc.json",
            "input_kind": InputKind.ALIGNED_BAM,
            "sniffles_policy": None,
            "alignment_policy": None,
            "quiet_seconds": 0.0,
        }
        values.update(overrides)
        return WatchSettings(**values)  # type: ignore[arg-type]

    def _drop(self, name: str, files: dict[str, bytes]) -> Path:
        directory = self.watch / name
        directory.mkdir(parents=True, exist_ok=True)
        stale = time.time() - 10_000
        for relative, payload in files.items():
            target = directory / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(payload)
            # Backdated so quiescence is satisfied against the real clock; pinning `now`
            # instead would date every freshly written fixture into the future.
            os.utime(target, (stale, stale))
        return directory

    def _aligned_drop(self, name: str = "SAMPLE_A") -> Path:
        return self._drop(name, {"reads.bam": b"BAM", "reads.bam.bai": b"BAI"})


class FindInputTests(WatchfolderCase):
    def test_an_aligned_bam_and_its_index_are_found(self) -> None:
        directory = self._aligned_drop()
        bam, index = find_input(directory, InputKind.ALIGNED_BAM)
        self.assertEqual(bam.name, "reads.bam")
        self.assertEqual(index.name, "reads.bam.bai")

    def test_a_sibling_style_index_is_accepted(self) -> None:
        directory = self._drop("SAMPLE_A", {"reads.bam": b"BAM", "reads.bai": b"BAI"})
        _, index = find_input(directory, InputKind.ALIGNED_BAM)
        self.assertEqual(index.name, "reads.bai")

    def test_a_missing_index_is_refused_rather_than_reinterpreted(self) -> None:
        """The dangerous alternative is inferring "no index means unaligned" and realigning."""
        directory = self._drop("SAMPLE_A", {"reads.bam": b"BAM"})
        with self.assertRaises(DropRejected) as raised:
            find_input(directory, InputKind.ALIGNED_BAM)
        self.assertIn("not treated as evidence", str(raised.exception))

    def test_two_bams_are_refused_rather_than_chosen_between(self) -> None:
        directory = self._drop("SAMPLE_A", {"a.bam": b"BAM", "b.bam": b"BAM"})
        with self.assertRaises(DropRejected) as raised:
            find_input(directory, InputKind.UNALIGNED_BAM)
        self.assertIn("exactly one", str(raised.exception))

    def test_no_matching_file_is_refused(self) -> None:
        directory = self._drop("SAMPLE_A", {"notes.txt": b"hello"})
        with self.assertRaises(DropRejected):
            find_input(directory, InputKind.UNALIGNED_BAM)

    def test_an_unaligned_bam_needs_no_index(self) -> None:
        directory = self._drop("SAMPLE_A", {"reads.bam": b"BAM"})
        bam, index = find_input(directory, InputKind.UNALIGNED_BAM)
        self.assertEqual(bam.name, "reads.bam")
        self.assertIsNone(index)

    def test_pod5_resolves_to_the_containing_directory(self) -> None:
        directory = self._drop("SAMPLE_A", {"pod5/a.pod5": b"signal", "pod5/b.pod5": b"signal"})
        found, index = find_input(directory, InputKind.POD5)
        self.assertEqual(found, directory / "pod5")
        self.assertIsNone(index)

    def test_pod5_split_across_directories_is_refused(self) -> None:
        directory = self._drop("SAMPLE_A", {"one/a.pod5": b"s", "two/b.pod5": b"s"})
        with self.assertRaises(DropRejected) as raised:
            find_input(directory, InputKind.POD5)
        self.assertIn("exactly one", str(raised.exception))

    def test_hidden_files_are_not_candidates(self) -> None:
        """Partial copies from rsync and friends appear as dotfiles."""
        directory = self._drop("SAMPLE_A", {"reads.bam": b"BAM", ".partial.bam": b"BAM"})
        bam, _ = find_input(directory, InputKind.UNALIGNED_BAM)
        self.assertEqual(bam.name, "reads.bam")


class ResolveTests(WatchfolderCase):
    def test_a_valid_configuration_resolves(self) -> None:
        resolved = resolve(self._settings())
        self.assertEqual(resolved.reference_lock.reference_id, "REF_V1")
        self.assertIsNone(resolved.sniffles_policy)

    def test_a_missing_watch_folder_is_a_configuration_error(self) -> None:
        with self.assertRaises(WatchConfigurationError):
            resolve(self._settings(watch_dir=self.root / "nowhere"))

    def test_a_broken_template_fails_the_watcher_not_the_samples(self) -> None:
        """One mistake in the template is one mistake, not a failure per sample."""
        self._write("manifest.json", {"analysis": {"profile": "p", "modules": []}})
        with self.assertRaises(WatchConfigurationError) as raised:
            resolve(self._settings())
        self.assertIn("does not produce a valid manifest", str(raised.exception))

    def test_a_missing_policy_file_is_a_configuration_error(self) -> None:
        with self.assertRaises(WatchConfigurationError):
            resolve(self._settings(qc_policy=self.config / "absent.json"))

    def test_an_aligning_kind_without_a_reference_fasta_is_refused_up_front(self) -> None:
        with self.assertRaises(WatchConfigurationError) as raised:
            resolve(self._settings(input_kind=InputKind.UNALIGNED_BAM))
        self.assertIn("--reference-fasta", str(raised.exception))

    def test_the_probe_identifier_never_reaches_a_real_manifest(self) -> None:
        """It exists only to validate the template, so it must not leak into settings."""
        resolved = resolve(self._settings())
        self.assertNotIn("TEMPLATE_PROBE", json.dumps(resolved.manifest_template))


class _StubReport:
    def __init__(self, passed: bool) -> None:
        self.passed = passed
        self.verdict_reason = "PASS - stub" if passed else "FAIL - stub"


class SweepTests(WatchfolderCase):
    def _sweep(self, **overrides: object):
        return sweep(self._settings(**overrides))

    def test_a_ready_sample_is_run_and_recorded(self) -> None:
        self._aligned_drop("SAMPLE_A")
        with mock.patch(
            "ontseq_platform.watchfolder.run_pipeline",
            return_value=(_StubReport(True), None),
        ) as runner:
            result = self._sweep()
        self.assertEqual(len(result.attempted), 1)
        self.assertIs(result.attempted[0].outcome, Outcome.COMPLETED)
        self.assertEqual(runner.call_count, 1)

    def test_the_manifest_handed_to_the_pipeline_carries_the_derived_identity(self) -> None:
        self._aligned_drop("AML_0031")
        with mock.patch(
            "ontseq_platform.watchfolder.run_pipeline",
            return_value=(_StubReport(True), None),
        ) as runner:
            self._sweep()
        manifest = runner.call_args.args[0].manifest
        self.assertEqual(manifest.sample_id, "AML_0031")
        self.assertEqual(manifest.assay.reference_id, "REF_V1")
        self.assertTrue(manifest.input.path.endswith("reads.bam"))

    def test_a_failing_run_is_recorded_as_failed(self) -> None:
        self._aligned_drop("SAMPLE_A")
        with mock.patch(
            "ontseq_platform.watchfolder.run_pipeline",
            return_value=(_StubReport(False), None),
        ):
            result = self._sweep()
        self.assertIs(result.attempted[0].outcome, Outcome.FAILED)
        self.assertEqual(len(result.failures), 1)

    def test_a_completed_sample_is_not_repeated_on_the_next_sweep(self) -> None:
        self._aligned_drop("SAMPLE_A")
        with mock.patch(
            "ontseq_platform.watchfolder.run_pipeline",
            return_value=(_StubReport(True), None),
        ) as runner:
            self._sweep()
            second = self._sweep()
        self.assertEqual(runner.call_count, 1)
        self.assertEqual(len(second.attempted), 0)
        self.assertIn("already completed", second.skipped[0][1])

    def test_a_failed_sample_is_not_repeated_without_retry_failed(self) -> None:
        self._aligned_drop("SAMPLE_A")
        with mock.patch(
            "ontseq_platform.watchfolder.run_pipeline",
            return_value=(_StubReport(False), None),
        ) as runner:
            self._sweep()
            self._sweep()
        self.assertEqual(runner.call_count, 1)

    def test_a_failed_sample_is_repeated_with_retry_failed(self) -> None:
        self._aligned_drop("SAMPLE_A")
        with mock.patch(
            "ontseq_platform.watchfolder.run_pipeline",
            return_value=(_StubReport(False), None),
        ) as runner:
            self._sweep()
            self._sweep(retry_failed=True)
        self.assertEqual(runner.call_count, 2)

    def test_a_directory_with_an_unusable_name_is_rejected_without_running(self) -> None:
        self._drop("not a sample id", {"reads.bam": b"BAM", "reads.bam.bai": b"BAI"})
        with mock.patch("ontseq_platform.watchfolder.run_pipeline") as runner:
            result = self._sweep()
        self.assertIs(result.attempted[0].outcome, Outcome.REJECTED)
        self.assertIn("will not invent", result.attempted[0].detail)
        runner.assert_not_called()

    def test_a_directory_missing_its_index_is_rejected_without_running(self) -> None:
        self._drop("SAMPLE_A", {"reads.bam": b"BAM"})
        with mock.patch("ontseq_platform.watchfolder.run_pipeline") as runner:
            result = self._sweep()
        self.assertIs(result.attempted[0].outcome, Outcome.REJECTED)
        runner.assert_not_called()

    def test_several_samples_are_processed_in_name_order(self) -> None:
        for name in ("SAMPLE_C", "SAMPLE_A", "SAMPLE_B"):
            self._aligned_drop(name)
        with mock.patch(
            "ontseq_platform.watchfolder.run_pipeline",
            return_value=(_StubReport(True), None),
        ):
            result = self._sweep()
        self.assertEqual(
            [item.name for item in result.attempted], ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"]
        )

    def test_one_rejected_sample_does_not_stop_the_others(self) -> None:
        self._drop("SAMPLE_A", {"notes.txt": b"x"})
        self._aligned_drop("SAMPLE_B")
        with mock.patch(
            "ontseq_platform.watchfolder.run_pipeline",
            return_value=(_StubReport(True), None),
        ):
            result = self._sweep()
        outcomes = {item.name: item.outcome for item in result.attempted}
        self.assertIs(outcomes["SAMPLE_A"], Outcome.REJECTED)
        self.assertIs(outcomes["SAMPLE_B"], Outcome.COMPLETED)

    def test_a_not_ready_directory_is_skipped_with_its_reason(self) -> None:
        self._aligned_drop("SAMPLE_A")
        with mock.patch("ontseq_platform.watchfolder.run_pipeline") as runner:
            result = self._sweep(ready_marker=".ready")
        self.assertEqual(len(result.attempted), 0)
        self.assertIn(".ready", result.skipped[0][1])
        runner.assert_not_called()

    def test_the_ledger_is_written_beside_the_output_not_into_the_drop_folder(self) -> None:
        """The drop folder may be owned by the instrument or mounted read-only."""
        self._aligned_drop("SAMPLE_A")
        with mock.patch(
            "ontseq_platform.watchfolder.run_pipeline",
            return_value=(_StubReport(True), None),
        ):
            self._sweep()
        self.assertTrue((self.root / "out" / ".ontseq-watch.json").is_file())
        self.assertEqual(list(self.watch.glob(".ontseq*")), [])


if __name__ == "__main__":
    unittest.main()


class GridionPod5LayoutTests(WatchfolderCase):
    """MinKNOW splits a GridION run into pod5_pass and pod5_fail by qscore.

    Which of the two enters the analysis is a scientific decision, not a filesystem
    detail: failed reads change the depth distribution that depth-based copy-number
    methods assume. So it is declared, and an undeclared split is refused.
    """

    def _gridion_run(self, *, fail_reads: bool = True) -> Path:
        sample = self.watch / "AML_00123"
        run = sample / "20260818_1030_X1_FAV12345_abcdef12"
        (run / "pod5_pass").mkdir(parents=True)
        (run / "pod5_pass" / "reads_0.pod5").write_bytes(b"POD5")
        if fail_reads:
            (run / "pod5_fail").mkdir(parents=True)
            (run / "pod5_fail" / "reads_0.pod5").write_bytes(b"POD5")
        return sample

    def test_an_undeclared_split_is_refused_with_the_directories_named(self) -> None:
        sample = self._gridion_run()
        with self.assertRaises(DropRejected) as caught:
            find_input(sample, InputKind.POD5)
        message = str(caught.exception)
        self.assertIn("pod5_fail", message)
        self.assertIn("pod5_pass", message)
        self.assertIn("--pod5-subdir", message)

    def test_declaring_the_read_class_resolves_it(self) -> None:
        sample = self._gridion_run()
        found, index = find_input(sample, InputKind.POD5, pod5_subdirectory="pod5_pass")
        self.assertEqual(found.name, "pod5_pass")
        self.assertIsNone(index)

    def test_the_other_read_class_is_equally_available(self) -> None:
        """The point is that the choice is stated, not that one of them is correct."""
        sample = self._gridion_run()
        found, _ = find_input(sample, InputKind.POD5, pod5_subdirectory="pod5_fail")
        self.assertEqual(found.name, "pod5_fail")

    def test_naming_a_directory_that_does_not_exist_is_refused(self) -> None:
        sample = self._gridion_run()
        with self.assertRaises(DropRejected) as caught:
            find_input(sample, InputKind.POD5, pod5_subdirectory="pod5_skipped")
        self.assertIn("pod5_skipped", str(caught.exception))

    def test_an_unsplit_run_needs_no_declaration(self) -> None:
        """A run written without qscore splitting is unambiguous on its own."""
        sample = self._gridion_run(fail_reads=False)
        found, _ = find_input(sample, InputKind.POD5)
        self.assertEqual(found.name, "pod5_pass")

    def test_a_declaration_is_honoured_even_when_nothing_is_split(self) -> None:
        sample = self._gridion_run(fail_reads=False)
        found, _ = find_input(sample, InputKind.POD5, pod5_subdirectory="pod5_pass")
        self.assertEqual(found.name, "pod5_pass")

    def test_the_sample_identifier_comes_from_the_directory_minknow_named(self) -> None:
        """MinKNOW nests <experiment>/<sample>/<run>/, so the watch dir is the experiment."""
        self._gridion_run()
        self.assertEqual(sample_id_from_directory("AML_00123"), "AML_00123")

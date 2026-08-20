"""Tests for the Dorado adapter, which no CI job can execute against the real binary.

That is exactly why these exist. "Never run against Dorado" and "logic never checked" are
different statements, and only the first one is true here. Everything below drives the
adapter through a fake runner: the argument vector it builds, every gate it fails closed
on, and the honesty markers it is required to carry. What remains unverified afterwards is
Dorado's own behaviour — not this module's.
"""

from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.basecall import (
    ADAPTER_VERIFICATION,
    UNVERIFIED_NOTICE,
    BasecallInputs,
    BasecallPolicy,
    build_basecaller_argv,
    model_signature,
    run_basecalling,
)
from ontseq_platform.execution import CommandResult


def _policy(**overrides: object) -> BasecallPolicy:
    values: dict[str, object] = {
        "profile_id": "test_profile",
        "status": "technical_defaults_only",
        "expected_version": "0.9.1",
        "model": "dna_r10.4.1_e8.2_400bps_sup@v5.0.0",
        "note": "test",
    }
    values.update(overrides)
    return BasecallPolicy(**values)  # type: ignore[arg-type]


class _FakeRunner:
    """Records commands, answers the version probe, and writes whatever it is told to."""

    def __init__(
        self,
        *,
        version_output: str = "dorado 0.9.1+abc123",
        version_returncode: int = 0,
        basecall_returncode: int = 0,
        payload: bytes = b"BAM\x01synthetic",
    ) -> None:
        self.version_output = version_output
        self.version_returncode = version_returncode
        self.basecall_returncode = basecall_returncode
        self.payload = payload
        self.commands: list[list[str]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        self.commands.append(list(argv))
        return CommandResult(
            argv=tuple(argv),
            returncode=self.version_returncode,
            stdout=self.version_output,
            stderr="",
        )

    def run_to_file(
        self, argv: Sequence[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult:
        self.commands.append(list(argv))
        if self.basecall_returncode == 0:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_bytes(self.payload)
        return CommandResult(
            argv=tuple(argv),
            returncode=self.basecall_returncode,
            stdout="",
            stderr="dorado: something went wrong",
        )


class PolicyTests(unittest.TestCase):
    def test_the_policy_cannot_claim_validation(self) -> None:
        """The whole point of the marker is that it cannot be set aside by configuration."""
        with self.assertRaises(ValueError) as raised:
            _policy(status="validated")
        self.assertIn("must not be marked validated", str(raised.exception))

    def test_duplicate_modified_base_models_are_refused(self) -> None:
        with self.assertRaises(ValueError):
            _policy(modified_bases=["5mCG_5hmCG", "5mCG_5hmCG"])

    def test_the_version_lock_must_be_a_full_version(self) -> None:
        """A loose lock cannot be enforced, and basecalling is the worst place to guess."""
        with self.assertRaises(ValueError):
            _policy(expected_version="0.9")

    def test_a_model_checksum_must_look_like_one(self) -> None:
        with self.assertRaises(ValueError):
            _policy(model_sha256="not-a-checksum")


class ArgvTests(unittest.TestCase):
    def test_model_and_input_directory_are_positional_and_ordered(self) -> None:
        argv = build_basecaller_argv(
            dorado="dorado", policy=_policy(), pod5_directory=Path("/pod5"), threads=None
        )
        self.assertEqual(argv[:2], ["dorado", "basecaller"])
        self.assertEqual(argv[2], "dna_r10.4.1_e8.2_400bps_sup@v5.0.0")
        self.assertEqual(argv[3], "/pod5")

    def test_the_device_is_always_explicit(self) -> None:
        """Left to a default, a CPU fallback would silently take days instead of failing."""
        argv = build_basecaller_argv(
            dorado="dorado", policy=_policy(), pod5_directory=Path("/pod5"), threads=None
        )
        self.assertIn("--device", argv)
        self.assertEqual(argv[argv.index("--device") + 1], "cuda:all")

    def test_each_modified_base_model_gets_its_own_flag(self) -> None:
        argv = build_basecaller_argv(
            dorado="dorado",
            policy=_policy(modified_bases=["5mCG_5hmCG", "6mA"]),
            pod5_directory=Path("/pod5"),
            threads=None,
        )
        self.assertEqual(argv.count("--modified-bases"), 2)
        self.assertIn("5mCG_5hmCG", argv)
        self.assertIn("6mA", argv)

    def test_optional_flags_are_absent_when_not_configured(self) -> None:
        argv = build_basecaller_argv(
            dorado="dorado", policy=_policy(), pod5_directory=Path("/pod5"), threads=None
        )
        for flag in ("--modified-bases", "--min-qscore", "--emit-moves", "--threads"):
            self.assertNotIn(flag, argv)

    def test_configured_optional_flags_are_present(self) -> None:
        argv = build_basecaller_argv(
            dorado="dorado",
            policy=_policy(minimum_qscore=10, emit_moves=True),
            pod5_directory=Path("/pod5"),
            threads=8,
        )
        self.assertEqual(argv[argv.index("--min-qscore") + 1], "10")
        self.assertIn("--emit-moves", argv)
        self.assertEqual(argv[argv.index("--threads") + 1], "8")

    def test_no_output_path_appears_in_the_command(self) -> None:
        """Dorado writes BAM to stdout; an -o here would silently produce an empty file."""
        argv = build_basecaller_argv(
            dorado="dorado", policy=_policy(), pod5_directory=Path("/pod5"), threads=None
        )
        self.assertNotIn("-o", argv)
        self.assertNotIn("--output-dir", argv)


class ModelSignatureTests(unittest.TestCase):
    def test_a_named_model_has_no_signature_rather_than_a_made_up_one(self) -> None:
        self.assertIsNone(model_signature("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"))

    def test_a_directory_is_fingerprinted_reproducibly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "weights.bin").write_bytes(b"weights")
            (root / "config.toml").write_text("model = 1\n", encoding="utf-8")
            self.assertEqual(model_signature(str(root)), model_signature(str(root)))

    def test_changing_a_file_changes_the_signature(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weights = root / "weights.bin"
            weights.write_bytes(b"weights")
            before = model_signature(str(root))
            weights.write_bytes(b"tampered")
            self.assertNotEqual(before, model_signature(str(root)))

    def test_renaming_a_file_changes_the_signature(self) -> None:
        """Content alone is not enough: a model is its layout as well as its bytes."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "a.bin").write_bytes(b"same")
            before = model_signature(str(root))
            (root / "a.bin").rename(root / "b.bin")
            self.assertNotEqual(before, model_signature(str(root)))


class _RunCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.root = Path(self._temporary.name)
        self.pod5 = self.root / "pod5"
        self.pod5.mkdir()
        (self.pod5 / "batch_0.pod5").write_bytes(b"signal")
        (self.pod5 / "batch_1.pod5").write_bytes(b"signal")
        self.output = self.root / "sample.unaligned.bam"

    def _run(self, policy: BasecallPolicy | None = None, **runner_kwargs: object):
        runner = _FakeRunner(**runner_kwargs)  # type: ignore[arg-type]
        report = run_basecalling(
            BasecallInputs(pod5_directory=self.pod5),
            policy or _policy(),
            sample_id="SYNTHETIC_BC_001",
            output_bam=self.output,
            runner=runner,
        )
        return report, runner


class FailClosedTests(_RunCase):
    def test_a_missing_pod5_directory_is_refused(self) -> None:
        with self.assertRaises(ValueError):
            run_basecalling(
                BasecallInputs(pod5_directory=self.root / "nowhere"),
                _policy(),
                sample_id="S",
                output_bam=self.output,
                runner=_FakeRunner(),
            )

    def test_an_empty_pod5_directory_is_refused_before_any_command_runs(self) -> None:
        empty = self.root / "empty"
        empty.mkdir()
        runner = _FakeRunner()
        with self.assertRaises(ValueError):
            run_basecalling(
                BasecallInputs(pod5_directory=empty),
                _policy(),
                sample_id="S",
                output_bam=self.output,
                runner=runner,
            )
        self.assertEqual(runner.commands, [])

    def test_an_existing_output_is_never_overwritten(self) -> None:
        self.output.write_bytes(b"previous run")
        with self.assertRaises(ValueError):
            self._run()
        self.assertEqual(self.output.read_bytes(), b"previous run")

    def test_a_version_mismatch_stops_the_run(self) -> None:
        with self.assertRaises(ValueError) as raised:
            self._run(version_output="dorado 0.8.0")
        self.assertIn("0.8.0", str(raised.exception))
        self.assertFalse(self.output.exists())

    def test_an_unparseable_version_stops_the_run(self) -> None:
        with self.assertRaises(ValueError):
            self._run(version_output="dorado, the friendly basecaller")

    def test_a_failed_version_probe_stops_the_run(self) -> None:
        with self.assertRaises(ValueError):
            self._run(version_returncode=127)

    def test_a_locked_model_checksum_without_a_local_model_is_refused(self) -> None:
        """A lock that cannot be checked is worse than no lock; it reads as reassurance."""
        with self.assertRaises(ValueError) as raised:
            self._run(_policy(model_sha256="b" * 64))
        self.assertIn("not a local directory", str(raised.exception))

    def test_a_model_that_does_not_match_its_lock_is_refused(self) -> None:
        model_dir = self.root / "model"
        model_dir.mkdir()
        (model_dir / "weights.bin").write_bytes(b"weights")
        with self.assertRaises(ValueError) as raised:
            self._run(_policy(model=str(model_dir), model_sha256="c" * 64))
        self.assertIn("does not match the policy lock", str(raised.exception))

    def test_a_matching_model_lock_is_accepted(self) -> None:
        model_dir = self.root / "model"
        model_dir.mkdir()
        (model_dir / "weights.bin").write_bytes(b"weights")
        signature = model_signature(str(model_dir))
        assert signature is not None
        report, _ = self._run(_policy(model=str(model_dir), model_sha256=signature))
        self.assertEqual(report.tool.parameters["model_sha256"], signature)

    def test_a_failing_basecaller_surfaces_its_diagnostic(self) -> None:
        with self.assertRaises(ValueError) as raised:
            self._run(basecall_returncode=1)
        self.assertIn("something went wrong", str(raised.exception))

    def test_an_empty_output_is_treated_as_a_failure(self) -> None:
        """Exit code zero and no reads is a failure mode, not an empty-but-valid result."""
        with self.assertRaises(ValueError) as raised:
            self._run(payload=b"")
        self.assertIn("produced no output", str(raised.exception))


class ReportTests(_RunCase):
    def test_the_report_is_marked_unverified(self) -> None:
        report, _ = self._run()
        self.assertEqual(report.adapter_verification, ADAPTER_VERIFICATION)
        self.assertIn(UNVERIFIED_NOTICE, report.warnings)
        self.assertIn(UNVERIFIED_NOTICE, report.limitations)

    def test_the_tool_record_carries_the_verification_marker(self) -> None:
        """Provenance readers look at tools, not only at the surrounding report."""
        report, _ = self._run()
        self.assertEqual(report.tool.parameters["adapter_verification"], ADAPTER_VERIFICATION)

    def test_every_pod5_file_is_counted(self) -> None:
        report, _ = self._run()
        self.assertEqual(report.pod5_file_count, 2)

    def test_only_the_file_name_is_recorded_not_its_location(self) -> None:
        report, _ = self._run()
        self.assertEqual(report.unaligned_bam_relative_path, "sample.unaligned.bam")
        self.assertNotIn(str(self.root), report.model_dump_json())

    def test_the_output_is_fingerprinted(self) -> None:
        report, _ = self._run()
        self.assertEqual(report.unaligned_bam_fingerprint.size_bytes, self.output.stat().st_size)

    def test_omitting_modified_bases_is_warned_because_the_loss_is_silent(self) -> None:
        report, _ = self._run()
        self.assertTrue(any("MM/ML" in warning for warning in report.warnings))

    def test_requesting_modified_bases_removes_that_warning(self) -> None:
        report, _ = self._run(_policy(modified_bases=["5mCG_5hmCG"]))
        self.assertFalse(any("MM/ML" in warning for warning in report.warnings))
        self.assertEqual(report.modified_bases_requested, ["5mCG_5hmCG"])

    def test_a_named_model_warns_that_provenance_is_incomplete(self) -> None:
        report, _ = self._run()
        self.assertTrue(any("model checksum" in warning for warning in report.warnings))

    def test_a_fingerprinted_model_does_not_warn(self) -> None:
        model_dir = self.root / "model"
        model_dir.mkdir()
        (model_dir / "weights.bin").write_bytes(b"weights")
        report, _ = self._run(_policy(model=str(model_dir)))
        self.assertFalse(any("model checksum" in warning for warning in report.warnings))

    def test_the_version_probe_precedes_the_basecall(self) -> None:
        _, runner = self._run()
        self.assertEqual(runner.commands[0], ["dorado", "--version"])
        self.assertEqual(runner.commands[1][:2], ["dorado", "basecaller"])


if __name__ == "__main__":
    unittest.main()

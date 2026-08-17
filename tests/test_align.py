from __future__ import annotations

import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.align import (
    AlignmentInputs,
    AlignmentPolicy,
    build_minimap2_argv,
    header_with_read_groups,
    read_group_lines,
    run_alignment,
)
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import GenomeBuild


def _policy(**overrides: object) -> AlignmentPolicy:
    values: dict[str, object] = {
        "profile_id": "test_profile",
        "status": "technical_defaults_only",
        "note": "test",
    }
    values.update(overrides)
    return AlignmentPolicy(**values)  # type: ignore[arg-type]


class PolicyTests(unittest.TestCase):
    def test_hard_clipping_is_refused(self) -> None:
        """SV callers need the clipped sequence, so the policy cannot opt out."""
        with self.assertRaises(ValueError):
            _policy(soft_clip_supplementary=False)

    def test_the_version_lock_must_look_like_a_version(self) -> None:
        """A free-text lock would make the fail-closed version gate unenforceable."""
        with self.assertRaises(ValueError):
            _policy(expected_minimap2_version="latest")


class ArgvTests(unittest.TestCase):
    def _argv(self, policy: AlignmentPolicy) -> list[str]:
        return build_minimap2_argv(
            minimap2="minimap2",
            policy=policy,
            reference_fasta=Path("/ref.fa"),
            reads_fastq=Path("/reads.fastq"),
            output_sam=Path("/out.sam"),
            threads=3,
        )

    def test_the_ont_preset_and_thread_count_are_explicit(self) -> None:
        argv = self._argv(_policy())
        self.assertEqual(argv[:6], ["minimap2", "-a", "-x", "map-ont", "-t", "3"])

    def test_tag_carrying_flags_are_present_by_default(self) -> None:
        argv = self._argv(_policy())
        for flag in ("--MD", "-Y", "-y"):
            self.assertIn(flag, argv)

    def test_dropping_modified_bases_drops_only_that_flag(self) -> None:
        argv = self._argv(_policy(preserve_modified_base_tags=False))
        self.assertNotIn("-y", argv)
        self.assertIn("-Y", argv)

    def test_reference_precedes_reads(self) -> None:
        """minimap2 is positional; swapping these silently aligns the reference."""
        argv = self._argv(_policy())
        self.assertEqual(argv[-2:], ["/ref.fa", "/reads.fastq"])


class ReadGroupHeaderTests(unittest.TestCase):
    HEADER = (
        "@HD\tVN:1.6\tSO:unknown\n"
        "@RG\tID:RG1\tSM:S1\tPL:ONT\n"
        "@RG\tID:RG2\tSM:S1\tPL:ONT\n"
        "@PG\tID:dorado\tPN:dorado\n"
    )

    def test_every_read_group_is_extracted_in_order(self) -> None:
        self.assertEqual(
            read_group_lines(self.HEADER),
            ["@RG\tID:RG1\tSM:S1\tPL:ONT", "@RG\tID:RG2\tSM:S1\tPL:ONT"],
        )

    def test_a_header_without_read_groups_yields_nothing(self) -> None:
        self.assertEqual(read_group_lines("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n"), [])

    def test_a_read_name_starting_with_rg_is_not_mistaken_for_a_header(self) -> None:
        self.assertEqual(read_group_lines("@RGX\tID:no\n"), [])

    def test_read_groups_land_after_the_sequence_dictionary(self) -> None:
        aligned = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n@PG\tID:minimap2\n"
        merged = header_with_read_groups(aligned, read_group_lines(self.HEADER)).splitlines()
        self.assertEqual(
            merged,
            [
                "@HD\tVN:1.6\tSO:coordinate",
                "@SQ\tSN:chr1\tLN:10",
                "@RG\tID:RG1\tSM:S1\tPL:ONT",
                "@RG\tID:RG2\tSM:S1\tPL:ONT",
                "@PG\tID:minimap2",
            ],
        )

    def test_all_sequence_records_stay_above_the_read_groups(self) -> None:
        aligned = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n@SQ\tSN:chr2\tLN:20\n@PG\tID:minimap2\n"
        merged = header_with_read_groups(aligned, ["@RG\tID:RG1"]).splitlines()
        self.assertEqual(merged.index("@RG\tID:RG1"), 3)

    def test_reapplying_the_same_read_groups_is_idempotent(self) -> None:
        aligned = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n@PG\tID:minimap2\n"
        groups = read_group_lines(self.HEADER)
        once = header_with_read_groups(aligned, groups)
        self.assertEqual(header_with_read_groups(once, groups), once)

    def test_a_header_with_no_sequence_records_still_places_read_groups(self) -> None:
        merged = header_with_read_groups("@PG\tID:minimap2\n", ["@RG\tID:RG1"]).splitlines()
        self.assertEqual(merged, ["@RG\tID:RG1", "@PG\tID:minimap2"])


class _FakeSamtools:
    """Stands in for samtools and minimap2, dispatching on the argument vector.

    Every filesystem side effect a real run would have is reproduced, because the adapter
    checks for them: sort writes a BAM, reheader writes a BAM, index writes an index. What
    is *not* reproduced is any alignment — these tests pin the command sequence and the
    branching, not minimap2's behaviour, which CI exercises for real.
    """

    def __init__(
        self,
        *,
        source_header: str,
        minimap2_version: str = "2.28-r1209",
        samtools_version: str = "samtools 1.24\nUsing htslib 1.24",
        failing_label: str | None = None,
    ) -> None:
        self.source_header = source_header
        self.minimap2_version = minimap2_version
        self.samtools_version = samtools_version
        self.failing_label = failing_label
        self.commands: list[list[str]] = []

    def _label(self, argv: list[str]) -> str:
        """Name a command precisely enough that a version probe is not mistaken for work."""
        if argv[0] == "minimap2":
            return "minimap2 --version" if "--version" in argv else "minimap2 align"
        return f"samtools {argv[1]}" if len(argv) > 1 else "samtools"

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        argv = list(argv)
        self.commands.append(argv)
        label = self._label(argv)
        returncode = 1 if label == self.failing_label else 0
        stdout = ""
        if argv[:2] == ["minimap2", "--version"]:
            stdout = self.minimap2_version
        elif argv[:2] == ["samtools", "--version"]:
            stdout = self.samtools_version
        elif argv[:3] == ["samtools", "view", "-H"]:
            target = Path(argv[3])
            stdout = (
                self.source_header
                if target.name.endswith("unaligned.bam")
                else "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@PG\tID:minimap2\n"
            )
        elif argv[1] == "fastq":
            Path(argv[argv.index("-0") + 1]).write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
        elif argv[1] == "sort" and returncode == 0:
            Path(argv[argv.index("-o") + 1]).write_bytes(b"BAM\x01sorted")
        elif argv[1] == "index" and returncode == 0:
            Path(argv[-1]).write_bytes(b"BAI")
        return CommandResult(argv=tuple(argv), returncode=returncode, stdout=stdout, stderr="")

    def run_to_file(
        self, argv: Sequence[str], output_path: Path, *, timeout_seconds: int = 300
    ) -> CommandResult:
        argv = list(argv)
        self.commands.append(argv)
        returncode = 1 if self._label(argv) == self.failing_label else 0
        if returncode == 0:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_bytes(b"BAM\x01reheadered")
        return CommandResult(argv=tuple(argv), returncode=returncode, stdout="", stderr="failed")

    def labels(self) -> list[str]:
        return [self._label(argv) for argv in self.commands]


ONE_READ_GROUP = (
    "@HD\tVN:1.6\tSO:unknown\n@RG\tID:RG1\tSM:SAMPLE_A\tPL:ONT\n@PG\tID:dorado\tPN:dorado\n"
)
NO_READ_GROUP = "@HD\tVN:1.6\tSO:unknown\n@PG\tID:dorado\tPN:dorado\n"


class RunAlignmentCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.root = Path(self._temporary.name)
        self.unaligned = self.root / "sample.unaligned.bam"
        self.unaligned.write_bytes(b"BAM\x01unaligned")
        self.reference = self.root / "reference.fasta"
        self.reference.write_text(">chr1\nACGT\n", encoding="utf-8")
        self.scratch = self.root / "work"
        self.output = self.root / "alignment" / "sample.bam"

    def _run(self, runner: _FakeSamtools, policy: AlignmentPolicy | None = None):
        return run_alignment(
            AlignmentInputs(unaligned_bam=self.unaligned, reference_fasta=self.reference),
            policy or _policy(),
            sample_id="SAMPLE_A",
            genome_build=GenomeBuild.GRCH38,
            reference_id="REF_V1",
            scratch_dir=self.scratch,
            output_bam=self.output,
            runner=runner,
            threads=2,
        )


class ReadGroupPreservationTests(RunAlignmentCase):
    def test_read_groups_are_carried_and_reattached(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP)
        report = self._run(runner)
        self.assertEqual(report.read_group_count, 1)
        self.assertIn("samtools reheader", runner.labels())

    def test_the_read_group_tag_rides_the_fastq_comment(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP)
        self._run(runner)
        fastq = next(argv for argv in runner.commands if argv[1] == "fastq")
        self.assertEqual(fastq[fastq.index("-T") + 1], "MM,ML,RG")

    def test_without_read_groups_nothing_is_reattached(self) -> None:
        runner = _FakeSamtools(source_header=NO_READ_GROUP)
        report = self._run(runner)
        self.assertEqual(report.read_group_count, 0)
        self.assertNotIn("samtools reheader", runner.labels())

    def test_without_read_groups_only_modified_base_tags_are_carried(self) -> None:
        runner = _FakeSamtools(source_header=NO_READ_GROUP)
        self._run(runner)
        fastq = next(argv for argv in runner.commands if argv[1] == "fastq")
        self.assertEqual(fastq[fastq.index("-T") + 1], "MM,ML")

    def test_missing_read_groups_are_warned_rather_than_passed_over(self) -> None:
        runner = _FakeSamtools(source_header=NO_READ_GROUP)
        report = self._run(runner)
        self.assertTrue(any("@RG" in warning for warning in report.warnings))

    def test_minimap2_is_never_asked_to_stamp_a_read_group(self) -> None:
        """-R would collapse distinct read groups into one and contradict the per-read tags."""
        runner = _FakeSamtools(source_header=ONE_READ_GROUP)
        self._run(runner)
        minimap2 = next(argv for argv in runner.commands if argv[0] == "minimap2" and "-a" in argv)
        self.assertNotIn("-R", minimap2)

    def test_dropping_modified_base_tags_still_carries_read_groups(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP)
        self._run(runner, _policy(preserve_modified_base_tags=False))
        fastq = next(argv for argv in runner.commands if argv[1] == "fastq")
        self.assertEqual(fastq[fastq.index("-T") + 1], "RG")


class CommandSequenceTests(RunAlignmentCase):
    def test_the_full_sequence_is_probed_then_executed_in_order(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP)
        self._run(runner)
        self.assertEqual(
            runner.labels(),
            [
                "minimap2 --version",
                "samtools --version",
                "samtools view",
                "samtools fastq",
                "minimap2 align",
                "samtools sort",
                "samtools view",
                "samtools reheader",
                "samtools index",
            ],
        )

    def test_scratch_files_do_not_survive_a_successful_run(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP)
        self._run(runner)
        self.assertEqual(sorted(path.name for path in self.scratch.iterdir()), [])

    def test_the_final_bam_and_index_exist(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP)
        report = self._run(runner)
        self.assertTrue(self.output.is_file())
        self.assertTrue(Path(f"{self.output}.bai").is_file())
        self.assertEqual(report.aligned_bam_relative_path, "sample.bam")

    def test_only_file_names_are_recorded_not_locations(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP)
        report = self._run(runner)
        self.assertNotIn(str(self.root), report.model_dump_json())


class FailClosedTests(RunAlignmentCase):
    def test_a_minimap2_version_mismatch_stops_the_run(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP, minimap2_version="2.26-r1175")
        with self.assertRaises(ValueError) as raised:
            self._run(runner)
        self.assertIn("2.26", str(raised.exception))
        self.assertFalse(self.output.exists())

    def test_a_samtools_version_mismatch_stops_the_run(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP, samtools_version="samtools 1.20")
        with self.assertRaises(ValueError):
            self._run(runner)

    def test_an_existing_output_is_never_overwritten(self) -> None:
        self.output.parent.mkdir(parents=True, exist_ok=True)
        self.output.write_bytes(b"previous run")
        with self.assertRaises(ValueError):
            self._run(_FakeSamtools(source_header=ONE_READ_GROUP))
        self.assertEqual(self.output.read_bytes(), b"previous run")

    def test_a_missing_unaligned_bam_is_refused(self) -> None:
        self.unaligned.unlink()
        with self.assertRaises(ValueError):
            self._run(_FakeSamtools(source_header=ONE_READ_GROUP))

    def test_a_missing_reference_is_refused(self) -> None:
        self.reference.unlink()
        with self.assertRaises(ValueError):
            self._run(_FakeSamtools(source_header=ONE_READ_GROUP))

    def test_a_failing_reheader_leaves_no_final_bam(self) -> None:
        """The BAM must not appear unless every step that shapes it succeeded."""
        runner = _FakeSamtools(source_header=ONE_READ_GROUP, failing_label="samtools reheader")
        with self.assertRaises(ValueError):
            self._run(runner)
        self.assertFalse(self.output.exists())

    def test_a_failing_alignment_stops_before_sorting(self) -> None:
        runner = _FakeSamtools(source_header=ONE_READ_GROUP, failing_label="minimap2 align")
        with self.assertRaises(ValueError):
            self._run(runner)
        self.assertNotIn("samtools sort", runner.labels())

    def test_threads_must_be_positive(self) -> None:
        with self.assertRaises(ValueError):
            run_alignment(
                AlignmentInputs(unaligned_bam=self.unaligned, reference_fasta=self.reference),
                _policy(),
                sample_id="SAMPLE_A",
                genome_build=GenomeBuild.GRCH38,
                reference_id="REF_V1",
                scratch_dir=self.scratch,
                output_bam=self.output,
                runner=_FakeSamtools(source_header=ONE_READ_GROUP),
                threads=0,
            )


if __name__ == "__main__":
    unittest.main()

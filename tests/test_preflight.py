from __future__ import annotations

import hashlib
import json
import os
import socket
import subprocess
import sys
import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.align import AlignmentPolicy
from ontseq_platform.basecall import BasecallPolicy, model_signature
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import ReferenceLock, SampleManifest
from ontseq_platform.pipeline.checks import Check, CheckStatus
from ontseq_platform.pipeline.lock import LOCK_FILENAME
from ontseq_platform.pipeline.stages import StageId, VerificationStatus
from ontseq_platform.preflight import PreflightRequest, preflight
from ontseq_platform.target_coverage import TargetCoveragePolicy

ALIGNMENT_POLICY = AlignmentPolicy(
    profile_id="test",
    status="technical_defaults_only",
    expected_minimap2_version="2.28",
    expected_samtools_version="1.24",
    note="test",
)

TARGET_COVERAGE_POLICY = TargetCoveragePolicy(
    profile_id="test",
    status="technical_defaults_only",
    expected_version="0.3.14",
    note="test",
)

BASECALL_POLICY = BasecallPolicy(
    profile_id="test",
    status="technical_defaults_only",
    expected_version="0.9.0",
    model="dna_r10.4.1_e8.2_400bps_hac@v5.0.0",
    modified_bases=["5mCG_5hmCG"],
    note="test",
)

DEFAULT_VERSIONS = {
    "minimap2": "2.28-r1209",
    "samtools": "samtools 1.24\nUsing htslib 1.24",
    "cramino": "cramino 0.14.5",
    "sniffles": "Sniffles2, Version 2.8.0",
    "dorado": "dorado 0.9.0+abcdef",
    "mosdepth": "mosdepth 0.3.14",
}


class FakeRunner:
    """Answers ``--version`` for the tools a test says exist, and nothing else."""

    def __init__(self, versions: dict[str, str]) -> None:
        self.versions = versions
        self.calls: list[list[str]] = []

    def run(self, argv: Sequence[str], *, timeout_seconds: int = 300) -> CommandResult:
        self.calls.append(list(argv))
        name = Path(argv[0]).name
        if name not in self.versions:
            return CommandResult(argv=tuple(argv), returncode=127, stdout="", stderr="not found")
        return CommandResult(argv=tuple(argv), returncode=0, stdout=self.versions[name], stderr="")


class PreflightCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.root = Path(self._temporary.name)
        self.output = self.root / "runs"

        # Every tool resolves to a real file, so shutil.which is never consulted and the
        # tests do not depend on what happens to be installed on the machine running them.
        self.bin = self.root / "bin"
        self.bin.mkdir()
        for name in DEFAULT_VERSIONS:
            executable = self.bin / name
            executable.write_text("#!/bin/sh\n", encoding="utf-8")
            executable.chmod(0o755)

        self.bam = self.root / "sample.bam"
        self.bam.write_bytes(b"BAM\x01payload")
        self.bai = self.root / "sample.bam.bai"
        self.bai.write_bytes(b"BAI")

        self.fasta = self.root / "reference.fa"
        self.fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
        self.fai = self.root / "reference.fa.fai"
        self.fai.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
        self.fai_sha256 = hashlib.sha256(self.fai.read_bytes()).hexdigest()

    def manifest(
        self,
        kind: str = "aligned_bam",
        *,
        assay: dict[str, object] | None = None,
        modules: list[str] | None = None,
    ) -> SampleManifest:
        payload: dict[str, object] = {
            "schema_version": "0.1.0",
            "sample_id": "SAMPLE_A",
            "run_id": "RUN_001",
            "input": {"kind": kind, "path": str(self.bam)},
            "assay": assay or {"mode": "lcwgs", "genome_build": "GRCh38", "reference_id": "REF_V1"},
            "analysis": {"profile": "lcwgs", "modules": modules or ["qc"]},
        }
        if kind == "aligned_bam":
            payload["input"] = {"kind": kind, "path": str(self.bam), "index_path": str(self.bai)}
        return SampleManifest.model_validate(payload)

    def lock(self, *, fai_sha256: str | None = None, reference_id: str = "REF_V1") -> ReferenceLock:
        return ReferenceLock.model_validate(
            {
                "reference_id": reference_id,
                "genome_build": "GRCh38",
                "contigs": [{"name": "chr1", "length": 4}],
                "source_fai_sha256": fai_sha256 or self.fai_sha256,
            }
        )

    def request(self, **overrides: object) -> PreflightRequest:
        values: dict[str, object] = {
            "manifest": self.manifest(),
            "reference_lock": self.lock(),
            "output_base": self.output,
            "run_id": "RUN_001",
            "executables": {name: str(self.bin / name) for name in DEFAULT_VERSIONS},
            "reference_fasta": self.fasta,
            "alignment_policy": ALIGNMENT_POLICY,
            "sniffles_policy": None,
            "basecall_policy": BASECALL_POLICY,
        }
        values.update(overrides)
        return PreflightRequest(**values)  # type: ignore[arg-type]

    def results(
        self,
        request: PreflightRequest | None = None,
        *,
        versions: dict[str, str] | None = None,
    ) -> dict[str, Check]:
        runner = FakeRunner(versions if versions is not None else dict(DEFAULT_VERSIONS))
        return {check.name: check for check in preflight(request or self.request(), runner=runner)}

    def pod5_request(self, **overrides: object) -> PreflightRequest:
        directory = self.root / "pod5"
        directory.mkdir(exist_ok=True)
        (directory / "a.pod5").write_bytes(b"POD5")
        values: dict[str, object] = {
            "manifest": self.manifest("pod5"),
            "pod5_directory": directory,
        }
        values.update(overrides)
        return self.request(**values)


class InputTests(PreflightCase):
    def test_a_complete_aligned_bam_run_passes_its_input_checks(self) -> None:
        found = self.results()
        self.assertIs(found["input.bam"].status, CheckStatus.OK)
        self.assertIs(found["input.bam.index"].status, CheckStatus.OK)

    def test_a_missing_input_fails(self) -> None:
        self.bam.unlink()
        self.assertIs(self.results()["input.bam"].status, CheckStatus.FAILED)

    def test_an_empty_input_fails(self) -> None:
        self.bam.write_bytes(b"")
        self.assertIs(self.results()["input.bam"].status, CheckStatus.FAILED)

    def test_a_missing_index_fails_with_the_command_that_makes_one(self) -> None:
        self.bai.unlink()
        check = self.results()["input.bam.index"]
        self.assertIs(check.status, CheckStatus.FAILED)
        self.assertIn("samtools index", check.remedy)

    def test_an_index_older_than_its_bam_warns_rather_than_blocks(self) -> None:
        """It may still be correct, so refusing to start would be presumptuous."""
        os.utime(self.bai, (1_000_000, 1_000_000))
        self.assertIs(self.results()["input.bam.index"].status, CheckStatus.WARNING)

    def test_an_unaligned_bam_run_needs_no_index(self) -> None:
        request = self.request(manifest=self.manifest("unaligned_bam"))
        self.assertIs(self.results(request)["input.bam.index"].status, CheckStatus.SKIPPED)

    def test_a_pod5_run_without_a_directory_fails(self) -> None:
        request = self.request(manifest=self.manifest("pod5"))
        check = self.results(request)["input.pod5"]
        self.assertIs(check.status, CheckStatus.FAILED)
        self.assertIn("--pod5-dir", check.remedy)

    def test_a_pod5_directory_with_no_signal_files_fails(self) -> None:
        empty = self.root / "empty-pod5"
        empty.mkdir()
        request = self.request(manifest=self.manifest("pod5"), pod5_directory=empty)
        self.assertIs(self.results(request)["input.pod5"].status, CheckStatus.FAILED)

    def test_pod5_files_are_found_below_the_top_level(self) -> None:
        directory = self.root / "pod5"
        (directory / "nested").mkdir(parents=True)
        (directory / "nested" / "a.pod5").write_bytes(b"POD5")
        request = self.request(manifest=self.manifest("pod5"), pod5_directory=directory)
        self.assertIs(self.results(request)["input.pod5"].status, CheckStatus.OK)


class ReferenceTests(PreflightCase):
    def test_a_matching_reference_passes(self) -> None:
        found = self.results(self.request(manifest=self.manifest("unaligned_bam")))
        self.assertIs(found["reference.fasta"].status, CheckStatus.OK)
        self.assertIs(found["reference.fai"].status, CheckStatus.OK)

    def test_a_reference_that_is_not_the_locked_one_fails(self) -> None:
        """The classic silent disaster: right chromosome names, wrong assembly."""
        request = self.request(
            manifest=self.manifest("unaligned_bam"),
            reference_lock=self.lock(fai_sha256="b" * 64),
        )
        check = self.results(request)["reference.fai"]
        self.assertIs(check.status, CheckStatus.FAILED)
        self.assertIn("lock", check.detail)

    def test_a_missing_fai_fails_with_the_command_that_makes_one(self) -> None:
        self.fai.unlink()
        request = self.request(manifest=self.manifest("unaligned_bam"))
        self.assertIn("samtools faidx", self.results(request)["reference.fai"].remedy)

    def test_an_aligned_bam_run_needs_no_reference_fasta(self) -> None:
        found = self.results()
        self.assertIs(found["reference.fasta"].status, CheckStatus.SKIPPED)
        self.assertIs(found["reference.fai"].status, CheckStatus.SKIPPED)

    def test_an_aligning_run_without_a_fasta_fails(self) -> None:
        request = self.request(manifest=self.manifest("unaligned_bam"), reference_fasta=None)
        self.assertIs(self.results(request)["reference.fasta"].status, CheckStatus.FAILED)

    def test_a_manifest_and_lock_that_disagree_on_the_reference_fail(self) -> None:
        request = self.request(reference_lock=self.lock(reference_id="SOMETHING_ELSE"))
        self.assertIs(self.results(request)["reference.id"].status, CheckStatus.FAILED)

    def test_the_reference_checks_are_always_present(self) -> None:
        """A varying set of check names would break a scheduler reading the JSON."""
        for kind in ("aligned_bam", "unaligned_bam"):
            found = self.results(self.request(manifest=self.manifest(kind)))
            for name in ("reference.build", "reference.id", "reference.fasta", "reference.fai"):
                self.assertIn(name, found, f"{name} missing for {kind}")


class ToolTests(PreflightCase):
    def test_locked_versions_that_match_pass(self) -> None:
        found = self.results()
        self.assertIs(found["tool.samtools"].status, CheckStatus.OK)
        self.assertIs(found["tool.cramino"].status, CheckStatus.OK)

    def test_a_version_that_does_not_match_the_lock_fails(self) -> None:
        versions = {**DEFAULT_VERSIONS, "samtools": "samtools 1.19\nUsing htslib 1.19"}
        request = self.request(manifest=self.manifest("unaligned_bam"))
        check = self.results(request, versions=versions)["tool.samtools"]
        self.assertIs(check.status, CheckStatus.FAILED)

    def test_a_lock_no_planned_stage_enforces_is_not_enforced_here(self) -> None:
        """The alignment policy locks samtools, but an aligned-BAM run never aligns.

        Enforcing it anyway would refuse runs that `ontseq run` completes, which is the
        one failure mode a preflight must not have.
        """
        versions = {**DEFAULT_VERSIONS, "samtools": "samtools 1.19\nUsing htslib 1.19"}
        check = self.results(versions=versions)["tool.samtools"]
        self.assertIs(check.status, CheckStatus.OK)
        self.assertIn("no policy lock", check.detail)

    def test_a_missing_required_tool_fails(self) -> None:
        (self.bin / "samtools").unlink()
        self.assertIs(self.results()["tool.samtools"].status, CheckStatus.FAILED)

    def test_a_missing_optional_tool_only_warns(self) -> None:
        """Sniffles serves the optional SV stage, so the run completes without it."""
        (self.bin / "sniffles").unlink()
        request = self.request(manifest=self.manifest(modules=["qc", "sv"]))
        check = self.results(request)["tool.sniffles"]
        self.assertIs(check.status, CheckStatus.WARNING)
        self.assertIn("NOT_RUN", check.detail)

    def test_a_run_that_did_not_ask_for_sv_is_not_told_about_its_callers(self) -> None:
        """A CNV-only run needs neither the SV binaries nor advice about them.

        Reporting a missing sniffles to a run whose manifest never requested structural
        variants trains an operator to ignore the tool section, which is the one section
        that has to keep meaning something.
        """
        (self.bin / "sniffles").unlink()
        results = self.results()
        self.assertNotIn("tool.sniffles", results)
        self.assertNotIn("tool.cutesv", results)
        self.assertIs(results["sv.callers"].status, CheckStatus.SKIPPED)

    def test_a_tool_whose_probe_fails_is_reported_as_unidentifiable(self) -> None:
        versions = {name: text for name, text in DEFAULT_VERSIONS.items() if name != "minimap2"}
        request = self.request(manifest=self.manifest("unaligned_bam"))
        check = self.results(request, versions=versions)["tool.minimap2"]
        self.assertIs(check.status, CheckStatus.FAILED)
        self.assertIn("could not be identified", check.detail)

    def test_an_aligned_bam_run_never_probes_the_basecaller(self) -> None:
        runner = FakeRunner(dict(DEFAULT_VERSIONS))
        preflight(self.request(), runner=runner)
        self.assertNotIn("dorado", [Path(call[0]).name for call in runner.calls])

    def test_a_tool_with_no_policy_loaded_passes_on_being_identifiable(self) -> None:
        request = self.request(manifest=self.manifest("unaligned_bam"), alignment_policy=None)
        check = self.results(request)["tool.minimap2"]
        self.assertIs(check.status, CheckStatus.OK)
        self.assertIn("no policy lock", check.detail)

    def test_cramino_is_accepted_exactly_as_the_qc_adapter_accepts_it(self) -> None:
        """No policy locks Cramino, and its adapter tolerates unparseable output; so must this.

        The invariant is agreement, not strictness: a preflight that rejected what the run
        accepts would block runs that were always going to work.
        """
        versions = {**DEFAULT_VERSIONS, "cramino": ""}
        self.assertIs(self.results(versions=versions)["tool.cramino"].status, CheckStatus.OK)


class BasecallTests(PreflightCase):
    def test_an_aligned_bam_run_skips_every_basecalling_check(self) -> None:
        found = self.results()
        self.assertIs(found["basecall.model"].status, CheckStatus.SKIPPED)
        self.assertIs(found["basecall.modified_bases"].status, CheckStatus.SKIPPED)

    def test_a_named_model_warns_that_provenance_cannot_record_the_weights(self) -> None:
        check = self.results(self.pod5_request())["basecall.model"]
        self.assertIs(check.status, CheckStatus.WARNING)
        self.assertIn("provenance", check.detail)

    def test_a_model_directory_that_matches_its_lock_passes(self) -> None:
        model = self.root / "model"
        model.mkdir()
        (model / "weights.bin").write_bytes(b"weights")
        policy = BASECALL_POLICY.model_copy(
            update={"model": str(model), "model_sha256": model_signature(str(model))}
        )
        self.assertIs(
            self.results(self.pod5_request(basecall_policy=policy))["basecall.model"].status,
            CheckStatus.OK,
        )

    def test_a_model_directory_that_does_not_match_its_lock_fails(self) -> None:
        """Hours of basecalling with the wrong weights is unrecoverable and undetectable."""
        model = self.root / "model"
        model.mkdir()
        (model / "weights.bin").write_bytes(b"weights")
        policy = BASECALL_POLICY.model_copy(update={"model": str(model), "model_sha256": "c" * 64})
        self.assertIs(
            self.results(self.pod5_request(basecall_policy=policy))["basecall.model"].status,
            CheckStatus.FAILED,
        )

    def test_a_locked_checksum_on_a_named_model_fails(self) -> None:
        policy = BASECALL_POLICY.model_copy(update={"model_sha256": "c" * 64})
        self.assertIs(
            self.results(self.pod5_request(basecall_policy=policy))["basecall.model"].status,
            CheckStatus.FAILED,
        )

    def test_no_modified_base_model_warns_before_the_run_rather_than_after(self) -> None:
        policy = BASECALL_POLICY.model_copy(update={"modified_bases": []})
        check = self.results(self.pod5_request(basecall_policy=policy))["basecall.modified_bases"]
        self.assertIs(check.status, CheckStatus.WARNING)
        self.assertIn("methylation", check.detail)

    def test_a_pod5_run_without_a_basecalling_policy_fails(self) -> None:
        request = self.pod5_request(basecall_policy=None)
        self.assertIs(self.results(request)["basecall.model"].status, CheckStatus.FAILED)

    def test_a_dorado_version_that_does_not_match_the_lock_fails(self) -> None:
        versions = {**DEFAULT_VERSIONS, "dorado": "dorado 0.8.3"}
        check = self.results(self.pod5_request(), versions=versions)["tool.dorado"]
        self.assertIs(check.status, CheckStatus.FAILED)


class EnvelopeTests(PreflightCase):
    def _lock_file(self, *, pid: int, hostname: str | None = None) -> None:
        root = self.output / "RUN_001" / "SAMPLE_A"
        root.mkdir(parents=True, exist_ok=True)
        (root / LOCK_FILENAME).write_text(
            json.dumps(
                {
                    "pid": pid,
                    "hostname": hostname or socket.gethostname(),
                    "acquired_at": "2026-08-17T09:00:00+00:00",
                    "run_id": "RUN_001",
                    "sample_id": "SAMPLE_A",
                    "pipeline_version": "0.0.0-test",
                }
            ),
            encoding="utf-8",
        )

    def test_a_free_envelope_passes(self) -> None:
        self.assertIs(self.results()["envelope.lock"].status, CheckStatus.OK)

    def test_a_live_lock_blocks_because_the_run_would_exit_four(self) -> None:
        self._lock_file(pid=os.getpid())
        self.assertIs(self.results()["envelope.lock"].status, CheckStatus.FAILED)

    def test_a_dead_lock_only_warns_because_the_run_reclaims_it(self) -> None:
        process = subprocess.Popen([sys.executable, "-c", ""])
        process.wait()
        self._lock_file(pid=process.pid)
        check = self.results()["envelope.lock"]
        self.assertIs(check.status, CheckStatus.WARNING)
        self.assertIn("resumes", check.detail)

    def test_preflight_leaves_nothing_behind(self) -> None:
        """A run that was never started must not leave an output tree suggesting it was."""
        self.results()
        self.assertFalse(self.output.exists())

    def test_writability_is_checked_on_the_nearest_existing_ancestor(self) -> None:
        request = self.request(output_base=self.root / "deep" / "nested" / "runs")
        check = self.results(request)["output.writable"]
        self.assertIs(check.status, CheckStatus.OK)
        self.assertFalse((self.root / "deep").exists())

    @unittest.skipIf(getattr(os, "geteuid", lambda: -1)() == 0, "root ignores permissions")
    def test_an_output_directory_that_cannot_be_written_fails(self) -> None:
        blocked = self.root / "blocked"
        blocked.mkdir()
        blocked.chmod(0o500)
        self.addCleanup(blocked.chmod, 0o700)
        request = self.request(output_base=blocked / "runs")
        self.assertIs(self.results(request)["output.writable"].status, CheckStatus.FAILED)


class DiskTests(PreflightCase):
    def test_free_space_is_reported_not_judged_without_a_requirement(self) -> None:
        """No measured size model exists here, and an invented one would look validated."""
        check = self.results()["disk.free"]
        self.assertIs(check.status, CheckStatus.UNKNOWN)
        self.assertIn("--require-free-gb", check.remedy)

    def test_a_stated_requirement_that_is_met_passes(self) -> None:
        request = self.request(require_free_gb=0.000001)
        self.assertIs(self.results(request)["disk.free"].status, CheckStatus.OK)

    def test_a_stated_requirement_that_is_not_met_fails(self) -> None:
        request = self.request(require_free_gb=1_000_000_000.0)
        self.assertIs(self.results(request)["disk.free"].status, CheckStatus.FAILED)


class AdapterTests(PreflightCase):
    def test_an_aligned_bam_run_runs_only_exercised_adapters(self) -> None:
        self.assertIs(self.results()["adapters.verification"].status, CheckStatus.OK)

    def test_a_pod5_run_is_told_basecalling_has_never_been_executed(self) -> None:
        check = self.results(self.pod5_request())["adapters.verification"]
        self.assertIs(check.status, CheckStatus.WARNING)
        self.assertIn("basecall", check.detail)

    def test_stages_with_no_adapter_are_reported_separately(self) -> None:
        """A stage with no adapter records NOT_RUN; that is not a claim about code quality.

        Collapsing the two would tell an operator that CNV rests on unexecuted code, when
        in fact no CNV caller is wired in at all — a materially different claim.
        """
        check = self.results()["stages.not_implemented"]
        self.assertIs(check.status, CheckStatus.WARNING)
        self.assertIn("cnv", check.detail)
        self.assertNotIn("target_coverage", check.detail)
        self.assertIn("not a negative biological finding", check.detail)

    def test_a_registered_adapter_is_reported_as_the_run_will_see_it(self) -> None:
        """Preflight must describe the run it is checking, not the bare declared graph.

        The CNV lane is installed by registration, so the graph alone calls it
        `not_implemented`. Preflight used not to register and reported exactly that — for
        runs that went on to execute a real QDNAseq/ACE analysis. Announcing "this stage
        has no adapter and will record NOT_RUN" about a stage that then runs against real
        R tooling is the one failure a preflight must not have.
        """
        request = self.request(
            stage_verification={StageId.CNV: VerificationStatus.VERIFIED_WITH_REAL_TOOL}
        )
        found = self.results(request)
        self.assertIs(found["stages.not_implemented"].status, CheckStatus.OK)
        self.assertNotIn("cnv", found["stages.not_implemented"].detail)
        self.assertIs(found["adapters.verification"].status, CheckStatus.OK)

    def test_an_unverified_registered_adapter_is_still_flagged(self) -> None:
        """Overriding the graph is not a way to silence the warning, only to correct it."""
        request = self.request(
            stage_verification={StageId.CNV: VerificationStatus.UNVERIFIED_ADAPTER}
        )
        found = self.results(request)
        check = found["adapters.verification"]
        self.assertIs(check.status, CheckStatus.WARNING)
        self.assertIn("cnv", check.detail)
        self.assertNotIn("cnv", found["stages.not_implemented"].detail)

    def test_the_two_adapter_claims_never_name_the_same_stage(self) -> None:
        found = self.results(self.pod5_request())
        self.assertIn("basecall", found["adapters.verification"].detail)
        self.assertNotIn("basecall", found["stages.not_implemented"].detail)


class ReportingTests(PreflightCase):
    def test_every_problem_is_reported_not_just_the_first(self) -> None:
        """An operator fixing a broken setup needs the whole list, not one per attempt."""
        self.bam.unlink()
        (self.bin / "samtools").unlink()
        request = self.request(reference_lock=self.lock(fai_sha256="d" * 64))
        failed = [
            name
            for name, check in self.results(request).items()
            if check.status is CheckStatus.FAILED
        ]
        self.assertIn("input.bam", failed)
        self.assertIn("tool.samtools", failed)

    def test_a_clean_aligned_bam_run_blocks_on_nothing(self) -> None:
        blocking = [
            name for name, check in self.results().items() if check.status is CheckStatus.FAILED
        ]
        self.assertEqual(blocking, [])


class TargetCoverageTests(PreflightCase):
    """Adaptive sampling has preconditions the run fails closed on. Preflight must too.

    ``ontseq run`` refuses an adaptive-sampling run without a target-coverage policy or a
    readable target BED, and probes Mosdepth before the stage. All three refusals happen
    after the envelope exists and the lock is taken, and a FAILED target-coverage stage
    fails the whole run — so a preflight that stayed silent about them would clear a run
    that could not succeed.
    """

    def adaptive_manifest(self, target_bed: Path | None) -> SampleManifest:
        assay: dict[str, object] = {
            "mode": "adaptive_sampling",
            "genome_build": "GRCh38",
            "reference_id": "REF_V1",
            "target_bed": str(target_bed) if target_bed is not None else None,
            "target_bed_version": "ROI_V1",
        }
        return self.manifest(assay=assay)

    def target_bed(self, text: str = "chr1\t1000\t2000\tTARGET_A\n", *, name: str = "roi") -> Path:
        """Write a target BED under its own name, so two in one test cannot overwrite each other."""
        bed = self.root / f"{name}.bed"
        bed.write_text(text, encoding="utf-8")
        return bed

    def adaptive_request(self, **overrides: object) -> PreflightRequest:
        values: dict[str, object] = {"target_coverage_policy": TARGET_COVERAGE_POLICY}
        values.update(overrides)
        # Built only when the test did not supply one, so writing the default BED cannot
        # overwrite a BED the test wrote for the manifest it passed in.
        values.setdefault("manifest", self.adaptive_manifest(self.target_bed()))
        return self.request(**values)

    def test_a_complete_adaptive_run_passes_its_target_checks(self) -> None:
        found = self.results(self.adaptive_request())
        self.assertIs(found["target_coverage.policy"].status, CheckStatus.OK)
        self.assertIs(found["target_coverage.bed"].status, CheckStatus.OK)
        self.assertIs(found["tool.mosdepth"].status, CheckStatus.OK)

    def test_the_bed_check_reports_what_was_actually_read(self) -> None:
        bed = self.target_bed("chr1\t1000\t2000\tA\nchr2\t0\t500\tB\n", name="two")
        found = self.results(self.adaptive_request(manifest=self.adaptive_manifest(bed)))
        self.assertIn("2 target(s)", found["target_coverage.bed"].detail)
        self.assertIn("1500 bp", found["target_coverage.bed"].detail)

    def test_a_missing_target_coverage_policy_blocks_the_run(self) -> None:
        found = self.results(self.adaptive_request(target_coverage_policy=None))
        self.assertIs(found["target_coverage.policy"].status, CheckStatus.FAILED)
        self.assertTrue(found["target_coverage.policy"].remedy)

    def test_an_absent_target_bed_blocks_the_run(self) -> None:
        manifest = self.adaptive_manifest(self.root / "never-written.bed")
        found = self.results(self.adaptive_request(manifest=manifest))
        self.assertIs(found["target_coverage.bed"].status, CheckStatus.FAILED)

    def test_an_unparseable_target_bed_blocks_the_run(self) -> None:
        """The BED is parsed, not merely stat'ed: a truncated ROI fails the stage."""
        bed = self.target_bed("chr1\t1000\n", name="truncated")
        found = self.results(self.adaptive_request(manifest=self.adaptive_manifest(bed)))
        self.assertIs(found["target_coverage.bed"].status, CheckStatus.FAILED)

    def test_a_target_bed_on_a_non_canonical_contig_blocks_the_run(self) -> None:
        bed = self.target_bed("chrUn_GL000220v1\t100\t200\tA\n", name="noncanonical")
        found = self.results(self.adaptive_request(manifest=self.adaptive_manifest(bed)))
        self.assertIs(found["target_coverage.bed"].status, CheckStatus.FAILED)

    def test_missing_mosdepth_blocks_an_adaptive_run_rather_than_warning(self) -> None:
        """The stage is optional in the graph; for this assay its failure fails the run."""
        request = self.adaptive_request(
            executables={
                name: str(self.bin / name) for name in DEFAULT_VERSIONS if name != "mosdepth"
            }
            | {"mosdepth": str(self.root / "absent" / "mosdepth")}
        )
        self.assertIs(self.results(request)["tool.mosdepth"].status, CheckStatus.FAILED)

    def test_missing_mosdepth_only_warns_on_a_run_that_never_measures_targets(self) -> None:
        """An lcWGS run records targets as out of scope, so it must not be blocked."""
        request = self.request(
            executables={
                name: str(self.bin / name) for name in DEFAULT_VERSIONS if name != "mosdepth"
            }
            | {"mosdepth": str(self.root / "absent" / "mosdepth")}
        )
        self.assertIs(self.results(request)["tool.mosdepth"].status, CheckStatus.WARNING)

    def test_the_target_checks_are_skipped_rather_than_passed_for_lcwgs(self) -> None:
        """Not applicable and checked-and-fine are different claims."""
        found = self.results()
        self.assertIs(found["target_coverage.policy"].status, CheckStatus.SKIPPED)
        self.assertIs(found["target_coverage.bed"].status, CheckStatus.SKIPPED)
        self.assertIn("lcwgs", found["target_coverage.bed"].detail)

    def test_the_mosdepth_version_lock_is_enforced_for_an_adaptive_run(self) -> None:
        versions = dict(DEFAULT_VERSIONS) | {"mosdepth": "mosdepth 0.3.10"}
        found = self.results(self.adaptive_request(), versions=versions)
        self.assertIs(found["tool.mosdepth"].status, CheckStatus.FAILED)
        self.assertIn("0.3.14", found["tool.mosdepth"].detail)

    def test_the_mosdepth_version_lock_is_not_applied_to_an_lcwgs_run(self) -> None:
        """An lcWGS run never invokes Mosdepth, so its lock must not refuse the run."""
        versions = dict(DEFAULT_VERSIONS) | {"mosdepth": "mosdepth 0.3.10"}
        request = self.request(target_coverage_policy=TARGET_COVERAGE_POLICY)
        self.assertIs(
            self.results(request, versions=versions)["tool.mosdepth"].status, CheckStatus.OK
        )


if __name__ == "__main__":
    unittest.main()

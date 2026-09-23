"""Opt-in real-worker failure tests; all genomic input is generated synthetic data."""

import os
import random
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ontseq_platform.cutesv import run_cutesv
from ontseq_platform.cutesv_build import BUILD_ID, executable_identity, prepare_executable
from ontseq_platform.execution import SubprocessRunner
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisSpec,
    AssaySpec,
    CuteSvPolicy,
    InputSpec,
    SampleManifest,
)

# Change only the worker entry points. Scheduling, result retrieval and VCF output
# remain those of the exact prepared real executable, including multiprocessing.
FAULT_WRAPPERS = """
_original_multi_run_wrapper = multi_run_wrapper
_original_run_del = run_del
def multi_run_wrapper(args):
    if os.environ.get("ONTSEQ_TEST_FAULT") == "extraction" and args[6][0] == "chr2":
        logging.error("SYNTHETIC_WORKER_FAILURE extraction chr2")
        raise RuntimeError("SYNTHETIC extraction failure")
    return _original_multi_run_wrapper(args)
def run_del(args):
    if os.environ.get("ONTSEQ_TEST_FAULT") == "clustering" and args[1] == "chr2":
        from pathlib import Path
        import time
        deadline = time.monotonic() + 10
        while not Path(os.environ["ONTSEQ_TEST_SUCCESS"]).is_file():
            if time.monotonic() >= deadline:
                raise RuntimeError("Fixture never observed the successful worker")
            time.sleep(0.01)
        logging.error("SYNTHETIC_WORKER_FAILURE clustering chr2")
        raise RuntimeError("SYNTHETIC clustering failure")
    result = _original_run_del(args)
    if args[1] == "chr1":
        from pathlib import Path
        Path(os.environ["ONTSEQ_TEST_SUCCESS"]).write_text(str(len(result[1])))
    return result
"""


def synthetic_inputs(root):
    import pysam

    rng = random.Random(739)
    sequence = "".join(rng.choices("ACGT", k=200_000))
    fasta = root / "reference.fa"
    fasta.write_text(f">chr1\n{sequence}\n>chr2\n{sequence}\n")
    pysam.faidx(str(fasta))
    bam = root / "synthetic.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": 200_000} for c in ("chr1", "chr2")],
    }
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for chrom in range(2):
            for i in range(80):
                start = 97_000 + i * 10
                left = 100_000 - start
                deleted = i % 2 == 0
                read = pysam.AlignedSegment()
                read.query_name = f"SYNTHETIC_{chrom}_{i}"
                read.reference_id = chrom
                read.reference_start = start
                read.flag = 0
                read.mapping_quality = 60
                read.cigartuples = (
                    [(0, left), (2, 300), (0, 5000 - left)] if deleted else [(0, 5000)]
                )
                read.query_sequence = (
                    sequence[start:100_000] + sequence[100_300 : 100_300 + 5000 - left]
                    if deleted
                    else sequence[start : start + 5000]
                )
                read.query_qualities = pysam.qualitystring_to_array("I" * 5000)
                out.write(read)
    pysam.index(str(bam))
    manifest = SampleManifest(
        sample_id="SYNTHETIC_WORKERS",
        run_id="SYNTHETIC_RUN",
        input=InputSpec(kind="aligned_bam", path=str(bam), index_path=str(bam) + ".bai"),
        assay=AssaySpec(mode="lcwgs", genome_build="GRCh38", reference_id="SYNTHETIC_REF"),
        analysis=AnalysisSpec(profile="synthetic", modules=["sv"]),
    )
    intake = AlignedBamIntakeReport(
        sample_id=manifest.sample_id,
        reference_id="SYNTHETIC_REF",
        genome_build="GRCh38",
        checks=[],
        verdict="PASS",
    )
    return manifest, intake, fasta


@unittest.skipUnless(os.environ.get("ONTSEQ_CUTESV_REAL_TOOL") == "1", "Opt-in real cuteSV")
class CuteSvWorkerFailureTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.manifest, self.intake, self.fasta = synthetic_inputs(self.root)
        self.tool = shutil.which("cuteSV")
        self.assertIsNotNone(self.tool)
        self.identity = executable_identity(self.tool)
        self.assertEqual(self.identity["cutesv_build_id"], BUILD_ID)
        self.policy = CuteSvPolicy(
            profile_id="synthetic",
            status="technical_defaults_only",
            note="Worker failure qualification",
        )
        self.environment = patch.dict(
            os.environ, {"ONTSEQ_CUTESV_SCRATCH_ROOT": str(self.root / "scratch")}
        )
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def test_successful_workers_preserve_both_known_deletions(self):
        report = run_cutesv(
            self.manifest,
            self.intake,
            self.policy,
            reference_fasta=self.fasta,
            output_vcf=self.root / "control.vcf",
            cutesv=self.tool,
            threads=2,
        )
        self.assertEqual(report.accepted_record_count, 2)
        records = [
            row.split("\t")
            for row in (self.root / "control.vcf").read_text().splitlines()
            if row and not row.startswith("#")
        ]
        self.assertEqual({row[0] for row in records}, {"chr1", "chr2"})
        for row in records:
            info = dict(item.split("=", 1) for item in row[7].split(";") if "=" in item)
            self.assertEqual(row[1], "100000")
            self.assertEqual(info["SVLEN"], "-300")
            self.assertEqual(info["RE"], "40")
        self.assertEqual(executable_identity(self.tool), self.identity)

    def test_legacy_integer_only_build_is_upgraded_and_qualified_copy_is_idempotent(self):
        from ontseq_platform.cutesv_build import FIXED_BODY_SHA256, INTEGER_BODY_SHA256

        stock = Path(self.tool).read_bytes()
        old = b"batch_size=local_ref_len/(int(i[1]/mapped_unit)+1)"
        integer = b"batch_size=max(1, local_ref_len//(int(i[1]/mapped_unit)+1))"
        self.assertEqual(stock.count(old), 1)
        legacy = self.root / "cuteSV-legacy"
        legacy.write_bytes(stock.replace(old, integer))
        legacy.chmod(0o700)
        before = legacy.read_bytes()
        identity = executable_identity(str(legacy))
        self.assertEqual(identity["cutesv_source_body_sha256"], INTEGER_BODY_SHA256)
        prepared = prepare_executable(str(legacy), self.root, identity)
        fixed = executable_identity(prepared)
        self.assertEqual(fixed["cutesv_source_body_sha256"], FIXED_BODY_SHA256)
        self.assertEqual(prepare_executable(prepared, self.root, fixed), prepared)
        for name, tool in (("legacy", str(legacy)), ("prepared", prepared)):
            report = run_cutesv(
                self.manifest,
                self.intake,
                self.policy,
                reference_fasta=self.fasta,
                output_vcf=self.root / f"{name}.vcf",
                cutesv=tool,
                threads=1,
            )
            self.assertEqual(report.accepted_record_count, 2)
            self.assertEqual(
                report.tool.parameters["cutesv_execution_body_sha256"], FIXED_BODY_SHA256
            )
        self.assertEqual(legacy.read_bytes(), before)

    def test_real_worker_failures_abort_without_promoting_partial_output(self):
        prepared = Path(prepare_executable(self.tool, self.root, self.identity))
        source = prepared.read_text()
        marker = "if __name__ == '__main__':"
        self.assertEqual(source.count(marker), 1)
        injected = self.root / "fault-cuteSV"
        injected.write_text(source.replace(marker, FAULT_WRAPPERS + "\n" + marker))
        injected.chmod(0o700)
        for stage in ("extraction", "clustering"):
            for workers in (1, 2):
                with self.subTest(stage=stage, workers=workers):
                    success = self.root / f"success-{stage}-{workers}"
                    output = self.root / f"{stage}-{workers}.vcf"
                    commands = []

                    class RecordingRunner:
                        def run(inner, argv, *, timeout_seconds=300):
                            result = SubprocessRunner().run(argv, timeout_seconds=60)
                            if "--version" not in argv:
                                inner.commands.append(result)
                            return result

                    recording_runner = RecordingRunner()
                    recording_runner.commands = commands
                    with (
                        patch(
                            "ontseq_platform.cutesv.prepare_executable",
                            return_value=str(injected),
                        ),
                        patch.dict(
                            os.environ,
                            {"ONTSEQ_TEST_FAULT": stage, "ONTSEQ_TEST_SUCCESS": str(success)},
                        ),
                        self.assertRaisesRegex(ValueError, "failed with exit code"),
                    ):
                        run_cutesv(
                            self.manifest,
                            self.intake,
                            self.policy,
                            reference_fasta=self.fasta,
                            output_vcf=output,
                            cutesv=self.tool,
                            threads=workers,
                            runner=recording_runner,
                        )
                    self.assertEqual(len(commands), 1)
                    self.assertNotEqual(commands[0].returncode, 0)
                    self.assertIn(f"SYNTHETIC_WORKER_FAILURE {stage}", commands[0].stderr)
                    if stage == "clustering":
                        self.assertEqual(success.read_text(), "1")
                    self.assertFalse(output.exists())
                    self.assertEqual(list((self.root / "scratch").iterdir()), [])
                    self.assertEqual(list(self.root.glob(".cutesv-*")), [])
        self.assertEqual(executable_identity(self.tool), self.identity)

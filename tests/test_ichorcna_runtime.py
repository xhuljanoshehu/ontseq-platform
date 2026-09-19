from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from pydantic import ValidationError

from ontseq_platform.cnv.ichorcna import (
    IchorCnaResourceArtifact,
    IchorCnaResourceBundle,
    IchorCnaResourceRole,
    IchorCnaUlpWgsPolicy,
    build_ichorcna_argv,
    run_ichorcna_ulp_wgs,
)
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import EventType, GenomeBuild, ModuleRunStatus
from ontseq_platform.multicaller_contracts import CallerAssayRegime


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _resources(
    *,
    coverage: Path,
    gc: Path,
    mappability: Path,
    centromere: Path,
    bin_size_bp: int = 1_000_000,
) -> IchorCnaResourceBundle:
    return IchorCnaResourceBundle(
        sample_id="CFDNA_001",
        assay_regime=CallerAssayRegime.CFDNA_ULP_WGS,
        genome_build=GenomeBuild.GRCH38,
        bin_size_bp=bin_size_bp,
        resources=[
            IchorCnaResourceArtifact(
                resource_id="coverage",
                role=IchorCnaResourceRole.COVERAGE_WIG,
                genome_build=GenomeBuild.GRCH38,
                bin_size_bp=bin_size_bp,
                size_bytes=coverage.stat().st_size,
                sha256=_sha(coverage),
            ),
            IchorCnaResourceArtifact(
                resource_id="gc",
                role=IchorCnaResourceRole.GC_WIG,
                genome_build=GenomeBuild.GRCH38,
                bin_size_bp=bin_size_bp,
                size_bytes=gc.stat().st_size,
                sha256=_sha(gc),
            ),
            IchorCnaResourceArtifact(
                resource_id="map",
                role=IchorCnaResourceRole.MAPPABILITY_WIG,
                genome_build=GenomeBuild.GRCH38,
                bin_size_bp=bin_size_bp,
                size_bytes=mappability.stat().st_size,
                sha256=_sha(mappability),
            ),
            IchorCnaResourceArtifact(
                resource_id="centromere",
                role=IchorCnaResourceRole.CENTROMERE,
                genome_build=GenomeBuild.GRCH38,
                bin_size_bp=None,
                size_bytes=centromere.stat().st_size,
                sha256=_sha(centromere),
            ),
        ],
    )


def _policy(script: Path, *, bin_size_bp: int = 1_000_000) -> IchorCnaUlpWgsPolicy:
    return IchorCnaUlpWgsPolicy(
        profile_id="ichorcna-ulp-test",
        expected_version="0.5.1",
        assay_regime=CallerAssayRegime.CFDNA_ULP_WGS,
        genome_build=GenomeBuild.GRCH38,
        genome_build_label="hg38",
        genome_style="UCSC",
        bin_size_bp=bin_size_bp,
        script_sha256=_sha(script),
        ploidy_candidates=[2, 3],
        normal_fraction_candidates=[0.5, 0.6, 0.7, 0.8, 0.9, 0.95],
        max_copy_number=5,
        minimum_mappability=0.75,
        estimate_normal=True,
        estimate_ploidy=True,
        estimate_subclone_prevalence=True,
        include_homozygous_deletion=False,
        transaction_e=0.9999,
        transaction_strength=10_000,
        timeout_seconds=300,
        note="Synthetic ichorCNA ULP-WGS cfDNA contract.",
    )


class FakeIchorRunner:
    def __init__(self, *, version: str = "0.5.1") -> None:
        self.version = version
        self.calls: list[tuple[str, ...]] = []

    def run(self, argv, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        del timeout_seconds
        args = tuple(str(item) for item in argv)
        self.calls.append(args)
        if "-e" in args:
            return CommandResult(args, 0, f"{self.version}\n", "")

        out_dir = Path(args[args.index("--outDir") + 1])
        sample = args[args.index("--id") + 1]
        out_dir.mkdir(parents=True, exist_ok=True)
        _write(
            out_dir / f"{sample}.params.txt",
            "Gender: unknown\n"
            "Tumor Fraction: 0.12\n"
            "Ploidy: 2.4\n"
            "Coverage: 0.18\n"
            "\n"
            "init n_est phi_est BIC Frac_genome_subclonal Frac_CNA_subclonal loglik\n"
            "n0.9-p2 0.88 2.1 120.0 0.05 0.10 4100\n"
            "n0.8-p3 0.79 2.8 125.0 0.08 0.12 4080\n",
        )
        _write(
            out_dir / f"{sample}.seg",
            "ID\tchrom\tstart\tend\tnum.mark\tseg.median.logR\tcopy.number\tcall\t"
            "subclone.status\tlogR_Copy_Number\tCorrected_Copy_Number\tCorrected_Call\n"
            "x1\tchr7\t1000001\t5000000\t4\t-0.35\t1\tHETD\tFALSE\t0.9\t1\tHETD\n"
            "x1\tchr7\t5000001\t9000000\t4\t0.40\t3\tGAIN\tFALSE\t3.1\t3\tGAIN\n",
        )
        _write(
            out_dir / f"{sample}.cna.seg",
            "chr\tstart\tend\tx1.copy.number\tx1.event\nchr7\t1000001\t2000000\t1\tHETD\n",
        )
        _write(
            out_dir / f"{sample}.correctedDepth.txt",
            "chr\tstart\tend\tlog2_TNratio_corrected\nchr7\t1\t1000000\t0.01\n",
        )
        (out_dir / f"{sample}.RData").write_bytes(b"synthetic-rdata")
        return CommandResult(args, 0, "", "")


class IchorCnaRuntimeTests(unittest.TestCase):
    def test_policy_pins_active_upstream_version_0_5_1(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / "runIchorCNA.R"
            _write(script, "# synthetic runIchorCNA.R\n")
            payload = _policy(script).model_dump()
            payload["expected_version"] = "0.5.1"

            upgraded = IchorCnaUlpWgsPolicy.model_validate(payload)

            self.assertEqual(upgraded.expected_version, "0.5.1")
            payload["expected_version"] = "0.2.0"
            with self.assertRaisesRegex(ValidationError, "0.5.1"):
                IchorCnaUlpWgsPolicy.model_validate(payload)

    def test_command_uses_read_depth_resources_and_never_a_bam(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            script = root / "runIchorCNA.R"
            _write(script, "# synthetic runIchorCNA.R\n")
            policy = _policy(script)
            argv = build_ichorcna_argv(
                rscript="Rscript",
                ichorcna_script=script,
                coverage_wig=Path("/tmp/sample.wig"),
                gc_wig=Path("/tmp/gc.wig"),
                mappability_wig=Path("/tmp/map.wig"),
                centromere=Path("/tmp/centromere.txt"),
                sample_id="CFDNA_001",
                output_dir=Path("/tmp/ichor"),
                policy=policy,
            )

        self.assertEqual(argv[:2], ("Rscript", str(script)))
        self.assertIn("--WIG", argv)
        self.assertIn("--gcWig", argv)
        self.assertIn("--mapWig", argv)
        self.assertIn("--centromere", argv)
        self.assertIn("--genomeBuild", argv)
        self.assertNotIn("--bam", argv)
        self.assertNotIn("--target-bam", argv)

    def test_runtime_rejects_adaptive_sampling_before_version_probe(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coverage = root / "sample.wig"
            gc = root / "gc.wig"
            mappability = root / "map.wig"
            centromere = root / "centromere.txt"
            script = root / "runIchorCNA.R"
            for path in (coverage, gc, mappability, centromere, script):
                _write(path, "synthetic\n")
            runner = FakeIchorRunner()

            with self.assertRaisesRegex(ValueError, "cfDNA|ULP-WGS|regime"):
                run_ichorcna_ulp_wgs(
                    coverage_wig=coverage,
                    gc_wig=gc,
                    mappability_wig=mappability,
                    centromere=centromere,
                    ichorcna_script=script,
                    sample_id="CFDNA_001",
                    output_dir=root / "out",
                    resources=_resources(
                        coverage=coverage,
                        gc=gc,
                        mappability=mappability,
                        centromere=centromere,
                    ),
                    policy=_policy(script),
                    assay_regime=CallerAssayRegime.ADAPTIVE_SAMPLING,
                    runner=runner,
                )

        self.assertEqual(runner.calls, [])

    def test_resource_bundle_rejects_bin_size_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coverage = root / "sample.wig"
            gc = root / "gc.wig"
            mappability = root / "map.wig"
            centromere = root / "centromere.txt"
            for path in (coverage, gc, mappability, centromere):
                _write(path, "synthetic\n")

            resources = _resources(
                coverage=coverage,
                gc=gc,
                mappability=mappability,
                centromere=centromere,
            )
            payload = resources.model_dump()
            payload["resources"][1]["bin_size_bp"] = 500_000
            with self.assertRaisesRegex(ValidationError, "bin size"):
                IchorCnaResourceBundle.model_validate(payload)

    def test_successful_run_parses_selected_fit_solutions_segments_and_native_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coverage = root / "sample.wig"
            gc = root / "gc.wig"
            mappability = root / "map.wig"
            centromere = root / "centromere.txt"
            script = root / "runIchorCNA.R"
            _write(coverage, "fixedStep chrom=chr7 start=1 step=1000000 span=1000000\n10\n")
            _write(gc, "fixedStep chrom=chr7 start=1 step=1000000 span=1000000\n0.42\n")
            _write(mappability, "fixedStep chrom=chr7 start=1 step=1000000 span=1000000\n0.98\n")
            _write(centromere, "chr\tstart\tend\nchr7\t60000000\t62000000\n")
            _write(script, "# synthetic runIchorCNA.R\n")
            runner = FakeIchorRunner()

            report = run_ichorcna_ulp_wgs(
                coverage_wig=coverage,
                gc_wig=gc,
                mappability_wig=mappability,
                centromere=centromere,
                ichorcna_script=script,
                sample_id="CFDNA_001",
                output_dir=root / "out",
                resources=_resources(
                    coverage=coverage,
                    gc=gc,
                    mappability=mappability,
                    centromere=centromere,
                ),
                policy=_policy(script),
                assay_regime=CallerAssayRegime.CFDNA_ULP_WGS,
                runner=runner,
            )

        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.tumor_fraction, 0.12)
        self.assertEqual(report.ploidy, 2.4)
        self.assertEqual(report.coverage_x, 0.18)
        self.assertEqual(len(report.solutions), 2)
        self.assertEqual(report.solutions[0].normal_fraction, 0.88)
        self.assertEqual(report.solutions[0].ploidy, 2.1)
        self.assertEqual(len(report.segments), 2)
        self.assertEqual(report.segments[0].start, 1_000_000)
        self.assertEqual(report.segments[0].end, 5_000_000)
        self.assertEqual(report.segments[0].copy_number, 1.0)
        self.assertEqual(report.segments[0].event_type, EventType.DELETION)
        self.assertEqual(report.segments[1].event_type, EventType.DUPLICATION)
        self.assertFalse(report.segments[0].reportable)
        native_roles = {item.role for item in report.native_artifacts}
        self.assertIn("parameters", native_roles)
        self.assertIn("segments", native_roles)
        self.assertIn("bin_level_cna", native_roles)
        self.assertIn("corrected_depth", native_roles)
        self.assertIn("r_workspace", native_roles)
        self.assertTrue(all(item.fingerprint.sha256 for item in report.native_artifacts))

    def test_resource_checksum_mismatch_fails_before_version_probe(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coverage = root / "sample.wig"
            gc = root / "gc.wig"
            mappability = root / "map.wig"
            centromere = root / "centromere.txt"
            script = root / "runIchorCNA.R"
            for path in (coverage, gc, mappability, centromere, script):
                _write(path, "synthetic\n")
            resources = _resources(
                coverage=coverage,
                gc=gc,
                mappability=mappability,
                centromere=centromere,
            )
            _write(gc, "tampered\n")
            runner = FakeIchorRunner()

            with self.assertRaisesRegex(ValueError, "SHA-256"):
                run_ichorcna_ulp_wgs(
                    coverage_wig=coverage,
                    gc_wig=gc,
                    mappability_wig=mappability,
                    centromere=centromere,
                    ichorcna_script=script,
                    sample_id="CFDNA_001",
                    output_dir=root / "out",
                    resources=resources,
                    policy=_policy(script),
                    assay_regime=CallerAssayRegime.CFDNA_ULP_WGS,
                    runner=runner,
                )

        self.assertEqual(runner.calls, [])

    def test_version_drift_fails_before_ichorcna_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coverage = root / "sample.wig"
            gc = root / "gc.wig"
            mappability = root / "map.wig"
            centromere = root / "centromere.txt"
            script = root / "runIchorCNA.R"
            for path in (coverage, gc, mappability, centromere, script):
                _write(path, "synthetic\n")
            runner = FakeIchorRunner(version="0.5.0")

            with self.assertRaisesRegex(ValueError, "version"):
                run_ichorcna_ulp_wgs(
                    coverage_wig=coverage,
                    gc_wig=gc,
                    mappability_wig=mappability,
                    centromere=centromere,
                    ichorcna_script=script,
                    sample_id="CFDNA_001",
                    output_dir=root / "out",
                    resources=_resources(
                        coverage=coverage,
                        gc=gc,
                        mappability=mappability,
                        centromere=centromere,
                    ),
                    policy=_policy(script),
                    assay_regime=CallerAssayRegime.CFDNA_ULP_WGS,
                    runner=runner,
                )

        self.assertEqual(len(runner.calls), 1)
        self.assertIn("-e", runner.calls[0])


if __name__ == "__main__":
    unittest.main()

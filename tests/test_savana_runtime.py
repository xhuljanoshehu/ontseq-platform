from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.execution import CommandResult
from ontseq_platform.models import EventType, GenomeBuild, ModuleRunStatus
from ontseq_platform.multicaller_contracts import CallerInputRole
from ontseq_platform.tumor.savana import (
    SavanaPairedPolicy,
    SavanaTumorOnlyPolicy,
    build_savana_paired_argv,
    build_savana_tumor_only_argv,
    run_savana_paired,
    run_savana_tumor_only,
)
from ontseq_platform.tumor_inputs import (
    TumorAuxiliaryInputArtifact,
    TumorInputBundle,
)


def _sha_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _write_bam(path: Path, payload: bytes) -> None:
    path.write_bytes(payload)
    Path(f"{path}.bai").write_bytes(b"index:" + payload)


def _write_reference(path: Path) -> str:
    path.write_text(">chr7\n" + "A" * 2000 + "\n", encoding="utf-8")
    Path(f"{path}.fai").write_text("chr7\t2000\t6\t2000\t2001\n", encoding="utf-8")
    return _sha_bytes(path.read_bytes())


def _write_snp_vcf(path: Path) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr7\t500\trs1\tA\tG\t60\tPASS\t.\n",
        encoding="utf-8",
    )


def _bundle(
    *,
    tumor_bam: Path,
    snp_vcf: Path,
    reference_sha256: str,
    normal_bam: Path | None = None,
) -> TumorInputBundle:
    artifacts = [
        TumorAuxiliaryInputArtifact(
            artifact_id="tumor-bam",
            role=CallerInputRole.TUMOR_BAM,
            analysis_sample_id="TUMOR_001",
            source_sample_id="TUMOR_001",
            genome_build=GenomeBuild.GRCH38,
            reference_id="GRCh38-test",
            reference_sha256=reference_sha256,
            sha256=_sha_bytes(tumor_bam.read_bytes()),
        ),
        TumorAuxiliaryInputArtifact(
            artifact_id="snp-vcf",
            role=CallerInputRole.SNP_VCF,
            analysis_sample_id="TUMOR_001",
            source_sample_id="NORMAL_001" if normal_bam is not None else "POPULATION_001",
            genome_build=GenomeBuild.GRCH38,
            reference_id="GRCh38-test",
            reference_sha256=reference_sha256,
            sha256=_sha_bytes(snp_vcf.read_bytes()),
        ),
    ]
    if normal_bam is not None:
        artifacts.append(
            TumorAuxiliaryInputArtifact(
                artifact_id="normal-bam",
                role=CallerInputRole.MATCHED_NORMAL_BAM,
                analysis_sample_id="TUMOR_001",
                source_sample_id="NORMAL_001",
                genome_build=GenomeBuild.GRCH38,
                reference_id="GRCh38-test",
                reference_sha256=reference_sha256,
                sha256=_sha_bytes(normal_bam.read_bytes()),
            )
        )
    return TumorInputBundle(
        analysis_sample_id="TUMOR_001",
        genome_build=GenomeBuild.GRCH38,
        reference_id="GRCh38-test",
        reference_sha256=reference_sha256,
        artifacts=artifacts,
    )


def _paired_policy(reference_sha256: str) -> SavanaPairedPolicy:
    return SavanaPairedPolicy(
        profile_id="savana-paired-test",
        expected_version="1.3.8",
        genome_build=GenomeBuild.GRCH38,
        reference_id="GRCh38-test",
        reference_sha256=reference_sha256,
        minimum_sv_length=30,
        minimum_mapping_quality=5,
        minimum_support=3,
        minimum_allele_fraction=0.01,
        cn_bin_size_kbp=10,
        timeout_seconds=300,
        note="Synthetic SAVANA paired SV and CNA adapter contract.",
    )


def _tumor_only_policy(reference_sha256: str) -> SavanaTumorOnlyPolicy:
    return SavanaTumorOnlyPolicy(
        profile_id="savana-tumor-only-test",
        expected_version="1.3.8",
        genome_build=GenomeBuild.GRCH38,
        reference_id="GRCh38-test",
        reference_sha256=reference_sha256,
        minimum_sv_length=30,
        minimum_mapping_quality=5,
        minimum_support=3,
        minimum_allele_fraction=0.01,
        cn_bin_size_kbp=10,
        timeout_seconds=300,
        note="Synthetic SAVANA tumor-only SV and CNA adapter contract.",
    )


class FakeSavanaRunner:
    def __init__(
        self,
        *,
        version: str = "1.3.8",
        empty_sv: bool = False,
        malformed_cna: bool = False,
        ranked_solutions_text: str | None = None,
        fitted_purity_ploidy_text: str | None = None,
    ) -> None:
        self.version = version
        self.empty_sv = empty_sv
        self.malformed_cna = malformed_cna
        self.ranked_solutions_text = ranked_solutions_text
        self.fitted_purity_ploidy_text = fitted_purity_ploidy_text
        self.calls: list[tuple[str, ...]] = []

    def run(self, argv, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        del timeout_seconds
        args = tuple(str(item) for item in argv)
        self.calls.append(args)
        if "--version" in args:
            return CommandResult(args, 0, f"SAVANA {self.version}\n", "")

        outdir = Path(args[args.index("--outdir") + 1])
        sample = args[args.index("--sample") + 1]
        outdir.mkdir(parents=True, exist_ok=True)
        records = (
            ""
            if self.empty_sv
            else (
                "chr7\t1001\tsavana_1\tN\t<INS>\t60\tPASS\t"
                "SVTYPE=INS;SVLEN=100;CLASS=SOMATIC;"
                "TUMOUR_READ_SUPPORT=5;TUMOUR_AF=0.25\n"
            )
        )
        (outdir / f"{sample}.classified.somatic.vcf").write_text(
            "##fileformat=VCFv4.2\n"
            "##source=SAVANA\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + records,
            encoding="utf-8",
        )
        (outdir / f"{sample}.classified.somatic.bedpe").write_text(
            "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\n",
            encoding="utf-8",
        )
        (outdir / f"{sample}_sv_breakpoints_read_support.tsv").write_text(
            "variant_id\ttumour_read_ids\tnormal_read_ids\nsavana_1\treadA,readB\t\n",
            encoding="utf-8",
        )
        (outdir / f"{sample}_ranked_solutions.tsv").write_text(
            self.ranked_solutions_text
            if self.ranked_solutions_text is not None
            else "purity\tploidy\tdistance\trank\n0.45\t2.1\t0.12\t1\n0.35\t2.6\t0.14\t2\n",
            encoding="utf-8",
        )
        (outdir / f"{sample}_fitted_purity_ploidy.tsv").write_text(
            self.fitted_purity_ploidy_text
            if self.fitted_purity_ploidy_text is not None
            else "purity\tploidy\tdistance\trank\n0.45\t2.1\t0.12\t1\n",
            encoding="utf-8",
        )
        copy_number = "not-a-number" if self.malformed_cna else "1.2"
        (outdir / f"{sample}_segmented_absolute_copy_number.tsv").write_text(
            "chromosome\tstart\tend\tsegment_id\tbin_count\t"
            "sum_of_bin_lengths\tweight\tcopyNumber\tminorAlleleCopyNumber\t"
            "meanBAF\tno_hetSNPs\n"
            f"chr7\t0\t1000000\tseg1\t100\t1000000\t1.0\t{copy_number}\t0.4\t0.45\t20\n",
            encoding="utf-8",
        )
        (outdir / f"{sample}_10_raw_read_counts.tsv").write_text(
            "chromosome\tstart\tend\ttumour\tnormal\nchr7\t0\t10000\t10\t8\n",
            encoding="utf-8",
        )
        return CommandResult(args, 0, "", "")


class SavanaRuntimeTests(unittest.TestCase):
    def test_paired_and_tumor_only_commands_are_explicitly_distinct(self) -> None:
        paired = build_savana_paired_argv(
            savana="savana",
            tumor_bam=Path("/tmp/tumor.bam"),
            normal_bam=Path("/tmp/normal.bam"),
            snp_vcf=Path("/tmp/snps.vcf"),
            reference_fasta=Path("/tmp/ref.fa"),
            sample_id="TUMOR_001",
            output_dir=Path("/tmp/savana-paired"),
            policy=_paired_policy("a" * 64),
            threads=4,
        )
        tumor_only = build_savana_tumor_only_argv(
            savana="savana",
            tumor_bam=Path("/tmp/tumor.bam"),
            snp_vcf=Path("/tmp/snps.vcf"),
            reference_fasta=Path("/tmp/ref.fa"),
            sample_id="TUMOR_001",
            output_dir=Path("/tmp/savana-to"),
            policy=_tumor_only_policy("a" * 64),
            threads=4,
        )

        self.assertEqual(paired[0], "savana")
        self.assertNotEqual(paired[1], "to")
        self.assertIn("--normal", paired)
        self.assertIn("--snp_vcf", paired)
        self.assertEqual(tumor_only[:2], ("savana", "to"))
        self.assertNotIn("--normal", tumor_only)
        self.assertIn("--snp_vcf", tumor_only)

    def test_paired_runtime_never_falls_back_when_matched_normal_is_missing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            snps = root / "snps.vcf"
            ref = root / "ref.fa"
            _write_bam(tumor, b"tumor")
            _write_snp_vcf(snps)
            ref_sha = _write_reference(ref)

            with self.assertRaisesRegex(ValueError, "matched normal"):
                run_savana_paired(
                    tumor_bam=tumor,
                    normal_bam=root / "normal.bam",
                    snp_vcf=snps,
                    reference_fasta=ref,
                    sample_id="TUMOR_001",
                    normal_sample_id="NORMAL_001",
                    output_dir=root / "paired",
                    inputs=_bundle(
                        tumor_bam=tumor,
                        snp_vcf=snps,
                        reference_sha256=ref_sha,
                    ),
                    policy=_paired_policy(ref_sha),
                    runner=FakeSavanaRunner(),
                )

    def test_paired_runtime_retains_sv_cna_fits_and_sensitive_native_support(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            snps = root / "snps.vcf"
            ref = root / "ref.fa"
            _write_bam(tumor, b"tumor")
            _write_bam(normal, b"normal")
            _write_snp_vcf(snps)
            ref_sha = _write_reference(ref)

            report = run_savana_paired(
                tumor_bam=tumor,
                normal_bam=normal,
                snp_vcf=snps,
                reference_fasta=ref,
                sample_id="TUMOR_001",
                normal_sample_id="NORMAL_001",
                output_dir=root / "paired",
                inputs=_bundle(
                    tumor_bam=tumor,
                    normal_bam=normal,
                    snp_vcf=snps,
                    reference_sha256=ref_sha,
                ),
                policy=_paired_policy(ref_sha),
                runner=FakeSavanaRunner(),
                threads=2,
            )

        self.assertEqual(report.mode, "paired")
        self.assertEqual(report.sv_status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.cna_status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertEqual(len(report.events), 1)
        self.assertEqual(report.events[0].event_type, EventType.INSERTION)
        self.assertFalse(report.events[0].reportable)
        self.assertEqual(report.selected_fit.purity, 0.45)
        self.assertEqual(report.selected_fit.ploidy, 2.1)
        self.assertEqual([item.rank for item in report.ranked_fits], [1, 2])
        self.assertEqual(report.copy_number_segments[0].copy_number, 1.2)
        self.assertEqual(report.copy_number_segments[0].minor_allele_copy_number, 0.4)

        native = {item.relative_path: item for item in report.native_artifacts}
        support = native["TUMOR_001_sv_breakpoints_read_support.tsv"]
        self.assertTrue(support.sensitive_output)
        self.assertFalse(support.exportable)
        self.assertNotIn("readA", report.model_dump_json())

    def test_malformed_existing_cna_output_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            snps = root / "snps.vcf"
            ref = root / "ref.fa"
            _write_bam(tumor, b"tumor")
            _write_bam(normal, b"normal")
            _write_snp_vcf(snps)
            ref_sha = _write_reference(ref)

            with self.assertRaisesRegex(ValueError, "copyNumber"):
                run_savana_paired(
                    tumor_bam=tumor,
                    normal_bam=normal,
                    snp_vcf=snps,
                    reference_fasta=ref,
                    sample_id="TUMOR_001",
                    normal_sample_id="NORMAL_001",
                    output_dir=root / "paired",
                    inputs=_bundle(
                        tumor_bam=tumor,
                        normal_bam=normal,
                        snp_vcf=snps,
                        reference_sha256=ref_sha,
                    ),
                    policy=_paired_policy(ref_sha),
                    runner=FakeSavanaRunner(malformed_cna=True),
                )

    def test_selected_fit_disagreeing_with_rank_one_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            snps = root / "snps.vcf"
            ref = root / "ref.fa"
            _write_bam(tumor, b"tumor")
            _write_bam(normal, b"normal")
            _write_snp_vcf(snps)
            ref_sha = _write_reference(ref)

            with self.assertRaisesRegex(ValueError, "rank-1"):
                run_savana_paired(
                    tumor_bam=tumor,
                    normal_bam=normal,
                    snp_vcf=snps,
                    reference_fasta=ref,
                    sample_id="TUMOR_001",
                    normal_sample_id="NORMAL_001",
                    output_dir=root / "paired",
                    inputs=_bundle(
                        tumor_bam=tumor,
                        normal_bam=normal,
                        snp_vcf=snps,
                        reference_sha256=ref_sha,
                    ),
                    policy=_paired_policy(ref_sha),
                    runner=FakeSavanaRunner(
                        ranked_solutions_text=(
                            "purity\tploidy\tdistance\trank\n"
                            "0.45\t2.1\t0.12\t1\n"
                            "0.35\t2.6\t0.14\t2\n"
                        ),
                        fitted_purity_ploidy_text=(
                            "purity\tploidy\tdistance\trank\n0.35\t2.6\t0.14\t2\n"
                        ),
                    ),
                )

    def test_duplicate_rank_one_solutions_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            snps = root / "snps.vcf"
            ref = root / "ref.fa"
            _write_bam(tumor, b"tumor")
            _write_bam(normal, b"normal")
            _write_snp_vcf(snps)
            ref_sha = _write_reference(ref)

            with self.assertRaisesRegex(ValueError, "duplicate ranks"):
                run_savana_paired(
                    tumor_bam=tumor,
                    normal_bam=normal,
                    snp_vcf=snps,
                    reference_fasta=ref,
                    sample_id="TUMOR_001",
                    normal_sample_id="NORMAL_001",
                    output_dir=root / "paired",
                    inputs=_bundle(
                        tumor_bam=tumor,
                        normal_bam=normal,
                        snp_vcf=snps,
                        reference_sha256=ref_sha,
                    ),
                    policy=_paired_policy(ref_sha),
                    runner=FakeSavanaRunner(
                        ranked_solutions_text=(
                            "purity\tploidy\tdistance\trank\n"
                            "0.45\t2.1\t0.12\t1\n"
                            "0.35\t2.6\t0.14\t1\n"
                        ),
                        fitted_purity_ploidy_text=(
                            "purity\tploidy\tdistance\trank\n0.45\t2.1\t0.12\t1\n"
                        ),
                    ),
                )

    def test_ranked_solutions_missing_rank_one_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            snps = root / "snps.vcf"
            ref = root / "ref.fa"
            _write_bam(tumor, b"tumor")
            _write_bam(normal, b"normal")
            _write_snp_vcf(snps)
            ref_sha = _write_reference(ref)

            with self.assertRaisesRegex(ValueError, "missing rank 1"):
                run_savana_paired(
                    tumor_bam=tumor,
                    normal_bam=normal,
                    snp_vcf=snps,
                    reference_fasta=ref,
                    sample_id="TUMOR_001",
                    normal_sample_id="NORMAL_001",
                    output_dir=root / "paired",
                    inputs=_bundle(
                        tumor_bam=tumor,
                        normal_bam=normal,
                        snp_vcf=snps,
                        reference_sha256=ref_sha,
                    ),
                    policy=_paired_policy(ref_sha),
                    runner=FakeSavanaRunner(
                        ranked_solutions_text=(
                            "purity\tploidy\tdistance\trank\n"
                            "0.45\t2.1\t0.12\t2\n"
                            "0.35\t2.6\t0.14\t3\n"
                        ),
                        fitted_purity_ploidy_text=(
                            "purity\tploidy\tdistance\trank\n0.45\t2.1\t0.12\t2\n"
                        ),
                    ),
                )

    def test_tumor_only_runtime_stays_distinct_and_does_not_use_normal(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            snps = root / "snps.vcf"
            ref = root / "ref.fa"
            _write_bam(tumor, b"tumor")
            _write_snp_vcf(snps)
            ref_sha = _write_reference(ref)
            runner = FakeSavanaRunner()

            report = run_savana_tumor_only(
                tumor_bam=tumor,
                snp_vcf=snps,
                reference_fasta=ref,
                sample_id="TUMOR_001",
                output_dir=root / "tumor-only",
                inputs=_bundle(
                    tumor_bam=tumor,
                    snp_vcf=snps,
                    reference_sha256=ref_sha,
                ),
                policy=_tumor_only_policy(ref_sha),
                runner=runner,
            )

        self.assertEqual(report.mode, "tumor_only")
        self.assertEqual(report.policy.mode, "tumor_only")
        run_argv = runner.calls[-1]
        self.assertEqual(run_argv[:2], ("savana", "to"))
        self.assertNotIn("--normal", run_argv)

    def test_empty_sv_callset_keeps_successful_cna_separate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            snps = root / "snps.vcf"
            ref = root / "ref.fa"
            _write_bam(tumor, b"tumor")
            _write_bam(normal, b"normal")
            _write_snp_vcf(snps)
            ref_sha = _write_reference(ref)

            report = run_savana_paired(
                tumor_bam=tumor,
                normal_bam=normal,
                snp_vcf=snps,
                reference_fasta=ref,
                sample_id="TUMOR_001",
                normal_sample_id="NORMAL_001",
                output_dir=root / "paired",
                inputs=_bundle(
                    tumor_bam=tumor,
                    normal_bam=normal,
                    snp_vcf=snps,
                    reference_sha256=ref_sha,
                ),
                policy=_paired_policy(ref_sha),
                runner=FakeSavanaRunner(empty_sv=True),
            )

        self.assertEqual(report.sv_status, ModuleRunStatus.NO_CALL)
        self.assertEqual(report.cna_status, ModuleRunStatus.COMPLETED)
        self.assertEqual(report.status, ModuleRunStatus.COMPLETED)
        self.assertTrue(any("not a biological negative" in item for item in report.warnings))

    def test_version_drift_fails_before_savana_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tumor = root / "tumor.bam"
            normal = root / "normal.bam"
            snps = root / "snps.vcf"
            ref = root / "ref.fa"
            _write_bam(tumor, b"tumor")
            _write_bam(normal, b"normal")
            _write_snp_vcf(snps)
            ref_sha = _write_reference(ref)
            runner = FakeSavanaRunner(version="1.3.7")

            with self.assertRaisesRegex(ValueError, "version"):
                run_savana_paired(
                    tumor_bam=tumor,
                    normal_bam=normal,
                    snp_vcf=snps,
                    reference_fasta=ref,
                    sample_id="TUMOR_001",
                    normal_sample_id="NORMAL_001",
                    output_dir=root / "paired",
                    inputs=_bundle(
                        tumor_bam=tumor,
                        normal_bam=normal,
                        snp_vcf=snps,
                        reference_sha256=ref_sha,
                    ),
                    policy=_paired_policy(ref_sha),
                    runner=runner,
                )

        self.assertEqual(len(runner.calls), 1)
        self.assertIn("--version", runner.calls[0])


if __name__ == "__main__":
    unittest.main()

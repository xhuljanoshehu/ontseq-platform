from __future__ import annotations

import re
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

from pydantic import Field

from .cutesv import CuteSVCallReport, CuteSVPolicy, normalize_cutesv_vcf
from .execution import CommandRunner, SubprocessRunner
from .models import (
    AlignedBamIntakeReport,
    InputKind,
    SampleManifest,
    StrictModel,
    ToolRecord,
    Verdict,
)
from .reference import contig_signature, reference_lock_from_fai

_VERSION = re.compile(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)")


class CuteSVExecutionPolicy(StrictModel):
    """Local cuteSV execution settings; these are not clinical thresholds."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    expected_version: str = Field(default="2.1.4", pattern=r"^\d+\.\d+\.\d+$")
    min_support: int = Field(default=5, ge=1)
    min_sv_length: int = Field(default=50, ge=1)
    max_size: Literal[-1] = -1
    min_mapq: int = Field(default=20, ge=0, le=60)
    max_cluster_bias_ins: int = Field(default=100, ge=0)
    diff_ratio_merging_ins: float = Field(default=0.3, ge=0, le=1)
    max_cluster_bias_del: int = Field(default=100, ge=0)
    diff_ratio_merging_del: float = Field(default=0.3, ge=0, le=1)
    max_cluster_bias_tra: int = Field(default=50, ge=0)
    diff_ratio_filtering_tra: float = Field(default=0.6, ge=0, le=1)
    report_read_ids: Literal[False] = False
    ignore_sequence: Literal[True] = True
    retain_work_dir: Literal[False] = False
    write_old_sigs: Literal[False] = False
    genotype: Literal[False] = False
    platform_profile: Literal["ont"] = "ont"
    clinically_validated: Literal[False] = False
    note: str = (
        "cuteSV execution parameters are software-engineering defaults for ONT comparison only; "
        "they are not assay-specific analytical or clinical thresholds."
    )

    def normalization_policy(self) -> CuteSVPolicy:
        return CuteSVPolicy(
            expected_version=self.expected_version,
            min_support=self.min_support,
            min_sv_length=self.min_sv_length,
            report_read_ids=False,
        )


def _cutesv_version(text: str) -> str:
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    first_line = text.splitlines()[0].strip() if text.strip() else "unknown"
    return first_line[:80]


def _reference_fai_path(reference_fasta: Path, reference_fai: Path | None) -> Path:
    return reference_fai if reference_fai is not None else Path(f"{reference_fasta}.fai")


def _validate_reference_compatibility(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    *,
    reference_fasta: Path,
    reference_fai: Path,
) -> tuple[str, str]:
    if not reference_fasta.is_file():
        raise ValueError("cuteSV reference FASTA is missing or unreadable")
    if not reference_fai.is_file():
        raise ValueError("cuteSV reference FASTA index is missing or unreadable")
    if intake.header is None:
        raise ValueError("cuteSV requires BAM header metadata for reference compatibility checking")

    reference_lock = reference_lock_from_fai(
        reference_fai,
        reference_id=manifest.assay.reference_id,
        genome_build=manifest.assay.genome_build,
    )
    reference_signature = contig_signature(
        (contig.name, contig.length) for contig in reference_lock.contigs
    )
    if reference_signature != intake.header.contig_signature_sha256:
        raise ValueError("cuteSV reference contig signature does not match the aligned BAM header")
    return reference_lock.source_fai_sha256, reference_signature


def run_cutesv(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    policy: CuteSVExecutionPolicy,
    *,
    reference_fasta: Path,
    output_vcf: Path,
    reference_fai: Path | None = None,
    runner: CommandRunner | None = None,
    cutesv: str = "cuteSV",
    threads: int = 4,
    timeout_seconds: int = 7200,
) -> CuteSVCallReport:
    """Run cuteSV locally without a shell and normalize only privacy-safe candidate evidence."""

    if manifest.input.kind != InputKind.ALIGNED_BAM:
        raise ValueError("cuteSV requires input.kind=aligned_bam")
    if manifest.sample_id != intake.sample_id:
        raise ValueError("Manifest and intake artifact must refer to the same sample")
    if manifest.assay.genome_build != intake.genome_build:
        raise ValueError("Manifest and intake artifact use different genome builds")
    if manifest.assay.reference_id != intake.reference_id:
        raise ValueError("Manifest and intake artifact use different reference IDs")
    if intake.verdict == Verdict.FAIL:
        raise ValueError("cuteSV cannot run after a failed aligned-BAM intake gate")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if timeout_seconds < 1:
        raise ValueError("timeout_seconds must be at least 1")
    if output_vcf.exists():
        raise ValueError("Refusing to overwrite an existing cuteSV VCF")

    resolved_fai = _reference_fai_path(reference_fasta, reference_fai)
    fai_sha256, reference_signature = _validate_reference_compatibility(
        manifest,
        intake,
        reference_fasta=reference_fasta,
        reference_fai=resolved_fai,
    )

    command_runner = runner or SubprocessRunner()
    version_result = command_runner.run([cutesv, "--version"], timeout_seconds=30)
    if version_result.returncode != 0:
        raise ValueError("cuteSV version probe returned a non-zero exit code")
    version = _cutesv_version(f"{version_result.stdout}\n{version_result.stderr}")
    if version != policy.expected_version:
        raise ValueError(
            f"cuteSV version {version!r} does not match policy lock {policy.expected_version!r}"
        )

    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    parameters: dict[str, str | int | float | bool] = {
        "threads": threads,
        "min_support": policy.min_support,
        "min_size": policy.min_sv_length,
        "max_size": policy.max_size,
        "min_mapq": policy.min_mapq,
        "max_cluster_bias_INS": policy.max_cluster_bias_ins,
        "diff_ratio_merging_INS": policy.diff_ratio_merging_ins,
        "max_cluster_bias_DEL": policy.max_cluster_bias_del,
        "diff_ratio_merging_DEL": policy.diff_ratio_merging_del,
        "max_cluster_bias_TRA": policy.max_cluster_bias_tra,
        "diff_ratio_filtering_TRA": policy.diff_ratio_filtering_tra,
        "report_read_ids": False,
        "ignore_sequence": True,
        "retain_work_dir": False,
        "write_old_sigs": False,
        "genotype": False,
        "platform_profile": policy.platform_profile,
        "clinically_validated": False,
        "expected_version": policy.expected_version,
        "reference_fai_sha256": fai_sha256,
        "reference_contig_signature_sha256": reference_signature,
    }

    with TemporaryDirectory(prefix="ontseq-cutesv-", dir=output_vcf.parent) as work_dir:
        argv = [
            cutesv,
            manifest.input.path,
            str(reference_fasta),
            str(output_vcf),
            work_dir,
            "--threads",
            str(threads),
            "--sample",
            manifest.sample_id,
            "--min_mapq",
            str(policy.min_mapq),
            "--min_support",
            str(policy.min_support),
            "--min_size",
            str(policy.min_sv_length),
            "--max_size",
            str(policy.max_size),
            "--max_cluster_bias_INS",
            str(policy.max_cluster_bias_ins),
            "--diff_ratio_merging_INS",
            str(policy.diff_ratio_merging_ins),
            "--max_cluster_bias_DEL",
            str(policy.max_cluster_bias_del),
            "--diff_ratio_merging_DEL",
            str(policy.diff_ratio_merging_del),
            "--max_cluster_bias_TRA",
            str(policy.max_cluster_bias_tra),
            "--diff_ratio_filtering_TRA",
            str(policy.diff_ratio_filtering_tra),
            "--ignore_sequence",
        ]
        result = command_runner.run(argv, timeout_seconds=timeout_seconds)
        if result.returncode != 0:
            output_vcf.unlink(missing_ok=True)
            raise ValueError("cuteSV returned a non-zero exit code")

    if not output_vcf.is_file():
        raise ValueError("cuteSV completed without producing the expected VCF")

    tool = ToolRecord(name="cuteSV", version=version, parameters=parameters)
    return normalize_cutesv_vcf(
        output_vcf,
        sample_id=manifest.sample_id,
        genome_build=manifest.assay.genome_build,
        policy=policy.normalization_policy(),
        tool=tool,
    )

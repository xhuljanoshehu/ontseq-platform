from __future__ import annotations

from pathlib import Path

from .bam_intake import AlignedBamInspector
from .execution import CommandRunner, SubprocessRunner
from .io import write_json
from .models import (
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CheckStatus,
    GenomeBuild,
    InputKind,
    InputSpec,
    LocalSmokeReport,
    QCPolicy,
    SampleManifest,
    SnifflesPolicy,
    ValidationCheck,
    Verdict,
)
from .mvp import assemble_aligned_bam_mvp
from .qc import run_cramino_qc
from .reference import reference_lock_from_fai
from .report import render_html
from .sniffles import run_sniffles
from .workbook import render_workbook

SYNTHETIC_SAMPLE_ID = "SYNTHETIC_SMOKE_001"
SYNTHETIC_REFERENCE_ID = "SYNTHETIC_NOT_REAL_GRCH38_V1"


def _sam_record(
    name: str,
    *,
    position: int,
    cigar: str,
    query_length: int,
    edit_distance: int,
) -> str:
    sequence = "A" * query_length
    quality = "I" * query_length
    fields = [
        name,
        "0",
        "chr1",
        str(position),
        "60",
        cigar,
        "*",
        "0",
        "0",
        sequence,
        quality,
        f"NM:i:{edit_distance}",
        "RG:Z:SYNTHETIC_RG",
    ]
    return "\t".join(fields)


def synthetic_sam_text() -> str:
    """Return an identifier-free long-read fixture with a supported 200 bp deletion."""
    lines = [
        "@HD\tVN:1.6\tSO:unsorted",
        "@SQ\tSN:chr1\tLN:100000",
        "@SQ\tSN:chr2\tLN:100000",
        "@RG\tID:SYNTHETIC_RG\tSM:SYNTHETIC_SMOKE_001\tPL:ONT",
        "@PG\tID:ontseq-smoke\tPN:ontseq-smoke\tVN:0.1",
    ]
    for index in range(12):
        lines.append(
            _sam_record(
                f"SYNTH_DEL_{index + 1:03d}",
                position=5001 + (index % 3),
                cigar="5000M200D5000M",
                query_length=10000,
                edit_distance=200,
            )
        )
    for index in range(12):
        lines.append(
            _sam_record(
                f"SYNTH_REF_{index + 1:03d}",
                position=5001 + (index % 3),
                cigar="10200M",
                query_length=10200,
                edit_distance=0,
            )
        )
    return "\n".join(lines) + "\n"


def _run_checked(
    runner: CommandRunner,
    argv: list[str],
    *,
    label: str,
    timeout_seconds: int,
) -> None:
    result = runner.run(argv, timeout_seconds=timeout_seconds)
    if result.returncode != 0:
        raise ValueError(f"{label} failed with exit code {result.returncode}")


def _assert_targets_absent(paths: list[Path]) -> None:
    existing = [path.name for path in paths if path.exists()]
    if existing:
        raise ValueError(
            "Refusing to overwrite existing local-smoke artifacts: " + ", ".join(existing)
        )


def run_local_smoke(
    output_dir: Path,
    qc_policy: QCPolicy,
    sniffles_policy: SnifflesPolicy,
    *,
    runner: CommandRunner | None = None,
    samtools: str = "samtools",
    cramino: str = "cramino",
    sniffles: str = "sniffles",
    threads: int = 2,
    pipeline_version: str = "UNKNOWN",
    git_commit: str = "LOCAL_SMOKE",
) -> LocalSmokeReport:
    if threads < 1:
        raise ValueError("threads must be at least 1")
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    sam_path = output_dir / "synthetic.input.sam"
    unsorted_bam = output_dir / "synthetic.unsorted.bam"
    bam_path = output_dir / "synthetic.sorted.bam"
    bai_path = output_dir / "synthetic.sorted.bam.bai"
    fai_path = output_dir / "synthetic.reference.fai"
    manifest_path = output_dir / "synthetic.manifest.json"
    reference_lock_path = output_dir / "synthetic.reference-lock.json"
    intake_path = output_dir / "synthetic.intake.json"
    qc_path = output_dir / "synthetic.qc.json"
    vcf_path = output_dir / "synthetic.sniffles.vcf"
    sniffles_report_path = output_dir / "synthetic.sniffles.json"
    report_path = output_dir / "local-smoke.report.json"
    result_path = output_dir / f"{SYNTHETIC_SAMPLE_ID}.result.json"
    html_path = output_dir / f"{SYNTHETIC_SAMPLE_ID}.report.html"
    workbook_path = output_dir / f"{SYNTHETIC_SAMPLE_ID}.results.xlsx"
    targets = [
        sam_path,
        unsorted_bam,
        bam_path,
        bai_path,
        fai_path,
        manifest_path,
        reference_lock_path,
        intake_path,
        qc_path,
        vcf_path,
        sniffles_report_path,
        report_path,
        result_path,
        html_path,
        workbook_path,
    ]
    _assert_targets_absent(targets)

    command_runner = runner or SubprocessRunner()
    sam_path.write_text(synthetic_sam_text(), encoding="utf-8")
    fai_path.write_text(
        "chr1\t100000\t0\t0\t0\nchr2\t100000\t0\t0\t0\n",
        encoding="utf-8",
    )
    _run_checked(
        command_runner,
        [samtools, "view", "-b", "-o", str(unsorted_bam), str(sam_path)],
        label="samtools BAM conversion",
        timeout_seconds=120,
    )
    _run_checked(
        command_runner,
        [samtools, "sort", "-@", str(threads), "-o", str(bam_path), str(unsorted_bam)],
        label="samtools coordinate sort",
        timeout_seconds=300,
    )
    _run_checked(
        command_runner,
        [samtools, "index", "-@", str(threads), str(bam_path), str(bai_path)],
        label="samtools index",
        timeout_seconds=300,
    )
    sam_path.unlink()
    unsorted_bam.unlink()

    manifest = SampleManifest(
        sample_id=SYNTHETIC_SAMPLE_ID,
        run_id="SYNTHETIC_SMOKE_RUN_001",
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path=str(bam_path),
            index_path=str(bai_path),
        ),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id=SYNTHETIC_REFERENCE_ID,
        ),
        analysis=AnalysisSpec(
            profile="synthetic-local-tool-smoke",
            modules=[AnalysisModule.QC, AnalysisModule.SV, AnalysisModule.REPORT],
        ),
    )
    reference_lock = reference_lock_from_fai(
        fai_path,
        reference_id=SYNTHETIC_REFERENCE_ID,
        genome_build=GenomeBuild.GRCH38,
    )
    write_json(manifest, manifest_path)
    write_json(reference_lock, reference_lock_path)

    intake = AlignedBamInspector(runner=command_runner, samtools=samtools).inspect(
        manifest,
        reference_lock,
        include_checksums=True,
    )
    write_json(intake, intake_path)
    if intake.verdict != Verdict.PASS:
        raise ValueError("Synthetic BAM did not pass every aligned-BAM intake check")

    qc_report = run_cramino_qc(
        manifest,
        qc_policy,
        runner=command_runner,
        cramino=cramino,
        threads=threads,
    )
    write_json(qc_report, qc_path)
    if qc_report.qc.verdict == Verdict.FAIL:
        raise ValueError("Synthetic BAM failed the configured Cramino QC gate")

    sniffles_report = run_sniffles(
        manifest,
        intake,
        sniffles_policy,
        output_vcf=vcf_path,
        runner=command_runner,
        sniffles=sniffles,
        threads=threads,
    )
    write_json(sniffles_report, sniffles_report_path)

    expected_deletions = [
        event
        for event in sniffles_report.events
        if event.event_type.value == "deletion"
        and event.primary.chromosome == "chr1"
        and event.primary.start < 10200
        and event.primary.end > 10000
        and event.length_bp is not None
        and 150 <= event.length_bp <= 250
        and event.evidence
        and (event.evidence[0].support_reads or 0) >= sniffles_policy.min_support
    ]
    deletion_detected = bool(expected_deletions)
    checks = [
        ValidationCheck(
            name="aligned_bam_intake",
            status=CheckStatus.PASS,
            message="Synthetic BAM passed the real samtools intake gate.",
        ),
        ValidationCheck(
            name="cramino_execution",
            status=CheckStatus.PASS,
            message="Cramino executed and returned a normalized QC artifact.",
        ),
        ValidationCheck(
            name="expected_synthetic_deletion",
            status=CheckStatus.PASS if deletion_detected else CheckStatus.FAIL,
            message=(
                "Sniffles2 recovered the expected synthetic deletion candidate."
                if deletion_detected
                else "Sniffles2 did not recover the expected synthetic deletion candidate."
            ),
            details={
                "accepted_sv_candidates": sniffles_report.accepted_record_count,
                "matching_expected_deletions": len(expected_deletions),
            },
        ),
        ValidationCheck(
            name="privacy_preserving_sniffles_mode",
            status=(
                CheckStatus.PASS
                if sniffles_report.tool.parameters.get("output_read_names") is False
                else CheckStatus.FAIL
            ),
            message="Read-name export is disabled in the Sniffles2 adapter.",
        ),
    ]
    verdict = (
        Verdict.FAIL if any(item.status == CheckStatus.FAIL for item in checks) else Verdict.PASS
    )
    report = LocalSmokeReport(
        sample_id=SYNTHETIC_SAMPLE_ID,
        verdict=verdict,
        intake=intake,
        qc=qc_report,
        sniffles=sniffles_report,
        checks=checks,
        limitations=[
            "This is a deterministic engineering smoke test using fully synthetic alignments.",
            "The synthetic reference is not GRCh38 despite exercising the GRCh38 namespace "
            "contract.",
            "A passing smoke test proves execution and normalization, not clinical performance.",
        ],
    )
    write_json(report, report_path)
    if report.verdict == Verdict.FAIL:
        failures = [
            f"{item.name}: {item.message} ({item.details})"
            for item in report.checks
            if item.status == CheckStatus.FAIL
        ]
        raise ValueError("Local real-tool smoke test failed: " + "; ".join(failures))
    pipeline_result = assemble_aligned_bam_mvp(
        manifest,
        intake,
        qc_report,
        pipeline_version=pipeline_version,
        git_commit=git_commit,
        sniffles_report=sniffles_report,
    )
    write_json(pipeline_result, result_path)
    render_html(pipeline_result, html_path)
    render_workbook(pipeline_result, workbook_path)
    return report

from __future__ import annotations

import argparse
from pathlib import Path

from pydantic import ValidationError

from . import __version__
from .bam_intake import AlignedBamInspector
from .benchmark import benchmark_case
from .demo import build_demo_result
from .dilution import (
    DilutionPolicy,
    DilutionSeriesPlan,
    LodPolicy,
    count_reads,
    evaluate_lod,
    execute_dilution_series,
    plan_dilution_series,
)
from .execution import SubprocessRunner, ToolExecutionError
from .io import load_model, write_json
from .methylation import MethylationPolicy, run_methylation
from .models import (
    AlignedBamIntakeReport,
    BenchmarkCase,
    BenchmarkReport,
    CraminoQCReport,
    GenomeBuild,
    PipelineResult,
    QCPolicy,
    ReferenceLock,
    SampleManifest,
    SnifflesCallReport,
    SnifflesPolicy,
    Verdict,
)
from .mvp import assemble_aligned_bam_mvp
from .qc import run_cramino_qc
from .reference import reference_lock_from_fai, validate_canonical_reference
from .report import render_html
from .smoke import run_local_smoke
from .sniffles import run_sniffles
from .target_coverage import TargetCoveragePolicy, run_target_coverage
from .workbook import render_workbook


def _render(result: PipelineResult, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = result.manifest.sample_id
    outputs = [
        write_json(result, output_dir / f"{stem}.result.json"),
        render_html(result, output_dir / f"{stem}.report.html"),
        render_workbook(result, output_dir / f"{stem}.results.xlsx"),
    ]
    return outputs


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ontseq", description="ONTSeq Platform core CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    demo = subparsers.add_parser("demo", help="Generate a synthetic HTML/Excel/JSON report")
    demo.add_argument("--output-dir", type=Path, default=Path("results/demo"))

    validate_manifest = subparsers.add_parser("validate-manifest")
    validate_manifest.add_argument("path", type=Path)

    validate_result = subparsers.add_parser("validate-result")
    validate_result.add_argument("path", type=Path)

    render = subparsers.add_parser("render", help="Render HTML and Excel from a result JSON")
    render.add_argument("result", type=Path)
    render.add_argument("--output-dir", type=Path, required=True)

    reference_lock = subparsers.add_parser(
        "reference-lock", help="Create a versioned reference lock from a FASTA .fai index"
    )
    reference_lock.add_argument("--fai", type=Path, required=True)
    reference_lock.add_argument("--reference-id", required=True)
    reference_lock.add_argument(
        "--genome-build", choices=[item.value for item in GenomeBuild], required=True
    )
    reference_lock.add_argument("--allow-extra-contigs", action="store_true")
    reference_lock.add_argument(
        "--require-canonical-assembly",
        action="store_true",
        help="Require complete canonical chromosomes 1-22, X and Y for the named build",
    )
    reference_lock.add_argument("--output", type=Path, required=True)

    validate_reference = subparsers.add_parser(
        "validate-reference", help="Validate and summarize an existing reference lock"
    )
    validate_reference.add_argument("path", type=Path)
    validate_reference.add_argument(
        "--expected-genome-build", choices=[item.value for item in GenomeBuild]
    )
    validate_reference.add_argument("--require-canonical-assembly", action="store_true")

    inspect_bam = subparsers.add_parser(
        "inspect-bam", help="Run the aligned-BAM integrity and reference gate"
    )
    inspect_bam.add_argument("manifest", type=Path)
    inspect_bam.add_argument("--reference-lock", type=Path, required=True)
    inspect_bam.add_argument("--samtools", default="samtools")
    inspect_bam.add_argument("--checksum", action="store_true")
    inspect_bam.add_argument("--output", type=Path, required=True)

    cramino_qc = subparsers.add_parser(
        "qc-cramino", help="Run Cramino and normalize descriptive BAM QC metrics"
    )
    cramino_qc.add_argument("manifest", type=Path)
    cramino_qc.add_argument("--policy", type=Path, required=True)
    cramino_qc.add_argument("--cramino", default="cramino")
    cramino_qc.add_argument("--threads", type=int, default=4)
    cramino_qc.add_argument("--output", type=Path, required=True)

    target_coverage = subparsers.add_parser(
        "qc-target-coverage",
        help="Run Mosdepth and normalize Adaptive Sampling target-region coverage",
    )
    target_coverage.add_argument("manifest", type=Path)
    target_coverage.add_argument("--intake", type=Path, required=True)
    target_coverage.add_argument("--policy", type=Path, required=True)
    target_coverage.add_argument("--mosdepth", default="mosdepth")
    target_coverage.add_argument("--threads", type=int, default=4)
    target_coverage.add_argument("--output-dir", type=Path, required=True)
    target_coverage.add_argument("--output", type=Path, required=True)

    call_methylation = subparsers.add_parser(
        "call-methylation",
        help="Run modkit pileup and normalize region-aggregated modified-base fractions",
    )
    call_methylation.add_argument("manifest", type=Path)
    call_methylation.add_argument("--intake", type=Path, required=True)
    call_methylation.add_argument("--policy", type=Path, required=True)
    call_methylation.add_argument(
        "--reference-fasta",
        type=Path,
        help="Required when the policy restricts the pileup to CpG sites",
    )
    call_methylation.add_argument("--modkit", default="modkit")
    call_methylation.add_argument("--samtools", default="samtools")
    call_methylation.add_argument("--threads", type=int, default=4)
    call_methylation.add_argument("--output-dir", type=Path, required=True)
    call_methylation.add_argument("--output", type=Path, required=True)

    call_sniffles = subparsers.add_parser(
        "call-sniffles", help="Run Sniffles2 and normalize conservative candidate SV evidence"
    )
    call_sniffles.add_argument("manifest", type=Path)
    call_sniffles.add_argument("--intake", type=Path, required=True)
    call_sniffles.add_argument("--policy", type=Path, required=True)
    call_sniffles.add_argument("--sniffles", default="sniffles")
    call_sniffles.add_argument("--threads", type=int, default=4)
    call_sniffles.add_argument("--vcf", type=Path, required=True)
    call_sniffles.add_argument("--output", type=Path, required=True)

    local_smoke = subparsers.add_parser(
        "local-smoke",
        help="Exercise samtools, Cramino and Sniffles2 with generated synthetic alignments",
    )
    local_smoke.add_argument("--output-dir", type=Path, default=Path("results/local-smoke"))
    local_smoke.add_argument("--qc-policy", type=Path, default=Path("configs/qc/defaults.yaml"))
    local_smoke.add_argument(
        "--sniffles-policy",
        type=Path,
        default=Path("configs/sv/sniffles2.conservative.technical.yaml"),
    )
    local_smoke.add_argument("--samtools", default="samtools")
    local_smoke.add_argument("--cramino", default="cramino")
    local_smoke.add_argument("--sniffles", default="sniffles")
    local_smoke.add_argument("--threads", type=int, default=2)
    local_smoke.add_argument("--git-commit", default="LOCAL_SMOKE")

    system_smoke = subparsers.add_parser(
        "system-smoke",
        help=(
            "Exercise the installed samtools/Cramino/Sniffles path plus canonical "
            "QDNAseq+ACE CNV, reporting, release checksums and resume"
        ),
    )
    system_smoke.add_argument("--output-dir", type=Path, default=Path("results/system-smoke"))
    system_smoke.add_argument("--qc-policy", type=Path, required=True)
    system_smoke.add_argument("--sniffles-policy", type=Path, required=True)
    system_smoke.add_argument("--cnv-policy", type=Path, required=True)
    system_smoke.add_argument("--qdnaseq-rscript", default="Rscript")
    system_smoke.add_argument("--qdnaseq-script", type=Path, required=True)
    system_smoke.add_argument("--samtools", default="samtools")
    system_smoke.add_argument("--cramino", default="cramino")
    system_smoke.add_argument("--sniffles", default="sniffles")
    system_smoke.add_argument("--threads", type=int, default=2)
    system_smoke.add_argument("--git-commit", default="SYSTEM_SMOKE")

    benchmark = subparsers.add_parser(
        "benchmark", help="Benchmark normalized CNV or SV events against a locked truth case"
    )
    benchmark.add_argument("case", type=Path)
    benchmark.add_argument("--output", type=Path, required=True)

    dilution_plan = subparsers.add_parser(
        "dilution-plan",
        help="Lay out an in-silico tumour dilution series from two source BAMs",
    )
    dilution_plan.add_argument("--policy", type=Path, required=True)
    dilution_plan.add_argument("--series-id", required=True)
    dilution_plan.add_argument("--tumor-bam", type=Path, required=True)
    dilution_plan.add_argument("--normal-bam", type=Path, required=True)
    dilution_plan.add_argument("--tumor-sample-id", required=True)
    dilution_plan.add_argument("--normal-sample-id", required=True)
    dilution_plan.add_argument(
        "--genome-build", choices=[item.value for item in GenomeBuild], required=True
    )
    dilution_plan.add_argument("--samtools", default="samtools")
    dilution_plan.add_argument("--threads", type=int, default=4)
    dilution_plan.add_argument("--output", type=Path, required=True)

    dilution_mix = subparsers.add_parser(
        "dilution-mix", help="Materialize the mixed BAMs of a planned dilution series"
    )
    dilution_mix.add_argument("plan", type=Path)
    dilution_mix.add_argument("--tumor-bam", type=Path, required=True)
    dilution_mix.add_argument("--normal-bam", type=Path, required=True)
    dilution_mix.add_argument("--samtools", default="samtools")
    dilution_mix.add_argument("--threads", type=int, default=4)
    dilution_mix.add_argument("--output-dir", type=Path, required=True)
    dilution_mix.add_argument("--output", type=Path, required=True)

    lod = subparsers.add_parser(
        "lod",
        help="Derive a technical detection limit from the benchmark reports of a series",
    )
    lod.add_argument("reports", type=Path, nargs="+")
    lod.add_argument("--policy", type=Path, required=True)
    lod.add_argument("--series-id", required=True)
    lod.add_argument("--output", type=Path, required=True)

    assemble = subparsers.add_parser(
        "assemble-aligned-mvp",
        help="Assemble intake, QC and optional candidate SV evidence into one result",
    )
    assemble.add_argument("manifest", type=Path)
    assemble.add_argument("--intake", type=Path, required=True)
    assemble.add_argument("--qc", type=Path, required=True)
    assemble.add_argument("--sniffles", type=Path)
    assemble.add_argument("--git-commit", default="UNKNOWN")
    assemble.add_argument("--output", type=Path, required=True)
    return parser


def main() -> None:
    args = _parser().parse_args()
    try:
        if args.command == "demo":
            for path in _render(build_demo_result(), args.output_dir):
                print(path)
        elif args.command == "validate-manifest":
            manifest = load_model(args.path, SampleManifest)
            print(f"VALID manifest: {manifest.sample_id}")
        elif args.command == "validate-result":
            result = load_model(args.path, PipelineResult)
            print(f"VALID result: {result.manifest.sample_id}")
        elif args.command == "render":
            result = load_model(args.result, PipelineResult)
            for path in _render(result, args.output_dir):
                print(path)
        elif args.command == "reference-lock":
            lock = reference_lock_from_fai(
                args.fai,
                reference_id=args.reference_id,
                genome_build=GenomeBuild(args.genome_build),
                allow_extra_contigs=args.allow_extra_contigs,
                require_canonical_assembly=args.require_canonical_assembly,
            )
            print(write_json(lock, args.output))
        elif args.command == "validate-reference":
            lock = load_model(args.path, ReferenceLock)
            if args.expected_genome_build and lock.genome_build != GenomeBuild(
                args.expected_genome_build
            ):
                raise ValueError(
                    f"reference lock is {lock.genome_build.value}, not {args.expected_genome_build}"
                )
            naming_style = "not checked"
            if args.require_canonical_assembly:
                summary = validate_canonical_reference(
                    ((item.name, item.length) for item in lock.contigs), lock.genome_build
                )
                naming_style = summary.naming_style
            total_bases = sum(item.length for item in lock.contigs)
            print(
                f"{lock.genome_build.value} \u00b7 {len(lock.contigs)} contigs \u00b7 "
                f"{total_bases} bp \u00b7 {lock.reference_id} \u00b7 {naming_style}"
            )
        elif args.command == "inspect-bam":
            manifest = load_model(args.manifest, SampleManifest)
            lock = load_model(args.reference_lock, ReferenceLock)
            intake_report = AlignedBamInspector(samtools=args.samtools).inspect(
                manifest, lock, include_checksums=args.checksum
            )
            print(write_json(intake_report, args.output))
            if intake_report.verdict == Verdict.FAIL:
                raise SystemExit(2)
        elif args.command == "qc-cramino":
            manifest = load_model(args.manifest, SampleManifest)
            qc_policy = load_model(args.policy, QCPolicy)
            cramino_report = run_cramino_qc(
                manifest,
                qc_policy,
                cramino=args.cramino,
                threads=args.threads,
            )
            print(write_json(cramino_report, args.output))
            if cramino_report.qc.verdict == Verdict.FAIL:
                raise SystemExit(2)
        elif args.command == "qc-target-coverage":
            manifest = load_model(args.manifest, SampleManifest)
            intake = load_model(args.intake, AlignedBamIntakeReport)
            coverage_policy = load_model(args.policy, TargetCoveragePolicy)
            coverage_report = run_target_coverage(
                manifest,
                intake,
                coverage_policy,
                output_dir=args.output_dir,
                mosdepth=args.mosdepth,
                threads=args.threads,
            )
            print(write_json(coverage_report, args.output))
        elif args.command == "call-methylation":
            manifest = load_model(args.manifest, SampleManifest)
            intake = load_model(args.intake, AlignedBamIntakeReport)
            methylation_policy = load_model(args.policy, MethylationPolicy)
            methylation_report = run_methylation(
                manifest,
                intake,
                methylation_policy,
                output_dir=args.output_dir,
                reference_fasta=args.reference_fasta,
                modkit=args.modkit,
                samtools=args.samtools,
                threads=args.threads,
            )
            print(write_json(methylation_report, args.output))
        elif args.command == "call-sniffles":
            manifest = load_model(args.manifest, SampleManifest)
            intake = load_model(args.intake, AlignedBamIntakeReport)
            sniffles_policy = load_model(args.policy, SnifflesPolicy)
            sniffles_call_report = run_sniffles(
                manifest,
                intake,
                sniffles_policy,
                output_vcf=args.vcf,
                sniffles=args.sniffles,
                threads=args.threads,
            )
            print(write_json(sniffles_call_report, args.output))
        elif args.command == "local-smoke":
            qc_policy = load_model(args.qc_policy, QCPolicy)
            sniffles_policy = load_model(args.sniffles_policy, SnifflesPolicy)
            local_smoke_report = run_local_smoke(
                args.output_dir,
                qc_policy,
                sniffles_policy,
                samtools=args.samtools,
                cramino=args.cramino,
                sniffles=args.sniffles,
                threads=args.threads,
                pipeline_version=__version__,
                git_commit=args.git_commit,
            )
            print(args.output_dir / "local-smoke.report.json")
            print(f"PASS: {local_smoke_report.sniffles.accepted_record_count} SV candidate(s)")
        elif args.command == "system-smoke":
            from .cnv.qdnaseq import QDNAseqPolicy
            from .system_smoke import run_system_smoke

            system_smoke_report = run_system_smoke(
                args.output_dir,
                load_model(args.qc_policy, QCPolicy),
                load_model(args.sniffles_policy, SnifflesPolicy),
                load_model(args.cnv_policy, QDNAseqPolicy),
                qdnaseq_script=args.qdnaseq_script,
                samtools=args.samtools,
                cramino=args.cramino,
                sniffles=args.sniffles,
                rscript=args.qdnaseq_rscript,
                threads=args.threads,
                pipeline_version=__version__,
                git_commit=args.git_commit,
            )
            print(args.output_dir / "system-smoke.report.json")
            print(f"PASS: {len(system_smoke_report.checks)} system checks")
        elif args.command == "benchmark":
            case = load_model(args.case, BenchmarkCase)
            print(write_json(benchmark_case(case), args.output))
        elif args.command == "dilution-plan":
            dilution_policy = load_model(args.policy, DilutionPolicy)
            command_runner = SubprocessRunner()
            plan = plan_dilution_series(
                dilution_policy,
                series_id=args.series_id,
                tumor_sample_id=args.tumor_sample_id,
                normal_sample_id=args.normal_sample_id,
                genome_build=GenomeBuild(args.genome_build),
                tumor_read_count=count_reads(
                    args.tumor_bam,
                    runner=command_runner,
                    samtools=args.samtools,
                    threads=args.threads,
                ),
                normal_read_count=count_reads(
                    args.normal_bam,
                    runner=command_runner,
                    samtools=args.samtools,
                    threads=args.threads,
                ),
            )
            print(write_json(plan, args.output))
            for warning in plan.warnings:
                print(f"WARNING: {warning}")
        elif args.command == "dilution-mix":
            series_plan = load_model(args.plan, DilutionSeriesPlan)
            series_report = execute_dilution_series(
                series_plan,
                tumor_bam=args.tumor_bam,
                normal_bam=args.normal_bam,
                output_dir=args.output_dir,
                samtools=args.samtools,
                threads=args.threads,
            )
            print(write_json(series_report, args.output))
        elif args.command == "lod":
            lod_policy = load_model(args.policy, LodPolicy)
            lod_report = evaluate_lod(
                [load_model(path, BenchmarkReport) for path in args.reports],
                lod_policy,
                series_id=args.series_id,
            )
            print(write_json(lod_report, args.output))
            limit = lod_report.detection_limit_fraction
            print(
                f"detection limit: {limit if limit is not None else 'not established'} "
                f"(bracketed: {lod_report.bracketed})"
            )
            for warning in lod_report.warnings:
                print(f"WARNING: {warning}")
        elif args.command == "assemble-aligned-mvp":
            manifest = load_model(args.manifest, SampleManifest)
            intake = load_model(args.intake, AlignedBamIntakeReport)
            qc_report = load_model(args.qc, CraminoQCReport)
            optional_sniffles_report = (
                load_model(args.sniffles, SnifflesCallReport) if args.sniffles else None
            )
            result = assemble_aligned_bam_mvp(
                manifest,
                intake,
                qc_report,
                pipeline_version=__version__,
                git_commit=args.git_commit,
                sniffles_report=optional_sniffles_report,
            )
            print(write_json(result, args.output))
    except (OSError, ValueError, ValidationError, ToolExecutionError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

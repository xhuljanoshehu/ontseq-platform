from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pydantic import ValidationError

from . import __version__
from .align import AlignmentPolicy
from .align_fixture import build_alignment_fixture
from .bam_intake import AlignedBamInspector
from .basecall import BasecallPolicy
from .benchmark import benchmark_case
from .cnv.cytobands import load_cytoband_file
from .cnv.demo import summarize_comparison, summarize_demo, write_demo_benchmark
from .cnv.evaluate import evaluate_case
from .cnv.models import CnvBenchmarkCase, CnvEvaluationReport
from .cnv.strata import aggregate, paired_detection_comparison
from .cnv.truth import truth_from_karyotype
from .demo import build_demo_result
from .execution import ToolExecutionError
from .io import load_model, write_json
from .models import (
    AlignedBamIntakeReport,
    BenchmarkCase,
    CraminoQCReport,
    GenomeBuild,
    InputKind,
    PipelineResult,
    QCPolicy,
    ReferenceLock,
    SampleManifest,
    SnifflesCallReport,
    SnifflesPolicy,
    Verdict,
)
from .mvp import assemble_aligned_bam_mvp
from .pipeline.checks import exit_code as check_exit_code
from .pipeline.checks import render_json as render_checks_json
from .pipeline.checks import render_text as render_checks_text
from .pipeline.lock import RunAlreadyRunning
from .pipeline.runner import RunConfiguration, run_pipeline
from .preflight import PreflightRequest, preflight
from .qc import run_cramino_qc
from .reference import reference_lock_from_fai
from .report import render_html
from .smoke import run_local_smoke
from .sniffles import run_sniffles
from .status import exit_code, render_json, render_ledger, render_text, scan
from .watchfolder import (
    PassResult,
    WatchConfigurationError,
    WatchSettings,
    watch,
)
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


def _print_pass(result: PassResult) -> None:
    """Report one sweep: what ran, and why everything else did not."""
    for attempt in result.attempted:
        print(f"  {attempt.name:<28} {attempt.outcome.value.upper():<10} {attempt.detail}")
    for name, reason in result.skipped:
        print(f"  {name:<28} {'skipped':<10} {reason}")
    if not result.attempted and not result.skipped:
        print("  nothing to do")


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
    reference_lock.add_argument("--output", type=Path, required=True)

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

    benchmark = subparsers.add_parser(
        "benchmark", help="Benchmark normalized CNV or SV events against a locked truth case"
    )
    benchmark.add_argument("case", type=Path)
    benchmark.add_argument("--output", type=Path, required=True)

    run = subparsers.add_parser(
        "run",
        help="Execute the whole pipeline for one sample into a resumable run envelope",
    )
    run.add_argument("manifest", type=Path)
    run.add_argument("--reference-lock", type=Path, required=True)
    run.add_argument("--qc-policy", type=Path, default=Path("configs/qc/defaults.yaml"))
    run.add_argument(
        "--sniffles-policy",
        type=Path,
        default=Path("configs/sv/sniffles2.conservative.technical.yaml"),
    )
    run.add_argument(
        "--alignment-policy",
        type=Path,
        default=Path("configs/alignment/minimap2.ont.technical.yaml"),
    )
    run.add_argument(
        "--basecall-policy",
        type=Path,
        default=Path("configs/basecalling/dorado.technical.yaml"),
    )
    run.add_argument("--reference-fasta", type=Path, help="Required when aligning")
    run.add_argument("--pod5-dir", type=Path, help="Required when starting from POD5")
    run.add_argument("--output-dir", type=Path, default=Path("results/runs"))
    run.add_argument("--run-id", required=True)
    run.add_argument("--threads", type=int, default=4)
    run.add_argument("--git-commit", default="UNKNOWN")
    run.add_argument("--samtools", default="samtools")
    run.add_argument("--cramino", default="cramino")
    run.add_argument("--sniffles", default="sniffles")
    run.add_argument("--minimap2", default="minimap2")
    run.add_argument("--dorado", default="dorado")
    run.add_argument(
        "--force",
        action="store_true",
        help="Re-run every stage instead of resuming unchanged ones",
    )

    preflight_parser = subparsers.add_parser(
        "preflight",
        help="Check a run's preconditions without starting it, creating nothing",
    )
    preflight_parser.add_argument("manifest", type=Path)
    preflight_parser.add_argument("--reference-lock", type=Path, required=True)
    preflight_parser.add_argument("--run-id", required=True)
    preflight_parser.add_argument("--output-dir", type=Path, default=Path("results/runs"))
    preflight_parser.add_argument(
        "--sniffles-policy",
        type=Path,
        default=Path("configs/sv/sniffles2.conservative.technical.yaml"),
    )
    preflight_parser.add_argument(
        "--alignment-policy",
        type=Path,
        default=Path("configs/alignment/minimap2.ont.technical.yaml"),
    )
    preflight_parser.add_argument(
        "--basecall-policy",
        type=Path,
        default=Path("configs/basecalling/dorado.technical.yaml"),
    )
    preflight_parser.add_argument("--reference-fasta", type=Path, help="Required when aligning")
    preflight_parser.add_argument("--pod5-dir", type=Path, help="Required when starting from POD5")
    preflight_parser.add_argument("--samtools", default="samtools")
    preflight_parser.add_argument("--cramino", default="cramino")
    preflight_parser.add_argument("--sniffles", default="sniffles")
    preflight_parser.add_argument("--minimap2", default="minimap2")
    preflight_parser.add_argument("--dorado", default="dorado")
    preflight_parser.add_argument(
        "--require-free-gb",
        type=float,
        help=(
            "Free space this run needs. Without it free space is reported, not judged: no "
            "measured size model for this lab's data exists in this repository"
        ),
    )
    preflight_parser.add_argument(
        "--verbose", action="store_true", help="Also list checks that do not apply"
    )
    preflight_parser.add_argument(
        "--json", action="store_true", dest="as_json", help="Emit JSON for a scheduler"
    )

    status_parser = subparsers.add_parser(
        "status",
        help="Report the state of every run envelope beneath an output directory",
    )
    status_parser.add_argument("output_dir", type=Path)
    status_parser.add_argument("--run-id", help="Restrict the report to one run")
    status_parser.add_argument(
        "--verbose", action="store_true", help="List every stage of every run"
    )
    status_parser.add_argument(
        "--json", action="store_true", dest="as_json", help="Emit JSON for a monitoring check"
    )

    watch_parser = subparsers.add_parser(
        "watch",
        help="Process every ready sample directory in a drop folder, once or continuously",
    )
    watch_parser.add_argument("watch_dir", type=Path)
    watch_parser.add_argument("--output-dir", type=Path, default=Path("results/runs"))
    watch_parser.add_argument(
        "--manifest-template",
        type=Path,
        required=True,
        help="Manifest with the assay-level constants; sample_id and input are filled in",
    )
    watch_parser.add_argument("--reference-lock", type=Path, required=True)
    watch_parser.add_argument(
        "--input-kind",
        required=True,
        choices=[item.value for item in InputKind],
        help="Declared, not sniffed: a drop folder does not alternate kind per sample",
    )
    watch_parser.add_argument("--qc-policy", type=Path, default=Path("configs/qc/defaults.yaml"))
    watch_parser.add_argument(
        "--sniffles-policy",
        type=Path,
        default=Path("configs/sv/sniffles2.conservative.technical.yaml"),
    )
    watch_parser.add_argument(
        "--alignment-policy",
        type=Path,
        default=Path("configs/alignment/minimap2.ont.technical.yaml"),
    )
    watch_parser.add_argument("--reference-fasta", type=Path, help="Required when aligning")
    watch_parser.add_argument(
        "--ready-marker",
        help=(
            "Glob for the file the producer writes when a run is complete, matched anywhere "
            "beneath the sample directory. Authoritative when set. For a GridION: "
            "'final_summary_*.txt'"
        ),
    )
    watch_parser.add_argument(
        "--pod5-subdir",
        help=(
            "Which POD5 directory to basecall when the instrument wrote several, e.g. "
            "'pod5_pass'. Declared rather than guessed: including failed reads changes the "
            "depth distribution that depth-based copy-number methods assume"
        ),
    )
    watch_parser.add_argument(
        "--quiet-seconds",
        type=float,
        default=300.0,
        help="Without a marker, how long a directory must be unmodified. A heuristic",
    )
    watch_parser.add_argument("--run-id-prefix", default="")
    watch_parser.add_argument("--threads", type=int, default=4)
    watch_parser.add_argument("--git-commit", default="UNKNOWN")
    watch_parser.add_argument(
        "--retry-failed",
        action="store_true",
        help="Re-attempt samples that failed before. Use once the cause is understood",
    )
    watch_parser.add_argument(
        "--once", action="store_true", help="Make a single pass and exit, for cron"
    )
    watch_parser.add_argument("--poll-seconds", type=float, default=60.0)

    align_fixture = subparsers.add_parser(
        "align-fixture",
        help=(
            "Write a synthetic unaligned BAM and reference so the alignment lane can be "
            "exercised with real tools"
        ),
    )
    align_fixture.add_argument("--output-dir", type=Path, default=Path("results/align-fixture"))
    align_fixture.add_argument("--samtools", default="samtools")

    cnv_evaluate = subparsers.add_parser(
        "cnv-evaluate",
        help="Score a CNV call set against a truth set over an explicit evaluable genome",
    )
    cnv_evaluate.add_argument("case", type=Path)
    cnv_evaluate.add_argument("--evaluation-id")
    cnv_evaluate.add_argument("--output", type=Path, required=True)

    cnv_aggregate = subparsers.add_parser(
        "cnv-aggregate",
        help="Pool CNV evaluations of one method into a stratified benchmark summary",
    )
    cnv_aggregate.add_argument("reports", type=Path, nargs="+")
    cnv_aggregate.add_argument("--aggregate-id", required=True)
    cnv_aggregate.add_argument("--target-detection-rate", type=float, default=0.95)
    cnv_aggregate.add_argument("--output", type=Path, required=True)

    cnv_compare = subparsers.add_parser(
        "cnv-compare-methods",
        help="Compare two CNV methods pairwise on the truth events both could assess",
    )
    cnv_compare.add_argument("--method-a", type=Path, nargs="+", required=True)
    cnv_compare.add_argument("--method-b", type=Path, nargs="+", required=True)
    cnv_compare.add_argument("--output", type=Path, required=True)

    cnv_karyotype = subparsers.add_parser(
        "cnv-karyotype-truth",
        help="Convert an ISCN karyotype into a band-resolved CNV truth set",
    )
    cnv_karyotype.add_argument("--karyotype", required=True)
    cnv_karyotype.add_argument("--cytobands", type=Path, required=True)
    cnv_karyotype.add_argument("--cytoband-resource-id", required=True)
    cnv_karyotype.add_argument(
        "--genome-build", choices=[item.value for item in GenomeBuild], required=True
    )
    cnv_karyotype.add_argument("--truth-id", required=True)
    cnv_karyotype.add_argument("--sample-id", required=True)
    cnv_karyotype.add_argument("--source-version", required=True)
    cnv_karyotype.add_argument("--resolution-bp", type=int, default=10_000_000)
    cnv_karyotype.add_argument("--output", type=Path, required=True)

    cnv_demo = subparsers.add_parser(
        "cnv-demo-benchmark",
        help="Run the fully synthetic CNV dilution and coverage benchmark end to end",
    )
    cnv_demo.add_argument("--output-dir", type=Path, default=Path("results/cnv-demo"))
    cnv_demo.add_argument("--replicates", type=int, default=3)
    cnv_demo.add_argument("--bin-size-bp", type=int, default=1_000_000)
    cnv_demo.add_argument("--seed", type=int, default=20260816)

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
            )
            print(write_json(lock, args.output))
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
            smoke_report = run_local_smoke(
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
            print(f"PASS: {smoke_report.sniffles.accepted_record_count} SV candidate(s)")
        elif args.command == "benchmark":
            case = load_model(args.case, BenchmarkCase)
            print(write_json(benchmark_case(case), args.output))
        elif args.command == "run":
            run_manifest = load_model(args.manifest, SampleManifest)
            configuration = RunConfiguration(
                manifest=run_manifest,
                reference_lock=load_model(args.reference_lock, ReferenceLock),
                output_base=args.output_dir,
                run_id=args.run_id,
                pipeline_version=__version__,
                git_commit=args.git_commit,
                qc_policy=load_model(args.qc_policy, QCPolicy),
                sniffles_policy=(
                    load_model(args.sniffles_policy, SnifflesPolicy)
                    if args.sniffles_policy.is_file()
                    else None
                ),
                alignment_policy=(
                    load_model(args.alignment_policy, AlignmentPolicy)
                    if args.alignment_policy.is_file()
                    else None
                ),
                basecall_policy=(
                    load_model(args.basecall_policy, BasecallPolicy)
                    if args.basecall_policy.is_file()
                    else None
                ),
                reference_fasta=args.reference_fasta,
                pod5_directory=args.pod5_dir,
                threads=args.threads,
                executables={
                    "samtools": args.samtools,
                    "cramino": args.cramino,
                    "sniffles": args.sniffles,
                    "minimap2": args.minimap2,
                    "dorado": args.dorado,
                },
                force=args.force,
            )
            run_report, release_bundle = run_pipeline(configuration)
            for stage_record in run_report.stages:
                marker = "resumed" if stage_record.resumed else stage_record.status.value
                print(f"  {stage_record.stage.value:<16} {marker:<10} {stage_record.reason}")
            outcome = "PASS" if run_report.passed else "FAIL"
            print(f"verdict: {outcome} - {run_report.verdict_reason}")
            if run_report.unverified_stages:
                names = ", ".join(item.value for item in run_report.unverified_stages)
                print(f"UNVERIFIED ADAPTERS COMPLETED: {names}")
            if release_bundle is not None:
                print(
                    f"release bundle: {len(release_bundle.artifacts)} artifact(s), "
                    f"{len(release_bundle.withheld_artifact_paths)} withheld, unsigned"
                )
            if not run_report.passed:
                raise SystemExit(2)
        elif args.command == "preflight":
            request = PreflightRequest(
                manifest=load_model(args.manifest, SampleManifest),
                reference_lock=load_model(args.reference_lock, ReferenceLock),
                output_base=args.output_dir,
                run_id=args.run_id,
                executables={
                    "samtools": args.samtools,
                    "cramino": args.cramino,
                    "sniffles": args.sniffles,
                    "minimap2": args.minimap2,
                    "dorado": args.dorado,
                },
                reference_fasta=args.reference_fasta,
                pod5_directory=args.pod5_dir,
                alignment_policy=(
                    load_model(args.alignment_policy, AlignmentPolicy)
                    if args.alignment_policy.is_file()
                    else None
                ),
                basecall_policy=(
                    load_model(args.basecall_policy, BasecallPolicy)
                    if args.basecall_policy.is_file()
                    else None
                ),
                sniffles_policy=(
                    load_model(args.sniffles_policy, SnifflesPolicy)
                    if args.sniffles_policy.is_file()
                    else None
                ),
                require_free_gb=args.require_free_gb,
            )
            checks = preflight(request)
            if args.as_json:
                print(render_checks_json(checks), end="")
            else:
                print(render_checks_text(checks, verbose=args.verbose))
            # 2 when at least one precondition makes the run impossible; a warning or an
            # unanswerable question does not block, or the command would be unusable on
            # exactly the machines it exists to help.
            code = check_exit_code(checks)
            if code:
                raise SystemExit(code)
        elif args.command == "status":
            statuses = scan(args.output_dir, run_id=args.run_id)
            if args.as_json:
                print(render_json(statuses), end="")
            else:
                print(render_text(statuses, verbose=args.verbose))
                ledger = render_ledger(args.output_dir)
                if ledger:
                    print(ledger)
            # 0 nothing wrong, 2 a run failed or its report is unreadable, 6 a run was
            # interrupted or never reached a verdict. A run merely in progress is not a
            # problem, and a check that fires during normal work teaches people to ignore it.
            code = exit_code(statuses)
            if code:
                raise SystemExit(code)
        elif args.command == "watch":
            settings = WatchSettings(
                watch_dir=args.watch_dir,
                output_dir=args.output_dir,
                manifest_template=args.manifest_template,
                reference_lock=args.reference_lock,
                qc_policy=args.qc_policy,
                input_kind=InputKind(args.input_kind),
                sniffles_policy=args.sniffles_policy,
                alignment_policy=args.alignment_policy,
                reference_fasta=args.reference_fasta,
                run_id_prefix=args.run_id_prefix,
                ready_marker=args.ready_marker,
                pod5_subdirectory=args.pod5_subdir,
                quiet_seconds=args.quiet_seconds,
                threads=args.threads,
                git_commit=args.git_commit,
                retry_failed=args.retry_failed,
            )
            passes = watch(
                settings,
                once=args.once,
                poll_seconds=args.poll_seconds,
                report=_print_pass,
            )
            # A failed sample makes the whole sweep non-zero so a cron job is noticed; a
            # sample that is merely not ready yet, or already done, does not.
            if any(item.failures for item in passes):
                raise SystemExit(2)
        elif args.command == "align-fixture":
            fixture = build_alignment_fixture(args.output_dir, samtools=args.samtools)
            for path in (
                fixture.reference_fasta,
                fixture.reference_fai,
                fixture.unaligned_bam,
                fixture.manifest,
                fixture.reference_lock,
            ):
                print(path)
        elif args.command == "cnv-evaluate":
            cnv_case = load_model(args.case, CnvBenchmarkCase)
            print(
                write_json(evaluate_case(cnv_case, evaluation_id=args.evaluation_id), args.output)
            )
        elif args.command == "cnv-aggregate":
            evaluations = [load_model(path, CnvEvaluationReport) for path in args.reports]
            summary = aggregate(
                evaluations,
                aggregate_id=args.aggregate_id,
                target_detection_rate=args.target_detection_rate,
            )
            print(write_json(summary, args.output))
            for line in summarize_demo(summary):
                print(line)
        elif args.command == "cnv-compare-methods":
            comparison = paired_detection_comparison(
                [load_model(path, CnvEvaluationReport) for path in args.method_a],
                [load_model(path, CnvEvaluationReport) for path in args.method_b],
            )
            print(write_json(comparison, args.output))
            print(
                f"{comparison.method_a} vs {comparison.method_b}: "
                f"{comparison.paired_events} paired event(s), "
                f"only-A={comparison.only_a_detected} only-B={comparison.only_b_detected}, "
                f"p={'undefined' if comparison.p_value is None else f'{comparison.p_value:.4f}'}"
            )
            print(comparison.note)
        elif args.command == "cnv-karyotype-truth":
            table = load_cytoband_file(
                args.cytobands,
                genome_build=GenomeBuild(args.genome_build),
                resource_id=args.cytoband_resource_id,
            )
            karyotype_truth, conversion = truth_from_karyotype(
                truth_id=args.truth_id,
                sample_id=args.sample_id,
                karyotype=args.karyotype,
                cytobands=table,
                source_version=args.source_version,
                resolution_bp=args.resolution_bp,
            )
            print(write_json(karyotype_truth, args.output))
            for construct in conversion.unsupported:
                print(f"UNSUPPORTED {construct.token}: {construct.reason}")
            for balanced in conversion.balanced_constructs:
                print(f"BALANCED {balanced}: asserts no copy-number change")
            if conversion.unsupported:
                # A partially converted karyotype is not a usable truth set. Exiting
                # non-zero stops a pipeline from scoring against an incomplete truth.
                raise SystemExit(3)
        elif args.command == "cnv-demo-benchmark":
            outputs = write_demo_benchmark(
                args.output_dir,
                replicates=args.replicates,
                bin_size_bp=args.bin_size_bp,
                seed=args.seed,
            )
            print(outputs.truth_path)
            for path in outputs.aggregate_paths:
                print(path)
            print(outputs.comparison_path)
            for summary in outputs.aggregates:
                for line in summarize_demo(summary):
                    print(line)
            for line in summarize_comparison(outputs.comparison):
                print(line)
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
    except WatchConfigurationError as exc:
        # 5, not 3: exit 3 already means a partially converted karyotype elsewhere. Run
        # outcomes keep their own codes — 2 a failed run, 4 a locked envelope, 5 a
        # configuration nothing can run under.
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(5) from exc
    except RunAlreadyRunning as exc:
        # Its own exit code, so a watcher or scheduler can tell "someone else is already
        # on this sample" apart from "this run failed" and simply move on.
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(4) from exc
    except (OSError, ValueError, ValidationError, ToolExecutionError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

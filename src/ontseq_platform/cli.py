from __future__ import annotations

import argparse
from pathlib import Path

from pydantic import ValidationError

from . import __version__
from .annotation import annotate_result, load_clinvar
from .annotation import describe as describe_annotation
from .bam_intake import AlignedBamInspector
from .benchmark import benchmark_case
from .cnv.adapters import call_set_from_qdnaseq_report
from .cnv.cytobands import load_cytoband_file
from .cnv.demo import summarize_comparison, summarize_demo, write_demo_benchmark
from .cnv.evaluate import evaluate_case
from .cnv.models import CnvBenchmarkCase, CnvDataBasis, CnvEvaluationReport
from .cnv.qdnaseq import QDNAseqCallReport
from .cnv.strata import aggregate, paired_detection_comparison
from .cnv.truth import truth_from_karyotype
from .demo import build_demo_result
from .execution import ToolExecutionError
from .io import load_model, write_json
from .knowledge.annotate import DEFAULT_EXACT_TOLERANCE_BP, DEFAULT_MINIMUM_OVERLAP
from .models import (
    AlignedBamIntakeReport,
    BenchmarkCase,
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
from .reference import reference_lock_from_fai
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

    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Attach knowledge-base records to a result's events, without classifying them",
    )
    annotate_parser.add_argument("result", type=Path, help="A validated result JSON")
    annotate_parser.add_argument(
        "--clinvar",
        type=Path,
        required=True,
        help="NCBI variant_summary.txt for the result's genome build",
    )
    annotate_parser.add_argument(
        "--release",
        required=True,
        help="The publisher's release identifier, e.g. 2026-08-01. Recorded with every "
        "annotation so a report can be reproduced",
    )
    annotate_parser.add_argument("--output", type=Path, required=True)
    annotate_parser.add_argument(
        "--minimum-overlap",
        type=float,
        default=DEFAULT_MINIMUM_OVERLAP,
        help="Reciprocal overlap for a partial match. An engineering default, not a "
        "validated concordance criterion",
    )
    annotate_parser.add_argument(
        "--exact-tolerance-bp",
        type=int,
        default=DEFAULT_EXACT_TOLERANCE_BP,
        help="Breakpoint slack within which a match counts as exact",
    )

    cnv_evaluate = subparsers.add_parser(
        "cnv-evaluate",
        help="Score a CNV call set against a truth set over an explicit evaluable genome",
    )
    cnv_evaluate.add_argument("case", type=Path)
    cnv_evaluate.add_argument("--evaluation-id")
    cnv_evaluate.add_argument("--output", type=Path, required=True)

    cnv_from_qdnaseq = subparsers.add_parser(
        "cnv-callset-from-qdnaseq",
        help="Normalize a QDNAseq/ACE run into a scoreable, non-reportable CNV call set",
    )
    cnv_from_qdnaseq.add_argument("report", type=Path, help="A run's *.qdnaseq.json report")
    cnv_from_qdnaseq.add_argument("--call-set-id", required=True)
    cnv_from_qdnaseq.add_argument(
        "--data-basis",
        required=True,
        choices=[item.value for item in CnvDataBasis],
        help="Which read population the estimate came from. Stated, never inferred: an "
        "adaptive-sampling run holds two populations whose depth behaviour differs",
    )
    cnv_from_qdnaseq.add_argument(
        "--segment-dir",
        type=Path,
        help="Directory holding the run's segment tables. Defaults to the qdnaseq/ "
        "directory beside the report",
    )
    cnv_from_qdnaseq.add_argument(
        "--reference-lock",
        type=Path,
        help="Contig lengths, so everything the segmentation does not cover is declared a "
        "no-call instead of being scored as agreement",
    )
    cnv_from_qdnaseq.add_argument(
        "--bin-size-kbp",
        type=int,
        help="Score this resolution instead of the run's primary fit",
    )
    cnv_from_qdnaseq.add_argument("--mean-coverage-x", type=float)
    cnv_from_qdnaseq.add_argument("--output", type=Path, required=True)

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
        elif args.command == "annotate":
            annotate_result_input = load_model(args.result, PipelineResult)
            clinvar_records, clinvar_lock = load_clinvar(
                args.clinvar,
                genome_build=annotate_result_input.manifest.assay.genome_build,
                release=args.release,
            )
            annotation_outcome = annotate_result(
                annotate_result_input,
                clinvar_records,
                lock=clinvar_lock,
                minimum_reciprocal_overlap=args.minimum_overlap,
                exact_tolerance_bp=args.exact_tolerance_bp,
            )
            print(write_json(annotation_outcome.result, args.output))
            for line in describe_annotation(annotation_outcome):
                print(line)
        elif args.command == "cnv-evaluate":
            cnv_case = load_model(args.case, CnvBenchmarkCase)
            print(
                write_json(evaluate_case(cnv_case, evaluation_id=args.evaluation_id), args.output)
            )
        elif args.command == "cnv-callset-from-qdnaseq":
            qdnaseq_report = load_model(args.report, QDNAseqCallReport)
            segment_dir = args.segment_dir or args.report.parent / "qdnaseq"
            contig_lengths = None
            if args.reference_lock:
                contig_lengths = {
                    contig.name: contig.length
                    for contig in load_model(args.reference_lock, ReferenceLock).contigs
                }
            qdnaseq_call_set = call_set_from_qdnaseq_report(
                qdnaseq_report,
                call_set_id=args.call_set_id,
                data_basis=CnvDataBasis(args.data_basis),
                output_dir=segment_dir,
                contig_lengths=contig_lengths,
                bin_size_kbp=args.bin_size_kbp,
                mean_coverage_x=args.mean_coverage_x,
            )
            print(write_json(qdnaseq_call_set, args.output))
            print(
                f"{qdnaseq_call_set.method_version}: {len(qdnaseq_call_set.segments)} segment(s), "
                f"{len(qdnaseq_call_set.no_call_regions)} declared no-call region(s), "
                f"status {qdnaseq_call_set.status.value}"
            )
            if contig_lengths is None:
                print(
                    "NOTE no reference lock was supplied, so uncovered regions are not "
                    "declared. Score against a mask that excludes them."
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
    except (OSError, ValueError, ValidationError, ToolExecutionError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

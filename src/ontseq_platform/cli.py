from __future__ import annotations

import argparse
from pathlib import Path

from pydantic import ValidationError

from . import __version__
from .bam_intake import AlignedBamInspector
from .benchmark import benchmark_case
from .demo import build_demo_result
from .execution import ToolExecutionError
from .io import load_model, write_json
from .models import (
    AlignedBamIntakeReport,
    BenchmarkCase,
    CraminoQCReport,
    GenomeBuild,
    PipelineResult,
    QCPolicy,
    ReferenceLock,
    SampleManifest,
    Verdict,
)
from .mvp import assemble_aligned_bam_mvp
from .qc import run_cramino_qc
from .reference import reference_lock_from_fai
from .report import render_html
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

    benchmark = subparsers.add_parser(
        "benchmark", help="Benchmark normalized CNV or SV events against a locked truth case"
    )
    benchmark.add_argument("case", type=Path)
    benchmark.add_argument("--output", type=Path, required=True)

    assemble = subparsers.add_parser(
        "assemble-aligned-mvp", help="Assemble intake and Cramino QC into a reportable MVP result"
    )
    assemble.add_argument("manifest", type=Path)
    assemble.add_argument("--intake", type=Path, required=True)
    assemble.add_argument("--qc", type=Path, required=True)
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
            policy = load_model(args.policy, QCPolicy)
            cramino_report = run_cramino_qc(
                manifest,
                policy,
                cramino=args.cramino,
                threads=args.threads,
            )
            print(write_json(cramino_report, args.output))
            if cramino_report.qc.verdict == Verdict.FAIL:
                raise SystemExit(2)
        elif args.command == "benchmark":
            case = load_model(args.case, BenchmarkCase)
            print(write_json(benchmark_case(case), args.output))
        elif args.command == "assemble-aligned-mvp":
            manifest = load_model(args.manifest, SampleManifest)
            intake = load_model(args.intake, AlignedBamIntakeReport)
            qc_report = load_model(args.qc, CraminoQCReport)
            result = assemble_aligned_bam_mvp(
                manifest,
                intake,
                qc_report,
                pipeline_version=__version__,
                git_commit=args.git_commit,
            )
            print(write_json(result, args.output))
    except (OSError, ValueError, ValidationError, ToolExecutionError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

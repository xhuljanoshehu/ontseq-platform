from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pydantic import ValidationError

from . import __version__
from .align import AlignmentPolicy
from .align_fixture import build_alignment_fixture
from .basecall import BasecallPolicy
from .execution import ToolExecutionError
from .io import load_model
from .model_lock import ModelLockError
from .model_lock import exit_code as model_lock_exit_code
from .model_lock import fingerprint as model_fingerprint
from .model_lock import render as render_model_lock
from .models import InputKind, QCPolicy, ReferenceLock, SampleManifest, SnifflesPolicy
from .pipeline.checks import exit_code as check_exit_code
from .pipeline.checks import render_json as render_checks_json
from .pipeline.checks import render_text as render_checks_text
from .pipeline.lock import RunAlreadyRunning
from .pipeline.review import Decision, ReviewError
from .pipeline.review import exit_code as review_exit_code
from .pipeline.runner import EnvelopeAlreadyReviewed, RunConfiguration, run_pipeline
from .preflight import PreflightRequest, preflight
from .review import inspect as inspect_review
from .review import record as record_review
from .review import render_json as render_review_json
from .review import render_text as render_review_text
from .service.app import ServiceConfig, serve
from .status import exit_code as status_exit_code
from .status import render_json as render_status_json
from .status import render_ledger, render_text as render_status_text, scan
from .watchfolder import PassResult, WatchConfigurationError, WatchSettings, watch

RUNTIME_COMMANDS = frozenset(
    {"run", "preflight", "model-lock", "serve", "review", "status", "watch", "align-fixture"}
)


def _alignment_policy(path: Path) -> AlignmentPolicy | None:
    return load_model(path, AlignmentPolicy) if path.is_file() else None


def _basecall_policy(path: Path) -> BasecallPolicy | None:
    return load_model(path, BasecallPolicy) if path.is_file() else None


def _sniffles_policy(path: Path) -> SnifflesPolicy | None:
    return load_model(path, SnifflesPolicy) if path.is_file() else None


def _add_execution_options(parser: argparse.ArgumentParser, *, include_qc: bool) -> None:
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--reference-lock", type=Path, required=True)
    if include_qc:
        parser.add_argument("--qc-policy", type=Path, default=Path("configs/qc/defaults.yaml"))
    parser.add_argument(
        "--sniffles-policy",
        type=Path,
        default=Path("configs/sv/sniffles2.conservative.technical.yaml"),
    )
    parser.add_argument(
        "--alignment-policy",
        type=Path,
        default=Path("configs/alignment/minimap2.ont.technical.yaml"),
    )
    parser.add_argument(
        "--basecall-policy",
        type=Path,
        default=Path("configs/basecalling/dorado.technical.yaml"),
    )
    parser.add_argument("--reference-fasta", type=Path)
    parser.add_argument("--pod5-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("results/runs"))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--cramino", default="cramino")
    parser.add_argument("--sniffles", default="sniffles")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--dorado", default="dorado")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ontseq",
        description="ONTSeq end-to-end execution and operational commands",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Execute one sample into a resumable run envelope")
    _add_execution_options(run, include_qc=True)
    run.add_argument("--threads", type=int, default=4)
    run.add_argument("--git-commit", default="UNKNOWN")
    run.add_argument("--force", action="store_true")

    pf = sub.add_parser("preflight", help="Check run preconditions without creating output")
    _add_execution_options(pf, include_qc=False)
    pf.add_argument("--require-free-gb", type=float)
    pf.add_argument("--verbose", action="store_true")
    pf.add_argument("--json", action="store_true", dest="as_json")

    lock = sub.add_parser("model-lock", help="Fingerprint a downloaded Dorado model directory")
    lock.add_argument("model", type=Path)
    lock.add_argument("--list-files", action="store_true")
    lock.add_argument("--json", action="store_true", dest="as_json")

    srv = sub.add_parser("serve", help="Run the loopback-only local operator service")
    srv.add_argument("--reference-lock", type=Path, required=True)
    srv.add_argument("--allow-root", type=Path, action="append", required=True, dest="allow_roots")
    srv.add_argument("--output-dir", type=Path, default=Path("results/runs"))
    srv.add_argument("--qc-policy", type=Path, default=Path("configs/qc/defaults.yaml"))
    srv.add_argument(
        "--sniffles-policy",
        type=Path,
        default=Path("configs/sv/sniffles2.conservative.technical.yaml"),
    )
    srv.add_argument("--port", type=int, default=8765)
    srv.add_argument("--threads", type=int, default=4)
    srv.add_argument("--no-browser", action="store_true")

    review = sub.add_parser("review", help="Record or inspect review state")
    review_sub = review.add_subparsers(dest="review_command", required=True)
    record = review_sub.add_parser("record")
    record.add_argument("envelope", type=Path)
    record.add_argument("--decision", required=True, choices=[item.value for item in Decision])
    record.add_argument("--reviewer", required=True)
    record.add_argument("--note", default="")
    review_status = review_sub.add_parser("status")
    review_status.add_argument("envelope", type=Path)
    review_status.add_argument("--verbose", action="store_true")
    review_status.add_argument("--json", action="store_true", dest="as_json")
    review_status.add_argument("--require-reviewers", type=int, default=0)

    status = sub.add_parser("status", help="Summarize run envelopes")
    status.add_argument("output_dir", type=Path)
    status.add_argument("--run-id")
    status.add_argument("--verbose", action="store_true")
    status.add_argument("--json", action="store_true", dest="as_json")

    watcher = sub.add_parser("watch", help="Process ready sample directories in a drop folder")
    watcher.add_argument("watch_dir", type=Path)
    watcher.add_argument("--output-dir", type=Path, default=Path("results/runs"))
    watcher.add_argument("--manifest-template", type=Path, required=True)
    watcher.add_argument("--reference-lock", type=Path, required=True)
    watcher.add_argument("--input-kind", required=True, choices=[item.value for item in InputKind])
    watcher.add_argument("--qc-policy", type=Path, default=Path("configs/qc/defaults.yaml"))
    watcher.add_argument(
        "--sniffles-policy",
        type=Path,
        default=Path("configs/sv/sniffles2.conservative.technical.yaml"),
    )
    watcher.add_argument(
        "--alignment-policy",
        type=Path,
        default=Path("configs/alignment/minimap2.ont.technical.yaml"),
    )
    watcher.add_argument("--reference-fasta", type=Path)
    watcher.add_argument("--ready-marker")
    watcher.add_argument("--pod5-subdir")
    watcher.add_argument("--quiet-seconds", type=float, default=300.0)
    watcher.add_argument("--run-id-prefix", default="")
    watcher.add_argument("--threads", type=int, default=4)
    watcher.add_argument("--git-commit", default="UNKNOWN")
    watcher.add_argument("--retry-failed", action="store_true")
    watcher.add_argument("--once", action="store_true")
    watcher.add_argument("--poll-seconds", type=float, default=60.0)

    fixture = sub.add_parser("align-fixture", help="Generate a synthetic real-alignment fixture")
    fixture.add_argument("--output-dir", type=Path, default=Path("results/align-fixture"))
    fixture.add_argument("--samtools", default="samtools")
    return parser


def _executables(args: argparse.Namespace) -> dict[str, str]:
    return {
        "samtools": args.samtools,
        "cramino": args.cramino,
        "sniffles": args.sniffles,
        "minimap2": args.minimap2,
        "dorado": args.dorado,
    }


def _print_pass(result: PassResult) -> None:
    for attempt in result.attempted:
        print(f"  {attempt.name:<28} {attempt.outcome.value.upper():<10} {attempt.detail}")
    for name, reason in result.skipped:
        print(f"  {name:<28} {'skipped':<10} {reason}")


def main() -> None:
    args = _parser().parse_args()
    try:
        if args.command == "run":
            config = RunConfiguration(
                manifest=load_model(args.manifest, SampleManifest),
                reference_lock=load_model(args.reference_lock, ReferenceLock),
                output_base=args.output_dir,
                run_id=args.run_id,
                pipeline_version=__version__,
                git_commit=args.git_commit,
                qc_policy=load_model(args.qc_policy, QCPolicy),
                sniffles_policy=_sniffles_policy(args.sniffles_policy),
                alignment_policy=_alignment_policy(args.alignment_policy),
                basecall_policy=_basecall_policy(args.basecall_policy),
                reference_fasta=args.reference_fasta,
                pod5_directory=args.pod5_dir,
                threads=args.threads,
                executables=_executables(args),
                force=args.force,
            )
            report, release = run_pipeline(config)
            for stage in report.stages:
                marker = "resumed" if stage.resumed else stage.status.value
                print(f"  {stage.stage.value:<16} {marker:<10} {stage.reason}")
            print(f"verdict: {'PASS' if report.passed else 'FAIL'} - {report.verdict_reason}")
            if report.unverified_stages:
                names = ", ".join(item.value for item in report.unverified_stages)
                print(f"UNVERIFIED ADAPTERS COMPLETED: {names}")
            if release is not None:
                print(f"release bundle: {len(release.artifacts)} artifact(s), unsigned")
            if not report.passed:
                raise SystemExit(2)

        elif args.command == "preflight":
            request = PreflightRequest(
                manifest=load_model(args.manifest, SampleManifest),
                reference_lock=load_model(args.reference_lock, ReferenceLock),
                output_base=args.output_dir,
                run_id=args.run_id,
                executables=_executables(args),
                reference_fasta=args.reference_fasta,
                pod5_directory=args.pod5_dir,
                alignment_policy=_alignment_policy(args.alignment_policy),
                basecall_policy=_basecall_policy(args.basecall_policy),
                sniffles_policy=_sniffles_policy(args.sniffles_policy),
                require_free_gb=args.require_free_gb,
            )
            checks = preflight(request)
            if args.as_json:
                print(render_checks_json(checks), end="")
            else:
                print(render_checks_text(checks, verbose=args.verbose))
            code = check_exit_code(checks)
            if code:
                raise SystemExit(code)

        elif args.command == "model-lock":
            model = model_fingerprint(args.model)
            if args.as_json:
                print(
                    json.dumps(
                        {
                            "path": str(model.path),
                            "sha256": model.signature,
                            "file_count": model.file_count,
                            "total_bytes": model.total_bytes,
                            "concerns": list(model.concerns),
                        },
                        indent=2,
                    )
                )
            else:
                print(render_model_lock(model, list_files=args.list_files))
            raise SystemExit(model_lock_exit_code(model))

        elif args.command == "serve":
            serve(
                ServiceConfig(
                    reference_lock=args.reference_lock,
                    output_dir=args.output_dir,
                    allowed_roots=list(args.allow_roots),
                    qc_policy=args.qc_policy,
                    sniffles_policy=args.sniffles_policy,
                    port=args.port,
                    threads=args.threads,
                ),
                open_browser=not args.no_browser,
            )

        elif args.command == "review":
            if args.review_command == "record":
                entry = record_review(
                    args.envelope,
                    decision=Decision(args.decision),
                    reviewer=args.reviewer,
                    note=args.note,
                )
                print(entry.describe())
                print(f"bound to release bundle {entry.release_sha256}")
            else:
                report = inspect_review(args.envelope)
                if args.as_json:
                    print(render_review_json(report), end="")
                else:
                    print(render_review_text(report, verbose=args.verbose))
                code = review_exit_code(
                    report.state,
                    reviewers=len(report.reviewers),
                    required_reviewers=args.require_reviewers,
                )
                if code:
                    raise SystemExit(code)

        elif args.command == "status":
            statuses = scan(args.output_dir, run_id=args.run_id)
            if args.as_json:
                print(render_status_json(statuses), end="")
            else:
                print(render_status_text(statuses, verbose=args.verbose))
                ledger = render_ledger(args.output_dir)
                if ledger:
                    print(ledger)
            code = status_exit_code(statuses)
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

    except WatchConfigurationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(5) from exc
    except EnvelopeAlreadyReviewed as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(7) from exc
    except ReviewError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
    except ModelLockError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
    except RunAlreadyRunning as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(4) from exc
    except (OSError, ValueError, ValidationError, ToolExecutionError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

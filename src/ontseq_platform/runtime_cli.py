from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
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
from .models import (
    AmlKnowledgeLock,
    CuteSvPolicy,
    InputKind,
    IntervalResourceLock,
    QCPolicy,
    ReferenceLock,
    SampleManifest,
    SnifflesPolicy,
    SvConsensusPolicy,
    SvEvidencePolicy,
)
from .pipeline.checks import exit_code as check_exit_code
from .pipeline.checks import render_json as render_checks_json
from .pipeline.checks import render_text as render_checks_text
from .pipeline.components import RunComponents
from .pipeline.lock import RunAlreadyRunning
from .pipeline.review import Decision, ReviewError
from .pipeline.review import exit_code as review_exit_code
from .pipeline.runner import EnvelopeAlreadyReviewed, RunConfiguration, run_pipeline
from .pipeline.stages import StageId
from .preflight import PreflightRequest, preflight
from .profile_analysis import AnalyzeSettings, build_profile_run_configuration
from .resource_commands import add_references_parser, handle_references_command
from .review import inspect as inspect_review
from .review import record as record_review
from .review import render_json as render_review_json
from .review import render_text as render_review_text
from .service.app import ServiceConfig, serve
from .status import exit_code as status_exit_code
from .status import render_json as render_status_json
from .status import render_ledger, scan
from .status import render_text as render_status_text
from .target_coverage import TargetCoveragePolicy
from .watchfolder import PassResult, WatchConfigurationError, WatchSettings, watch

RUNTIME_COMMANDS = frozenset(
    {
        "run",
        "analyze",
        "references",
        "preflight",
        "model-lock",
        "serve",
        "review",
        "status",
        "watch",
        "align-fixture",
    }
)


SELECTABLE_STAGES = (
    StageId.BASECALL,
    StageId.ALIGN,
    StageId.QC,
    StageId.TARGET_COVERAGE,
    StageId.CNV,
    StageId.SV,
)


def _selected_policy(selection: RunComponents | None, stage: StageId, fallback: Path) -> Path:
    """A selection may name the policy file, so one document configures the whole run."""
    choice = selection.choice_for(stage) if selection is not None else None
    if choice is not None and choice.policy:
        return Path(choice.policy)
    return fallback


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _alignment_policy(path: Path) -> AlignmentPolicy | None:
    return load_model(path, AlignmentPolicy) if path.is_file() else None


def _basecall_policy(path: Path) -> BasecallPolicy | None:
    return load_model(path, BasecallPolicy) if path.is_file() else None


def _sniffles_policy(path: Path) -> SnifflesPolicy | None:
    return load_model(path, SnifflesPolicy) if path.is_file() else None


def _cutesv_policy(path: Path | None) -> CuteSvPolicy | None:
    return load_model(path, CuteSvPolicy) if path is not None and path.is_file() else None


def _cutesv_policy_for_run(
    path: Path | None,
    reference_fasta: Path | None,
) -> CuteSvPolicy | None:
    """Enable the optional second SV caller only when its required FASTA is present.

    The advanced ``ontseq run`` command predates profile-managed references. Its default
    policy path must not turn an aligned-BAM CNV-only run into a failed cuteSV invocation
    when the operator did not provide ``--reference-fasta``. Profile runs always resolve
    and validate the FASTA before constructing their configuration, so they continue to
    run the pinned Sniffles+cuteSV lane.
    """

    if reference_fasta is None:
        return None
    return _cutesv_policy(path)


def _sv_consensus_policy(path: Path | None) -> SvConsensusPolicy | None:
    return load_model(path, SvConsensusPolicy) if path is not None and path.is_file() else None


def _sv_evidence_policy(path: Path | None) -> SvEvidencePolicy | None:
    return load_model(path, SvEvidencePolicy) if path is not None and path.is_file() else None


def _interval_resource(
    path: Path | None, lock_path: Path | None
) -> tuple[Path, IntervalResourceLock] | None:
    if path is None and lock_path is None:
        return None
    if path is None or lock_path is None:
        raise SystemExit("an interval resource requires both data and lock paths")
    return path, load_model(lock_path, IntervalResourceLock)


def _interval_resources(
    pairs: Sequence[Sequence[Path]],
) -> tuple[tuple[Path, IntervalResourceLock], ...]:
    resources: list[tuple[Path, IntervalResourceLock]] = []
    for pair in pairs:
        if len(pair) != 2:
            raise SystemExit("each interval resource requires DATA and LOCK")
        resolved = _interval_resource(pair[0], pair[1])
        assert resolved is not None
        resources.append(resolved)
    return tuple(resources)


def _aml_knowledge(
    path: Path | None, lock_path: Path | None
) -> tuple[Path, AmlKnowledgeLock] | None:
    if path is None and lock_path is None:
        return None
    if path is None or lock_path is None:
        raise SystemExit("AML knowledge requires both resource and lock paths")
    return path, load_model(lock_path, AmlKnowledgeLock)


def _target_coverage_policy(path: Path) -> TargetCoveragePolicy | None:
    return load_model(path, TargetCoveragePolicy) if path.is_file() else None


def _resolve_component_policies(selection: RunComponents, source: Path) -> RunComponents:
    """Resolve selection policy paths against the selection's packaged/repository root.

    Component documents deliberately use repository-root-relative paths such as
    ``configs/qc/...``. A packed Desktop runtime may be launched from any Windows working
    directory, so leaving those paths relative would make the selected policy depend on the
    operator's current directory. ``configs/components/default.yaml`` lives two directory
    levels below the root both in the repository and in ``share/ontseq`` after packaging.
    """
    root = source.resolve().parents[2]
    updated = dict(selection.components)
    for stage, choice in selection.components.items():
        if not choice.policy:
            continue
        policy = Path(choice.policy)
        if policy.is_absolute():
            continue
        updated[stage] = choice.model_copy(update={"policy": str((root / policy).resolve())})
    return selection.model_copy(update={"components": updated})


def _components(args: argparse.Namespace) -> RunComponents | None:
    """Resolve the component selection for this run, if the operator asked for one.

    ``--without`` is applied on top of the file rather than instead of it, so switching a
    stage off for one run does not require editing, copying or forking a selection.
    """
    selection: RunComponents | None = None
    path: Path | None = getattr(args, "components", None)
    if path is not None:
        if not path.is_file():
            raise SystemExit(f"component selection not found: {path}")
        selection = _resolve_component_policies(load_model(path, RunComponents), path)
    without = [StageId(name) for name in getattr(args, "without", []) or []]
    if without:
        base = selection or RunComponents(
            selection_id="command-line-only",
            status="technical_defaults_only",
            note="Created implicitly by --without; no versions are pinned.",
        )
        selection = base.without(*without)
    return selection


def _add_cnv_options(parser: argparse.ArgumentParser) -> None:
    root = _repo_root()
    parser.add_argument(
        "--cnv-policy",
        type=Path,
        default=root / "configs/cnv/qdnaseq_ace.technical.yaml",
    )
    parser.add_argument("--qdnaseq-rscript", default="Rscript")
    parser.add_argument(
        "--qdnaseq-script",
        type=Path,
        default=root / "scripts/run_qdnaseq_ace.R",
    )


def _register_cnv(args: argparse.Namespace, selection: RunComponents | None) -> None:
    """Install the QDNAseq/ACE lane unless this run deselected it.

    The lane still arrives by registration rather than as a first-class member of the
    graph, which remains the outstanding architectural debt. Gating it on the selection at
    least means a run that switched CNV off does not silently get it anyway.
    """
    from .cnv.extension import QDNAseqExtensionSettings, register_qdnaseq_extension
    from .cnv.qdnaseq import QDNAseqPolicy

    choice = selection.choice_for(StageId.CNV) if selection is not None else None
    if choice is not None and not choice.enabled:
        return
    if choice is not None and choice.policy:
        args.cnv_policy = Path(choice.policy)

    if args.cnv_policy.is_file():
        policy = load_model(args.cnv_policy, QDNAseqPolicy)
    else:
        policy = QDNAseqPolicy(
            profile_id="qdnaseq-ace-multibin-v1",
            cytoband_affected_fraction=0.66,
            note="Built-in fallback matching configs/cnv/qdnaseq_ace.technical.yaml",
        )
    register_qdnaseq_extension(
        QDNAseqExtensionSettings(
            policy=policy,
            rscript=args.qdnaseq_rscript,
            script=args.qdnaseq_script,
        )
    )


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
        "--cutesv-policy",
        type=Path,
        default=Path("configs/sv/cutesv.conservative.technical.yaml"),
    )
    parser.add_argument(
        "--sv-consensus-policy",
        type=Path,
        default=Path("configs/sv/sniffles2_cutesv.consensus.technical.yaml"),
    )
    parser.add_argument(
        "--sv-evidence-policy",
        type=Path,
        default=Path("configs/sv/evidence-priority.technical.yaml"),
    )
    parser.add_argument("--gene-annotation", type=Path)
    parser.add_argument("--gene-annotation-lock", type=Path)
    parser.add_argument("--cytoband-annotation", type=Path)
    parser.add_argument("--cytoband-annotation-lock", type=Path)
    parser.add_argument(
        "--sv-context-resource",
        type=Path,
        nargs=2,
        action="append",
        metavar=("DATA", "LOCK"),
        default=[],
    )
    parser.add_argument(
        "--aml-knowledge",
        type=Path,
        default=Path("configs/knowledge/aml_rearrangements.v0.1.json"),
    )
    parser.add_argument(
        "--aml-knowledge-lock",
        type=Path,
        default=Path("configs/knowledge/aml_rearrangements.v0.1.lock.json"),
    )
    parser.add_argument(
        "--sv-minimum-mean-depth",
        type=float,
        default=10.0,
        help="Unvalidated technical depth floor used only for observability labels",
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
    parser.add_argument(
        "--target-coverage-policy",
        type=Path,
        default=Path("configs/qc/adaptive_target_coverage.technical.yaml"),
    )
    parser.add_argument(
        "--components",
        type=Path,
        help="Component selection for this run: which provider and version runs each stage",
    )
    parser.add_argument(
        "--without",
        action="append",
        choices=sorted(item.value for item in SELECTABLE_STAGES),
        default=[],
        help="Switch a stage off for this run; repeatable",
    )
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--cramino", default="cramino")
    parser.add_argument("--sniffles", default="sniffles")
    parser.add_argument("--cutesv", default="cuteSV")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--mosdepth", default="mosdepth")
    parser.add_argument("--dorado", default="dorado")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ontseq",
        description="ONTSeq end-to-end execution and operational commands",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Execute one sample into a resumable run envelope")
    _add_execution_options(run, include_qc=True)
    _add_cnv_options(run)
    run.add_argument("--threads", type=int, default=4)
    run.add_argument("--git-commit", default="UNKNOWN")
    run.add_argument("--force", action="store_true")

    analyze = sub.add_parser(
        "analyze", help="Analyze one indexed GRCh38 BAM using an installed profile"
    )
    analyze.add_argument("bam", type=Path)
    analyze.add_argument(
        "--profile",
        required=True,
        choices=("AML_LCWGS_GRCh38", "AML_AS_111_GRCh38"),
    )
    analyze.add_argument("--resource-root", type=Path)
    analyze.add_argument("--config-root", type=Path)
    analyze.add_argument("--output-dir", type=Path, default=Path("results/runs"))
    analyze.add_argument("--sample-id")
    analyze.add_argument("--run-id")
    analyze.add_argument("--threads", type=int, default=4)
    analyze.add_argument("--git-commit", default="UNKNOWN")
    analyze.add_argument("--force", action="store_true")
    analyze.add_argument("--samtools", default="samtools")
    analyze.add_argument("--cramino", default="cramino")
    analyze.add_argument("--sniffles", default="sniffles")
    analyze.add_argument("--cutesv", default="cuteSV")
    analyze.add_argument("--minimap2", default="minimap2")
    analyze.add_argument("--mosdepth", default="mosdepth")
    analyze.add_argument("--dorado", default="dorado")
    _add_cnv_options(analyze)

    add_references_parser(sub)

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
    srv.add_argument("--reference-lock", type=Path)
    srv.add_argument(
        "--resource-root",
        type=Path,
        help="Manifested GRCh38 resource root used by profile-based Desktop runs",
    )
    srv.add_argument("--allow-root", type=Path, action="append", required=True, dest="allow_roots")
    srv.add_argument("--output-dir", type=Path, default=Path("results/runs"))
    srv.add_argument("--qc-policy", type=Path, default=Path("configs/qc/defaults.yaml"))
    srv.add_argument(
        "--sniffles-policy",
        type=Path,
        default=Path("configs/sv/sniffles2.conservative.technical.yaml"),
    )
    srv.add_argument(
        "--cutesv-policy",
        type=Path,
        default=Path("configs/sv/cutesv.conservative.technical.yaml"),
    )
    srv.add_argument(
        "--sv-consensus-policy",
        type=Path,
        default=Path("configs/sv/sniffles2_cutesv.consensus.technical.yaml"),
    )
    srv.add_argument(
        "--sv-evidence-policy",
        type=Path,
        default=Path("configs/sv/evidence-priority.technical.yaml"),
    )
    srv.add_argument("--reference-fasta", type=Path)
    srv.add_argument("--gene-annotation", type=Path)
    srv.add_argument("--gene-annotation-lock", type=Path)
    srv.add_argument("--cytoband-annotation", type=Path)
    srv.add_argument("--cytoband-annotation-lock", type=Path)
    srv.add_argument(
        "--sv-context-resource",
        type=Path,
        nargs=2,
        action="append",
        metavar=("DATA", "LOCK"),
        default=[],
    )
    srv.add_argument(
        "--aml-knowledge",
        type=Path,
        default=Path("configs/knowledge/aml_rearrangements.v0.1.json"),
    )
    srv.add_argument(
        "--aml-knowledge-lock",
        type=Path,
        default=Path("configs/knowledge/aml_rearrangements.v0.1.lock.json"),
    )
    srv.add_argument("--sv-minimum-mean-depth", type=float, default=10.0)
    srv.add_argument("--cutesv", default="cuteSV")
    srv.add_argument(
        "--target-coverage-policy",
        type=Path,
        default=Path("configs/qc/adaptive_target_coverage.technical.yaml"),
    )
    srv.add_argument(
        "--components",
        type=Path,
        help="Component selection applied to every run started by this local service",
    )
    _add_cnv_options(srv)
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
        "--cutesv-policy",
        type=Path,
        default=Path("configs/sv/cutesv.conservative.technical.yaml"),
    )
    watcher.add_argument(
        "--sv-consensus-policy",
        type=Path,
        default=Path("configs/sv/sniffles2_cutesv.consensus.technical.yaml"),
    )
    watcher.add_argument(
        "--sv-evidence-policy",
        type=Path,
        default=Path("configs/sv/evidence-priority.technical.yaml"),
    )
    watcher.add_argument(
        "--target-coverage-policy",
        type=Path,
        default=Path("configs/qc/adaptive_target_coverage.technical.yaml"),
    )
    watcher.add_argument("--gene-annotation", type=Path)
    watcher.add_argument("--gene-annotation-lock", type=Path)
    watcher.add_argument("--cytoband-annotation", type=Path)
    watcher.add_argument("--cytoband-annotation-lock", type=Path)
    watcher.add_argument(
        "--sv-context-resource",
        type=Path,
        nargs=2,
        action="append",
        metavar=("DATA", "LOCK"),
        default=[],
    )
    watcher.add_argument(
        "--aml-knowledge",
        type=Path,
        default=Path("configs/knowledge/aml_rearrangements.v0.1.json"),
    )
    watcher.add_argument(
        "--aml-knowledge-lock",
        type=Path,
        default=Path("configs/knowledge/aml_rearrangements.v0.1.lock.json"),
    )
    watcher.add_argument("--sv-minimum-mean-depth", type=float, default=10.0)
    watcher.add_argument(
        "--alignment-policy",
        type=Path,
        default=Path("configs/alignment/minimap2.ont.technical.yaml"),
    )
    _add_cnv_options(watcher)
    watcher.add_argument("--reference-fasta", type=Path)
    watcher.add_argument("--cutesv", default="cuteSV")
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
        "cutesv": args.cutesv,
        "minimap2": args.minimap2,
        "mosdepth": getattr(args, "mosdepth", "mosdepth"),
        "dorado": args.dorado,
    }


def _print_pass(result: PassResult) -> None:
    for attempt in result.attempted:
        print(f"  {attempt.name:<28} {attempt.outcome.value.upper():<10} {attempt.detail}")
    for name, reason in result.skipped:
        print(f"  {name:<28} {'skipped':<10} {reason}")


def main() -> None:
    args = _parser().parse_args()
    # Preflight resolves the selection too. Its whole value is that it agrees with the
    # run: checking the default policies while `ontseq run` would use the ones a component
    # selection names is how a preflight clears a run that then fails on what it checked.
    selection = _components(args) if args.command in {"run", "serve", "preflight"} else None
    if args.command in {"run", "analyze", "serve", "watch"}:
        _register_cnv(args, selection)
    try:
        if handle_references_command(args):
            return

        if args.command == "analyze":
            config = build_profile_run_configuration(
                AnalyzeSettings(
                    bam=args.bam,
                    profile_id=args.profile,
                    resource_root=args.resource_root,
                    output_dir=args.output_dir,
                    configuration_root=args.config_root,
                    sample_id=args.sample_id,
                    run_id=args.run_id,
                    pipeline_version=__version__,
                    git_commit=args.git_commit,
                    threads=args.threads,
                    force=args.force,
                    executables=_executables(args),
                )
            )
            print(
                f"profile: {args.profile}; detected build: "
                f"{config.manifest.assay.genome_build.value}"
            )
            run_report, release = run_pipeline(config)
            for stage in run_report.stages:
                marker = "resumed" if stage.resumed else stage.status.value
                print(f"  {stage.stage.value:<16} {marker:<10} {stage.reason}")
            outcome = "PASS" if run_report.passed else "FAIL"
            print(f"verdict: {outcome} - {run_report.verdict_reason}")
            if run_report.unverified_stages:
                names = ", ".join(item.value for item in run_report.unverified_stages)
                print(f"UNVERIFIED ADAPTERS COMPLETED: {names}")
            if release is not None:
                print(f"release bundle: {len(release.artifacts)} artifact(s), unsigned")
            if not run_report.passed:
                raise SystemExit(2)

        elif args.command == "run":
            config = RunConfiguration(
                manifest=load_model(args.manifest, SampleManifest),
                reference_lock=load_model(args.reference_lock, ReferenceLock),
                output_base=args.output_dir,
                run_id=args.run_id,
                pipeline_version=__version__,
                git_commit=args.git_commit,
                qc_policy=load_model(args.qc_policy, QCPolicy),
                sniffles_policy=_sniffles_policy(
                    _selected_policy(selection, StageId.SV, args.sniffles_policy)
                ),
                cutesv_policy=_cutesv_policy_for_run(
                    args.cutesv_policy,
                    args.reference_fasta,
                ),
                sv_consensus_policy=_sv_consensus_policy(args.sv_consensus_policy),
                sv_evidence_policy=_sv_evidence_policy(args.sv_evidence_policy),
                gene_annotation=_interval_resource(args.gene_annotation, args.gene_annotation_lock),
                cytoband_annotation=_interval_resource(
                    args.cytoband_annotation, args.cytoband_annotation_lock
                ),
                sv_context_resources=_interval_resources(args.sv_context_resource),
                aml_knowledge=_aml_knowledge(args.aml_knowledge, args.aml_knowledge_lock),
                sv_minimum_mean_depth=args.sv_minimum_mean_depth,
                target_coverage_policy=_target_coverage_policy(
                    _selected_policy(
                        selection, StageId.TARGET_COVERAGE, args.target_coverage_policy
                    )
                ),
                alignment_policy=_alignment_policy(
                    _selected_policy(selection, StageId.ALIGN, args.alignment_policy)
                ),
                basecall_policy=_basecall_policy(
                    _selected_policy(selection, StageId.BASECALL, args.basecall_policy)
                ),
                components=selection,
                reference_fasta=args.reference_fasta,
                pod5_directory=args.pod5_dir,
                threads=args.threads,
                executables=_executables(args),
                force=args.force,
            )
            if selection is not None:
                print(f"component selection: {selection.selection_id}")
                for line in selection.summary():
                    print(f"  {line}")
            run_report, release = run_pipeline(config)
            for stage in run_report.stages:
                marker = "resumed" if stage.resumed else stage.status.value
                print(f"  {stage.stage.value:<16} {marker:<10} {stage.reason}")
            outcome = "PASS" if run_report.passed else "FAIL"
            print(f"verdict: {outcome} - {run_report.verdict_reason}")
            if run_report.unverified_stages:
                names = ", ".join(item.value for item in run_report.unverified_stages)
                print(f"UNVERIFIED ADAPTERS COMPLETED: {names}")
            if release is not None:
                print(f"release bundle: {len(release.artifacts)} artifact(s), unsigned")
            if not run_report.passed:
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
                alignment_policy=_alignment_policy(
                    _selected_policy(selection, StageId.ALIGN, args.alignment_policy)
                ),
                basecall_policy=_basecall_policy(
                    _selected_policy(selection, StageId.BASECALL, args.basecall_policy)
                ),
                sniffles_policy=_sniffles_policy(
                    _selected_policy(selection, StageId.SV, args.sniffles_policy)
                ),
                cutesv_policy=_cutesv_policy_for_run(
                    args.cutesv_policy,
                    args.reference_fasta,
                ),
                target_coverage_policy=_target_coverage_policy(
                    _selected_policy(
                        selection, StageId.TARGET_COVERAGE, args.target_coverage_policy
                    )
                ),
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
            if args.reference_lock is None and args.resource_root is None:
                raise ValueError("serve requires --resource-root or legacy --reference-lock")
            serve(
                ServiceConfig(
                    reference_lock=args.reference_lock,
                    output_dir=args.output_dir,
                    allowed_roots=list(args.allow_roots),
                    qc_policy=args.qc_policy,
                    sniffles_policy=_selected_policy(selection, StageId.SV, args.sniffles_policy),
                    target_coverage_policy=_selected_policy(
                        selection, StageId.TARGET_COVERAGE, args.target_coverage_policy
                    ),
                    components=selection,
                    cutesv_policy=args.cutesv_policy,
                    sv_consensus_policy=args.sv_consensus_policy,
                    sv_evidence_policy=args.sv_evidence_policy,
                    reference_fasta=args.reference_fasta,
                    gene_annotation=_interval_resource(
                        args.gene_annotation, args.gene_annotation_lock
                    ),
                    cytoband_annotation=_interval_resource(
                        args.cytoband_annotation, args.cytoband_annotation_lock
                    ),
                    sv_context_resources=_interval_resources(args.sv_context_resource),
                    aml_knowledge=(
                        _aml_knowledge(args.aml_knowledge, args.aml_knowledge_lock)
                        if args.reference_lock is not None
                        else None
                    ),
                    sv_minimum_mean_depth=args.sv_minimum_mean_depth,
                    cutesv_executable=args.cutesv,
                    port=args.port,
                    threads=args.threads,
                    resource_root=args.resource_root,
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
                review_report = inspect_review(args.envelope)
                if args.as_json:
                    print(render_review_json(review_report), end="")
                else:
                    print(render_review_text(review_report, verbose=args.verbose))
                code = review_exit_code(
                    review_report.state,
                    reviewers=len(review_report.reviewers),
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
                cutesv_policy=args.cutesv_policy,
                sv_consensus_policy=args.sv_consensus_policy,
                sv_evidence_policy=args.sv_evidence_policy,
                target_coverage_policy=args.target_coverage_policy,
                gene_annotation=_interval_resource(args.gene_annotation, args.gene_annotation_lock),
                cytoband_annotation=_interval_resource(
                    args.cytoband_annotation, args.cytoband_annotation_lock
                ),
                sv_context_resources=_interval_resources(args.sv_context_resource),
                aml_knowledge=_aml_knowledge(args.aml_knowledge, args.aml_knowledge_lock),
                sv_minimum_mean_depth=args.sv_minimum_mean_depth,
                alignment_policy=args.alignment_policy,
                reference_fasta=args.reference_fasta,
                run_id_prefix=args.run_id_prefix,
                ready_marker=args.ready_marker,
                pod5_subdirectory=args.pod5_subdir,
                quiet_seconds=args.quiet_seconds,
                threads=args.threads,
                git_commit=args.git_commit,
                retry_failed=args.retry_failed,
                executables={"cutesv": args.cutesv},
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
    except (
        KeyError,
        OSError,
        RuntimeError,
        ValueError,
        ValidationError,
        ToolExecutionError,
    ) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

"""End-to-end execution: plan the stages, run them, record everything, resume safely.

The runner owns three responsibilities and nothing else:

**Planning.** Which stages apply follows from the declared input kind
(:mod:`ontseq_platform.pipeline.stages`), never from what happens to exist on disk.

**Execution with honest bookkeeping.** Every stage produces a record: what ran, under
which tool versions, producing which checksummed artifacts, and why it did not run when it
did not. A stage that raises becomes ``FAILED`` with the reason attached; its dependents
become ``NOT_RUN``. The run report is rewritten after every stage, so an interrupted run
leaves a truthful partial record rather than nothing.

**Content-addressed resume.** A completed stage is reused only when its input signature is
unchanged *and* every artifact it claimed still verifies byte for byte. Anything else
re-runs. Resume is an optimisation; it must never be the reason two incompatible results
end up in one envelope.

Adapters are called, never reimplemented. The scientific behaviour lives in
``bam_intake``, ``qc``, ``sniffles``, ``align``, ``basecall`` and ``mvp``; this module only
sequences them.
"""

from __future__ import annotations

import traceback
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..align import AlignmentInputs, AlignmentPolicy, run_alignment
from ..bam_intake import AlignedBamInspector
from ..basecall import BasecallInputs, BasecallPolicy, run_basecalling
from ..execution import StreamingCommandRunner, SubprocessRunner, ToolExecutionError
from ..models import (
    AlignedBamIntakeReport,
    CheckStatus,
    CraminoQCReport,
    InputKind,
    InputSpec,
    ModuleRunStatus,
    QCPolicy,
    ReferenceLock,
    SampleManifest,
    SnifflesCallReport,
    SnifflesPolicy,
    ToolRecord,
    ValidationCheck,
    Verdict,
)
from ..mvp import assemble_aligned_bam_mvp
from ..qc import run_cramino_qc
from ..reference import contig_signature, reference_lock_signature
from ..report import render_html
from ..sniffles import run_sniffles
from ..workbook import render_workbook
from .envelope import Artifact, RunEnvelope, sha256_file, stage_signature
from .lock import run_lock
from .review import RELEASE_RELATIVE, REVIEW_LOG, ReviewError, ReviewState
from .review import current_state as review_state
from .review import read_log as read_review_log
from .stages import (
    SPEC_BY_STAGE,
    InputKindName,
    StageId,
    StageOutcome,
    blocking_dependency,
    planned_stages,
    summarize,
)
from .state import ArtifactRecord, ReleaseBundle, RunReport, StageRecord

#: Where each stage's primary artifacts live inside the envelope.
BASECALL_BAM = "alignment/{sample}.unaligned.bam"
BASECALL_REPORT = "provenance/basecall.json"
ALIGNED_BAM = "alignment/{sample}.bam"
ALIGNED_BAI = "alignment/{sample}.bam.bai"
ALIGN_REPORT = "provenance/alignment.json"
INTAKE_REPORT = "manifest/intake.json"
QC_REPORT = "qc/cramino.json"
SV_VCF = "evidence/sv/{sample}.sniffles.vcf"
SV_REPORT = "evidence/sv/{sample}.sniffles.json"
RESULT_JSON = "normalized/{sample}.result.json"
REPORT_HTML = "reports/{sample}.report.html"
REPORT_XLSX = "reports/{sample}.results.xlsx"
RUN_REPORT = "provenance/run.json"
RELEASE_JSON = "release/release.json"
RELEASE_CHECKSUMS = "release/checksums.sha256"

STANDING_LIMITATIONS = (
    "Research use only. No output of this run may inform diagnosis or treatment.",
    "A stage recorded as NOT_RUN, FAILED or NO_CALL is not a negative biological finding.",
    "Tool versions are locked for reproducibility; none of the thresholds involved is a "
    "validated clinical limit.",
)


class StageFailure(RuntimeError):
    """Raised inside a stage to fail it with a specific, readable reason."""


@dataclass
class RunConfiguration:
    """Everything one run needs, resolved before execution starts."""

    manifest: SampleManifest
    reference_lock: ReferenceLock
    output_base: Path
    run_id: str
    pipeline_version: str
    git_commit: str
    qc_policy: QCPolicy
    sniffles_policy: SnifflesPolicy | None = None
    alignment_policy: AlignmentPolicy | None = None
    basecall_policy: BasecallPolicy | None = None
    reference_fasta: Path | None = None
    pod5_directory: Path | None = None
    threads: int = 4
    executables: Mapping[str, str] = field(
        default_factory=lambda: {
            "samtools": "samtools",
            "cramino": "cramino",
            "sniffles": "sniffles",
            "minimap2": "minimap2",
            "dorado": "dorado",
        }
    )
    #: Ignore any previous run state and execute every stage again.
    force: bool = False

    def executable(self, name: str) -> str:
        return self.executables.get(name, name)


@dataclass
class RunContext:
    """Mutable state threaded through the stages of one run."""

    config: RunConfiguration
    envelope: RunEnvelope
    runner: StreamingCommandRunner
    #: The manifest as it currently stands. Alignment rewrites its input to the BAM the
    #: pipeline just produced, so downstream adapters need no special casing.
    manifest: SampleManifest
    artifacts: dict[StageId, list[Artifact]] = field(default_factory=dict)

    @property
    def sample_id(self) -> str:
        return self.config.manifest.sample_id

    def path(self, template: str) -> str:
        return template.format(sample=self.sample_id)

    def upstream(self, stage: StageId, input_kind: InputKindName) -> list[Artifact]:
        collected: list[Artifact] = []
        for dependency in SPEC_BY_STAGE[stage].depends_on:
            collected.extend(self.artifacts.get(dependency, []))
        return collected


@dataclass(frozen=True)
class StagePlan:
    """What a stage will do, resolved before the resume decision."""

    parameters: dict[str, object]
    tool_versions: dict[str, str]
    external_inputs: tuple[tuple[str, str], ...] = ()


@dataclass(frozen=True)
class StageResult:
    """What a stage did."""

    status: ModuleRunStatus
    reason: str
    outputs: list[Artifact] = field(default_factory=list)
    tools: list[ToolRecord] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class StageImplementation:
    plan: Callable[[RunContext], StagePlan]
    execute: Callable[[RunContext, StagePlan], StageResult]
    #: Re-point the context at what the stage produced. Runs after a stage completes *and*
    #: after it resumes, because a resumed stage produced its artifacts just as surely as
    #: one that just ran. Putting this inside ``execute`` would mean a resumed alignment
    #: left the manifest pointing at the unaligned input, and every downstream stage would
    #: then either re-run against the wrong file or fail outright. It receives the recorded
    #: artifacts rather than re-reading the envelope, so adopting a multi-gigabyte BAM does
    #: not cost a second checksum pass over it.
    settle: Callable[[RunContext, Sequence[Artifact]], None] | None = None


def _probe(
    runner: StreamingCommandRunner, executable: str, argv: Sequence[str], *, tool: str
) -> str:
    import re

    result = runner.run(list(argv), timeout_seconds=60)
    if result.returncode != 0:
        raise StageFailure(f"{tool} version probe returned exit code {result.returncode}")
    match = re.search(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)", f"{result.stdout}\n{result.stderr}")
    if not match:
        raise StageFailure(f"could not determine the {tool} version")
    return match.group(1)


def _stable_digest(path: Path) -> tuple[str, bool]:
    """Hash a regular file and report whether its size/mtime stayed fixed while reading."""

    if not path.is_file():
        raise StageFailure("required external input is missing")
    before = path.stat()
    digest = sha256_file(path)
    after = path.stat()
    stable = (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    return digest, stable


def _external_fingerprint(path: Path, *, label: str | None = None) -> tuple[str, str]:
    """Fingerprint an input from outside the envelope by name and content."""
    digest, stable = _stable_digest(path)
    if not stable:
        raise StageFailure(
            f"{label or 'required external input'} changed while it was being fingerprinted"
        )
    return (label or path.name, digest)


# --------------------------------------------------------------------------------------
# Stage implementations
# --------------------------------------------------------------------------------------


def _basecall_plan(ctx: RunContext) -> StagePlan:
    policy = ctx.config.basecall_policy
    if policy is None:
        raise StageFailure("a POD5 run requires a basecalling policy")
    if ctx.config.pod5_directory is None:
        raise StageFailure("a POD5 run requires --pod5-dir")
    version = _probe(
        ctx.runner,
        ctx.config.executable("dorado"),
        [ctx.config.executable("dorado"), "--version"],
        tool="dorado",
    )
    return StagePlan(
        parameters={
            "model": policy.model,
            "modified_bases": list(policy.modified_bases),
            "device": policy.device,
            "minimum_qscore": policy.minimum_qscore,
        },
        tool_versions={"dorado": version},
        external_inputs=(("pod5_directory", str(ctx.config.pod5_directory.name)),),
    )


def _recorded(outputs: Sequence[Artifact], relative_path: str, *, stage: str) -> Artifact:
    for artifact in outputs:
        if artifact.relative_path == relative_path:
            return artifact
    raise StageFailure(f"the {stage} stage recorded no artifact at {relative_path}")


def _basecall_settle(ctx: RunContext, outputs: Sequence[Artifact]) -> None:
    """Point the manifest at the unaligned BAM basecalling produced."""
    relative = ctx.path(BASECALL_BAM)
    bam = _recorded(outputs, relative, stage="basecall")
    ctx.manifest = ctx.manifest.model_copy(
        update={
            "input": InputSpec(
                kind=InputKind.UNALIGNED_BAM,
                path=str(ctx.envelope.path(relative)),
                sha256=bam.sha256,
            )
        }
    )


def _basecall_execute(ctx: RunContext, plan: StagePlan) -> StageResult:
    policy = ctx.config.basecall_policy
    assert policy is not None and ctx.config.pod5_directory is not None
    target = ctx.envelope.path(ctx.path(BASECALL_BAM))
    target.unlink(missing_ok=True)
    report = run_basecalling(
        BasecallInputs(pod5_directory=ctx.config.pod5_directory),
        policy,
        sample_id=ctx.sample_id,
        output_bam=target,
        runner=ctx.runner,
        dorado=ctx.config.executable("dorado"),
    )
    bam = ctx.envelope.fingerprint(ctx.path(BASECALL_BAM))
    record = ctx.envelope.atomic_write_text(
        BASECALL_REPORT, report.model_dump_json(indent=2) + "\n"
    )
    return StageResult(
        status=ModuleRunStatus.COMPLETED,
        reason=f"Basecalled {report.pod5_file_count} POD5 file(s) into an unaligned BAM.",
        outputs=[bam, record],
        tools=[report.tool],
        warnings=report.warnings,
        limitations=report.limitations,
    )


def _align_plan(ctx: RunContext) -> StagePlan:
    policy = ctx.config.alignment_policy
    if policy is None:
        raise StageFailure("an unaligned input requires an alignment policy")
    if ctx.config.reference_fasta is None:
        raise StageFailure("alignment requires --reference-fasta")
    reference = _external_fingerprint(ctx.config.reference_fasta)
    minimap2 = ctx.config.executable("minimap2")
    samtools = ctx.config.executable("samtools")
    return StagePlan(
        parameters={
            "preset": policy.preset,
            "md_tag": policy.emit_md_tag,
            "soft_clip_supplementary": policy.soft_clip_supplementary,
            "modified_base_tags": policy.preserve_modified_base_tags,
            "threads": ctx.config.threads,
        },
        tool_versions={
            "minimap2": _probe(ctx.runner, minimap2, [minimap2, "--version"], tool="minimap2"),
            "samtools": _probe(ctx.runner, samtools, [samtools, "--version"], tool="samtools"),
        },
        external_inputs=(reference, _external_fingerprint(Path(ctx.manifest.input.path))),
    )


def _align_settle(ctx: RunContext, outputs: Sequence[Artifact]) -> None:
    """Point the manifest at the aligned BAM, so downstream adapters need no special case."""
    bam = _recorded(outputs, ctx.path(ALIGNED_BAM), stage="align")
    _recorded(outputs, ctx.path(ALIGNED_BAI), stage="align")
    bam_path = ctx.envelope.path(bam.relative_path)
    ctx.manifest = ctx.manifest.model_copy(
        update={
            "input": InputSpec(
                kind=InputKind.ALIGNED_BAM,
                path=str(bam_path),
                index_path=f"{bam_path}.bai",
                sha256=bam.sha256,
            )
        }
    )


def _align_execute(ctx: RunContext, plan: StagePlan) -> StageResult:
    policy = ctx.config.alignment_policy
    assert policy is not None and ctx.config.reference_fasta is not None
    bam_path = ctx.envelope.path(ctx.path(ALIGNED_BAM))
    for stale in (bam_path, Path(f"{bam_path}.bai")):
        stale.unlink(missing_ok=True)
    report = run_alignment(
        AlignmentInputs(
            unaligned_bam=Path(ctx.manifest.input.path),
            reference_fasta=ctx.config.reference_fasta,
        ),
        policy,
        sample_id=ctx.sample_id,
        genome_build=ctx.manifest.assay.genome_build,
        reference_id=ctx.manifest.assay.reference_id,
        scratch_dir=ctx.envelope.path("work"),
        output_bam=bam_path,
        runner=ctx.runner,
        minimap2=ctx.config.executable("minimap2"),
        samtools=ctx.config.executable("samtools"),
        threads=ctx.config.threads,
    )
    bam = ctx.envelope.fingerprint(ctx.path(ALIGNED_BAM))
    bai = ctx.envelope.fingerprint(ctx.path(ALIGNED_BAI))
    record = ctx.envelope.atomic_write_text(ALIGN_REPORT, report.model_dump_json(indent=2) + "\n")
    return StageResult(
        status=ModuleRunStatus.COMPLETED,
        reason="Reads aligned to the locked reference and coordinate sorted.",
        outputs=[bam, bai, record],
        tools=report.tools,
        warnings=report.warnings,
        limitations=report.limitations,
    )


def _intake_plan(ctx: RunContext) -> StagePlan:
    samtools = ctx.config.executable("samtools")
    reference = ctx.config.reference_lock
    index_path = ctx.manifest.input.index_path
    if index_path is None:
        raise StageFailure("aligned-BAM intake requires a BAM index")
    return StagePlan(
        parameters={
            "manifest_reference_id": ctx.manifest.assay.reference_id,
            "manifest_genome_build": ctx.manifest.assay.genome_build.value,
            "manifest_input_sha256": ctx.manifest.input.sha256,
            "bam_extension": Path(ctx.manifest.input.path).suffix.lower(),
            "reference_id": reference.reference_id,
            "reference_genome_build": reference.genome_build.value,
            "reference_source_fai_sha256": reference.source_fai_sha256,
            "reference_dictionary_sha256": contig_signature(
                (item.name, item.length) for item in reference.contigs
            ),
            "reference_allow_extra_contigs": reference.allow_extra_contigs,
            "reference_lock_sha256": reference_lock_signature(reference),
        },
        tool_versions={
            "samtools": _probe(ctx.runner, samtools, [samtools, "--version"], tool="samtools")
        },
        external_inputs=(
            _external_fingerprint(Path(ctx.manifest.input.path), label="aligned_bam"),
            _external_fingerprint(Path(index_path), label="bam_index"),
        ),
    )


def _intake_execute(ctx: RunContext, plan: StagePlan) -> StageResult:
    report = AlignedBamInspector(
        runner=ctx.runner, samtools=ctx.config.executable("samtools")
    ).inspect(ctx.manifest, ctx.config.reference_lock, include_checksums=False)

    planned = dict(plan.external_inputs)

    def final_digest(path: Path) -> tuple[str | None, bool]:
        try:
            digest, stable = _stable_digest(path)
        except (OSError, StageFailure):
            return None, False
        return digest, stable

    index_path = ctx.manifest.input.index_path
    assert index_path is not None
    bam_digest, bam_stable = final_digest(Path(ctx.manifest.input.path))
    index_digest, index_stable = final_digest(Path(index_path))
    bam_matches_plan = bam_stable and bam_digest == planned.get("aligned_bam")
    index_matches_plan = index_stable and index_digest == planned.get("bam_index")
    inputs_stable = bam_matches_plan and index_matches_plan

    input_fingerprint = report.input_fingerprint
    if input_fingerprint is not None and bam_digest is not None:
        input_fingerprint = input_fingerprint.model_copy(update={"sha256": bam_digest})
    index_fingerprint = report.index_fingerprint
    if index_fingerprint is not None and index_digest is not None:
        index_fingerprint = index_fingerprint.model_copy(update={"sha256": index_digest})
    stability_check = ValidationCheck(
        name="input_stability",
        status=CheckStatus.PASS if inputs_stable else CheckStatus.FAIL,
        message=(
            "BAM and index stayed identical to the planned inputs"
            if inputs_stable
            else "BAM or index changed between intake planning and verification"
        ),
        details={
            "bam_matches_plan": bam_matches_plan,
            "index_matches_plan": index_matches_plan,
        },
    )
    report = report.model_copy(
        update={
            "input_fingerprint": input_fingerprint,
            "index_fingerprint": index_fingerprint,
            "checks": [*report.checks, stability_check],
            "verdict": report.verdict if inputs_stable else Verdict.FAIL,
        }
    )
    artifact = ctx.envelope.atomic_write_text(
        INTAKE_REPORT, report.model_dump_json(indent=2) + "\n"
    )
    if report.verdict == Verdict.FAIL:
        failed = [item for item in report.checks if item.status.value == "FAIL"]
        detail = "; ".join(f"{item.name}: {item.message}" for item in failed) or "unspecified check"
        return StageResult(
            status=ModuleRunStatus.FAILED,
            reason=(f"Aligned-BAM intake failed: {detail}. Full details: {INTAKE_REPORT}"),
            tools=[report.tool] if report.tool else [],
            limitations=[
                *report.limitations,
                f"Failure diagnostic: {INTAKE_REPORT} sha256:{artifact.sha256}",
            ],
        )
    return StageResult(
        status=ModuleRunStatus.COMPLETED,
        reason=f"Aligned-BAM intake gate returned {report.verdict.value}.",
        outputs=[artifact],
        tools=[report.tool] if report.tool else [],
        limitations=report.limitations,
    )


def _qc_plan(ctx: RunContext) -> StagePlan:
    cramino = ctx.config.executable("cramino")
    return StagePlan(
        parameters={"threads": ctx.config.threads, "policy": ctx.config.qc_policy.status},
        tool_versions={
            "cramino": _probe(ctx.runner, cramino, [cramino, "--version"], tool="cramino")
        },
        external_inputs=(_external_fingerprint(Path(ctx.manifest.input.path)),),
    )


def _qc_execute(ctx: RunContext, plan: StagePlan) -> StageResult:
    report = run_cramino_qc(
        ctx.manifest,
        ctx.config.qc_policy,
        runner=ctx.runner,
        cramino=ctx.config.executable("cramino"),
        threads=ctx.config.threads,
    )
    artifact = ctx.envelope.atomic_write_text(QC_REPORT, report.model_dump_json(indent=2) + "\n")
    if report.qc.verdict == Verdict.FAIL:
        raise StageFailure(
            "QC gate failed: " + ", ".join(report.qc.failed_gates or ["unspecified gate"])
        )
    return StageResult(
        status=ModuleRunStatus.COMPLETED,
        reason=f"Descriptive read QC returned {report.qc.verdict.value}.",
        outputs=[artifact],
        tools=[report.tool],
        warnings=report.qc.warnings,
        limitations=report.limitations,
    )


def _sv_plan(ctx: RunContext) -> StagePlan:
    policy = ctx.config.sniffles_policy
    if policy is None:
        raise StageFailure("structural-variant calling requires a Sniffles policy")
    sniffles = ctx.config.executable("sniffles")
    return StagePlan(
        parameters={
            "profile": policy.profile_id,
            "min_support": policy.min_support,
            "min_sv_length": policy.min_sv_length,
            "mapq": policy.mapq,
            "threads": ctx.config.threads,
        },
        tool_versions={
            "sniffles": _probe(ctx.runner, sniffles, [sniffles, "--version"], tool="sniffles")
        },
        external_inputs=(_external_fingerprint(Path(ctx.manifest.input.path)),),
    )


def _sv_execute(ctx: RunContext, plan: StagePlan) -> StageResult:
    policy = ctx.config.sniffles_policy
    assert policy is not None
    intake = AlignedBamIntakeReport.model_validate_json(
        ctx.envelope.path(INTAKE_REPORT).read_text(encoding="utf-8")
    )
    vcf_path = ctx.envelope.path(ctx.path(SV_VCF))
    vcf_path.unlink(missing_ok=True)
    report = run_sniffles(
        ctx.manifest,
        intake,
        policy,
        output_vcf=vcf_path,
        runner=ctx.runner,
        sniffles=ctx.config.executable("sniffles"),
        threads=ctx.config.threads,
    )
    outputs = [
        ctx.envelope.fingerprint(ctx.path(SV_VCF)),
        ctx.envelope.atomic_write_text(
            ctx.path(SV_REPORT), report.model_dump_json(indent=2) + "\n"
        ),
    ]
    reason = (
        f"Normalized {report.accepted_record_count} candidate SV event(s) from "
        f"{report.raw_record_count} record(s)."
        if report.status == ModuleRunStatus.COMPLETED
        else "No record passed the technical policy; this NO_CALL is not a biological negative."
    )
    return StageResult(
        status=report.status,
        reason=reason,
        outputs=outputs,
        tools=[report.tool],
        warnings=report.warnings,
        limitations=report.limitations,
    )


def _assemble_plan(ctx: RunContext) -> StagePlan:
    return StagePlan(
        parameters={
            "pipeline_version": ctx.config.pipeline_version,
            "git_commit": ctx.config.git_commit,
        },
        tool_versions={},
    )


def _assemble_execute(ctx: RunContext, plan: StagePlan) -> StageResult:
    intake = AlignedBamIntakeReport.model_validate_json(
        ctx.envelope.path(INTAKE_REPORT).read_text(encoding="utf-8")
    )
    qc = CraminoQCReport.model_validate_json(
        ctx.envelope.path(QC_REPORT).read_text(encoding="utf-8")
    )
    sv_path = ctx.envelope.path(ctx.path(SV_REPORT))
    sniffles = (
        SnifflesCallReport.model_validate_json(sv_path.read_text(encoding="utf-8"))
        if sv_path.is_file()
        else None
    )
    result = assemble_aligned_bam_mvp(
        ctx.manifest,
        intake,
        qc,
        pipeline_version=ctx.config.pipeline_version,
        git_commit=ctx.config.git_commit,
        sniffles_report=sniffles,
    )
    artifact = ctx.envelope.atomic_write_text(
        ctx.path(RESULT_JSON), result.model_dump_json(indent=2) + "\n"
    )
    return StageResult(
        status=ModuleRunStatus.COMPLETED,
        reason="Module outcomes assembled into the validated result contract.",
        outputs=[artifact],
        warnings=["Structural-variant evidence was omitted from the result."]
        if sniffles is None
        else [],
    )


def _report_plan(ctx: RunContext) -> StagePlan:
    return StagePlan(parameters={"formats": ["json", "html", "xlsx"]}, tool_versions={})


def _report_execute(ctx: RunContext, plan: StagePlan) -> StageResult:
    from ..models import PipelineResult

    result = PipelineResult.model_validate_json(
        ctx.envelope.path(ctx.path(RESULT_JSON)).read_text(encoding="utf-8")
    )
    render_html(result, ctx.envelope.path(ctx.path(REPORT_HTML)))
    render_workbook(result, ctx.envelope.path(ctx.path(REPORT_XLSX)))
    return StageResult(
        status=ModuleRunStatus.COMPLETED,
        reason="Reviewer artifacts rendered as HTML and Excel.",
        outputs=[
            ctx.envelope.fingerprint(ctx.path(REPORT_HTML)),
            ctx.envelope.fingerprint(ctx.path(REPORT_XLSX)),
        ],
    )


def _release_plan(ctx: RunContext) -> StagePlan:
    return StagePlan(parameters={"bundle": "checksummed"}, tool_versions={})


def _release_execute(ctx: RunContext, plan: StagePlan) -> StageResult:
    # The release stage is written by build_release_bundle after the run report exists,
    # because the bundle must checksum the run report itself.
    return StageResult(
        status=ModuleRunStatus.COMPLETED,
        reason="Release bundle prepared; artifacts are checksummed after the run report.",
    )


IMPLEMENTATIONS: dict[StageId, StageImplementation] = {
    StageId.BASECALL: StageImplementation(_basecall_plan, _basecall_execute, _basecall_settle),
    StageId.ALIGN: StageImplementation(_align_plan, _align_execute, _align_settle),
    StageId.INTAKE: StageImplementation(_intake_plan, _intake_execute),
    StageId.QC: StageImplementation(_qc_plan, _qc_execute),
    StageId.SV: StageImplementation(_sv_plan, _sv_execute),
    StageId.ASSEMBLE: StageImplementation(_assemble_plan, _assemble_execute),
    StageId.REPORT: StageImplementation(_report_plan, _report_execute),
    StageId.RELEASE: StageImplementation(_release_plan, _release_execute),
}


# --------------------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------------------


def _load_previous(envelope: RunEnvelope) -> RunReport | None:
    path = envelope.path(RUN_REPORT)
    if not path.is_file():
        return None
    try:
        return RunReport.model_validate_json(path.read_text(encoding="utf-8"))
    except ValueError:
        # A corrupt or schema-incompatible state file must not silently disable resume
        # checks; discarding it means the run simply repeats every stage.
        return None


def _build_report(
    config: RunConfiguration,
    records: Sequence[StageRecord],
    started_at: datetime,
    run_warnings: Sequence[str] = (),
) -> RunReport:
    outcomes = {item.stage: StageOutcome(item.status.value) for item in records}
    kind = InputKindName(config.manifest.input.kind.value)
    verdict = summarize(kind, outcomes)
    unverified = [
        item.stage
        for item in records
        if item.status in {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL}
        and item.verification.value in {"unverified_adapter", "not_implemented"}
    ]
    warnings: list[str] = list(run_warnings)
    if unverified:
        names = ", ".join(item.value for item in unverified)
        warnings.append(
            f"Stage(s) {names} completed on an adapter that has never been executed "
            "against the real tool. Their output is unverified."
        )
    return RunReport(
        run_id=config.run_id,
        sample_id=config.manifest.sample_id,
        input_kind=config.manifest.input.kind,
        genome_build=config.manifest.assay.genome_build,
        manifest=config.manifest,
        passed=verdict.passed,
        verdict_reason=verdict.describe(),
        stages=list(records),
        pipeline_version=config.pipeline_version,
        git_commit=config.git_commit,
        started_at=started_at,
        finished_at=datetime.now(UTC),
        unverified_stages=unverified,
        warnings=warnings,
        limitations=list(STANDING_LIMITATIONS),
    )


def _persist(envelope: RunEnvelope, report: RunReport) -> None:
    envelope.atomic_write_text(RUN_REPORT, report.model_dump_json(indent=2) + "\n")


def build_release_bundle(
    envelope: RunEnvelope, report: RunReport, config: RunConfiguration
) -> ReleaseBundle:
    """Checksum every exportable artifact of a run, including the run report itself."""
    _persist(envelope, report)
    run_report_artifact = envelope.fingerprint(RUN_REPORT)

    exportable: list[ArtifactRecord] = [ArtifactRecord.of(run_report_artifact)]
    withheld: list[str] = []
    seen = {run_report_artifact.relative_path}
    for record in report.stages:
        for artifact in record.outputs:
            if artifact.relative_path in seen:
                continue
            seen.add(artifact.relative_path)
            if artifact.exportable:
                exportable.append(artifact)
            else:
                withheld.append(artifact.relative_path)

    bundle = ReleaseBundle(
        run_id=report.run_id,
        sample_id=report.sample_id,
        pipeline_version=config.pipeline_version,
        git_commit=config.git_commit,
        run_report_sha256=run_report_artifact.sha256,
        artifacts=sorted(exportable, key=lambda item: item.relative_path),
        withheld_artifact_paths=sorted(withheld),
        total_bytes=sum(item.size_bytes for item in exportable),
        warnings=list(report.warnings),
        limitations=[
            *STANDING_LIMITATIONS,
            "Raw genomic artifacts are listed but deliberately excluded from this bundle.",
        ],
    )
    envelope.atomic_write_text(RELEASE_JSON, bundle.model_dump_json(indent=2) + "\n")
    envelope.atomic_write_text(RELEASE_CHECKSUMS, bundle.checksum_manifest())
    return bundle


class EnvelopeAlreadyReviewed(RuntimeError):
    """Raised when a run would modify an envelope somebody has already signed off.

    Nothing else in the design catches this. The lock stops two runs colliding *now*;
    content-addressed resume stops a stale artifact being accepted. Neither notices that a
    human accepted this envelope yesterday and a resumed run is about to rewrite what they
    accepted — which would leave the review pointing at content nobody reviewed.

    Deliberately not overridable by a flag. A flag would be used, and the situation it
    covers has a correct answer that costs nothing: use a new run id. The old envelope then
    keeps its review, and the new one gets its own.
    """


def _refuse_if_reviewed(envelope_root: Path) -> None:
    """Refuse to write into an envelope whose latest review accepts its current content.

    A rejected or stale review does not block: a rejection is often precisely why somebody
    re-runs, and a stale review already says it no longer describes what is on disk.
    """
    release = envelope_root / RELEASE_RELATIVE
    if not release.is_file():
        return
    digest = sha256_file(release)
    try:
        entries = read_review_log(envelope_root / REVIEW_LOG)
    except ReviewError:
        # An unreadable trail is reported by `ontseq review`, whose job that is. Blocking
        # the run here as well would make a corrupt log impossible to move past.
        return
    state, detail = review_state(entries, digest)
    if state is ReviewState.ACCEPTED:
        raise EnvelopeAlreadyReviewed(
            f"{envelope_root} carries an accepted review ({detail}). Running again would "
            "rewrite the content that was signed off. Use a different --run-id; the "
            "reviewed envelope then keeps its review and this run gets its own."
        )


def run_pipeline(
    config: RunConfiguration, *, runner: StreamingCommandRunner | None = None
) -> tuple[RunReport, ReleaseBundle | None]:
    """Execute every applicable stage and return the run report and release bundle."""
    started_at = datetime.now(UTC)
    envelope = RunEnvelope.create(
        config.output_base, run_id=config.run_id, sample_id=config.manifest.sample_id
    )
    _refuse_if_reviewed(envelope.root)
    with run_lock(
        envelope.root,
        run_id=config.run_id,
        sample_id=config.manifest.sample_id,
        pipeline_version=config.pipeline_version,
    ) as lock_warnings:
        return _run_locked(config, envelope, started_at, runner, lock_warnings)


def _run_locked(
    config: RunConfiguration,
    envelope: RunEnvelope,
    started_at: datetime,
    runner: StreamingCommandRunner | None,
    run_warnings: list[str],
) -> tuple[RunReport, ReleaseBundle | None]:
    """Execute the run. Split out so the lock covers every write, including the first."""
    envelope.atomic_write_text(
        "manifest/sample.manifest.json", config.manifest.model_dump_json(indent=2) + "\n"
    )
    envelope.atomic_write_text(
        "manifest/reference.lock.json", config.reference_lock.model_dump_json(indent=2) + "\n"
    )

    previous = None if config.force else _load_previous(envelope)
    context = RunContext(
        config=config,
        envelope=envelope,
        runner=runner or SubprocessRunner(),
        manifest=config.manifest,
    )
    kind = InputKindName(config.manifest.input.kind.value)
    outcomes: dict[StageId, StageOutcome] = {}
    records: list[StageRecord] = []

    for stage in planned_stages(kind):
        spec = SPEC_BY_STAGE[stage]
        record = _execute_stage(stage, context, kind, outcomes, previous)
        records.append(record)
        outcomes[stage] = StageOutcome(record.status.value)
        context.artifacts[stage] = [item.to_artifact() for item in record.outputs]
        if record.status == ModuleRunStatus.FAILED and spec.required:
            # Continue the loop so every remaining stage is recorded as NOT_RUN with a
            # reason, rather than vanishing from the report.
            pass
        _persist(envelope, _build_report(config, records, started_at, run_warnings))

    report = _build_report(config, records, started_at, run_warnings)
    bundle: ReleaseBundle | None = None
    if outcomes.get(StageId.RELEASE) == StageOutcome.COMPLETED:
        bundle = build_release_bundle(envelope, report, config)
    else:
        _persist(envelope, report)
    return report, bundle


def _execute_stage(
    stage: StageId,
    context: RunContext,
    kind: InputKindName,
    outcomes: Mapping[StageId, StageOutcome],
    previous: RunReport | None,
) -> StageRecord:
    spec = SPEC_BY_STAGE[stage]
    base = {
        "stage": stage,
        "title": spec.title,
        "verification": spec.verification,
        "required": spec.required,
    }

    blocked = blocking_dependency(stage, kind, outcomes)
    if blocked is not None:
        return StageRecord(**base, status=ModuleRunStatus.NOT_RUN, reason=blocked.describe())

    implementation = IMPLEMENTATIONS.get(stage)
    if implementation is None:
        return StageRecord(
            **base,
            status=ModuleRunStatus.NOT_RUN,
            reason=f"No adapter is wired in for this stage. {spec.purpose}",
        )

    started = datetime.now(UTC)
    try:
        plan = implementation.plan(context)
    except (StageFailure, ToolExecutionError, ValueError, OSError) as error:
        return StageRecord(
            **base,
            status=ModuleRunStatus.FAILED,
            reason=f"Stage could not be prepared: {error}",
            started_at=started,
            finished_at=datetime.now(UTC),
        )

    signature = stage_signature(
        stage=stage.value,
        upstream=context.upstream(stage, kind),
        parameters={
            **plan.parameters,
            "_runner_pipeline_version": context.config.pipeline_version,
            "_runner_git_commit": context.config.git_commit,
        },
        tool_versions=plan.tool_versions,
        external_inputs=plan.external_inputs,
    )

    prior = previous.record_for(stage) if previous else None
    if (
        prior is not None
        and prior.status in StageRecord.CONCLUDED
        and prior.signature == signature
        and not context.envelope.verify([item.to_artifact() for item in prior.outputs])
    ):
        if implementation.settle is not None:
            try:
                implementation.settle(context, [item.to_artifact() for item in prior.outputs])
            except (StageFailure, ValueError, OSError) as error:
                # The artifacts verified, so this is a bug rather than stale state; fail
                # the stage instead of letting a half-settled context reach the next one.
                return StageRecord(
                    **base,
                    status=ModuleRunStatus.FAILED,
                    reason=f"Stage resumed but its outputs could not be adopted: {error}",
                    signature=signature,
                    started_at=started,
                    finished_at=datetime.now(UTC),
                )
        return prior.model_copy(
            update={
                "resumed": True,
                "reason": f"Resumed unchanged from a previous run. {prior.reason}",
            }
        )

    try:
        result = implementation.execute(context, plan)
    except (StageFailure, ToolExecutionError, ValueError, OSError) as error:
        return StageRecord(
            **base,
            status=ModuleRunStatus.FAILED,
            reason=str(error),
            signature=signature,
            started_at=started,
            finished_at=datetime.now(UTC),
            limitations=[f"Diagnostic: {traceback.format_exc(limit=1).strip().splitlines()[-1]}"],
        )

    if implementation.settle is not None and result.status in {
        ModuleRunStatus.COMPLETED,
        ModuleRunStatus.NO_CALL,
    }:
        try:
            implementation.settle(context, result.outputs)
        except (StageFailure, ValueError, OSError) as error:
            return StageRecord(
                **base,
                status=ModuleRunStatus.FAILED,
                reason=f"Stage ran but its outputs could not be adopted: {error}",
                signature=signature,
                started_at=started,
                finished_at=datetime.now(UTC),
            )

    finished = datetime.now(UTC)
    return StageRecord(
        **base,
        status=result.status,
        reason=result.reason,
        signature=signature,
        started_at=started,
        finished_at=finished,
        duration_seconds=(finished - started).total_seconds(),
        outputs=[ArtifactRecord.of(item) for item in result.outputs],
        tools=result.tools,
        warnings=result.warnings,
        limitations=result.limitations,
    )

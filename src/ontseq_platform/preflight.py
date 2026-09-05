"""Check a run's preconditions before it starts, instead of stage by stage while it runs.

Every gate this command applies already exists somewhere in the pipeline, and each one
fails closed correctly. The problem is *when* they fire. A POD5 run discovers a missing
Dorado model after the envelope has been created and the lock taken; an alignment run
discovers a missing reference index after intake; a run into a busy envelope discovers
that only when it tries to lock. Individually correct, collectively expensive: the
feedback arrives after the run has been queued, scheduled and partly executed.

So this asks the same questions up front, in a couple of seconds, and answers them without
side effects. Nothing here creates an envelope, takes a lock or writes an artifact.

Two properties matter more than the check list:

**Preflight must agree with the run.** A preflight that clears a run which then fails on
the very thing it checked is worse than no preflight, because it converts a fast, honest
failure into a slow, surprising one. So version strings are parsed by the adapters' own
parsers (:func:`ontseq_platform.align.probe_versions`,
:func:`ontseq_platform.qc.cramino_version`, :func:`ontseq_platform.sniffles.sniffles_version`,
:func:`ontseq_platform.basecall.dorado_version`) rather than re-implemented here, and the
reference is compared through the ``source_fai_sha256`` the lock already records.

**Not knowing is a distinct answer.** GPU availability on a compute node, and whether the
free disk is enough for a dataset whose size nobody has measured here, are not knowable
from this process. They are reported ``UNKNOWN`` with the raw figure attached, never
guessed into an ``OK``. Free space becomes a real check only when the caller states a
requirement with ``--require-free-gb``; inventing a threshold would produce a number that
looks validated and is not.
"""

from __future__ import annotations

import shutil
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from .align import AlignmentPolicy, parse_version
from .basecall import BasecallPolicy, dorado_version, model_signature
from .cutesv import cutesv_version
from .execution import CommandRunner, SubprocessRunner
from .methylation import MethylationPolicy, MethylationRegionSource, modkit_version
from .models import (
    AnalysisModule,
    AssayMode,
    CuteSvPolicy,
    InputKind,
    ReferenceLock,
    SampleManifest,
    SnifflesPolicy,
)
from .pipeline.checks import Check, CheckList, ToolRequirement, required_tools
from .pipeline.lock import LOCK_FILENAME, holder_is_running, read_holder
from .pipeline.stages import (
    SPEC_BY_STAGE,
    InputKindName,
    StageId,
    VerificationStatus,
    planned_stages,
    unverified_specs,
    verification_of,
)
from .qc import cramino_version
from .reference import sha256_file
from .sniffles import sniffles_version
from .target_coverage import TargetCoveragePolicy, load_target_bed, mosdepth_version

GIGABYTE = 1024**3


@dataclass(frozen=True)
class PreflightRequest:
    """Everything preflight needs, resolved the same way ``ontseq run`` resolves it."""

    manifest: SampleManifest
    reference_lock: ReferenceLock
    output_base: Path
    run_id: str
    executables: Mapping[str, str] = field(default_factory=dict)
    reference_fasta: Path | None = None
    pod5_directory: Path | None = None
    alignment_policy: AlignmentPolicy | None = None
    basecall_policy: BasecallPolicy | None = None
    sniffles_policy: SnifflesPolicy | None = None
    cutesv_policy: CuteSvPolicy | None = None
    target_coverage_policy: TargetCoveragePolicy | None = None
    methylation_policy: MethylationPolicy | None = None
    #: What each stage's *registered implementation* says it has been verified against,
    #: overriding the declared graph. Supplied by the caller because it is a property of
    #: what the run will install, not of the manifest. Without it a preflight reports the
    #: bare graph while the run it is checking installs adapters the graph knows nothing
    #: about — which is how CNV came to be announced as having no adapter for runs that
    #: then executed a real QDNAseq analysis.
    stage_verification: Mapping[StageId, VerificationStatus] = field(default_factory=dict)
    #: Free space the caller knows this run needs. Without it, space is reported, not judged.
    require_free_gb: float | None = None

    def executable(self, name: str) -> str:
        return self.executables.get(name, name)

    @property
    def input_kind(self) -> InputKindName:
        return InputKindName(self.manifest.input.kind.value)

    @property
    def envelope_root(self) -> Path:
        return self.output_base / self.run_id / self.manifest.sample_id


def _check_input(request: PreflightRequest, checks: CheckList) -> None:
    """The declared input exists and has the shape its kind promises."""
    kind = request.manifest.input.kind
    if kind == InputKind.POD5:
        directory = request.pod5_directory
        if directory is None:
            checks.failed(
                "input.pod5",
                "the manifest declares POD5 input but no POD5 directory was given",
                remedy="pass --pod5-dir",
                stage=StageId.BASECALL,
            )
            return
        if not directory.is_dir():
            checks.failed(
                "input.pod5",
                f"POD5 directory does not exist: {directory}",
                stage=StageId.BASECALL,
            )
            return
        pod5_files = sorted(directory.rglob("*.pod5"))
        if not pod5_files:
            checks.failed(
                "input.pod5",
                f"no .pod5 files beneath {directory}",
                remedy="check the sequencing run finished writing before starting",
                stage=StageId.BASECALL,
            )
            return
        total = sum(item.stat().st_size for item in pod5_files)
        checks.ok(
            "input.pod5",
            f"{len(pod5_files)} POD5 file(s), {total / GIGABYTE:.2f} GiB",
            stage=StageId.BASECALL,
        )
        return

    path = Path(request.manifest.input.path)
    if not path.is_file():
        checks.failed("input.bam", f"declared input does not exist: {path}")
        return
    if path.stat().st_size == 0:
        checks.failed("input.bam", f"declared input is empty: {path}")
        return
    checks.ok("input.bam", f"{path} ({path.stat().st_size / GIGABYTE:.2f} GiB)")

    if kind != InputKind.ALIGNED_BAM:
        checks.skipped(
            "input.bam.index",
            "an unaligned BAM carries no index; the pipeline produces one when it aligns",
        )
        return
    # The manifest model already requires index_path for an aligned BAM, so this is about
    # the file being present rather than the field being filled in.
    index = Path(request.manifest.input.index_path or "")
    if not index.is_file():
        checks.failed(
            "input.bam.index",
            f"BAM index does not exist: {index}",
            remedy=f"samtools index {path}",
            stage=StageId.INTAKE,
        )
        return
    if index.stat().st_mtime < path.stat().st_mtime:
        checks.warning(
            "input.bam.index",
            "the index is older than the BAM it indexes",
            remedy=f"samtools index {path}",
            stage=StageId.INTAKE,
        )
        return
    checks.ok("input.bam.index", str(index), stage=StageId.INTAKE)


def _check_reference(request: PreflightRequest, checks: CheckList) -> None:
    """The reference the run will align against is the one the lock was built from."""
    manifest_build = request.manifest.assay.genome_build
    if manifest_build != request.reference_lock.genome_build:
        checks.failed(
            "reference.build",
            f"manifest declares {manifest_build.value} but the lock is "
            f"{request.reference_lock.genome_build.value}",
        )
    else:
        checks.ok(
            "reference.build",
            f"{manifest_build.value}, {len(request.reference_lock.contigs)} contigs",
        )

    if request.manifest.assay.reference_id != request.reference_lock.reference_id:
        checks.failed(
            "reference.id",
            f"manifest declares reference {request.manifest.assay.reference_id!r} but the "
            f"lock is {request.reference_lock.reference_id!r}",
        )
    else:
        checks.ok("reference.id", request.reference_lock.reference_id)

    needs_fasta = (
        StageId.ALIGN in planned_stages(request.input_kind) or request.cutesv_policy is not None
    )
    reference_stage = (
        StageId.ALIGN if StageId.ALIGN in planned_stages(request.input_kind) else StageId.SV
    )
    if not needs_fasta:
        checks.skipped(
            "reference.fasta", "this run neither aligns nor runs cuteSV, so no FASTA is needed"
        )
        checks.skipped(
            "reference.fai", "this run neither aligns nor runs cuteSV, so no FASTA index is needed"
        )
        return

    fasta = request.reference_fasta
    if fasta is None:
        checks.failed(
            "reference.fasta",
            "this run aligns or runs cuteSV but no reference FASTA was given",
            remedy="pass --reference-fasta",
            stage=reference_stage,
        )
        return
    if not fasta.is_file():
        checks.failed(
            "reference.fasta", f"reference FASTA does not exist: {fasta}", stage=reference_stage
        )
        return
    checks.ok("reference.fasta", str(fasta), stage=reference_stage)

    fai = fasta.with_suffix(fasta.suffix + ".fai")
    if not fai.is_file():
        checks.failed(
            "reference.fai",
            f"reference index does not exist: {fai}",
            remedy=f"samtools faidx {fasta}",
            stage=reference_stage,
        )
        return
    # The lock records the checksum of the .fai it was generated from, so this is an exact
    # answer to "is this the reference the lock describes?" rather than a contig-by-contig
    # comparison that could pass on a different build with the same chromosome names.
    observed = sha256_file(fai)
    if observed != request.reference_lock.source_fai_sha256:
        checks.failed(
            "reference.fai",
            "the reference index does not match the one the lock was built from",
            remedy=(
                "point --reference-fasta at the locked reference, or regenerate the lock "
                "with `ontseq reference-lock` if the reference genuinely changed"
            ),
            stage=reference_stage,
        )
        return
    checks.ok("reference.fai", f"{fai} matches the lock checksum", stage=reference_stage)


def _measures_targets(request: PreflightRequest) -> bool:
    """Whether this run will actually measure per-target coverage.

    The stage applies to every input kind, but only an adaptive-sampling run does anything
    in it; any other mode records that targets are out of scope without touching Mosdepth.
    """
    return (
        StageId.TARGET_COVERAGE in planned_stages(request.input_kind)
        and request.manifest.assay.mode == AssayMode.ADAPTIVE_SAMPLING
    )


def _fatal_stages(request: PreflightRequest) -> frozenset[StageId]:
    """Stages this particular run cannot get away without, beyond the declared-required set.

    ``StageSpec.required`` is a property of the graph, not of a run. Target coverage is
    declared optional because an lcWGS run legitimately records it as out of scope — but an
    adaptive-sampling run neither skips it nor survives it failing: the runner refuses to
    continue without a policy and a target BED, and ``summarize`` fails a run on any FAILED
    stage whether or not the graph called it required. Preflight has to apply the same rule,
    or it clears a run that cannot succeed and reports the missing tool as a warning.
    """
    fatal = {stage for stage in planned_stages(request.input_kind) if SPEC_BY_STAGE[stage].required}
    if _measures_targets(request):
        fatal.add(StageId.TARGET_COVERAGE)
    return frozenset(fatal)


def _check_target_coverage(request: PreflightRequest, checks: CheckList) -> None:
    """The adaptive-sampling inputs exist and parse, before the envelope is created.

    Both of these fail the run closed inside the target-coverage stage, which is right —
    a run whose enrichment was never measured must not produce a report that looks
    complete. But they fail it after the envelope exists and the lock is taken, and both
    are answerable here in milliseconds.
    """
    if not _measures_targets(request):
        reason = (
            f"assay mode is {request.manifest.assay.mode.value}; per-target coverage does not apply"
        )
        checks.skipped("target_coverage.policy", reason)
        checks.skipped("target_coverage.bed", reason)
        return

    if request.target_coverage_policy is None:
        checks.failed(
            "target_coverage.policy",
            "adaptive sampling was selected but no target-coverage policy was supplied",
            remedy="pass --target-coverage-policy with the technical policy for this assay",
            stage=StageId.TARGET_COVERAGE,
        )
    else:
        checks.ok(
            "target_coverage.policy",
            f"{request.target_coverage_policy.profile_id} "
            f"({request.target_coverage_policy.status})",
            stage=StageId.TARGET_COVERAGE,
        )

    declared = request.manifest.assay.target_bed
    if not declared:
        checks.failed(
            "target_coverage.bed",
            "adaptive sampling requires assay.target_bed",
            remedy="name the controlled analysis ROI BED in the manifest",
            stage=StageId.TARGET_COVERAGE,
        )
        return
    bed = Path(declared)
    try:
        regions = load_target_bed(bed)
    except (ValueError, OSError) as error:
        checks.failed(
            "target_coverage.bed",
            f"{bed} could not be read as a target BED: {error}",
            remedy="point assay.target_bed at the controlled, readable ROI BED for this panel",
            stage=StageId.TARGET_COVERAGE,
        )
        return
    bases = sum(region.length for region in regions)
    checks.ok(
        "target_coverage.bed",
        f"{len(regions)} target(s) over {bases} bp "
        f"({request.manifest.assay.target_bed_role.value})",
        stage=StageId.TARGET_COVERAGE,
    )


def _expected_version(request: PreflightRequest, tool: str) -> str | None:
    """What the policies lock this tool to, when a *planned* stage will enforce that lock.

    Gated on the plan rather than on which policies happen to be loaded. The alignment
    policy locks a samtools version, but an aligned-BAM run never aligns and never applies
    that lock — enforcing it here would refuse runs that ``ontseq run`` completes happily,
    which is the one failure mode a preflight must not have.
    """
    planned = planned_stages(request.input_kind)
    if request.alignment_policy is not None and StageId.ALIGN in planned:
        if tool == "minimap2":
            return request.alignment_policy.expected_minimap2_version
        if tool == "samtools":
            return request.alignment_policy.expected_samtools_version
    if (
        tool == "sniffles"
        and request.sniffles_policy is not None
        and StageId.SV in planned
        and _analyses_structural_variants(request)
    ):
        return request.sniffles_policy.expected_version
    if (
        tool == "cutesv"
        and request.cutesv_policy is not None
        and StageId.SV in planned
        and _analyses_structural_variants(request)
    ):
        return request.cutesv_policy.expected_version
    if tool == "dorado" and request.basecall_policy is not None and StageId.BASECALL in planned:
        return request.basecall_policy.expected_version
    # Gated on the assay rather than only on the plan: the stage is planned for every input
    # kind, but an lcWGS run never invokes Mosdepth and so never applies its lock.
    if (
        tool == "mosdepth"
        and request.target_coverage_policy is not None
        and _measures_targets(request)
    ):
        return request.target_coverage_policy.expected_version
    if (
        tool == "modkit"
        and request.methylation_policy is not None
        and _analyses_methylation(request)
    ):
        return request.methylation_policy.expected_version
    return None


def _probe_version(runner: CommandRunner, tool: str, executable: str) -> str:
    """Ask one tool its version, parsing it exactly as the adapter that will run it does.

    One tool at a time on purpose. ``align.probe_versions`` covers minimap2 and samtools
    together and fails closed if either is unidentifiable, which is right for a run and
    wrong here: preflight has to report each binary separately, so a missing samtools does
    not come back as a minimap2 problem.
    """
    result = runner.run([executable, "--version"], timeout_seconds=60)
    if result.returncode != 0:
        raise ValueError(f"version probe returned exit code {result.returncode}")
    combined = f"{result.stdout}\n{result.stderr}"
    if tool in {"minimap2", "samtools"}:
        return parse_version(combined, tool=tool)
    if tool == "cramino":
        return cramino_version(result.stdout)
    if tool == "sniffles":
        return sniffles_version(combined)
    if tool == "cutesv":
        return cutesv_version(combined)
    if tool == "dorado":
        return dorado_version(combined)
    if tool == "mosdepth":
        return mosdepth_version(combined)
    if tool == "modkit":
        return modkit_version(combined)
    raise ValueError(f"no version parser for {tool!r}")


def _check_tools(request: PreflightRequest, runner: CommandRunner, checks: CheckList) -> None:
    """Every binary the planned stages will invoke is present, runnable and locked-version.

    A tool serving only optional stages is a warning when missing rather than a failure:
    the run completes without it, records the stage as ``NOT_RUN``, and that is a
    legitimate outcome the operator should be told about in advance, not blocked on.
    """
    fatal = _fatal_stages(request)
    requirements = list(required_tools(request.input_kind))
    if request.sniffles_policy is None and request.cutesv_policy is not None:
        requirements = [item for item in requirements if item.name != "sniffles"]
    if not _analyses_methylation(request):
        # The stage is planned for every input kind but invokes modkit only when the
        # manifest asks for it, so an absent binary is not news to a run that never
        # wanted one.
        requirements = [item for item in requirements if item.name != "modkit"]
    if not _analyses_structural_variants(request):
        requirements = [item for item in requirements if item.name not in {"sniffles", "cutesv"}]
    if request.cutesv_policy is not None and StageId.SV in planned_stages(request.input_kind):
        requirements.append(ToolRequirement(name="cutesv", stages=(StageId.SV,), required=False))
    for requirement in requirements:
        name = f"tool.{requirement.name}"
        executable = request.executable(requirement.name)
        stage = requirement.stages[0]
        # A tool is fatal when any stage needing it is one this run cannot do without.
        # Mosdepth is the case that matters: the graph calls target coverage optional, but
        # an adaptive-sampling run fails outright without it.
        blocking = requirement.required or any(item in fatal for item in requirement.stages)
        located = shutil.which(executable)
        if located is None and not Path(executable).is_file():
            detail = f"{executable} is not on PATH"
            remedy = f"install {requirement.name}, or pass --{requirement.name} with its path"
            if blocking:
                checks.failed(name, detail, remedy=remedy, stage=stage)
            else:
                stage_names = ", ".join(item.value for item in requirement.stages)
                checks.warning(
                    name,
                    f"{detail}; the optional stage(s) {stage_names} will record NOT_RUN",
                    remedy=remedy,
                    stage=stage,
                )
            continue

        try:
            observed = _probe_version(runner, requirement.name, executable)
        except (ValueError, OSError) as error:
            checks.failed(
                name,
                f"{executable} could not be identified: {error}",
                remedy="check the installation is complete and executable",
                stage=stage,
            )
            continue

        expected = _expected_version(request, requirement.name)
        if expected is None:
            checks.ok(name, f"{observed} (no policy lock declares a version)", stage=stage)
        elif observed != expected:
            checks.failed(
                name,
                f"version {observed!r} does not match the policy lock {expected!r}",
                remedy=(
                    f"install {requirement.name} {expected}, or change the policy lock "
                    "deliberately and record why"
                ),
                stage=stage,
            )
        else:
            checks.ok(name, f"{observed} matches the policy lock", stage=stage)


def _check_sv_configuration(request: PreflightRequest, checks: CheckList) -> None:
    if StageId.SV not in planned_stages(request.input_kind):
        return
    if not _analyses_structural_variants(request):
        checks.skipped("sv.callers", "the manifest does not request the structural-variant module")
        return
    if request.sniffles_policy is None and request.cutesv_policy is None:
        checks.warning(
            "sv.callers",
            "no caller policy was loaded; the optional SV stage will fail closed if executed",
            remedy="supply at least one version-locked caller policy",
            stage=StageId.SV,
        )
        return
    callers = []
    if request.sniffles_policy is not None:
        callers.append("Sniffles2")
    if request.cutesv_policy is not None:
        callers.append("cuteSV")
    checks.ok("sv.callers", " + ".join(callers), stage=StageId.SV)


def _check_basecalling(request: PreflightRequest, checks: CheckList) -> None:
    """The Dorado model exists and is the one the policy locks, before hours are spent."""
    if StageId.BASECALL not in planned_stages(request.input_kind):
        checks.skipped("basecall.model", "this run does not basecall")
        checks.skipped("basecall.modified_bases", "this run does not basecall")
        return
    policy = request.basecall_policy
    if policy is None:
        checks.failed(
            "basecall.model",
            "this run basecalls but no basecalling policy was loaded",
            remedy="pass --basecall-policy",
            stage=StageId.BASECALL,
        )
        checks.skipped("basecall.modified_bases", "no basecalling policy to read")
        return

    observed = model_signature(policy.model)
    if policy.model_sha256 is None:
        if observed is None:
            checks.warning(
                "basecall.model",
                f"model {policy.model!r} is a name, not a local directory, so which "
                "weights were used cannot be recorded in provenance",
                remedy=(
                    "point the policy at a downloaded model directory and record its "
                    "checksum in model_sha256"
                ),
                stage=StageId.BASECALL,
            )
        else:
            checks.warning(
                "basecall.model",
                f"model directory present but the policy locks no checksum "
                f"(observed {observed[:12]}...)",
                remedy=f"set model_sha256: {observed} in the basecalling policy",
                stage=StageId.BASECALL,
            )
    elif observed is None:
        checks.failed(
            "basecall.model",
            f"the policy locks a model checksum but {policy.model!r} is not a local "
            "directory that can be fingerprinted",
            stage=StageId.BASECALL,
        )
    elif observed != policy.model_sha256:
        checks.failed(
            "basecall.model",
            "the model on disk does not match the policy lock",
            remedy="restore the locked model, or re-lock deliberately and record why",
            stage=StageId.BASECALL,
        )
    else:
        checks.ok("basecall.model", f"{policy.model} matches the lock", stage=StageId.BASECALL)

    if policy.modified_bases:
        checks.ok(
            "basecall.modified_bases",
            ", ".join(policy.modified_bases),
            stage=StageId.BASECALL,
        )
    else:
        # Discovered after a twelve-hour basecall, this costs the whole run: the reads
        # carry no MM/ML tags and no later stage can recover them.
        checks.warning(
            "basecall.modified_bases",
            "no modified-base model requested, so the reads will carry no MM/ML tags and "
            "cannot support a later methylation analysis",
            remedy="set modified_bases in the basecalling policy if methylation is wanted",
            stage=StageId.BASECALL,
        )


def _analyses_structural_variants(request: PreflightRequest) -> bool:
    """Whether this run asked for structural variants at all.

    The stage is planned for every input kind but invokes its callers only when the
    manifest requests the module, so a CNV-only run needs neither the binaries nor the
    advice about them.
    """
    return AnalysisModule.SV in request.manifest.analysis.modules


def _analyses_methylation(request: PreflightRequest) -> bool:
    """Whether this run asked for modified bases at all.

    Unlike target coverage, the methylation lane is gated on the requested analysis rather
    than on the assay mode: an lcWGS and an Adaptive Sampling run can both carry MM/ML
    tags, and neither should be told about a missing modkit it will never invoke.
    """
    return AnalysisModule.METHYLATION in request.manifest.analysis.modules


def _check_methylation(request: PreflightRequest, checks: CheckList) -> None:
    """The methylation lane has everything it needs, before the envelope is created.

    The one precondition preflight cannot answer is the important one: whether the reads
    carry ``MM``/``ML`` tags at all. Reading that means scanning the BAM, which is the
    stage's job. It is stated as a warning here so an operator learns before the run that
    a BAM basecalled without a modified-base model will fail the stage rather than quietly
    produce an empty pileup.
    """
    if not _analyses_methylation(request):
        reason = "the manifest does not request the methylation module"
        checks.skipped("methylation.policy", reason)
        checks.skipped("methylation.reference", reason)
        return

    policy = request.methylation_policy
    if policy is None:
        checks.failed(
            "methylation.policy",
            "the methylation module was requested but no methylation policy was supplied",
            remedy="pass --methylation-policy with the technical policy for this assay",
            stage=StageId.METHYLATION,
        )
        return
    checks.ok(
        "methylation.policy",
        f"{policy.profile_id} ({policy.status})",
        stage=StageId.METHYLATION,
    )

    if policy.region_source == MethylationRegionSource.TARGET_BED and not (
        request.manifest.assay.target_bed
    ):
        checks.failed(
            "methylation.regions",
            "the policy aggregates over the target design but the manifest declares no target BED",
            remedy="name the target BED in the manifest, or set region_source: chromosome",
            stage=StageId.METHYLATION,
        )

    if not policy.cpg_only:
        checks.skipped(
            "methylation.reference",
            "the policy does not restrict the pileup to a reference motif",
        )
    elif request.reference_fasta is None:
        checks.failed(
            "methylation.reference",
            "the policy restricts the pileup to CpG sites, which is a property of the "
            "reference, but no reference FASTA was given",
            remedy="pass --reference-fasta",
            stage=StageId.METHYLATION,
        )
    elif not request.reference_fasta.is_file():
        checks.failed(
            "methylation.reference",
            f"reference FASTA does not exist: {request.reference_fasta}",
            stage=StageId.METHYLATION,
        )
    else:
        checks.ok("methylation.reference", str(request.reference_fasta), stage=StageId.METHYLATION)

    if request.manifest.input.kind == InputKind.ALIGNED_BAM:
        checks.warning(
            "methylation.modified_base_tags",
            "whether the aligned BAM carries MM/ML tags cannot be answered without reading "
            "it; the stage verifies this and fails closed rather than emitting an empty "
            "pileup, which would read as an unmethylated sample",
            remedy="confirm the BAM was basecalled with a modified-base model",
            stage=StageId.METHYLATION,
        )


def _check_envelope(request: PreflightRequest, checks: CheckList) -> None:
    """The output location is writable and nobody else is working in this envelope."""
    root = request.envelope_root
    lock_path = root / LOCK_FILENAME
    if lock_path.exists():
        holder = read_holder(lock_path)
        alive = holder_is_running(holder) if holder is not None else None
        described = holder.describe() if holder is not None else "an unreadable lock file"
        if alive is False:
            checks.warning(
                "envelope.lock",
                f"a previous run died here ({described}); this run reclaims the lock and resumes",
            )
        else:
            checks.failed(
                "envelope.lock",
                f"the envelope is in use by {described}",
                remedy="wait for that run to finish, or use a different --run-id",
            )
    else:
        checks.ok("envelope.lock", f"{root} is free")

    # Probed on the nearest existing ancestor rather than by creating the output directory,
    # because preflight promises to leave nothing behind: a run that is never started must
    # not have an output tree sitting there afterwards suggesting it was.
    existing = request.output_base.resolve()
    while not existing.exists() and existing != existing.parent:
        existing = existing.parent
    try:
        probe = existing / ".ontseq-preflight-probe"
        probe.write_text("", encoding="utf-8")
        probe.unlink()
    except OSError as error:
        checks.failed(
            "output.writable",
            f"cannot write beneath {existing}: {error}",
            remedy="check the path, its permissions and whether the filesystem is mounted",
        )
        return
    checks.ok("output.writable", str(existing))


def _check_disk(request: PreflightRequest, checks: CheckList) -> None:
    """Report free space; judge it only against a requirement the caller states.

    There is no measured relationship in this repository between an input's size and the
    space a run consumes — that depends on the lab's chemistry, depth and retention policy,
    and nobody has measured it here. A multiplier invented in this function would look like
    a validated figure and would not be one, so free space is reported as ``UNKNOWN`` until
    somebody who knows the answer supplies it.
    """
    existing = request.output_base.resolve()
    while not existing.exists() and existing != existing.parent:
        existing = existing.parent
    try:
        usage = shutil.disk_usage(existing)
    except OSError as error:
        checks.unknown("disk.free", f"free space could not be determined: {error}")
        return
    free_gb = usage.free / GIGABYTE
    if request.require_free_gb is None:
        checks.unknown(
            "disk.free",
            f"{free_gb:.1f} GiB free; no requirement was stated, so this is reported, not judged",
            remedy=(
                "pass --require-free-gb once the space a run of this kind consumes has "
                "been measured on real data"
            ),
        )
        return
    if free_gb < request.require_free_gb:
        checks.failed(
            "disk.free",
            f"{free_gb:.1f} GiB free, below the required {request.require_free_gb:.1f} GiB",
            remedy="free space or point --output-dir at a larger filesystem",
        )
        return
    checks.ok(
        "disk.free",
        f"{free_gb:.1f} GiB free, at or above the required {request.require_free_gb:.1f} GiB",
    )


def _check_adapters(request: PreflightRequest, checks: CheckList) -> None:
    """Say before the run what its stages will and will not actually produce.

    Two claims that ``unverified_specs`` groups together but a reader must not: a stage on
    an ``unverified_adapter`` *will run*, on code nobody has executed against the real tool,
    and its output is an assumption. A ``not_implemented`` stage has no adapter at all — it
    will record ``NOT_RUN``, which is not a statement about code quality and, crucially, not
    a negative biological finding either. Reporting the second as "an adapter that has never
    been executed" would be false.
    """
    stages = tuple(
        stage
        for stage in planned_stages(request.input_kind)
        # The methylation stage is planned for every input kind but invokes its adapter
        # only when the manifest asks for modified bases. Warning that an unexercised
        # adapter "will run" on a run that never calls it would train operators to ignore
        # this line, which is the one line that has to keep meaning something.
        if stage is not StageId.METHYLATION or _analyses_methylation(request)
    )
    resolved = unverified_specs(stages, request.stage_verification)
    unverified = [
        spec
        for spec in resolved
        if verification_of(spec.stage, request.stage_verification)
        is VerificationStatus.UNVERIFIED_ADAPTER
    ]
    missing = [
        spec
        for spec in resolved
        if verification_of(spec.stage, request.stage_verification)
        is VerificationStatus.NOT_IMPLEMENTED
    ]

    if unverified:
        names = ", ".join(spec.stage.value for spec in unverified)
        checks.warning(
            "adapters.verification",
            f"planned stage(s) {names} use adapters that have never been executed against "
            "the real tool; treat this run as a test of them",
        )
    else:
        checks.ok("adapters.verification", "every stage that will run has an exercised adapter")

    if missing:
        names = ", ".join(spec.stage.value for spec in missing)
        checks.warning(
            "stages.not_implemented",
            f"planned stage(s) {names} have no adapter wired in and will record NOT_RUN, "
            "which is not a negative biological finding",
        )
    else:
        checks.ok("stages.not_implemented", "every planned stage has an adapter wired in")


def preflight(request: PreflightRequest, *, runner: CommandRunner | None = None) -> list[Check]:
    """Answer every precondition question that can be answered without starting the run.

    Ordered the way a reader wants to read it: the input first, then the reference it is
    interpreted against, then the tools, then the place the output goes. Checks are never
    short-circuited on the first failure — an operator fixing a broken setup should get the
    whole list, not one problem per attempt.
    """
    command_runner = runner or SubprocessRunner()
    checks = CheckList()
    _check_input(request, checks)
    _check_reference(request, checks)
    _check_sv_configuration(request, checks)
    _check_tools(request, command_runner, checks)
    _check_basecalling(request, checks)
    _check_target_coverage(request, checks)
    _check_methylation(request, checks)
    _check_envelope(request, checks)
    _check_disk(request, checks)
    _check_adapters(request, checks)
    return checks.checks

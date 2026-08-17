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
from .execution import CommandRunner, SubprocessRunner
from .models import InputKind, ReferenceLock, SampleManifest, SnifflesPolicy
from .pipeline.checks import Check, CheckList, required_tools
from .pipeline.lock import LOCK_FILENAME, holder_is_running, read_holder
from .pipeline.stages import InputKindName, StageId, planned_stages, unverified_specs
from .qc import cramino_version
from .reference import sha256_file
from .sniffles import sniffles_version

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

    if StageId.ALIGN not in planned_stages(request.input_kind):
        checks.skipped("reference.fasta", "this run does not align, so no FASTA is needed")
        checks.skipped("reference.fai", "this run does not align, so no FASTA index is needed")
        return

    fasta = request.reference_fasta
    if fasta is None:
        checks.failed(
            "reference.fasta",
            "this run aligns but no reference FASTA was given",
            remedy="pass --reference-fasta",
            stage=StageId.ALIGN,
        )
        return
    if not fasta.is_file():
        checks.failed(
            "reference.fasta", f"reference FASTA does not exist: {fasta}", stage=StageId.ALIGN
        )
        return
    checks.ok("reference.fasta", str(fasta), stage=StageId.ALIGN)

    fai = fasta.with_suffix(fasta.suffix + ".fai")
    if not fai.is_file():
        checks.failed(
            "reference.fai",
            f"reference index does not exist: {fai}",
            remedy=f"samtools faidx {fasta}",
            stage=StageId.ALIGN,
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
            stage=StageId.ALIGN,
        )
        return
    checks.ok("reference.fai", f"{fai} matches the lock checksum", stage=StageId.ALIGN)


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
    if tool == "sniffles" and request.sniffles_policy is not None and StageId.SV in planned:
        return request.sniffles_policy.expected_version
    if tool == "dorado" and request.basecall_policy is not None and StageId.BASECALL in planned:
        return request.basecall_policy.expected_version
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
    if tool == "dorado":
        return dorado_version(combined)
    raise ValueError(f"no version parser for {tool!r}")


def _check_tools(request: PreflightRequest, runner: CommandRunner, checks: CheckList) -> None:
    """Every binary the planned stages will invoke is present, runnable and locked-version.

    A tool serving only optional stages is a warning when missing rather than a failure:
    the run completes without it, records the stage as ``NOT_RUN``, and that is a
    legitimate outcome the operator should be told about in advance, not blocked on.
    """
    for requirement in required_tools(request.input_kind):
        name = f"tool.{requirement.name}"
        executable = request.executable(requirement.name)
        stage = requirement.stages[0]
        located = shutil.which(executable)
        if located is None and not Path(executable).is_file():
            detail = f"{executable} is not on PATH"
            remedy = f"install {requirement.name}, or pass --{requirement.name} with its path"
            if requirement.required:
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
                f"a previous run died here ({described}); this run reclaims the lock and "
                "resumes",
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
            f"{free_gb:.1f} GiB free; no requirement was stated, so this is reported, "
            "not judged",
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
    """Say before the run which of its stages rest on code nobody has executed for real."""
    stages = planned_stages(request.input_kind)
    unverified = unverified_specs(stages)
    if not unverified:
        checks.ok("adapters.verification", "every planned stage has an exercised adapter")
        return
    names = ", ".join(spec.stage.value for spec in unverified)
    checks.warning(
        "adapters.verification",
        f"planned stage(s) {names} use adapters that have never been executed against the "
        "real tool; treat this run as a test of them",
    )


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
    _check_tools(request, command_runner, checks)
    _check_basecalling(request, checks)
    _check_envelope(request, checks)
    _check_disk(request, checks)
    _check_adapters(request, checks)
    return checks.checks

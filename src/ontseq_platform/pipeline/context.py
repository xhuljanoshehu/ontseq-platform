"""The stage contract: what one run is configured with and what each stage plans and returns.

These types are shared by the runner, which sequences stages, and by the analysis lanes that
implement them (:mod:`ontseq_platform.pipeline.marlin`, :mod:`ontseq_platform.cnv.lane`). They
live apart from the runner so that a lane can depend on the contract without importing the
orchestrator that in turn imports every lane.

Everything a run needs is part of :class:`RunConfiguration`. No stage implementation, stage
specification or lane setting is installed into process-global state; a long-running service
therefore executes each run exactly as that run's configuration says.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from ..align import AlignmentPolicy
from ..basecall import BasecallPolicy
from ..execution import StreamingCommandRunner
from ..methylation import MethylationPolicy
from ..models import (
    AmlKnowledgeLock,
    AnalysisModule,
    AssayMode,
    CuteSvPolicy,
    IntervalResourceLock,
    ModuleRunStatus,
    QCPolicy,
    ReferenceLock,
    ResolvedResourceContext,
    SampleManifest,
    SnifflesPolicy,
    SvConsensusPolicy,
    SvEvidencePolicy,
    ToolRecord,
)
from ..target_coverage import TargetCoveragePolicy
from .components import RunComponents
from .envelope import Artifact, RunEnvelope
from .input_digest import RunInputDigestCache
from .stages import SPEC_BY_STAGE, InputKindName, StageId
from .state import StageRecord

if TYPE_CHECKING:
    from ..cnv.lane import CnvLaneSettings


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
    cutesv_policy: CuteSvPolicy | None = None
    sv_consensus_policy: SvConsensusPolicy | None = None
    sv_evidence_policy: SvEvidencePolicy | None = None
    gene_annotation: tuple[Path, IntervalResourceLock] | None = None
    cytoband_annotation: tuple[Path, IntervalResourceLock] | None = None
    sv_context_resources: tuple[tuple[Path, IntervalResourceLock], ...] = ()
    aml_knowledge: tuple[Path, AmlKnowledgeLock] | None = None
    sv_minimum_mean_depth: float = 10.0
    target_coverage_policy: TargetCoveragePolicy | None = None
    methylation_policy: MethylationPolicy | None = None
    marlin_installation: Path | None = None
    #: The copy-number lane of this run (QDNAseq/ACE policy and R runner). ``None`` means the
    #: run has no CNV lane: the stage records that scope statement instead of calling copy
    #: number. Part of the configuration rather than of process state, so two runs in one
    #: service process can differ and neither inherits the other's lane.
    cnv_lane: CnvLaneSettings | None = None
    #: Which component runs each stage, and at which version. ``None`` keeps the built-in
    #: defaults and pins nothing, which is what every run did before selection existed.
    components: RunComponents | None = None
    alignment_policy: AlignmentPolicy | None = None
    basecall_policy: BasecallPolicy | None = None
    reference_fasta: Path | None = None
    pod5_directory: Path | None = None
    threads: int = 4
    # cuteSV rebuilds genome-wide signature lists in each worker. Keep its memory
    # concurrency independent of the other callers and record it in the SV plan.
    cutesv_threads: int = 1
    executables: Mapping[str, str] = field(
        default_factory=lambda: {
            "samtools": "samtools",
            "cramino": "cramino",
            "sniffles": "sniffles",
            "cutesv": "cuteSV",
            "minimap2": "minimap2",
            "mosdepth": "mosdepth",
            "modkit": "modkit",
            "dorado": "dorado",
        }
    )
    #: Ignore any previous run state and execute every stage again.
    force: bool = False
    resource_context: ResolvedResourceContext | None = None
    annotation_cache: Path | None = None
    selection_target_bed: Path | None = None
    context_resource_paths: Mapping[str, Path] = field(default_factory=dict)

    def __post_init__(self) -> None:
        modules = self.manifest.analysis.modules
        if AnalysisModule.METHYLATION in modules and AnalysisModule.MARLIN not in modules:
            self.manifest = self.manifest.model_copy(
                update={
                    "analysis": self.manifest.analysis.model_copy(
                        update={
                            "modules": [*modules, AnalysisModule.MARLIN],
                        }
                    ),
                }
            )
        if (
            isinstance(self.cutesv_threads, bool)
            or not isinstance(self.cutesv_threads, int)
            or self.cutesv_threads < 1
        ):
            raise ValueError("cutesv_threads must be a positive integer")

    def executable(self, name: str) -> str:
        return self.executables.get(name, name)


#: Stages whose outputs the reviewer report reads directly, besides the assembled result.
#: The report renders the CNV plots, the methylation heatmap and table and the MARLIN view
#: from those artifacts, so a change in them must re-render the report even when the
#: assembled result JSON happens to be byte-identical.
_REPORT_EVIDENCE_STAGES = (StageId.CNV, StageId.METHYLATION, StageId.MARLIN)


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
    stage_records: dict[StageId, StageRecord] = field(default_factory=dict)
    input_digests: RunInputDigestCache = field(default_factory=RunInputDigestCache, repr=False)

    def fingerprint_external_input(
        self, path: Path, *, label: str | None = None
    ) -> tuple[str, str]:
        if not path.is_file():
            raise StageFailure("required external input is missing")
        try:
            digest, stable = self.input_digests.digest(path)
        except OSError as exc:
            raise StageFailure("required external input could not be fingerprinted") from exc
        if not stable:
            raise StageFailure(
                f"{label or 'required external input'} changed while it was being fingerprinted"
            )
        return (label or path.name, digest)

    @property
    def sample_id(self) -> str:
        return self.config.manifest.sample_id

    def path(self, template: str) -> str:
        return template.format(sample=self.sample_id)

    def upstream(self, stage: StageId, input_kind: InputKindName) -> list[Artifact]:
        collected: list[Artifact] = []
        dependencies = list(SPEC_BY_STAGE[stage].depends_on)
        if (
            stage in {StageId.SV, StageId.REPORT}
            and self.manifest.assay.mode == AssayMode.ADAPTIVE_SAMPLING
        ):
            dependencies.append(StageId.TARGET_COVERAGE)
        if stage is StageId.ASSEMBLE:
            # Optional evidence affects the assembled result even though a missing lane
            # must not block assembly. Track current outputs in the resume signature;
            # checking old files on disk would retain an earlier opt-in after deselection.
            dependencies.extend((StageId.CNV, StageId.SV, StageId.METHYLATION, StageId.MARLIN))
        if stage is StageId.REPORT:
            dependencies.extend(_REPORT_EVIDENCE_STAGES)
        for dependency in dependencies:
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


def current_artifact(ctx: RunContext, stage: StageId, relative_path: str) -> Artifact | None:
    """Return the artifact ``stage`` recorded at ``relative_path`` for this run, if any.

    The only way a lane's output becomes evidence for a later stage. A file at the same
    path left by an earlier attempt, a failed re-run or a since-deselected stage is not a
    current artifact and is never returned; a recorded artifact that no longer verifies byte
    for byte fails the reading stage.
    """
    current = next(
        (
            artifact
            for artifact in ctx.artifacts.get(stage, [])
            if artifact.relative_path == relative_path
        ),
        None,
    )
    if current is None:
        return None
    if ctx.envelope.verify([current]):
        raise StageFailure(f"The current {stage.value} artifact failed its checksum verification")
    return current

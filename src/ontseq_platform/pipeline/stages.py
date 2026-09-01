"""The pipeline stage graph, its applicability rules and its failure propagation.

This module is the declarative heart of end-to-end execution. It answers three questions
without touching the filesystem, without running a tool and without pydantic, so that all
three can be unit tested in isolation:

1. **Which stages apply at all?** A run that starts from an aligned BAM must not report
   basecalling as failed; basecalling was never required. Applicability is a property of
   the declared input kind, not of the run's outcome.
2. **In what order may stages run?** Dependencies are declared, and the order is derived,
   so adding a stage cannot silently produce a cycle or an unmet dependency.
3. **What happens downstream of a problem?** A stage whose dependency did not complete is
   ``NOT_RUN``, never ``FAILED``. It did not fail; it never started. Collapsing the two
   would turn one real failure into a cascade of apparent failures and hide the cause.

Verification status is declared here as well, next to each stage, because it is a property
of the adapter rather than of a particular run: some adapters are exercised against the
real binary in continuous integration, and some cannot be.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import StrEnum


class StageId(StrEnum):
    """Every stage the pipeline knows about."""

    BASECALL = "basecall"
    ALIGN = "align"
    INTAKE = "intake"
    QC = "qc"
    TARGET_COVERAGE = "target_coverage"
    CNV = "cnv"
    SV = "sv"
    ASSEMBLE = "assemble"
    REPORT = "report"
    RELEASE = "release"


class InputKindName(StrEnum):
    """Mirror of :class:`ontseq_platform.models.InputKind`.

    Duplicated deliberately so this module carries no pydantic dependency. The two are
    kept in step by a test that compares the member values.
    """

    POD5 = "pod5"
    UNALIGNED_BAM = "unaligned_bam"
    ALIGNED_BAM = "aligned_bam"


class VerificationStatus(StrEnum):
    """How much confidence the *adapter* itself has earned.

    This is not about scientific validity, which no stage in this repository has. It is
    the narrower engineering question of whether the code has ever been executed against
    the real tool.
    """

    #: Continuous integration runs the real binary against synthetic input.
    VERIFIED_WITH_REAL_TOOL = "verified_with_real_tool"
    #: Pure-Python stage with no external tool, covered by unit tests.
    VERIFIED_PURE_PYTHON = "verified_pure_python"
    #: The adapter exists and is structured, but has never been run against the real
    #: binary here or in CI. Its behaviour on real input is an assumption.
    UNVERIFIED_ADAPTER = "unverified_adapter"
    #: Declared so the graph is complete, but no implementation is wired in yet.
    NOT_IMPLEMENTED = "not_implemented"


class StageOutcome(StrEnum):
    """Result of a single stage. Mirrors ``ModuleRunStatus`` plus nothing else.

    Deliberately the same vocabulary the repository already uses for modules, so a stage
    record and a module outcome cannot drift apart in meaning.
    """

    COMPLETED = "COMPLETED"
    NOT_RUN = "NOT_RUN"
    FAILED = "FAILED"
    NO_CALL = "NO_CALL"


#: Outcomes that allow a dependent stage to proceed.
SATISFYING_OUTCOMES: frozenset[StageOutcome] = frozenset({StageOutcome.COMPLETED})


@dataclass(frozen=True)
class StageSpec:
    """Declarative description of one stage."""

    stage: StageId
    title: str
    depends_on: tuple[StageId, ...]
    #: Input kinds for which this stage is meaningful at all.
    applicable_for: frozenset[InputKindName]
    verification: VerificationStatus
    #: A stage the run cannot be considered complete without.
    required: bool
    #: Why the stage exists, surfaced in the run report so a reader never has to guess.
    purpose: str
    #: Recorded when the stage is skipped because it does not apply.
    not_applicable_reason: str = ""


_ALL_KINDS: frozenset[InputKindName] = frozenset(InputKindName)
_FROM_UNALIGNED: frozenset[InputKindName] = frozenset(
    {InputKindName.POD5, InputKindName.UNALIGNED_BAM}
)

STAGE_SPECS: tuple[StageSpec, ...] = (
    StageSpec(
        stage=StageId.BASECALL,
        title="Dorado basecalling",
        depends_on=(),
        applicable_for=frozenset({InputKindName.POD5}),
        verification=VerificationStatus.UNVERIFIED_ADAPTER,
        required=True,
        purpose="Convert POD5 signal into an unaligned BAM with modified-base tags retained.",
        not_applicable_reason="The run does not start from POD5, so no basecalling is required.",
    ),
    StageSpec(
        stage=StageId.ALIGN,
        title="Minimap2 alignment",
        depends_on=(StageId.BASECALL,),
        applicable_for=_FROM_UNALIGNED,
        # CI executes this stage against real minimap2 and samtools on a synthetic
        # reference (see ``align_fixture`` and the alignment-lane steps in ci.yml),
        # asserting that reads map, that read groups survive the FASTQ round trip and
        # that modified-base tags are still present afterwards. That job is what earns
        # the claim below; without it this field would have stayed UNVERIFIED_ADAPTER.
        verification=VerificationStatus.VERIFIED_WITH_REAL_TOOL,
        required=True,
        purpose="Align reads to the locked reference and produce a sorted, indexed BAM.",
        not_applicable_reason="The run already starts from an aligned BAM.",
    ),
    StageSpec(
        stage=StageId.INTAKE,
        title="Aligned-BAM integrity gate",
        depends_on=(StageId.ALIGN,),
        applicable_for=_ALL_KINDS,
        verification=VerificationStatus.VERIFIED_WITH_REAL_TOOL,
        required=True,
        purpose="Fail closed on BAM/BAI, sort order, read-group and reference-build problems.",
    ),
    StageSpec(
        stage=StageId.QC,
        title="Cramino read QC",
        depends_on=(StageId.INTAKE,),
        applicable_for=_ALL_KINDS,
        verification=VerificationStatus.VERIFIED_WITH_REAL_TOOL,
        required=True,
        purpose="Normalize descriptive alignment QC without exporting read-level data.",
    ),
    StageSpec(
        stage=StageId.TARGET_COVERAGE,
        title="Adaptive-sampling target coverage",
        depends_on=(StageId.INTAKE,),
        applicable_for=_ALL_KINDS,
        verification=VerificationStatus.VERIFIED_WITH_REAL_TOOL,
        required=False,
        purpose=(
            "Per-target depth over the locked target design. Runs only for "
            "assay.mode=adaptive_sampling; for any other mode the stage records that it "
            "does not apply, which is a scope statement rather than a coverage result."
        ),
    ),
    StageSpec(
        stage=StageId.CNV,
        title="Copy-number evidence",
        depends_on=(StageId.QC,),
        applicable_for=_ALL_KINDS,
        verification=VerificationStatus.NOT_IMPLEMENTED,
        required=False,
        purpose=(
            "Copy-number calling. No production caller is selected; the benchmark "
            "subsystem exists to make that choice on evidence."
        ),
    ),
    StageSpec(
        stage=StageId.SV,
        title="Sniffles2 candidate structural variants",
        depends_on=(StageId.INTAKE,),
        applicable_for=_ALL_KINDS,
        verification=VerificationStatus.VERIFIED_WITH_REAL_TOOL,
        required=False,
        purpose="Conservative, non-reportable candidate SV evidence.",
    ),
    StageSpec(
        stage=StageId.ASSEMBLE,
        title="Normalized result assembly",
        depends_on=(StageId.QC,),
        applicable_for=_ALL_KINDS,
        verification=VerificationStatus.VERIFIED_PURE_PYTHON,
        required=True,
        purpose="Merge every module outcome into one validated result contract.",
    ),
    StageSpec(
        stage=StageId.REPORT,
        title="Reviewer artifacts",
        depends_on=(StageId.ASSEMBLE,),
        applicable_for=_ALL_KINDS,
        verification=VerificationStatus.VERIFIED_PURE_PYTHON,
        required=True,
        purpose="Render JSON, HTML and Excel for human review.",
    ),
    StageSpec(
        stage=StageId.RELEASE,
        title="Immutable release bundle",
        depends_on=(StageId.REPORT,),
        applicable_for=_ALL_KINDS,
        verification=VerificationStatus.VERIFIED_PURE_PYTHON,
        required=True,
        purpose=(
            "Checksum every artifact so the run can be archived and later verified "
            "byte for byte. Signing is a separate, authorised step."
        ),
    ),
)

SPEC_BY_STAGE: Mapping[StageId, StageSpec] = {spec.stage: spec for spec in STAGE_SPECS}


def _validate_graph() -> tuple[StageId, ...]:
    """Topologically order the graph, rejecting cycles and unknown dependencies."""
    known = set(SPEC_BY_STAGE)
    for spec in STAGE_SPECS:
        unknown = set(spec.depends_on) - known
        if unknown:
            raise ValueError(f"stage {spec.stage} depends on unknown stage(s): {sorted(unknown)}")

    ordered: list[StageId] = []
    permanent: set[StageId] = set()
    temporary: set[StageId] = set()

    def visit(stage: StageId) -> None:
        if stage in permanent:
            return
        if stage in temporary:
            raise ValueError(f"stage graph contains a cycle involving {stage}")
        temporary.add(stage)
        for dependency in SPEC_BY_STAGE[stage].depends_on:
            visit(dependency)
        temporary.discard(stage)
        permanent.add(stage)
        ordered.append(stage)

    for spec in STAGE_SPECS:
        visit(spec.stage)
    return tuple(ordered)


#: Stages in a valid execution order, computed once at import so a malformed graph is a
#: startup error rather than a runtime surprise halfway through a run.
EXECUTION_ORDER: tuple[StageId, ...] = _validate_graph()


def applies_to(stage: StageId, input_kind: InputKindName) -> bool:
    """Return whether a stage is meaningful for a run starting from ``input_kind``."""
    return input_kind in SPEC_BY_STAGE[stage].applicable_for


def effective_dependencies(stage: StageId, input_kind: InputKindName) -> tuple[StageId, ...]:
    """Resolve dependencies, transparently skipping stages that do not apply.

    An aligned-BAM run has no ``align`` stage, so ``intake`` must not wait for it. Rather
    than encoding one graph per input kind, a non-applicable dependency is replaced by
    *its* dependencies, recursively. The declared graph therefore stays single and
    readable while every input kind gets a correct one.
    """
    resolved: list[StageId] = []
    for dependency in SPEC_BY_STAGE[stage].depends_on:
        if applies_to(dependency, input_kind):
            resolved.append(dependency)
        else:
            resolved.extend(effective_dependencies(dependency, input_kind))
    # Preserve order while removing duplicates introduced by diamond dependencies.
    seen: set[StageId] = set()
    unique: list[StageId] = []
    for item in resolved:
        if item not in seen:
            seen.add(item)
            unique.append(item)
    return tuple(unique)


def planned_stages(input_kind: InputKindName) -> tuple[StageId, ...]:
    """Return the stages that apply to an input kind, in execution order."""
    return tuple(stage for stage in EXECUTION_ORDER if applies_to(stage, input_kind))


@dataclass(frozen=True)
class BlockedReason:
    """Why a stage cannot start."""

    blocking_stage: StageId
    blocking_outcome: StageOutcome

    def describe(self) -> str:
        return (
            f"dependency {self.blocking_stage.value} did not complete "
            f"(outcome {self.blocking_outcome.value}); this stage was never started"
        )


def blocking_dependency(
    stage: StageId,
    input_kind: InputKindName,
    outcomes: Mapping[StageId, StageOutcome],
) -> BlockedReason | None:
    """Return the first unmet dependency of a stage, or ``None`` when it may run.

    A stage with an unmet dependency becomes ``NOT_RUN``. It is never ``FAILED``: it did
    not fail, it never ran, and a reader tracing a problem needs to reach the one stage
    that actually broke rather than a cascade of look-alikes.
    """
    for dependency in effective_dependencies(stage, input_kind):
        outcome = outcomes.get(dependency)
        if outcome is None or outcome not in SATISFYING_OUTCOMES:
            return BlockedReason(
                blocking_stage=dependency,
                blocking_outcome=outcome if outcome is not None else StageOutcome.NOT_RUN,
            )
    return None


@dataclass(frozen=True)
class RunVerdict:
    """Overall assessment of a run, derived from its stage outcomes."""

    passed: bool
    failed_stages: tuple[StageId, ...]
    incomplete_required_stages: tuple[StageId, ...]
    skipped_optional_stages: tuple[StageId, ...]
    unverified_stages: tuple[StageId, ...] = field(default=())

    def describe(self) -> str:
        if self.failed_stages:
            names = ", ".join(item.value for item in self.failed_stages)
            return f"Run failed: {names} reported FAILED."
        if self.incomplete_required_stages:
            names = ", ".join(item.value for item in self.incomplete_required_stages)
            return f"Run incomplete: required stage(s) {names} did not complete."
        if self.skipped_optional_stages:
            names = ", ".join(item.value for item in self.skipped_optional_stages)
            return (
                f"Run completed with optional stage(s) not run: {names}. "
                "A stage that did not run is not a negative biological finding."
            )
        return "Run completed; every applicable stage finished."


def summarize(
    input_kind: InputKindName,
    outcomes: Mapping[StageId, StageOutcome],
) -> RunVerdict:
    """Derive the run verdict from stage outcomes.

    ``NO_CALL`` counts as a completed stage: the module ran, looked, and correctly
    declined to assert anything. Treating it as an incomplete run would push a legitimate
    scientific answer into an operational failure.
    """
    failed: list[StageId] = []
    incomplete: list[StageId] = []
    skipped: list[StageId] = []
    unverified: list[StageId] = []

    for stage in planned_stages(input_kind):
        spec = SPEC_BY_STAGE[stage]
        outcome = outcomes.get(stage, StageOutcome.NOT_RUN)
        if outcome == StageOutcome.FAILED:
            failed.append(stage)
        elif outcome in {StageOutcome.COMPLETED, StageOutcome.NO_CALL}:
            if spec.verification in {
                VerificationStatus.UNVERIFIED_ADAPTER,
                VerificationStatus.NOT_IMPLEMENTED,
            }:
                unverified.append(stage)
        elif spec.required:
            incomplete.append(stage)
        else:
            skipped.append(stage)

    return RunVerdict(
        passed=not failed and not incomplete,
        failed_stages=tuple(failed),
        incomplete_required_stages=tuple(incomplete),
        skipped_optional_stages=tuple(skipped),
        unverified_stages=tuple(unverified),
    )


def unverified_specs(stages: Iterable[StageId]) -> tuple[StageSpec, ...]:
    """Return the specs of stages whose adapter has never met the real tool."""
    return tuple(
        SPEC_BY_STAGE[stage]
        for stage in stages
        if SPEC_BY_STAGE[stage].verification
        in {VerificationStatus.UNVERIFIED_ADAPTER, VerificationStatus.NOT_IMPLEMENTED}
    )

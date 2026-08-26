"""The vocabulary of preflight checks, and which of them a given run needs.

A run fails stage by stage: a missing reference is discovered after intake, a missing
Dorado model after the envelope exists, a locked envelope only once the lock is attempted.
Each of those is correct behaviour and each of them is expensive to learn late — a POD5
run can spend hours before reaching the stage that was always going to fail.

This module answers "what must be true before the run starts?" as a property of the
declared input kind, in the same spirit as :mod:`ontseq_platform.pipeline.stages`: which
external tools are needed follows from which stages are planned, and whether a missing
tool is fatal follows from whether its stage is required. Nothing here touches the
filesystem or executes anything, so all of it is unit tested without pydantic, a binary or
a reference genome.

Five statuses, and the distinctions between them are the point:

``OK``
    Checked, and the precondition holds.
``FAILED``
    Checked, and the run cannot succeed. Preflight exits non-zero on this alone.
``WARNING``
    Checked, the run can proceed, and somebody should know anyway — an optional caller is
    missing, an unverified adapter is about to run, a previous run died here.
``UNKNOWN``
    Genuinely not determinable from here. A GPU on a compute node, free space against an
    unknown data size. Never collapsed into OK, because "we did not look" and "we looked
    and it was fine" are different claims.
``SKIPPED``
    Does not apply to this input kind. An aligned-BAM run needs no basecalling model.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum

from .stages import SPEC_BY_STAGE, InputKindName, StageId, planned_stages


class CheckStatus(StrEnum):
    """The outcome of a single precondition check."""

    OK = "ok"
    FAILED = "failed"
    WARNING = "warning"
    UNKNOWN = "unknown"
    SKIPPED = "skipped"


@dataclass(frozen=True)
class Check:
    """One precondition, its outcome and — when it went wrong — what to do about it."""

    name: str
    status: CheckStatus
    detail: str
    #: Concrete next step. Empty when there is nothing to do.
    remedy: str = ""
    #: The stage this check protects, when it protects exactly one.
    stage: StageId | None = None

    @property
    def blocking(self) -> bool:
        return self.status is CheckStatus.FAILED


@dataclass(frozen=True)
class ToolRequirement:
    """An external binary a planned run will invoke."""

    name: str
    #: Stages that need it, in execution order.
    stages: tuple[StageId, ...]
    #: False when every stage needing it is optional, so its absence is not fatal.
    required: bool
    #: Version the policies lock it to, when one is declared.
    expected_version: str | None = None

    def with_expected_version(self, version: str | None) -> ToolRequirement:
        return ToolRequirement(
            name=self.name,
            stages=self.stages,
            required=self.required,
            expected_version=version,
        )


#: Which binary each stage invokes. Stages absent from this map need no external tool.
TOOLS_BY_STAGE: Mapping[StageId, tuple[str, ...]] = {
    StageId.BASECALL: ("dorado",),
    StageId.ALIGN: ("minimap2", "samtools"),
    StageId.INTAKE: ("samtools",),
    StageId.QC: ("cramino",),
    StageId.TARGET_COVERAGE: ("mosdepth",),
    StageId.SV: ("sniffles",),
}


def required_tools(input_kind: InputKindName) -> tuple[ToolRequirement, ...]:
    """Derive the binaries a run of this input kind will invoke.

    A tool is *required* when at least one required stage needs it. ``sniffles`` serves
    only the optional SV stage, so a machine without it can still complete a run and
    preflight says so as a warning rather than refusing to start. Deriving this from
    ``StageSpec.required`` rather than listing it here means adding a stage cannot leave
    the two descriptions disagreeing.
    """
    stages_for: dict[str, list[StageId]] = {}
    required_for: dict[str, bool] = {}
    for stage in planned_stages(input_kind):
        for tool in TOOLS_BY_STAGE.get(stage, ()):
            stages_for.setdefault(tool, []).append(stage)
            required_for[tool] = required_for.get(tool, False) or SPEC_BY_STAGE[stage].required
    return tuple(
        ToolRequirement(name=name, stages=tuple(stages), required=required_for[name])
        for name, stages in sorted(stages_for.items())
    )


def worst(checks: Iterable[Check]) -> CheckStatus:
    """The most serious status present, in the order a reader cares about."""
    seen = {check.status for check in checks}
    for status in (
        CheckStatus.FAILED,
        CheckStatus.WARNING,
        CheckStatus.UNKNOWN,
        CheckStatus.OK,
    ):
        if status in seen:
            return status
    return CheckStatus.SKIPPED


def exit_code(checks: Sequence[Check]) -> int:
    """2 when the run cannot succeed, 0 otherwise.

    A warning or an unknown does not block: preflight advises, and refusing to start on
    "we could not check the GPU from here" would make the command unusable on exactly the
    machines it is meant to help.
    """
    return 2 if any(check.blocking for check in checks) else 0


@dataclass
class CheckList:
    """Accumulator, so a caller records checks in the order it makes them."""

    checks: list[Check] = field(default_factory=list)

    def add(
        self,
        name: str,
        status: CheckStatus,
        detail: str,
        *,
        remedy: str = "",
        stage: StageId | None = None,
    ) -> Check:
        check = Check(name=name, status=status, detail=detail, remedy=remedy, stage=stage)
        self.checks.append(check)
        return check

    def ok(self, name: str, detail: str, *, stage: StageId | None = None) -> Check:
        return self.add(name, CheckStatus.OK, detail, stage=stage)

    def failed(
        self, name: str, detail: str, *, remedy: str = "", stage: StageId | None = None
    ) -> Check:
        return self.add(name, CheckStatus.FAILED, detail, remedy=remedy, stage=stage)

    def warning(
        self, name: str, detail: str, *, remedy: str = "", stage: StageId | None = None
    ) -> Check:
        return self.add(name, CheckStatus.WARNING, detail, remedy=remedy, stage=stage)

    def unknown(
        self, name: str, detail: str, *, remedy: str = "", stage: StageId | None = None
    ) -> Check:
        return self.add(name, CheckStatus.UNKNOWN, detail, remedy=remedy, stage=stage)

    def skipped(self, name: str, detail: str, *, stage: StageId | None = None) -> Check:
        return self.add(name, CheckStatus.SKIPPED, detail, stage=stage)


_MARKERS: Mapping[CheckStatus, str] = {
    CheckStatus.OK: "ok  ",
    CheckStatus.FAILED: "FAIL",
    CheckStatus.WARNING: "warn",
    CheckStatus.UNKNOWN: "????",
    CheckStatus.SKIPPED: "--  ",
}


def render_text(checks: Sequence[Check], *, verbose: bool = False) -> str:
    """Render for a person about to start a run.

    Skipped checks are hidden unless asked for: a reader looking at an aligned-BAM run
    does not need six lines saying basecalling is not involved. Everything else is shown,
    including what passed, because a preflight that only prints problems leaves the reader
    unable to tell a clean run from an incomplete check.
    """
    shown = [check for check in checks if verbose or check.status is not CheckStatus.SKIPPED]
    if not shown:
        return "no checks apply"
    lines: list[str] = []
    for check in shown:
        lines.append(f"[{_MARKERS[check.status]}] {check.name:<28} {check.detail}")
        if check.remedy:
            lines.append(f"           -> {check.remedy}")
    counts: dict[str, int] = {}
    for check in checks:
        counts[check.status.value] = counts.get(check.status.value, 0) + 1
    summary = ", ".join(f"{count} {status}" for status, count in sorted(counts.items()))
    lines.append(f"preflight: {worst(checks).value.upper()} — {summary}")
    return "\n".join(lines)


def render_json(checks: Sequence[Check]) -> str:
    """Render the same information for a scheduler deciding whether to submit a job."""
    payload = {
        "verdict": worst(checks).value,
        "blocking": [check.name for check in checks if check.blocking],
        "checks": [
            {
                "name": check.name,
                "status": check.status.value,
                "detail": check.detail,
                "remedy": check.remedy,
                "stage": None if check.stage is None else check.stage.value,
            }
            for check in checks
        ],
    }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"

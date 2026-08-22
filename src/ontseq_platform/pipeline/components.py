"""Per-run component selection: which tool runs a stage, and at exactly which version.

Two things were previously impossible to express. A run could not say *which* caller it
wanted, because the answer was compiled in — the CNV lane was installed by mutating the
runner's implementation table at import time, so the executable graph depended on process
state rather than on the run. And a run could not say which *version* it wanted, because
versions lived inside individual policy files and were only enforced, if at all, by the
adapter that happened to read them.

A :class:`RunComponents` document fixes both. It names, per stage, the provider and the
exact tool version this run is asking for, plus the policy file that parameterises it.
Three consequences follow, and they are the point of the module:

**Selection is data.** Comparing Sniffles 2.4 against 2.8.0 is two files and two runs, not
a code change. The comparison is therefore something a reviewer can reproduce.

**Version drift fails closed.** Every stage already probes the tool it is about to run and
reports the result in its plan. Checking that probe against the declared version in one
place means a run cannot silently execute against whatever happens to be installed. A
mismatch fails the stage with both numbers named, before any output exists.

**Deselection is explicit and visible.** A stage switched off for this run is recorded as
``NOT_RUN`` with the selection named as the reason — never as an absent or negative result.

Stages that run no external tool (assemble, report, release) are deliberately not
selectable: they are the pipeline's own bookkeeping, and making them optional would let a
run produce evidence with no report or no checksum bundle.
"""

from __future__ import annotations

from typing import Literal

from pydantic import Field, model_validator

from ..models import StrictModel
from .stages import StageId

#: Providers this repository has an adapter for, per stage. A selection naming anything
#: else is refused when it is loaded rather than when the stage is reached, so a typo
#: cannot cost a run that has already spent hours on alignment.
SUPPORTED_PROVIDERS: dict[StageId, frozenset[str]] = {
    StageId.BASECALL: frozenset({"dorado"}),
    StageId.ALIGN: frozenset({"minimap2"}),
    StageId.INTAKE: frozenset({"samtools"}),
    StageId.QC: frozenset({"cramino"}),
    StageId.TARGET_COVERAGE: frozenset({"mosdepth"}),
    StageId.CNV: frozenset({"qdnaseq_ace"}),
    StageId.SV: frozenset({"sniffles2"}),
}

#: The key each provider reports its version under in :attr:`StagePlan.tool_versions`.
#: Kept separate from the provider name because a provider is a scientific choice while the
#: key is whatever the adapter happens to call its executable.
PROVIDER_VERSION_KEY: dict[str, str] = {
    "dorado": "dorado",
    "minimap2": "minimap2",
    "samtools": "samtools",
    "cramino": "cramino",
    "mosdepth": "mosdepth",
    "qdnaseq_ace": "QDNAseq",
    "sniffles2": "sniffles",
}


class ComponentVersionMismatch(RuntimeError):
    """The installed tool is not the one this run selected."""


class ComponentChoice(StrictModel):
    """One stage's component for one run."""

    provider: str = Field(min_length=1)
    #: The exact version this run asks for. ``None`` means the run declines to pin, which
    #: is allowed but recorded, because an unpinned run is not reproducible.
    version: str | None = Field(default=None, pattern=r"^\d+\.\d+(?:\.\d+)?$")
    enabled: bool = True
    #: Policy file that parameterises the provider. Relative paths resolve against the
    #: repository root, so a selection file stays portable between machines.
    policy: str | None = None
    note: str | None = None

    def verify(self, observed: dict[str, str], *, stage: StageId) -> None:
        """Raise unless the probed tool version is the one this run selected."""
        if self.version is None:
            return
        key = PROVIDER_VERSION_KEY.get(self.provider, self.provider)
        actual = observed.get(key)
        if actual is None:
            raise ComponentVersionMismatch(
                f"{stage.value}: the selection pins {self.provider} to {self.version}, but the "
                f"stage reported no version for {key!r}"
            )
        if actual != self.version:
            raise ComponentVersionMismatch(
                f"{stage.value}: the selection pins {self.provider} to {self.version}, but the "
                f"installed tool reports {actual}. Refusing to run a selected component at an "
                f"unselected version"
            )


class RunComponents(StrictModel):
    """The component selection for one run."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    selection_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    components: dict[StageId, ComponentChoice] = Field(default_factory=dict)
    note: str | None = None

    @model_validator(mode="after")
    def providers_are_supported(self) -> RunComponents:
        for stage, choice in self.components.items():
            supported = SUPPORTED_PROVIDERS.get(stage)
            if supported is None:
                raise ValueError(
                    f"{stage.value} runs no external tool and cannot be selected per run"
                )
            if choice.provider not in supported:
                allowed = ", ".join(sorted(supported))
                raise ValueError(
                    f"{stage.value}: provider {choice.provider!r} has no adapter in this "
                    f"repository (supported: {allowed})"
                )
        return self

    def choice_for(self, stage: StageId) -> ComponentChoice | None:
        return self.components.get(stage)

    def disabled_stages(self) -> tuple[StageId, ...]:
        return tuple(
            stage for stage, choice in sorted(self.components.items()) if not choice.enabled
        )

    def unpinned_stages(self) -> tuple[StageId, ...]:
        return tuple(
            stage
            for stage, choice in sorted(self.components.items())
            if choice.enabled and choice.version is None
        )

    def without(self, *stages: StageId) -> RunComponents:
        """Return a copy with the named stages switched off."""
        updated = dict(self.components)
        for stage in stages:
            existing = updated.get(stage)
            if existing is None:
                supported = SUPPORTED_PROVIDERS.get(stage)
                if supported is None:
                    raise ValueError(
                        f"{stage.value} runs no external tool and cannot be deselected"
                    )
                updated[stage] = ComponentChoice(
                    provider=sorted(supported)[0],
                    enabled=False,
                    note="Deselected on the command line without a prior selection entry.",
                )
            else:
                updated[stage] = existing.model_copy(update={"enabled": False})
        return self.model_copy(update={"components": updated})

    def summary(self) -> list[str]:
        """One human-readable line per selected stage, for the run log and provenance."""
        lines: list[str] = []
        for stage, choice in sorted(self.components.items()):
            state = "on " if choice.enabled else "off"
            version = choice.version or "unpinned"
            lines.append(f"{state} {stage.value:<16} {choice.provider} {version}")
        return lines

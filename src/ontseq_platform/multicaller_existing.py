from __future__ import annotations

from typing import Literal

from pydantic import Field, model_validator

from .models import StrictModel
from .multicaller_contracts import CallerMode


class ExistingCallerBridge(StrictModel):
    """Reference an already implemented ONTSeq adapter without wrapping or rewriting it."""

    mode_id: CallerMode
    execution_module: str = Field(min_length=1)
    execution_symbol: str | None = Field(default=None, min_length=1)
    evidence_module: str | None = Field(default=None, min_length=1)
    evidence_symbol: str | None = Field(default=None, min_length=1)
    runtime_qualified: bool
    qualification_reason: str = Field(min_length=12)
    retain_native_output: Literal[True] = True
    caller_agreement_is_truth: Literal[False] = False
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_bridge(self) -> ExistingCallerBridge:
        if self.runtime_qualified and self.execution_symbol is None:
            raise ValueError("Runtime-qualified existing bridge requires an execution symbol")
        if (self.evidence_module is None) != (self.evidence_symbol is None):
            raise ValueError("Evidence module and evidence symbol must be declared together")
        return self


EXISTING_CALLER_BRIDGES: dict[CallerMode, ExistingCallerBridge] = {
    CallerMode.QDNASEQ_ACE_MULTIBIN: ExistingCallerBridge(
        mode_id=CallerMode.QDNASEQ_ACE_MULTIBIN,
        execution_module="ontseq_platform.cnv.qdnaseq",
        execution_symbol="run_qdnaseq_ace",
        evidence_module="ontseq_platform.cnv_validation_qdnaseq",
        evidence_symbol="build_qdnaseq_ace_validation_manifest",
        runtime_qualified=True,
        qualification_reason=(
            "Existing QDNAseq+ACE runtime plus the verified Block-4 full-evidence adapter."
        ),
    ),
    CallerMode.SNIFFLES2_STANDARD: ExistingCallerBridge(
        mode_id=CallerMode.SNIFFLES2_STANDARD,
        execution_module="ontseq_platform.sniffles",
        execution_symbol="run_sniffles",
        runtime_qualified=True,
        qualification_reason="Existing standard Sniffles2 runtime and normalized SV adapter.",
    ),
    CallerMode.CUTESV_STANDARD: ExistingCallerBridge(
        mode_id=CallerMode.CUTESV_STANDARD,
        execution_module="ontseq_platform.cutesv",
        execution_symbol="run_cutesv",
        runtime_qualified=True,
        qualification_reason="Existing standard cuteSV runtime and normalized SV adapter.",
    ),
    CallerMode.SNIFFLES2_MOSAIC: ExistingCallerBridge(
        mode_id=CallerMode.SNIFFLES2_MOSAIC,
        execution_module="ontseq_platform.sniffles",
        execution_symbol=None,
        runtime_qualified=False,
        qualification_reason=(
            "Mosaic mode remains a separate catalog entry and needs separate runtime qualification."
        ),
    ),
}


def get_existing_caller_bridge(mode: CallerMode | str) -> ExistingCallerBridge:
    try:
        parsed_mode = mode if isinstance(mode, CallerMode) else CallerMode(mode)
    except ValueError as exc:
        raise ValueError(f"Unknown multi-caller provider mode: {mode!r}") from exc
    try:
        return EXISTING_CALLER_BRIDGES[parsed_mode]
    except KeyError as exc:
        raise ValueError(
            f"{parsed_mode.value} has no bridge to an already implemented ONTSeq adapter"
        ) from exc

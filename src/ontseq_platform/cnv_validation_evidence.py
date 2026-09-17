from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from enum import StrEnum
from pathlib import PurePosixPath, PureWindowsPath
from typing import Any, Literal

from pydantic import Field, field_validator, model_validator

from .cnv_validation_contracts import CnvDataBasis, CnvRepeatKind
from .models import EventType, GenomeBuild, Locus, StrictModel

SHA256 = r"^[0-9a-f]{64}$"
ID = r"^[A-Za-z0-9][A-Za-z0-9._-]{1,127}$"
ROLE = r"^[a-z][a-z0-9_]{0,63}$"


def _validate_json_tree(value: object, *, path: str = "$") -> None:
    """Reject non-JSON and non-finite caller metadata before it enters a content lock."""
    if value is None or isinstance(value, str | bool | int):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"Non-finite numeric value at {path}")
        return
    if isinstance(value, Mapping):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ValueError(f"Caller metadata keys must be strings at {path}")
            _validate_json_tree(item, path=f"{path}.{key}")
        return
    if isinstance(value, list | tuple):
        for index, item in enumerate(value):
            _validate_json_tree(item, path=f"{path}[{index}]")
        return
    raise ValueError(f"Unsupported caller metadata type at {path}: {type(value).__name__}")


def canonical_evidence_sha256(value: StrictModel | Mapping[str, Any]) -> str:
    """Hash canonical JSON while refusing NaN/Inf and arbitrary Python objects."""
    payload: object = value.model_dump(mode="json") if isinstance(value, StrictModel) else dict(value)
    _validate_json_tree(payload)
    rendered = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )
    return hashlib.sha256(rendered.encode("utf-8")).hexdigest()


def _safe_relative_path(value: str) -> str:
    if "\\" in value:
        raise ValueError("Evidence artifact paths must use portable POSIX separators")
    posix = PurePosixPath(value)
    windows = PureWindowsPath(value)
    if value in {"", "."} or posix.is_absolute() or windows.is_absolute() or windows.drive:
        raise ValueError("Evidence artifact path must be relative")
    if ".." in posix.parts or ".." in windows.parts:
        raise ValueError("Evidence artifact path cannot traverse outside its bundle")
    return posix.as_posix()


class CnvEvidenceRecordKind(StrEnum):
    RUN_SUMMARY = "run_summary"
    CALLER_FIT = "caller_fit"
    CALLER_BIN = "caller_bin"
    CALLER_SEGMENT = "caller_segment"
    CALLER_EVENT = "caller_event"
    CHROMOSOME_SUMMARY = "chromosome_summary"


class CnvContributionStatus(StrEnum):
    USED_FOR_PRIMARY_ANALYSIS = "USED_FOR_PRIMARY_ANALYSIS"
    SECONDARY = "SECONDARY"
    EXPLORATORY = "EXPLORATORY"
    EXCLUDED_FROM_PRIMARY_METRIC = "EXCLUDED_FROM_PRIMARY_METRIC"
    NOT_ASSESSABLE = "NOT_ASSESSABLE"
    NO_CALL = "NO_CALL"
    FAILED = "FAILED"


class CnvRunOutcomeState(StrEnum):
    OBSERVED = "OBSERVED"
    BIOLOGICAL_NEGATIVE_ORTHOGONAL_TRUTH = "BIOLOGICAL_NEGATIVE_ORTHOGONAL_TRUTH"
    NO_CALL = "NO_CALL"
    NOT_ASSESSABLE = "NOT_ASSESSABLE"
    FAILED = "FAILED"


class CnvNumericMeasurement(StrictModel):
    """Caller-native scalar retained without assigning biological meaning to it."""

    name: str = Field(pattern=ROLE)
    value: float
    unit: str | None = Field(default=None, min_length=1)
    source_field: str = Field(min_length=1)

    @model_validator(mode="after")
    def finite_value(self) -> CnvNumericMeasurement:
        if not math.isfinite(self.value):
            raise ValueError("CNV evidence measurements must be finite")
        return self


class CnvNativeArtifactReference(StrictModel):
    """Fingerprint-only reference to a caller-native file retained outside table rows."""

    artifact_id: str = Field(pattern=ID)
    role: str = Field(pattern=ROLE)
    relative_path: str = Field(min_length=1)
    sha256: str = Field(pattern=SHA256)
    size_bytes: int = Field(ge=0)
    media_type: str = Field(min_length=1)
    description: str = ""
    sensitive_output: Literal[True] = True
    research_only: Literal[True] = True

    @field_validator("relative_path")
    @classmethod
    def portable_relative_path(cls, value: str) -> str:
        return _safe_relative_path(value)


class CnvFullEvidenceRecord(StrictModel):
    """One immutable caller-native evidence row; exclusion never removes the row."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    record_id: str = Field(pattern=ID)
    registration_sha256: str = Field(pattern=SHA256)

    specimen_id: str = Field(pattern=ID)
    biological_specimen_id: str = Field(pattern=ID)
    caller_id: str = Field(pattern=r"^[a-z][a-z0-9_]{1,63}$")
    caller_version: str = Field(min_length=1)
    adapter_policy_sha256: str = Field(pattern=SHA256)
    execution_identity_sha256: str = Field(pattern=SHA256)

    genome_build: GenomeBuild
    data_basis: CnvDataBasis
    reference_id: str = Field(min_length=1)
    reference_sha256: str = Field(pattern=SHA256)
    input_sha256: str = Field(pattern=SHA256)

    coverage_x: float | None = Field(default=None, ge=0)
    coverage_definition: str | None = None
    tumor_fraction: float | None = Field(default=None, ge=0, le=1)
    tumor_fraction_method: str | None = None
    tumor_fraction_timepoint: str | None = None
    bin_size_kbp: int | None = Field(default=None, gt=0)

    repeat_kind: CnvRepeatKind
    replicate_id: str = Field(pattern=ID)
    repeat_group_id: str | None = Field(default=None, pattern=ID)

    record_kind: CnvEvidenceRecordKind
    native_record_id: str = Field(pattern=ID)
    run_outcome: CnvRunOutcomeState
    contribution_status: CnvContributionStatus
    outcome_reason: str | None = Field(default=None, min_length=3)
    contribution_reason: str | None = Field(default=None, min_length=3)
    truth_support_resource_ids: list[str] = Field(default_factory=list)

    caller_parameters: dict[str, Any]
    caller_parameters_sha256: str = Field(pattern=SHA256)
    dependency_versions: dict[str, str] = Field(default_factory=dict)
    native_artifact_ids: list[str] = Field(default_factory=list)

    fit_group_id: str | None = Field(default=None, pattern=ID)
    selected_fit: bool | None = None
    cellularity: float | None = Field(default=None, ge=0, le=1)
    ploidy: float | None = Field(default=None, gt=0)
    fit_error: float | None = Field(default=None, ge=0)

    event_type: EventType | None = None
    primary: Locus | None = None
    copy_number: float | None = Field(default=None, ge=0)
    raw_depth: float | None = Field(default=None, ge=0)
    corrected_depth: float | None = Field(default=None, ge=0)
    log2_ratio: float | None = None
    numeric_measurements: list[CnvNumericMeasurement] = Field(default_factory=list)

    sensitive_output: Literal[True] = True
    research_only: Literal[True] = True

    @property
    def address_key(self) -> str:
        """Deterministic analytical address; record_id is intentionally not the key."""
        payload = {
            "specimen_id": self.specimen_id,
            "caller_id": self.caller_id,
            "caller_version": self.caller_version,
            "genome_build": self.genome_build.value,
            "data_basis": self.data_basis.value,
            "coverage_x": self.coverage_x,
            "tumor_fraction": self.tumor_fraction,
            "bin_size_kbp": self.bin_size_kbp,
            "record_kind": self.record_kind.value,
            "native_record_id": self.native_record_id,
            "replicate_id": self.replicate_id,
        }
        return canonical_evidence_sha256(payload)

    @field_validator("truth_support_resource_ids", "native_artifact_ids")
    @classmethod
    def unique_links(cls, value: list[str]) -> list[str]:
        if len(value) != len(set(value)):
            raise ValueError("Evidence links must be unique")
        return value

    @field_validator("dependency_versions")
    @classmethod
    def dependency_versions_are_explicit(cls, value: dict[str, str]) -> dict[str, str]:
        if any(not key or not version for key, version in value.items()):
            raise ValueError("Dependency names and versions must be non-empty")
        return value

    @model_validator(mode="after")
    def coherent_evidence(self) -> CnvFullEvidenceRecord:
        if self.coverage_x is not None and not math.isfinite(self.coverage_x):
            raise ValueError("Coverage must be finite")
        if self.coverage_x is None and self.coverage_definition is not None:
            raise ValueError("Coverage definition cannot exist without measured coverage")
        if self.coverage_x is not None and not self.coverage_definition:
            raise ValueError("Measured coverage requires its definition")

        if self.tumor_fraction is not None and not math.isfinite(self.tumor_fraction):
            raise ValueError("Tumour fraction must be finite")
        fraction_metadata = [self.tumor_fraction_method, self.tumor_fraction_timepoint]
        if self.tumor_fraction is None and any(item is not None for item in fraction_metadata):
            raise ValueError("Tumour-fraction metadata cannot imply a missing fraction")
        if self.tumor_fraction is not None and any(not item for item in fraction_metadata):
            raise ValueError("Measured tumour fraction requires method and timepoint")

        for name, value in {
            "cellularity": self.cellularity,
            "ploidy": self.ploidy,
            "fit_error": self.fit_error,
            "copy_number": self.copy_number,
            "raw_depth": self.raw_depth,
            "corrected_depth": self.corrected_depth,
            "log2_ratio": self.log2_ratio,
        }.items():
            if value is not None and not math.isfinite(value):
                raise ValueError(f"{name} must be finite")

        _validate_json_tree(self.caller_parameters, path="$.caller_parameters")
        if canonical_evidence_sha256(self.caller_parameters) != self.caller_parameters_sha256:
            raise ValueError("Caller parameters do not match caller_parameters_sha256")

        if self.repeat_kind == CnvRepeatKind.INDEPENDENT and self.repeat_group_id is not None:
            raise ValueError("Independent evidence cannot declare a repeat group")
        if self.repeat_kind != CnvRepeatKind.INDEPENDENT and self.repeat_group_id is None:
            raise ValueError("Repeated evidence requires repeat_group_id")

        terminal_contribution = {
            CnvRunOutcomeState.FAILED: CnvContributionStatus.FAILED,
            CnvRunOutcomeState.NO_CALL: CnvContributionStatus.NO_CALL,
            CnvRunOutcomeState.NOT_ASSESSABLE: CnvContributionStatus.NOT_ASSESSABLE,
        }
        required_contribution = terminal_contribution.get(self.run_outcome)
        if required_contribution is not None:
            if self.contribution_status != required_contribution:
                raise ValueError("Terminal run outcome requires its matching contribution status")
            if not self.outcome_reason:
                raise ValueError("Terminal run outcome requires an explicit reason")
        elif self.contribution_status in set(terminal_contribution.values()):
            raise ValueError("Observed evidence cannot use a terminal contribution status")

        if self.run_outcome == CnvRunOutcomeState.BIOLOGICAL_NEGATIVE_ORTHOGONAL_TRUTH:
            if self.record_kind != CnvEvidenceRecordKind.RUN_SUMMARY:
                raise ValueError("Orthogonal biological negative must be a run-summary record")
            if not self.truth_support_resource_ids:
                raise ValueError("Orthogonal biological negative requires truth resource linkage")
            if not self.outcome_reason:
                raise ValueError("Orthogonal biological negative requires an explicit reason")

        if (
            self.contribution_status == CnvContributionStatus.EXCLUDED_FROM_PRIMARY_METRIC
            and not self.contribution_reason
        ):
            raise ValueError("Excluded evidence requires an explicit contribution reason")

        if (
            self.record_kind == CnvEvidenceRecordKind.CALLER_FIT
            and (self.fit_group_id is None or self.selected_fit is None)
        ):
            raise ValueError("Caller-fit evidence requires fit_group_id and selected_fit")

        return self

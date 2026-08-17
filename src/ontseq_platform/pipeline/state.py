"""Persisted run state: what ran, what it produced, and what remains unproven.

The run report is the operational counterpart to ``PipelineResult``. Where the result
contract describes *biology*, this describes *execution*: which stages ran, in what order,
against which tool versions, producing which checksummed artifacts, and which of them rest
on an adapter that has never met the real tool.

It is written after every stage rather than only at the end, so an interrupted run leaves
a truthful record of how far it got and can be resumed from it.
"""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Literal

from pydantic import Field, model_validator

from ..models import (
    GenomeBuild,
    InputKind,
    ModuleRunStatus,
    SampleManifest,
    StrictModel,
    ToolRecord,
)
from .envelope import Artifact
from .stages import StageId, VerificationStatus


class ArtifactRecord(StrictModel):
    """A checksummed artifact, addressed relative to the run envelope."""

    relative_path: str = Field(min_length=1)
    size_bytes: int = Field(ge=0)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    #: False for intermediate data that must not leave the execution system.
    exportable: bool

    @model_validator(mode="after")
    def path_is_relative(self) -> ArtifactRecord:
        if self.relative_path.startswith("/") or ".." in self.relative_path.split("/"):
            raise ValueError(
                "artifact paths must be relative to the run envelope; absolute source "
                "paths must never appear in reviewer artifacts"
            )
        return self

    @classmethod
    def of(cls, artifact: Artifact) -> ArtifactRecord:
        return cls(
            relative_path=artifact.relative_path,
            size_bytes=artifact.size_bytes,
            sha256=artifact.sha256,
            exportable=artifact.exportable,
        )

    def to_artifact(self) -> Artifact:
        return Artifact(
            relative_path=self.relative_path,
            size_bytes=self.size_bytes,
            sha256=self.sha256,
            exportable=self.exportable,
        )


class StageRecord(StrictModel):
    """Outcome of one stage."""

    stage: StageId
    title: str
    status: ModuleRunStatus
    verification: VerificationStatus
    required: bool
    reason: str = Field(min_length=1)
    #: Content-addressed signature of the stage's inputs; ``None`` when it never ran.
    signature: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    #: True when a previous run's output was reused unchanged.
    resumed: bool = False
    started_at: datetime | None = None
    finished_at: datetime | None = None
    duration_seconds: float | None = Field(default=None, ge=0)
    outputs: list[ArtifactRecord] = Field(default_factory=list)
    tools: list[ToolRecord] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def outcome_is_self_consistent(self) -> StageRecord:
        if self.status == ModuleRunStatus.COMPLETED and not self.signature:
            raise ValueError("a completed stage must record the signature it completed under")
        if self.status in {ModuleRunStatus.NOT_RUN, ModuleRunStatus.FAILED} and self.outputs:
            raise ValueError(f"a {self.status.value} stage must not claim outputs")
        if self.resumed and self.status != ModuleRunStatus.COMPLETED:
            raise ValueError("only a completed stage can be resumed")
        return self


class RunReport(StrictModel):
    """The complete execution record of one run."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    run_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    sample_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$")
    input_kind: InputKind
    genome_build: GenomeBuild
    manifest: SampleManifest
    passed: bool
    verdict_reason: str = Field(min_length=1)
    stages: list[StageRecord] = Field(min_length=1)
    pipeline_version: str
    git_commit: str
    started_at: datetime
    finished_at: datetime
    #: Stages that completed on an adapter never executed against the real tool.
    unverified_stages: list[StageId] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def stages_are_unique_and_consistent(self) -> RunReport:
        seen = [item.stage for item in self.stages]
        if len(seen) != len(set(seen)):
            raise ValueError("run report contains duplicate stage records")
        if self.finished_at < self.started_at:
            raise ValueError("run finished before it started")
        completed_unverified = {
            item.stage
            for item in self.stages
            if item.status in {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL}
            and item.verification
            in {VerificationStatus.UNVERIFIED_ADAPTER, VerificationStatus.NOT_IMPLEMENTED}
        }
        if set(self.unverified_stages) != completed_unverified:
            raise ValueError(
                "unverified_stages must list exactly the completed stages whose adapter "
                "has not been executed against the real tool"
            )
        return self

    def record_for(self, stage: StageId) -> StageRecord | None:
        for item in self.stages:
            if item.stage == stage:
                return item
        return None


class ReleaseBundle(StrictModel):
    """An immutable, verifiable inventory of a run's exportable artifacts.

    Non-exportable artifacts are listed by path and deliberately carry no checksum: the
    bundle documents that intermediate data existed and stayed behind, without becoming a
    vehicle for moving it.

    Signing is not performed here. An unsigned bundle that says so is honest; a bundle
    with a fabricated signature field would not be.
    """

    schema_version: Literal["0.1.0"] = "0.1.0"
    run_id: str
    sample_id: str
    pipeline_version: str
    git_commit: str
    run_report_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    artifacts: list[ArtifactRecord] = Field(min_length=1)
    withheld_artifact_paths: list[str] = Field(default_factory=list)
    total_bytes: int = Field(ge=0)
    signature_status: Literal["unsigned"] = "unsigned"
    signature_note: str = (
        "This bundle is unsigned. Electronic signature requires an authorised key, a "
        "named responsible person and institutional change control."
    )
    release_status: Literal["REVIEW_REQUIRED"] = "REVIEW_REQUIRED"
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def bundle_is_consistent(self) -> ReleaseBundle:
        if any(not item.exportable for item in self.artifacts):
            raise ValueError("a release bundle must not contain non-exportable artifacts")
        if self.total_bytes != sum(item.size_bytes for item in self.artifacts):
            raise ValueError("release bundle total_bytes does not match its artifacts")
        paths = [item.relative_path for item in self.artifacts]
        if len(paths) != len(set(paths)):
            raise ValueError("release bundle lists an artifact twice")
        return self

    def checksum_manifest(self) -> str:
        """Render a ``sha256sum``-compatible manifest for offline verification."""
        return (
            "\n".join(
                f"{item.sha256}  {item.relative_path}"
                for item in sorted(self.artifacts, key=lambda entry: entry.relative_path)
            )
            + "\n"
        )

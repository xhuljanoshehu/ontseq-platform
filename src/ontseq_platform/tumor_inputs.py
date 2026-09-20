from __future__ import annotations

from typing import Literal

from pydantic import Field, model_validator

from .models import GenomeBuild, StrictModel
from .multicaller_contracts import CallerInputRole

SHA256 = r"^[0-9a-f]{64}$"
ID = r"^[A-Za-z0-9][A-Za-z0-9._-]{1,127}$"

TUMOR_AUXILIARY_INPUT_ROLES = frozenset(
    {
        CallerInputRole.TUMOR_BAM,
        CallerInputRole.MATCHED_NORMAL_BAM,
        CallerInputRole.SNP_VCF,
        CallerInputRole.BAF_TABLE,
        CallerInputRole.PHASED_VARIANTS,
        CallerInputRole.PANEL_OF_NORMALS,
        CallerInputRole.SEVERUS_BREAKPOINTS,
    }
)


class TumorAuxiliaryInputArtifact(StrictModel):
    """Fingerprint-bound auxiliary input for tumor-aware multi-caller lanes."""

    artifact_id: str = Field(pattern=ID)
    role: CallerInputRole
    analysis_sample_id: str = Field(pattern=ID)
    source_sample_id: str = Field(pattern=ID)
    genome_build: GenomeBuild
    reference_id: str = Field(min_length=1)
    reference_sha256: str = Field(pattern=SHA256)
    sha256: str = Field(pattern=SHA256)
    producer_lane_id: str | None = Field(default=None, pattern=ID)
    producer_lane_sha256: str | None = Field(default=None, pattern=SHA256)
    description: str = ""
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_auxiliary_input(self) -> TumorAuxiliaryInputArtifact:
        if self.role not in TUMOR_AUXILIARY_INPUT_ROLES:
            raise ValueError(f"{self.role.value} is not a tumor auxiliary input role")
        if (self.producer_lane_id is None) != (self.producer_lane_sha256 is None):
            raise ValueError("Producer lane ID and producer lane SHA-256 must be declared together")
        return self


class TumorInputBundle(StrictModel):
    """One tumor-analysis identity with explicitly bound auxiliary artifacts."""

    analysis_sample_id: str = Field(pattern=ID)
    genome_build: GenomeBuild
    reference_id: str = Field(min_length=1)
    reference_sha256: str = Field(pattern=SHA256)
    artifacts: list[TumorAuxiliaryInputArtifact] = Field(min_length=1)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_bundle(self) -> TumorInputBundle:
        artifact_ids = [item.artifact_id for item in self.artifacts]
        roles = [item.role for item in self.artifacts]
        if len(artifact_ids) != len(set(artifact_ids)):
            raise ValueError("Tumor auxiliary input artifact IDs must be unique")
        if len(roles) != len(set(roles)):
            raise ValueError("Tumor auxiliary input roles must be unique")

        for artifact in self.artifacts:
            if artifact.analysis_sample_id != self.analysis_sample_id:
                raise ValueError("Tumor auxiliary input analysis sample does not match bundle")
            if artifact.genome_build != self.genome_build:
                raise ValueError("Tumor auxiliary input genome build does not match bundle")
            if artifact.reference_id != self.reference_id:
                raise ValueError("Tumor auxiliary input reference ID does not match bundle")
            if artifact.reference_sha256 != self.reference_sha256:
                raise ValueError("Tumor auxiliary input reference SHA-256 does not match bundle")
        return self

    def artifact_for(self, role: CallerInputRole) -> TumorAuxiliaryInputArtifact | None:
        return next((item for item in self.artifacts if item.role == role), None)

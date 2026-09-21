"""Isolated, versioned Mutect2 research contracts; no clinical-release authority."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, StrictBool, field_validator, model_validator

DIGEST = r"^[0-9a-f]{64}$"
IDENTIFIER = r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$"
# Explicit compatibility baseline, not a claim of latest release.
GATK_VERSION: Literal["4.6.2.0"] = "4.6.2.0"


class GatkModel(BaseModel):
    # Kept independent of ONTSeq's CNV/SV schema until the canonical integration is reviewed.
    model_config = ConfigDict(extra="forbid", frozen=True, allow_inf_nan=False)


class LockedFile(GatkModel):
    path: Path
    sha256: str = Field(pattern=DIGEST)

    @field_validator("path", mode="before")
    @classmethod
    def local_absolute_path(cls, value: object) -> str:
        if not isinstance(value, (str, Path)):
            raise ValueError("A local absolute path is required")
        text = str(value)
        if "://" in text or any(ord(c) < 32 for c in text):
            raise ValueError("URLs and control characters are not accepted in file paths")
        path = Path(text)
        if not path.is_absolute() or ".." in path.parts:
            raise ValueError("A local absolute path without parent traversal is required")
        # Keep JSON input a string; Pydantic constructs Path after this validator.
        return text


class ReferenceBundle(GatkModel):
    reference_id: str = Field(pattern=IDENTIFIER)
    genome_build: Literal["GRCh37", "GRCh38", "synthetic"]
    fasta: LockedFile
    fai: LockedFile
    dictionary: LockedFile


class BamSample(GatkModel):
    sample_id: str = Field(pattern=IDENTIFIER)
    bam: LockedFile
    bai: LockedFile
    reference_sha256: str = Field(pattern=DIGEST)


class VcfResource(GatkModel):
    name: str = Field(pattern=IDENTIFIER)
    vcf: LockedFile
    index: LockedFile
    reference_sha256: str = Field(pattern=DIGEST)
    kind: Literal["population", "ont_pon"]
    assay_id: str | None = Field(default=None, pattern=IDENTIFIER)


class RuntimeLock(GatkModel):
    gatk_jar: LockedFile
    expected_version: Literal["4.6.2.0"] = GATK_VERSION
    java_executable: str = "java"
    java_major: Literal[17] = 17
    heap_gb: int = Field(default=8, ge=1, le=128, strict=True)
    pair_hmm_threads: int = Field(default=2, ge=1, le=64, strict=True)

    @field_validator("java_executable")
    @classmethod
    def java_path_or_name(cls, value: str) -> str:
        if value != "java" and not Path(value).is_absolute():
            raise ValueError("java_executable must be java or an absolute executable path")
        if any(ord(c) < 32 for c in value):
            raise ValueError("Control characters are not allowed in executable paths")
        return value


class Mutect2Config(GatkModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    run_id: str = Field(pattern=IDENTIFIER)
    assay_id: str = Field(pattern=IDENTIFIER)
    mode: Literal["tumor_only", "tumor_normal"]
    enabled: StrictBool = False
    allow_experimental_ont: StrictBool = False
    reference: ReferenceBundle
    tumor: BamSample
    normal: BamSample | None = None
    germline_resource: VcfResource
    panel_of_normals: VcfResource | None = None
    estimate_contamination: StrictBool = False
    contamination_sites: VcfResource | None = None
    intervals: LockedFile | None = None
    runtime: RuntimeLock
    step_timeout_seconds: int = Field(default=14400, ge=1, le=604800, strict=True)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_inputs(self) -> Mutect2Config:
        if (self.mode == "tumor_normal") != (self.normal is not None):
            raise ValueError("tumor_normal requires a normal; tumor_only forbids a normal")
        if self.normal and (
            self.normal.sample_id == self.tumor.sample_id
            or self.normal.bam.path == self.tumor.bam.path
            or self.normal.bam.sha256 == self.tumor.bam.sha256
        ):
            raise ValueError("Tumor and normal must have distinct sample IDs and BAMs")
        reference_digest = self.reference.fasta.sha256
        subjects: list[BamSample | VcfResource] = [self.tumor, self.germline_resource]
        if self.normal:
            subjects.append(self.normal)
        if self.panel_of_normals:
            subjects.append(self.panel_of_normals)
        if self.contamination_sites:
            subjects.append(self.contamination_sites)
        if any(item.reference_sha256 != reference_digest for item in subjects):
            raise ValueError("All inputs must declare the same reference SHA-256")
        if self.germline_resource.kind != "population":
            raise ValueError("Germline resource must be an allele-frequency population resource")
        if self.panel_of_normals and (
            self.panel_of_normals.kind != "ont_pon"
            or self.panel_of_normals.assay_id != self.assay_id
        ):
            raise ValueError("Panel of normals must be an ONT PoN matched to this assay ID")
        if self.estimate_contamination != (self.contamination_sites is not None):
            raise ValueError("Contamination estimation and its sites resource must be set together")
        if self.contamination_sites and self.contamination_sites.kind != "population":
            raise ValueError("Contamination sites must be a population allele-frequency resource")
        return self


class SmallVariantCandidate(GatkModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    evidence_id: str = Field(pattern=DIGEST)
    run_id: str
    sample_id: str
    chromosome: str
    start0: int = Field(ge=0)
    end0: int = Field(ge=1)
    native_pos1: int = Field(ge=1)
    ref: str
    alt: str
    allele_index: int = Field(ge=1)
    variant_type: Literal["SNV", "MNV", "INS", "DEL", "COMPLEX"]
    tumor_af: float | None = Field(default=None, ge=0, le=1)
    normal_af: float | None = Field(default=None, ge=0, le=1)
    ref_depth: int | None = Field(default=None, ge=0)
    alt_depth: int | None = Field(default=None, ge=0)
    depth: int | None = Field(default=None, ge=0)
    native_site_filter: str
    native_allele_filter_status: str | None = None
    native_record_sha256: str = Field(pattern=DIGEST)
    source_vcf_sha256: str = Field(pattern=DIGEST)
    source_line_number: int = Field(ge=1)
    caller: Literal["GATK Mutect2"] = "GATK Mutect2"
    coordinate_system: Literal["zero_based_half_open"] = "zero_based_half_open"
    left_aligned: Literal[False] = False
    research_only: Literal[True] = True
    clinically_reportable: Literal[False] = False


class Mutect2Result(GatkModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    run_id: str
    sample_id: str
    mode: Literal["tumor_only", "tumor_normal"]
    status: Literal["NOT_RUN", "COMPLETED", "FAILED"]
    messages: tuple[str, ...] = ()
    candidates: tuple[SmallVariantCandidate, ...] = ()
    config_sha256: str = Field(pattern=DIGEST)
    input_sha256s: dict[str, str] = Field(default_factory=dict)
    output_sha256s: dict[str, str] = Field(default_factory=dict)
    step_exit_codes: dict[str, int] = Field(default_factory=dict)
    execution_evidence: Literal["NONE", "TEST_DOUBLE", "NATIVE_EXECUTION"] = "NONE"
    expected_gatk_version: Literal["4.6.2.0"] = GATK_VERSION
    biological_interpretation: Literal["NOT_ASSESSED"] = "NOT_ASSESSED"
    no_biological_negative_claim: Literal[True] = True
    caller_agreement_is_truth: Literal[False] = False
    analytically_validated: Literal[False] = False
    clinically_reportable: Literal[False] = False
    sensitive_local_only: Literal[True] = True
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def no_candidates_without_completion(self) -> Mutect2Result:
        if self.status != "COMPLETED" and self.candidates:
            raise ValueError("Only a completed lane can expose candidate evidence")
        return self

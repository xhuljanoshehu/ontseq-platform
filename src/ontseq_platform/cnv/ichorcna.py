from __future__ import annotations

import csv
import math
import os
import re
import shutil
import tempfile
from enum import StrEnum
from pathlib import Path, PurePosixPath
from typing import Literal

from pydantic import Field, field_validator, model_validator

from ..execution import CommandRunner, SubprocessRunner
from ..models import (
    EventType,
    FileFingerprint,
    GenomeBuild,
    ModuleRunStatus,
    StrictModel,
    ToolRecord,
)
from ..multicaller_contracts import CallerAssayRegime
from ..reference import sha256_file

_SHA256 = r"^[0-9a-f]{64}$"
_SAFE_ID = r"^[A-Za-z0-9][A-Za-z0-9._-]{1,127}$"
_VERSION = re.compile(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)")


class IchorCnaResourceRole(StrEnum):
    COVERAGE_WIG = "coverage_wig"
    GC_WIG = "gc_wig"
    MAPPABILITY_WIG = "mappability_wig"
    CENTROMERE = "centromere"
    NORMAL_PANEL = "normal_panel"


class IchorCnaResourceArtifact(StrictModel):
    resource_id: str = Field(pattern=_SAFE_ID)
    role: IchorCnaResourceRole
    genome_build: GenomeBuild
    bin_size_bp: int | None = Field(default=None, ge=1)
    size_bytes: int = Field(ge=0)
    sha256: str = Field(pattern=_SHA256)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def role_matches_bin_size_semantics(self) -> IchorCnaResourceArtifact:
        if self.role == IchorCnaResourceRole.CENTROMERE:
            if self.bin_size_bp is not None:
                raise ValueError("ichorCNA centromere resource must not declare a bin size")
        elif self.bin_size_bp is None:
            raise ValueError(f"ichorCNA {self.role.value} resource requires a bin size")
        return self


class IchorCnaResourceBundle(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str = Field(pattern=_SAFE_ID)
    assay_regime: CallerAssayRegime
    genome_build: GenomeBuild
    bin_size_bp: int = Field(ge=1)
    resources: list[IchorCnaResourceArtifact] = Field(min_length=4)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_resource_bundle(self) -> IchorCnaResourceBundle:
        if self.assay_regime != CallerAssayRegime.CFDNA_ULP_WGS:
            raise ValueError("ichorCNA resource bundle is restricted to cfDNA ULP-WGS")
        resource_ids = [item.resource_id for item in self.resources]
        roles = [item.role for item in self.resources]
        if len(resource_ids) != len(set(resource_ids)):
            raise ValueError("ichorCNA resource IDs must be unique")
        if len(roles) != len(set(roles)):
            raise ValueError("ichorCNA resource roles must be unique")

        required = {
            IchorCnaResourceRole.COVERAGE_WIG,
            IchorCnaResourceRole.GC_WIG,
            IchorCnaResourceRole.MAPPABILITY_WIG,
            IchorCnaResourceRole.CENTROMERE,
        }
        missing = required.difference(roles)
        if missing:
            rendered = ", ".join(sorted(item.value for item in missing))
            raise ValueError(f"ichorCNA resource bundle is missing required roles: {rendered}")

        for resource in self.resources:
            if resource.genome_build != self.genome_build:
                raise ValueError("ichorCNA resource genome build does not match bundle")
            if (
                resource.role != IchorCnaResourceRole.CENTROMERE
                and resource.bin_size_bp != self.bin_size_bp
            ):
                raise ValueError("ichorCNA resource bin size does not match bundle")
        return self

    def resource_for(self, role: IchorCnaResourceRole) -> IchorCnaResourceArtifact:
        for resource in self.resources:
            if resource.role == role:
                return resource
        raise ValueError(f"ichorCNA resource bundle does not contain {role.value}")


class IchorCnaUlpWgsPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["research_acceptance_candidate_unvalidated"] = (
        "research_acceptance_candidate_unvalidated"
    )
    expected_version: str = Field(pattern=r"^\d+\.\d+(?:\.\d+)?$")
    assay_regime: CallerAssayRegime
    genome_build: GenomeBuild
    genome_build_label: Literal["hg19", "hg38"]
    genome_style: Literal["NCBI", "UCSC"]
    bin_size_bp: int = Field(ge=1)
    script_sha256: str = Field(pattern=_SHA256)
    ploidy_candidates: list[int] = Field(min_length=1)
    normal_fraction_candidates: list[float] = Field(min_length=1)
    max_copy_number: int = Field(ge=2)
    minimum_mappability: float = Field(ge=0, le=1)
    estimate_normal: bool
    estimate_ploidy: bool
    estimate_subclone_prevalence: bool
    include_homozygous_deletion: bool
    transaction_e: float = Field(gt=0, lt=1)
    transaction_strength: int = Field(ge=1)
    timeout_seconds: int = Field(ge=60)
    real_tool_qualified: bool = False
    analytical_validation: Literal["not_validated"] = "not_validated"
    note: str = Field(min_length=12)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_policy(self) -> IchorCnaUlpWgsPolicy:
        if self.expected_version != "0.5.1":
            raise ValueError("ichorCNA expected_version must match the pinned 0.5.1 contract")
        if self.assay_regime != CallerAssayRegime.CFDNA_ULP_WGS:
            raise ValueError("ichorCNA policy is restricted to cfDNA ULP-WGS")
        expected_label = "hg38" if self.genome_build == GenomeBuild.GRCH38 else "hg19"
        if self.genome_build_label != expected_label:
            raise ValueError("ichorCNA genome build label does not match genome build")
        if len(self.ploidy_candidates) != len(set(self.ploidy_candidates)):
            raise ValueError("ichorCNA ploidy candidates must be unique")
        if any(value < 1 for value in self.ploidy_candidates):
            raise ValueError("ichorCNA ploidy candidates must be positive")
        if len(self.normal_fraction_candidates) != len(set(self.normal_fraction_candidates)):
            raise ValueError("ichorCNA normal fraction candidates must be unique")
        if any(not 0 <= value <= 1 for value in self.normal_fraction_candidates):
            raise ValueError("ichorCNA normal fraction candidates must be between zero and one")
        return self


class IchorCnaSolution(StrictModel):
    initialization: str = Field(min_length=1)
    normal_fraction: float = Field(ge=0, le=1)
    ploidy: float = Field(gt=0)
    bic: float | None = None
    fraction_genome_subclonal: float | None = Field(default=None, ge=0, le=1)
    fraction_cna_subclonal: float | None = Field(default=None, ge=0, le=1)
    log_likelihood: float | None = None


class IchorCnaSegment(StrictModel):
    chromosome: str = Field(pattern=r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    marker_count: int | None = Field(default=None, ge=0)
    median_log_r: float | None = None
    copy_number: float = Field(ge=0)
    native_call: str = Field(min_length=1)
    subclone_status: bool | None = None
    log_r_copy_number: float | None = None
    event_type: EventType | None = None
    reportable: Literal[False] = False
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def valid_interval(self) -> IchorCnaSegment:
        if self.end <= self.start:
            raise ValueError("ichorCNA segment end must be greater than start")
        return self


class IchorCnaNativeArtifact(StrictModel):
    relative_path: str = Field(min_length=1)
    role: str = Field(min_length=1)
    media_type: str = Field(min_length=1)
    fingerprint: FileFingerprint
    sensitive_output: Literal[True] = True
    exportable: Literal[False] = False
    research_only: Literal[True] = True

    @field_validator("relative_path")
    @classmethod
    def safe_relative_path(cls, value: str) -> str:
        if "\\" in value:
            raise ValueError("ichorCNA artifact paths must use POSIX separators")
        path = PurePosixPath(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("ichorCNA artifact path must remain inside output directory")
        return value


class IchorCnaReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str = Field(pattern=_SAFE_ID)
    assay_regime: CallerAssayRegime
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: IchorCnaUlpWgsPolicy
    tumor_fraction: float = Field(ge=0, le=1)
    ploidy: float = Field(gt=0)
    coverage_x: float | None = Field(default=None, ge=0)
    solutions: list[IchorCnaSolution] = Field(default_factory=list)
    segments: list[IchorCnaSegment] = Field(min_length=1)
    tool: ToolRecord
    script_fingerprint: FileFingerprint
    coverage_wig_fingerprint: FileFingerprint
    gc_wig_fingerprint: FileFingerprint
    mappability_wig_fingerprint: FileFingerprint
    centromere_fingerprint: FileFingerprint
    normal_panel_fingerprint: FileFingerprint | None = None
    native_artifacts: list[IchorCnaNativeArtifact] = Field(min_length=1)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def coherent_report(self) -> IchorCnaReport:
        if self.assay_regime != CallerAssayRegime.CFDNA_ULP_WGS:
            raise ValueError("ichorCNA report is restricted to cfDNA ULP-WGS")
        if self.assay_regime != self.policy.assay_regime:
            raise ValueError("ichorCNA report assay regime must match policy")
        if self.status != ModuleRunStatus.COMPLETED:
            raise ValueError("ichorCNA report currently represents completed executions only")
        return self


def _r_vector(values: list[int] | list[float]) -> str:
    return "c(" + ",".join(str(value) for value in values) + ")"


def _chromosome_expression(policy: IchorCnaUlpWgsPolicy, *, training: bool) -> str:
    if policy.genome_style == "UCSC":
        suffix = "c(1:22)" if training else 'c(1:22, "X")'
        return f"paste0('chr', {suffix})"
    return "c(1:22)" if training else 'c(1:22, "X")'


def build_ichorcna_argv(
    *,
    rscript: str,
    ichorcna_script: Path,
    coverage_wig: Path,
    gc_wig: Path,
    mappability_wig: Path,
    centromere: Path,
    sample_id: str,
    output_dir: Path,
    policy: IchorCnaUlpWgsPolicy,
    normal_panel: Path | None = None,
) -> tuple[str, ...]:
    argv = [
        rscript,
        str(ichorcna_script),
        "--id",
        sample_id,
        "--WIG",
        str(coverage_wig),
        "--ploidy",
        _r_vector(policy.ploidy_candidates),
        "--normal",
        _r_vector(policy.normal_fraction_candidates),
        "--maxCN",
        str(policy.max_copy_number),
        "--gcWig",
        str(gc_wig),
        "--mapWig",
        str(mappability_wig),
        "--centromere",
        str(centromere),
        "--minMapScore",
        str(policy.minimum_mappability),
        "--includeHOMD",
        str(policy.include_homozygous_deletion),
        "--chrs",
        _chromosome_expression(policy, training=False),
        "--chrTrain",
        _chromosome_expression(policy, training=True),
        "--chrNormalize",
        _chromosome_expression(policy, training=True),
        "--genomeBuild",
        policy.genome_build_label,
        "--genomeStyle",
        policy.genome_style,
        "--estimateNormal",
        str(policy.estimate_normal),
        "--estimatePloidy",
        str(policy.estimate_ploidy),
        "--estimateScPrevalence",
        str(policy.estimate_subclone_prevalence),
        "--txnE",
        str(policy.transaction_e),
        "--txnStrength",
        str(policy.transaction_strength),
        "--outDir",
        str(output_dir),
    ]
    if normal_panel is not None:
        argv.extend(["--normalPanel", str(normal_panel)])
    return tuple(argv)


def _fingerprint(path: Path) -> FileFingerprint:
    return FileFingerprint(size_bytes=path.stat().st_size, sha256=sha256_file(path))


def _validate_registered_path(
    *,
    path: Path,
    resource: IchorCnaResourceArtifact,
    label: str,
) -> FileFingerprint:
    if not path.is_file():
        raise ValueError(f"ichorCNA requires {label}: {path}")
    fingerprint = _fingerprint(path)
    if fingerprint.size_bytes != resource.size_bytes or fingerprint.sha256 != resource.sha256:
        raise ValueError(f"ichorCNA {label} SHA-256 or size does not match registered resource")
    return fingerprint


def _validate_resources(
    *,
    resources: IchorCnaResourceBundle,
    policy: IchorCnaUlpWgsPolicy,
    assay_regime: CallerAssayRegime,
    coverage_wig: Path,
    gc_wig: Path,
    mappability_wig: Path,
    centromere: Path,
    normal_panel: Path | None,
) -> tuple[
    FileFingerprint,
    FileFingerprint,
    FileFingerprint,
    FileFingerprint,
    FileFingerprint | None,
]:
    if assay_regime != CallerAssayRegime.CFDNA_ULP_WGS:
        raise ValueError("ichorCNA runtime is restricted to cfDNA ULP-WGS assay regime")
    if resources.assay_regime != assay_regime or policy.assay_regime != assay_regime:
        raise ValueError("ichorCNA assay regime does not match registered resources and policy")
    if resources.genome_build != policy.genome_build:
        raise ValueError("ichorCNA resource genome build does not match policy")
    if resources.bin_size_bp != policy.bin_size_bp:
        raise ValueError("ichorCNA resource bin size does not match policy")

    coverage_fp = _validate_registered_path(
        path=coverage_wig,
        resource=resources.resource_for(IchorCnaResourceRole.COVERAGE_WIG),
        label="coverage WIG",
    )
    gc_fp = _validate_registered_path(
        path=gc_wig,
        resource=resources.resource_for(IchorCnaResourceRole.GC_WIG),
        label="GC WIG",
    )
    map_fp = _validate_registered_path(
        path=mappability_wig,
        resource=resources.resource_for(IchorCnaResourceRole.MAPPABILITY_WIG),
        label="mappability WIG",
    )
    centromere_fp = _validate_registered_path(
        path=centromere,
        resource=resources.resource_for(IchorCnaResourceRole.CENTROMERE),
        label="centromere resource",
    )

    panel_resource = next(
        (item for item in resources.resources if item.role == IchorCnaResourceRole.NORMAL_PANEL),
        None,
    )
    if (normal_panel is None) != (panel_resource is None):
        raise ValueError(
            "ichorCNA normal panel path and registered resource must be declared together"
        )
    panel_fp = None
    if normal_panel is not None and panel_resource is not None:
        panel_fp = _validate_registered_path(
            path=normal_panel,
            resource=panel_resource,
            label="normal panel",
        )
    return coverage_fp, gc_fp, map_fp, centromere_fp, panel_fp


def _parse_optional_number(raw: str, *, label: str) -> float | None:
    value = raw.strip()
    if value.upper() in {"NA", "NAN", "NULL", ""}:
        return None
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"ichorCNA {label} is not numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"ichorCNA {label} must be finite")
    return number


def _required_number(raw: str, *, label: str) -> float:
    value = _parse_optional_number(raw, label=label)
    if value is None:
        raise ValueError(f"ichorCNA {label} is required")
    return value


def _parse_parameters(
    path: Path,
) -> tuple[float, float, float | None, list[IchorCnaSolution]]:
    if not path.is_file():
        raise ValueError("ichorCNA parameters output is missing")
    selected: dict[str, str] = {}
    table_header: list[str] | None = None
    table_rows: list[list[str]] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if table_header is not None:
            table_values = line.split()
            if len(table_values) == len(table_header):
                table_rows.append(table_values)
            continue
        if line.startswith("init "):
            table_header = line.split()
            continue
        key, separator, value = line.partition(":")
        if separator:
            selected[key.strip()] = value.strip()

    tumor_fraction = _required_number(
        selected.get("Tumor Fraction", ""),
        label="Tumor Fraction",
    )
    ploidy = _required_number(selected.get("Ploidy", ""), label="Ploidy")
    coverage = _parse_optional_number(selected.get("Coverage", ""), label="Coverage")
    if not 0 <= tumor_fraction <= 1:
        raise ValueError("ichorCNA Tumor Fraction must be between zero and one")
    if ploidy <= 0:
        raise ValueError("ichorCNA Ploidy must be positive")
    if coverage is not None and coverage < 0:
        raise ValueError("ichorCNA Coverage must be non-negative")

    solutions: list[IchorCnaSolution] = []
    if table_header is not None:
        for values in table_rows:
            solution_row = dict(zip(table_header, values, strict=True))
            normal_fraction = _required_number(
                solution_row["n_est"],
                label="solution n_est",
            )
            solution_ploidy = _required_number(
                solution_row["phi_est"],
                label="solution phi_est",
            )
            if not 0 <= normal_fraction <= 1:
                raise ValueError("ichorCNA solution normal fraction is outside zero to one")
            if solution_ploidy <= 0:
                raise ValueError("ichorCNA solution ploidy must be positive")
            solutions.append(
                IchorCnaSolution(
                    initialization=solution_row["init"],
                    normal_fraction=normal_fraction,
                    ploidy=solution_ploidy,
                    bic=_parse_optional_number(
                        solution_row["BIC"],
                        label="solution BIC",
                    ),
                    fraction_genome_subclonal=_parse_optional_number(
                        solution_row["Frac_genome_subclonal"],
                        label="solution fraction genome subclonal",
                    ),
                    fraction_cna_subclonal=_parse_optional_number(
                        solution_row["Frac_CNA_subclonal"],
                        label="solution fraction CNA subclonal",
                    ),
                    log_likelihood=_parse_optional_number(
                        solution_row["loglik"],
                        label="solution log likelihood",
                    ),
                )
            )
    return tumor_fraction, ploidy, coverage, solutions


def _event_type(call: str) -> EventType | None:
    normalized = call.upper()
    if normalized in {"HOMD", "HETD", "LOSS", "DEL"}:
        return EventType.DELETION
    if normalized in {"GAIN", "AMP", "HLAMP", "DUP"}:
        return EventType.DUPLICATION
    return None


def _parse_bool(raw: str) -> bool | None:
    normalized = raw.strip().upper()
    if normalized in {"TRUE", "1"}:
        return True
    if normalized in {"FALSE", "0"}:
        return False
    if normalized in {"", "NA"}:
        return None
    raise ValueError("ichorCNA subclone status is invalid")


def _parse_segments(path: Path) -> list[IchorCnaSegment]:
    if not path.is_file():
        raise ValueError("ichorCNA segment output is missing")
    result: list[IchorCnaSegment] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chromosome = row.get("chrom")
            if chromosome in {None, ""}:
                chromosome = row.get("chr")
            if chromosome in {None, ""}:
                raise ValueError("ichorCNA segment output is missing chrom/chr column")
            try:
                native_start = int(row["start"])
                native_end = int(row["end"])
            except (KeyError, ValueError) as exc:
                raise ValueError("ichorCNA segment coordinates are invalid") from exc
            if native_start < 1 or native_end < native_start:
                raise ValueError("ichorCNA segment coordinates are invalid")

            copy_raw = row.get("Corrected_Copy_Number")
            if copy_raw in {None, "", "NA"}:
                copy_raw = row.get("copy.number", "")
            copy_number = _required_number(copy_raw or "", label="segment copy number")
            if copy_number < 0:
                raise ValueError("ichorCNA segment copy number must be non-negative")

            call = row.get("Corrected_Call")
            if call in {None, "", "NA"}:
                call = row.get("call")
            if call in {None, ""}:
                raise ValueError("ichorCNA segment call is missing")

            marker_count = None
            raw_markers = row.get("num.mark")
            if raw_markers not in {None, "", "NA"}:
                try:
                    marker_count = int(raw_markers)
                except ValueError as exc:
                    raise ValueError("ichorCNA segment marker count is invalid") from exc
                if marker_count < 0:
                    raise ValueError("ichorCNA segment marker count must be non-negative")

            result.append(
                IchorCnaSegment(
                    chromosome=chromosome,
                    start=native_start - 1,
                    end=native_end,
                    marker_count=marker_count,
                    median_log_r=_parse_optional_number(
                        row.get("seg.median.logR", ""),
                        label="segment median logR",
                    ),
                    copy_number=copy_number,
                    native_call=call,
                    subclone_status=_parse_bool(row.get("subclone.status", "")),
                    log_r_copy_number=_parse_optional_number(
                        row.get("logR_Copy_Number", ""),
                        label="segment logR copy number",
                    ),
                    event_type=_event_type(call),
                )
            )
    if not result:
        raise ValueError("ichorCNA segment output is empty")
    return result


def _artifact_role(relative_path: str) -> str:
    name = PurePosixPath(relative_path).name
    if name.endswith(".params.txt"):
        return "parameters"
    if name.endswith(".cna.seg"):
        return "bin_level_cna"
    if name.endswith(".seg"):
        return "segments"
    if name.endswith(".correctedDepth.txt"):
        return "corrected_depth"
    if name.endswith(".RData"):
        return "r_workspace"
    if name.endswith((".pdf", ".png")):
        return "plot"
    return "caller_native"


def _media_type(relative_path: str) -> str:
    if relative_path.endswith((".txt", ".seg", ".wig")):
        return "text/plain"
    if relative_path.endswith(".pdf"):
        return "application/pdf"
    if relative_path.endswith(".png"):
        return "image/png"
    return "application/octet-stream"


def _native_artifacts(output_dir: Path) -> list[IchorCnaNativeArtifact]:
    artifacts: list[IchorCnaNativeArtifact] = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.is_symlink():
            raise ValueError("ichorCNA native output must not contain symbolic links")
        relative = path.relative_to(output_dir).as_posix()
        artifacts.append(
            IchorCnaNativeArtifact(
                relative_path=relative,
                role=_artifact_role(relative),
                media_type=_media_type(relative),
                fingerprint=_fingerprint(path),
            )
        )
    if not artifacts:
        raise ValueError("ichorCNA returned success without native artifacts")
    return artifacts


def _version(text: str) -> str:
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    return text.splitlines()[0].strip()[:80] if text.strip() else "unknown"


def _probe_version(
    runner: CommandRunner,
    *,
    rscript: str,
    expected_version: str,
) -> str:
    expression = 'cat(as.character(utils::packageVersion("ichorCNA")))'
    result = runner.run([rscript, "-e", expression], timeout_seconds=30)
    if result.returncode != 0:
        raise ValueError(f"ichorCNA version probe failed with exit code {result.returncode}")
    version = _version(f"{result.stdout}\n{result.stderr}")
    if version != expected_version:
        raise ValueError(
            f"ichorCNA version {version!r} does not match policy lock {expected_version!r}"
        )
    return version


def run_ichorcna_ulp_wgs(
    *,
    coverage_wig: Path,
    gc_wig: Path,
    mappability_wig: Path,
    centromere: Path,
    ichorcna_script: Path,
    sample_id: str,
    output_dir: Path,
    resources: IchorCnaResourceBundle,
    policy: IchorCnaUlpWgsPolicy,
    assay_regime: CallerAssayRegime,
    runner: CommandRunner | None = None,
    rscript: str = "Rscript",
    normal_panel: Path | None = None,
) -> IchorCnaReport:
    if assay_regime != CallerAssayRegime.CFDNA_ULP_WGS:
        raise ValueError("ichorCNA runtime is restricted to cfDNA ULP-WGS assay regime")
    if resources.sample_id != sample_id:
        raise ValueError("ichorCNA sample ID does not match registered resources")
    if output_dir.exists():
        raise ValueError("Refusing to overwrite an existing ichorCNA output directory")

    (
        coverage_fp,
        gc_fp,
        map_fp,
        centromere_fp,
        panel_fp,
    ) = _validate_resources(
        resources=resources,
        policy=policy,
        assay_regime=assay_regime,
        coverage_wig=coverage_wig,
        gc_wig=gc_wig,
        mappability_wig=mappability_wig,
        centromere=centromere,
        normal_panel=normal_panel,
    )
    if not ichorcna_script.is_file():
        raise ValueError("ichorCNA runtime script is missing")
    script_fp = _fingerprint(ichorcna_script)
    if script_fp.sha256 != policy.script_sha256:
        raise ValueError("ichorCNA runtime script SHA-256 does not match policy lock")

    command_runner = runner or SubprocessRunner()
    version = _probe_version(
        command_runner,
        rscript=rscript,
        expected_version=policy.expected_version,
    )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    promoted = False
    try:
        argv = build_ichorcna_argv(
            rscript=rscript,
            ichorcna_script=ichorcna_script,
            coverage_wig=coverage_wig,
            gc_wig=gc_wig,
            mappability_wig=mappability_wig,
            centromere=centromere,
            sample_id=sample_id,
            output_dir=staged,
            policy=policy,
            normal_panel=normal_panel,
        )
        result = command_runner.run(argv, timeout_seconds=policy.timeout_seconds)
        if result.returncode != 0:
            diagnostic = result.stderr.strip()[-4000:]
            suffix = f": {diagnostic}" if diagnostic else ""
            raise ValueError(f"ichorCNA failed with exit code {result.returncode}{suffix}")

        parameter_path = staged / f"{sample_id}.params.txt"
        segment_path = staged / f"{sample_id}.seg"
        tumor_fraction, ploidy, coverage, solutions = _parse_parameters(parameter_path)
        segments = _parse_segments(segment_path)
        report = IchorCnaReport(
            sample_id=sample_id,
            assay_regime=assay_regime,
            genome_build=policy.genome_build,
            status=ModuleRunStatus.COMPLETED,
            policy=policy,
            tumor_fraction=tumor_fraction,
            ploidy=ploidy,
            coverage_x=coverage,
            solutions=solutions,
            segments=segments,
            tool=ToolRecord(
                name="ichorCNA",
                version=version,
                parameters={
                    "assay_regime": assay_regime.value,
                    "bin_size_bp": policy.bin_size_bp,
                    "genome_build": policy.genome_build_label,
                    "genome_style": policy.genome_style,
                    "max_copy_number": policy.max_copy_number,
                    "minimum_mappability": policy.minimum_mappability,
                    "estimate_normal": policy.estimate_normal,
                    "estimate_ploidy": policy.estimate_ploidy,
                    "estimate_subclone_prevalence": policy.estimate_subclone_prevalence,
                    "include_homozygous_deletion": policy.include_homozygous_deletion,
                    "transaction_e": policy.transaction_e,
                    "transaction_strength": policy.transaction_strength,
                    "real_tool_qualified": policy.real_tool_qualified,
                    "analytical_validation": policy.analytical_validation,
                },
            ),
            script_fingerprint=script_fp,
            coverage_wig_fingerprint=coverage_fp,
            gc_wig_fingerprint=gc_fp,
            mappability_wig_fingerprint=map_fp,
            centromere_fingerprint=centromere_fp,
            normal_panel_fingerprint=panel_fp,
            native_artifacts=_native_artifacts(staged),
            warnings=[policy.note],
            limitations=[
                (
                    "This ichorCNA adapter is restricted to the registered cfDNA ULP-WGS "
                    "regime; ONT lcWGS and Adaptive Sampling require separate transfer studies."
                ),
                (
                    "Tumor fraction, ploidy and copy-number states are model-derived research "
                    "outputs and do not establish biological truth or clinical validity."
                ),
                "All ichorCNA outputs require human review before downstream interpretation.",
            ],
        )
        os.replace(staged, output_dir)
        promoted = True
        return report
    finally:
        if not promoted:
            shutil.rmtree(staged, ignore_errors=True)

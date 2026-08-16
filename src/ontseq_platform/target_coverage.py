from __future__ import annotations

import gzip
import re
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .execution import CommandRunner, SubprocessRunner
from .models import (
    AlignedBamIntakeReport,
    AssayMode,
    FileFingerprint,
    GenomeBuild,
    InputKind,
    ModuleRunStatus,
    SampleManifest,
    StrictModel,
    ToolRecord,
    Verdict,
)
from .reference import sha256_file

_VERSION = re.compile(r"(?<!\d)(\d+\.\d+\.\d+)(?!\d)")
_CANONICAL_CHROMOSOME = re.compile(r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")


class TargetCoveragePolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    expected_version: str = Field(default="0.3.14", pattern=r"^\d+\.\d+\.\d+$")
    thresholds: list[int] = Field(default_factory=lambda: [1, 10, 20, 30], min_length=1)
    mapq: int = Field(default=0, ge=0, le=60)
    exclude_flags: int = Field(default=1796, ge=0, le=65535)
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def thresholds_are_strictly_increasing(self) -> TargetCoveragePolicy:
        if any(value < 1 for value in self.thresholds):
            raise ValueError("Target coverage thresholds must be positive integers")
        if self.thresholds != sorted(set(self.thresholds)):
            raise ValueError("Target coverage thresholds must be unique and strictly increasing")
        return self


class TargetCoverageRegion(StrictModel):
    chromosome: str = Field(pattern=r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    region_id: str = Field(min_length=1)
    mean_depth: float = Field(ge=0)
    bases_at_threshold: dict[str, int]
    fraction_at_threshold: dict[str, float]

    @model_validator(mode="after")
    def validate_region(self) -> TargetCoverageRegion:
        if self.end <= self.start:
            raise ValueError("Target coverage region end must be greater than start")
        length = self.end - self.start
        if set(self.bases_at_threshold) != set(self.fraction_at_threshold):
            raise ValueError("Target coverage threshold keys are inconsistent")
        for key, bases in self.bases_at_threshold.items():
            if bases < 0 or bases > length:
                raise ValueError(f"Target coverage count {key!r} is outside the region length")
        for fraction in self.fraction_at_threshold.values():
            if not 0 <= fraction <= 1:
                raise ValueError("Target coverage fractions must be between 0 and 1")
        return self


class TargetCoverageReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    target_bed_version: str = Field(min_length=1)
    target_bed_role: Literal["analysis_roi_unbuffered"] = "analysis_roi_unbuffered"
    status: ModuleRunStatus
    policy: TargetCoveragePolicy
    summary_metrics: dict[str, float | int]
    regions: list[TargetCoverageRegion] = Field(min_length=1)
    target_bed_fingerprint: FileFingerprint
    tool: ToolRecord
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def report_is_consistent(self) -> TargetCoverageReport:
        if self.status != ModuleRunStatus.COMPLETED:
            raise ValueError("A normalized target coverage report must have COMPLETED status")
        if self.tool.version != self.policy.expected_version:
            raise ValueError("Mosdepth tool version does not match the target coverage policy lock")
        expected_keys = {f"{value}x" for value in self.policy.thresholds}
        for region in self.regions:
            if set(region.bases_at_threshold) != expected_keys:
                raise ValueError("Target coverage region thresholds do not match the policy")
        if self.summary_metrics.get("region_count") != len(self.regions):
            raise ValueError("Target coverage region_count is inconsistent")
        interval_bases = sum(region.end - region.start for region in self.regions)
        if self.summary_metrics.get("interval_bases") != interval_bases:
            raise ValueError("Target coverage interval_bases is inconsistent")
        return self


@dataclass(frozen=True)
class _BedRegion:
    chromosome: str
    start: int
    end: int
    region_id: str

    @property
    def key(self) -> tuple[str, int, int]:
        return (self.chromosome, self.start, self.end)

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class _RegionDepth:
    region: _BedRegion
    mean_depth: float


def _open_text(path: Path) -> list[str]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return handle.read().splitlines()
    return path.read_text(encoding="utf-8").splitlines()


def _parse_int(raw: str, *, field: str) -> int:
    try:
        return int(raw)
    except ValueError as exc:
        raise ValueError(f"Invalid integer in {field}") from exc


def _parse_float(raw: str, *, field: str) -> float:
    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError(f"Invalid numeric value in {field}") from exc
    if value < 0:
        raise ValueError(f"Negative value in {field}")
    return value


def load_target_bed(path: Path) -> list[_BedRegion]:
    if not path.is_file():
        raise ValueError("Target BED is missing or unreadable")
    regions: list[_BedRegion] = []
    seen: set[tuple[str, int, int]] = set()
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith(("#", "track ", "browser ")):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise ValueError(f"Target BED line {line_number} has fewer than three columns")
        chromosome = fields[0]
        if _CANONICAL_CHROMOSOME.fullmatch(chromosome) is None:
            raise ValueError(f"Target BED line {line_number} uses a non-canonical chromosome")
        start = _parse_int(fields[1], field=f"target BED line {line_number} start")
        end = _parse_int(fields[2], field=f"target BED line {line_number} end")
        if start < 0 or end <= start:
            raise ValueError(f"Target BED line {line_number} has invalid coordinates")
        key = (chromosome, start, end)
        if key in seen:
            raise ValueError("Target BED contains duplicate genomic intervals")
        seen.add(key)
        region_id = (
            fields[3].strip()
            if len(fields) >= 4 and fields[3].strip() not in {"", "."}
            else f"{chromosome}:{start}-{end}"
        )
        regions.append(
            _BedRegion(
                chromosome=chromosome,
                start=start,
                end=end,
                region_id=region_id,
            )
        )
    if not regions:
        raise ValueError("Target BED contains no usable regions")
    return regions


def _region_map(regions: list[_BedRegion]) -> dict[tuple[str, int, int], _BedRegion]:
    return {region.key: region for region in regions}


def parse_mosdepth_regions(path: Path, bed_regions: list[_BedRegion]) -> list[_RegionDepth]:
    if not path.is_file():
        raise ValueError("Mosdepth regions output is missing")
    expected = _region_map(bed_regions)
    parsed: dict[tuple[str, int, int], _RegionDepth] = {}
    for line_number, line in enumerate(_open_text(path), start=1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) not in {4, 5}:
            raise ValueError(f"Mosdepth regions line {line_number} has an unexpected column count")
        chromosome = fields[0]
        start = _parse_int(fields[1], field=f"mosdepth regions line {line_number} start")
        end = _parse_int(fields[2], field=f"mosdepth regions line {line_number} end")
        key = (chromosome, start, end)
        region = expected.get(key)
        if region is None:
            raise ValueError(
                "Mosdepth regions output contains an interval absent from the target BED"
            )
        if key in parsed:
            raise ValueError("Mosdepth regions output contains a duplicate interval")
        if len(fields) == 5:
            output_name = fields[3]
            if output_name not in {region.region_id, "unknown"}:
                raise ValueError("Mosdepth region name does not match the target BED")
            raw_depth = fields[4]
        else:
            raw_depth = fields[3]
        parsed[key] = _RegionDepth(
            region=region,
            mean_depth=_parse_float(raw_depth, field="mosdepth region mean depth"),
        )
    if set(parsed) != set(expected):
        raise ValueError("Mosdepth regions output does not cover every target BED interval")
    return [parsed[region.key] for region in bed_regions]


def parse_mosdepth_thresholds(
    path: Path,
    bed_regions: list[_BedRegion],
    thresholds: list[int],
) -> dict[tuple[str, int, int], dict[str, int]]:
    if not path.is_file():
        raise ValueError("Mosdepth thresholds output is missing")
    lines = _open_text(path)
    header = next((line for line in lines if line.startswith("#")), None)
    if header is None:
        raise ValueError("Mosdepth thresholds output has no header")
    header_fields = header.lstrip("#").split()
    if len(header_fields) < 5:
        raise ValueError("Mosdepth thresholds header is malformed")
    expected_labels = [f"{value}X" for value in thresholds]
    observed_labels = header_fields[-len(expected_labels) :]
    if observed_labels != expected_labels:
        raise ValueError("Mosdepth threshold columns do not match the configured policy")

    expected = _region_map(bed_regions)
    parsed: dict[tuple[str, int, int], dict[str, int]] = {}
    for line_number, line in enumerate(lines, start=1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 4 + len(thresholds):
            raise ValueError(
                f"Mosdepth thresholds line {line_number} has an unexpected column count"
            )
        chromosome = fields[0]
        start = _parse_int(fields[1], field=f"mosdepth thresholds line {line_number} start")
        end = _parse_int(fields[2], field=f"mosdepth thresholds line {line_number} end")
        key = (chromosome, start, end)
        region = expected.get(key)
        if region is None:
            raise ValueError(
                "Mosdepth thresholds output contains an interval absent from the target BED"
            )
        if key in parsed:
            raise ValueError("Mosdepth thresholds output contains a duplicate interval")
        output_name = fields[3]
        if output_name not in {region.region_id, "unknown"}:
            raise ValueError("Mosdepth threshold region name does not match the target BED")
        counts: dict[str, int] = {}
        for threshold, raw_count in zip(thresholds, fields[4:], strict=True):
            count = _parse_int(raw_count, field="mosdepth threshold base count")
            if count < 0 or count > region.length:
                raise ValueError("Mosdepth threshold count is outside the target interval length")
            counts[f"{threshold}x"] = count
        parsed[key] = counts
    if set(parsed) != set(expected):
        raise ValueError("Mosdepth thresholds output does not cover every target BED interval")
    return parsed


def _overlap_count(regions: list[_BedRegion]) -> int:
    count = 0
    by_chromosome: dict[str, list[_BedRegion]] = {}
    for region in regions:
        by_chromosome.setdefault(region.chromosome, []).append(region)
    for chromosome_regions in by_chromosome.values():
        ordered = sorted(chromosome_regions, key=lambda item: (item.start, item.end))
        maximum_end = -1
        for region in ordered:
            if region.start < maximum_end:
                count += 1
            maximum_end = max(maximum_end, region.end)
    return count


def normalize_target_coverage(
    *,
    sample_id: str,
    genome_build: GenomeBuild,
    target_bed: Path,
    target_bed_version: str,
    regions_path: Path,
    thresholds_path: Path,
    policy: TargetCoveragePolicy,
    tool: ToolRecord,
) -> TargetCoverageReport:
    if tool.version != policy.expected_version:
        raise ValueError(
            f"Mosdepth version {tool.version!r} does not match policy lock "
            f"{policy.expected_version!r}"
        )
    bed_regions = load_target_bed(target_bed)
    region_depths = parse_mosdepth_regions(regions_path, bed_regions)
    threshold_counts = parse_mosdepth_thresholds(
        thresholds_path,
        bed_regions,
        policy.thresholds,
    )

    regions: list[TargetCoverageRegion] = []
    for item in region_depths:
        counts = threshold_counts[item.region.key]
        fractions = {key: count / item.region.length for key, count in counts.items()}
        regions.append(
            TargetCoverageRegion(
                chromosome=item.region.chromosome,
                start=item.region.start,
                end=item.region.end,
                region_id=item.region.region_id,
                mean_depth=item.mean_depth,
                bases_at_threshold=counts,
                fraction_at_threshold=fractions,
            )
        )

    interval_bases = sum(region.end - region.start for region in regions)
    weighted_mean = sum(
        region.mean_depth * (region.end - region.start) for region in regions
    ) / interval_bases
    summary_metrics: dict[str, float | int] = {
        "region_count": len(regions),
        "interval_bases": interval_bases,
        "interval_weighted_mean_depth": weighted_mean,
        "minimum_region_mean_depth": min(region.mean_depth for region in regions),
        "median_region_mean_depth": statistics.median(region.mean_depth for region in regions),
        "maximum_region_mean_depth": max(region.mean_depth for region in regions),
    }
    for threshold in policy.thresholds:
        key = f"{threshold}x"
        summary_metrics[f"interval_bases_at_{key}_fraction"] = (
            sum(region.bases_at_threshold[key] for region in regions) / interval_bases
        )

    overlap_count = _overlap_count(bed_regions)
    summary_metrics["overlapping_interval_count"] = overlap_count
    warnings = [policy.note]
    if overlap_count:
        warnings.append(
            "Target BED intervals overlap; interval-weighted summaries count each BED interval "
            "independently and therefore must not be interpreted as unique genomic base coverage."
        )
    return TargetCoverageReport(
        sample_id=sample_id,
        genome_build=genome_build,
        target_bed_version=target_bed_version,
        status=ModuleRunStatus.COMPLETED,
        policy=policy,
        summary_metrics=summary_metrics,
        regions=regions,
        target_bed_fingerprint=FileFingerprint(
            size_bytes=target_bed.stat().st_size,
            sha256=sha256_file(target_bed),
        ),
        tool=tool,
        warnings=warnings,
        limitations=[
            "The target BED is interpreted as an unbuffered analysis ROI BED. The pipeline cannot "
            "infer whether a supplied BED contains Adaptive Sampling selection buffers.",
            "Coverage thresholds are descriptive technical bins and are not validated adequacy or "
            "reportability thresholds.",
            "Off-target enrichment and CNV inference are outside this adapter.",
            "The fourth BED column is retained as region_id when present; no gene semantics are "
            "inferred from region labels.",
            "No read names, source BAM path or per-read data are copied into this report.",
        ],
    )


def _mosdepth_version(text: str) -> str:
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    first_line = text.splitlines()[0].strip() if text.strip() else "unknown"
    return first_line[:80]


def run_target_coverage(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    policy: TargetCoveragePolicy,
    *,
    output_dir: Path,
    runner: CommandRunner | None = None,
    mosdepth: str = "mosdepth",
    threads: int = 4,
) -> TargetCoverageReport:
    if manifest.input.kind != InputKind.ALIGNED_BAM:
        raise ValueError("Target coverage requires input.kind=aligned_bam")
    if manifest.assay.mode != AssayMode.ADAPTIVE_SAMPLING:
        raise ValueError("Target coverage adapter requires assay.mode=adaptive_sampling")
    if manifest.sample_id != intake.sample_id:
        raise ValueError("Manifest and intake artifact must refer to the same sample")
    if manifest.assay.genome_build != intake.genome_build:
        raise ValueError("Manifest and intake artifact use different genome builds")
    if manifest.assay.reference_id != intake.reference_id:
        raise ValueError("Manifest and intake artifact use different reference IDs")
    if intake.verdict == Verdict.FAIL:
        raise ValueError("Target coverage cannot run after a failed aligned-BAM intake gate")
    if not manifest.assay.target_bed or not manifest.assay.target_bed_version:
        raise ValueError("Adaptive Sampling target coverage requires target BED metadata")
    if threads < 1:
        raise ValueError("threads must be at least 1")

    target_bed = Path(manifest.assay.target_bed)
    load_target_bed(target_bed)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / f"{manifest.sample_id}.target-coverage"
    regions_path = Path(f"{prefix}.regions.bed.gz")
    thresholds_path = Path(f"{prefix}.thresholds.bed.gz")
    protected_outputs = [
        regions_path,
        thresholds_path,
        Path(f"{prefix}.mosdepth.summary.txt"),
        Path(f"{prefix}.mosdepth.global.dist.txt"),
        Path(f"{prefix}.mosdepth.region.dist.txt"),
        Path(f"{prefix}.regions.bed.gz.csi"),
        Path(f"{prefix}.thresholds.bed.gz.csi"),
    ]
    if any(path.exists() for path in protected_outputs):
        raise ValueError("Refusing to overwrite existing Mosdepth target coverage outputs")

    command_runner = runner or SubprocessRunner()
    version_result = command_runner.run([mosdepth, "--version"], timeout_seconds=30)
    if version_result.returncode != 0:
        raise ValueError("Mosdepth version probe returned a non-zero exit code")
    version = _mosdepth_version(f"{version_result.stdout}\n{version_result.stderr}")
    if version != policy.expected_version:
        raise ValueError(
            f"Mosdepth version {version!r} does not match policy lock "
            f"{policy.expected_version!r}"
        )

    parameters: dict[str, object] = {
        "threads": threads,
        "no_per_base": True,
        "thresholds": policy.thresholds,
        "mapq": policy.mapq,
        "exclude_flags": policy.exclude_flags,
        "target_bed_role": "analysis_roi_unbuffered",
        "expected_version": policy.expected_version,
    }
    argv = [
        mosdepth,
        "--threads",
        str(threads),
        "--no-per-base",
        "--by",
        str(target_bed),
        "--thresholds",
        ",".join(str(value) for value in policy.thresholds),
        "--mapq",
        str(policy.mapq),
        "--flag",
        str(policy.exclude_flags),
        str(prefix),
        manifest.input.path,
    ]
    result = command_runner.run(argv, timeout_seconds=7200)
    if result.returncode != 0:
        raise ValueError(f"Mosdepth failed with exit code {result.returncode}")
    return normalize_target_coverage(
        sample_id=manifest.sample_id,
        genome_build=manifest.assay.genome_build,
        target_bed=target_bed,
        target_bed_version=manifest.assay.target_bed_version,
        regions_path=regions_path,
        thresholds_path=thresholds_path,
        policy=policy,
        tool=ToolRecord(name="mosdepth", version=version, parameters=parameters),
    )
from __future__ import annotations

import csv
import json
import os
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from ..execution import StreamingCommandRunner
from ..models import (
    EventType,
    Evidence,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    ReferenceLock,
    StrictModel,
    ToolRecord,
)


class QDNAseqPolicy(StrictModel):
    """Executable policy for the live QDNAseq + ACE CNV lane."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    bin_sizes_kbp: list[int] = Field(default_factory=lambda: [100, 500, 1000], min_length=1)
    primary_bin_size_kbp: int = 500
    ace_penalty: float = Field(default=0.6, ge=0, le=1)
    ploidy_min: float = Field(default=1.5, gt=0)
    ploidy_max: float = Field(default=4.5, gt=0)
    ploidy_step: float = Field(default=0.05, gt=0)
    minimum_segment_bins: int = Field(default=1, ge=1)
    whole_chromosome_fraction: float = Field(default=0.90, gt=0, le=1)
    timeout_seconds: int = Field(default=7200, ge=60)
    expected_qdnaseq_version: str | None = None
    expected_ace_version: str | None = None
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def validate_policy(self) -> QDNAseqPolicy:
        if len(self.bin_sizes_kbp) != len(set(self.bin_sizes_kbp)):
            raise ValueError("QDNAseq bin sizes must be unique")
        if any(item <= 0 for item in self.bin_sizes_kbp):
            raise ValueError("QDNAseq bin sizes must be positive")
        if self.primary_bin_size_kbp not in self.bin_sizes_kbp:
            raise ValueError("primary_bin_size_kbp must be one of bin_sizes_kbp")
        if self.ploidy_max <= self.ploidy_min:
            raise ValueError("ploidy_max must be greater than ploidy_min")
        return self


class CnvFit(StrictModel):
    bin_size_kbp: int = Field(gt=0)
    cellularity: float = Field(ge=0, le=1)
    ploidy: float = Field(gt=0)
    fit_error: float = Field(ge=0)
    candidate_count: int = Field(ge=1)
    segment_count: int = Field(ge=0)
    alternatives: list[dict[str, float]] = Field(default_factory=list)
    segment_file: str
    chromosome_file: str
    fit_plot: str
    copy_number_plot: str
    rds_file: str


class CnvChromosomeConsensus(StrictModel):
    chromosome: str = Field(pattern=r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    median_copy_number: float = Field(ge=0)
    rounded_copy_number: int = Field(ge=0)
    agreeing_bins: int = Field(ge=0)
    contributing_bins: int = Field(ge=1)
    min_copy_number: float = Field(ge=0)
    max_copy_number: float = Field(ge=0)


class QDNAseqCallReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    status: ModuleRunStatus
    primary_fit: CnvFit
    fits: list[CnvFit]
    chromosome_consensus: list[CnvChromosomeConsensus]
    events: list[GenomicEvent]
    tools: list[ToolRecord]
    output_files: list[str]
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def status_matches_events(self) -> QDNAseqCallReport:
        expected = ModuleRunStatus.COMPLETED if self.events else ModuleRunStatus.NO_CALL
        if self.status != expected:
            raise ValueError("CNV status is inconsistent with normalized events")
        if self.primary_fit not in self.fits:
            raise ValueError("primary CNV fit must also be present in fits")
        return self


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"QDNAseq output is missing or empty: {path.name}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"QDNAseq output has no header: {path.name}")
        return [dict(row) for row in reader]


def _float(row: Mapping[str, str], key: str) -> float:
    value = row.get(key, "")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid numeric value for {key}: {value!r}") from exc
    if parsed != parsed:
        raise ValueError(f"invalid NaN value for {key}")
    return parsed


def _int(row: Mapping[str, str], key: str) -> int:
    value = row.get(key, "")
    try:
        return int(float(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid integer value for {key}: {value!r}") from exc


def _as_float(value: object, key: str) -> float:
    if not isinstance(value, (int, float, str)):
        raise ValueError(f"invalid numeric value for {key}: {value!r}")
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"invalid numeric value for {key}: {value!r}") from exc
    if parsed != parsed:
        raise ValueError(f"invalid NaN value for {key}")
    return parsed


def _as_int(value: object, key: str) -> int:
    parsed = _as_float(value, key)
    if not parsed.is_integer():
        raise ValueError(f"invalid integer value for {key}: {value!r}")
    return int(parsed)


def _as_text(value: object, key: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"invalid text value for {key}: {value!r}")
    return value


def _contig_lengths(reference_lock: ReferenceLock) -> dict[str, int]:
    result: dict[str, int] = {}
    canonical = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    for item in reference_lock.contigs:
        name = item.name if item.name.startswith("chr") else f"chr{item.name}"
        if name in canonical:
            result[name] = item.length
    return result


def _tool_records(summary: Mapping[str, object], policy: QDNAseqPolicy) -> list[ToolRecord]:
    raw_versions = summary.get("package_versions")
    if not isinstance(raw_versions, dict):
        raise ValueError("QDNAseq summary is missing package_versions")
    versions = {str(key): str(value) for key, value in raw_versions.items() if value is not None}
    qdna = versions.get("QDNAseq", "UNKNOWN")
    ace = versions.get("ACE", "UNKNOWN")
    if policy.expected_qdnaseq_version and qdna != policy.expected_qdnaseq_version:
        raise ValueError(
            f"QDNAseq version mismatch: expected {policy.expected_qdnaseq_version}, observed {qdna}"
        )
    if policy.expected_ace_version and ace != policy.expected_ace_version:
        expected = policy.expected_ace_version
        raise ValueError(f"ACE version mismatch: expected {expected}, observed {ace}")
    shared = {
        "profile": policy.profile_id,
        "bin_sizes_kbp": policy.bin_sizes_kbp,
        "primary_bin_size_kbp": policy.primary_bin_size_kbp,
        "ace_penalty": policy.ace_penalty,
        "ploidy_min": policy.ploidy_min,
        "ploidy_max": policy.ploidy_max,
        "ploidy_step": policy.ploidy_step,
    }
    return [
        ToolRecord(name="QDNAseq", version=qdna, parameters=shared),
        ToolRecord(name="ACE", version=ace, parameters=shared),
        ToolRecord(name="R", version=versions.get("R", "UNKNOWN"), parameters={}),
    ]


def _parse_fit(raw: Mapping[str, object]) -> CnvFit:
    alternatives_raw = raw.get("alternatives", [])
    alternatives: list[dict[str, float]] = []
    if isinstance(alternatives_raw, list):
        for item in alternatives_raw:
            if not isinstance(item, dict):
                continue
            alternatives.append(
                {
                    "cellularity": _as_float(item.get("cellularity"), "alternative.cellularity"),
                    "ploidy": _as_float(item.get("ploidy"), "alternative.ploidy"),
                    "fit_error": _as_float(item.get("fit_error"), "alternative.fit_error"),
                }
            )
    return CnvFit(
        bin_size_kbp=_as_int(raw.get("bin_size_kbp"), "bin_size_kbp"),
        cellularity=_as_float(raw.get("cellularity"), "cellularity"),
        ploidy=_as_float(raw.get("ploidy"), "ploidy"),
        fit_error=_as_float(raw.get("fit_error"), "fit_error"),
        candidate_count=_as_int(raw.get("candidate_count"), "candidate_count"),
        segment_count=_as_int(raw.get("segment_count"), "segment_count"),
        alternatives=alternatives,
        segment_file=_as_text(raw.get("segment_file"), "segment_file"),
        chromosome_file=_as_text(raw.get("chromosome_file"), "chromosome_file"),
        fit_plot=_as_text(raw.get("fit_plot"), "fit_plot"),
        copy_number_plot=_as_text(raw.get("copy_number_plot"), "copy_number_plot"),
        rds_file=_as_text(raw.get("rds_file"), "rds_file"),
    )


def _parse_consensus(path: Path) -> list[CnvChromosomeConsensus]:
    rows = _read_tsv(path)
    result: list[CnvChromosomeConsensus] = []
    for row in rows:
        result.append(
            CnvChromosomeConsensus(
                chromosome=str(row["chromosome"]),
                median_copy_number=_float(row, "median_copy_number"),
                rounded_copy_number=_int(row, "rounded_copy_number"),
                agreeing_bins=_int(row, "agreeing_bins"),
                contributing_bins=_int(row, "contributing_bins"),
                min_copy_number=_float(row, "min_copy_number"),
                max_copy_number=_float(row, "max_copy_number"),
            )
        )
    return result


def _events_from_primary_segments(
    path: Path,
    *,
    sample_id: str,
    fit: CnvFit,
    tools: Sequence[ToolRecord],
    reference_lock: ReferenceLock,
    minimum_segment_bins: int,
    whole_chromosome_fraction: float,
    consensus: Mapping[str, CnvChromosomeConsensus],
) -> list[GenomicEvent]:
    rows = _read_tsv(path)
    lengths = _contig_lengths(reference_lock)
    events: list[GenomicEvent] = []
    caller_version = "/".join(item.version for item in tools if item.name in {"QDNAseq", "ACE"})
    baseline = int(round(fit.ploidy))
    serial = 0
    for row in rows:
        bins = _int(row, "bin_count")
        if bins < minimum_segment_bins:
            continue
        call = _float(row, "call")
        if abs(call) < 1e-12:
            continue
        chromosome = str(row["chromosome"])
        start = _int(row, "start")
        end = _int(row, "end")
        if end <= start:
            raise ValueError(f"QDNAseq segment end is not after start on {chromosome}")
        copy_number = max(0.0, _float(row, "absolute_copy_number"))
        contig_length = lengths.get(chromosome)
        fraction = (end - start) / contig_length if contig_length else 0.0
        direction_gain = copy_number > baseline
        if fraction >= whole_chromosome_fraction:
            event_type = EventType.CHROMOSOME_GAIN if direction_gain else EventType.CHROMOSOME_LOSS
        else:
            event_type = EventType.DUPLICATION if direction_gain else EventType.DELETION
        serial += 1
        agreement = consensus.get(chromosome)
        notes = [
            f"QDNAseq primary bin {fit.bin_size_kbp} kbp",
            (
                f"ACE cellularity={fit.cellularity:.3f}, ploidy={fit.ploidy:.3f}, "
                f"fit_error={fit.fit_error:.6g}"
            ),
            f"ACE call={call:.3g}; segment supported by {bins} QDNAseq bin(s)",
        ]
        if agreement is not None:
            notes.append(
                "Chromosome-level multi-bin agreement "
                f"{agreement.agreeing_bins}/{agreement.contributing_bins}; "
                f"median CN={agreement.median_copy_number:.3f}"
            )
        quality = None
        if row.get("qnorm_log10") not in {None, ""}:
            quality = abs(_float(row, "qnorm_log10"))
        events.append(
            GenomicEvent(
                event_id=f"CNV_{sample_id}_{serial:04d}",
                event_type=event_type,
                primary=Locus(chromosome=chromosome, start=start, end=end),
                length_bp=end - start,
                copy_number=copy_number,
                evidence=[
                    Evidence(
                        caller="QDNAseq+ACE",
                        caller_version=caller_version,
                        quality=quality,
                    )
                ],
                confidence="unclassified",
                reportable=False,
                notes=notes,
            )
        )
    return events


def _copy_outputs(staged: Path, final: Path) -> None:
    if final.exists():
        shutil.rmtree(final)
    os.replace(staged, final)


def run_qdnaseq_ace(
    *,
    bam: Path,
    sample_id: str,
    genome_build: GenomeBuild,
    reference_lock: ReferenceLock,
    policy: QDNAseqPolicy,
    output_dir: Path,
    script: Path,
    runner: StreamingCommandRunner,
    rscript: str = "Rscript",
    threads: int = 4,
) -> QDNAseqCallReport:
    """Run QDNAseq at multiple resolutions, fit ACE purity/ploidy and normalize events.

    The R process writes into a staging directory. Only a fully parseable result directory
    is promoted into ``output_dir`` so an interrupted R process cannot leave a plausible
    final CNV result behind.
    """
    if not bam.is_file():
        raise ValueError(f"CNV BAM does not exist: {bam}")
    if not script.is_file():
        raise ValueError(f"QDNAseq runner script does not exist: {script}")
    if reference_lock.genome_build != genome_build:
        raise ValueError("CNV reference lock and manifest use different genome builds")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    promoted = False
    try:
        argv = [
            rscript,
            str(script),
            "--bam",
            str(bam),
            "--output-dir",
            str(staged),
            "--sample-id",
            sample_id,
            "--genome-build",
            genome_build.value,
            "--bin-sizes-kbp",
            ",".join(str(item) for item in policy.bin_sizes_kbp),
            "--primary-bin-kbp",
            str(policy.primary_bin_size_kbp),
            "--ace-penalty",
            str(policy.ace_penalty),
            "--ploidy-min",
            str(policy.ploidy_min),
            "--ploidy-max",
            str(policy.ploidy_max),
            "--ploidy-step",
            str(policy.ploidy_step),
            "--threads",
            str(max(1, threads)),
        ]
        completed = runner.run(argv, timeout_seconds=policy.timeout_seconds)
        if completed.returncode != 0:
            tail = completed.stderr.strip()[-4000:]
            raise ValueError(
                f"QDNAseq/ACE exited with code {completed.returncode}"
                + (f": {tail}" if tail else "")
            )

        summary_path = staged / f"{sample_id}.qdnaseq-ace.summary.json"
        if not summary_path.is_file():
            raise ValueError("QDNAseq/ACE completed without a summary JSON")
        raw_summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if not isinstance(raw_summary, dict):
            raise ValueError("QDNAseq summary must be a JSON object")
        summary: Mapping[str, object] = raw_summary
        if summary.get("sample_id") != sample_id:
            raise ValueError("QDNAseq summary sample ID does not match the manifest")
        if summary.get("genome_build") != genome_build.value:
            raise ValueError("QDNAseq summary genome build does not match the manifest")
        observed_primary = _as_int(summary.get("primary_bin_size_kbp"), "primary_bin_size_kbp")
        if observed_primary != policy.primary_bin_size_kbp:
            raise ValueError("QDNAseq summary primary bin does not match policy")

        runs_raw = summary.get("runs")
        if not isinstance(runs_raw, list) or not runs_raw:
            raise ValueError("QDNAseq summary contains no bin-size runs")
        fits = [_parse_fit(item) for item in runs_raw if isinstance(item, dict)]
        if {item.bin_size_kbp for item in fits} != set(policy.bin_sizes_kbp):
            raise ValueError("QDNAseq summary does not contain exactly the configured bin sizes")
        primary = next(item for item in fits if item.bin_size_kbp == policy.primary_bin_size_kbp)

        consensus_name = _as_text(summary.get("consensus_file"), "consensus_file")
        if Path(consensus_name).name != consensus_name:
            raise ValueError("QDNAseq summary contains an invalid consensus filename")
        consensus = _parse_consensus(staged / consensus_name)
        consensus_by_chr = {item.chromosome: item for item in consensus}
        tools = _tool_records(summary, policy)

        events = _events_from_primary_segments(
            staged / primary.segment_file,
            sample_id=sample_id,
            fit=primary,
            tools=tools,
            reference_lock=reference_lock,
            minimum_segment_bins=policy.minimum_segment_bins,
            whole_chromosome_fraction=policy.whole_chromosome_fraction,
            consensus=consensus_by_chr,
        )

        expected_names = {summary_path.name, consensus_name}
        for fit in fits:
            expected_names.update(
                {
                    fit.segment_file,
                    fit.chromosome_file,
                    fit.fit_plot,
                    fit.copy_number_plot,
                    fit.rds_file,
                }
            )
        missing = sorted(name for name in expected_names if not (staged / name).is_file())
        if missing:
            raise ValueError("QDNAseq result is incomplete; missing: " + ", ".join(missing))

        warnings: list[str] = []
        for chromosome in consensus:
            if chromosome.agreeing_bins < chromosome.contributing_bins:
                disagreements = chromosome.contributing_bins - chromosome.agreeing_bins
                warnings.append(
                    f"{chromosome.chromosome}: rounded chromosome CN disagrees across "
                    f"{disagreements}/{chromosome.contributing_bins} bin sizes"
                )
        cellularities = [item.cellularity for item in fits]
        ploidies = [item.ploidy for item in fits]
        if max(cellularities) - min(cellularities) >= 0.20:
            warnings.append("ACE cellularity estimates differ by >=0.20 across bin sizes")
        if max(ploidies) - min(ploidies) >= 0.50:
            warnings.append("ACE ploidy estimates differ by >=0.50 across bin sizes")

        _copy_outputs(staged, output_dir)
        promoted = True
        return QDNAseqCallReport(
            sample_id=sample_id,
            genome_build=genome_build,
            status=ModuleRunStatus.COMPLETED if events else ModuleRunStatus.NO_CALL,
            primary_fit=primary,
            fits=fits,
            chromosome_consensus=consensus,
            events=events,
            tools=tools,
            output_files=sorted(expected_names),
            warnings=warnings,
            limitations=[
                (
                    "Segment-level event normalization uses the configured primary QDNAseq "
                    "bin size; chromosome-level agreement is retained across all configured "
                    "bin sizes."
                ),
                (
                    "ACE purity/ploidy is selected deterministically from squaremodel minima; "
                    "alternative low-error fits are retained in the CNV JSON for inspection."
                ),
            ],
        )
    finally:
        if not promoted and staged.exists():
            shutil.rmtree(staged, ignore_errors=True)

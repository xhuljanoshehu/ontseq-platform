from __future__ import annotations

import gzip
import os
import re
import shutil
import tempfile
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

from pydantic import ValidationError

from .execution import CommandRunner, SubprocessRunner
from .models import (
    AlignedBamIntakeReport,
    CuteSvCallReport,
    CuteSvPolicy,
    EventType,
    Evidence,
    FileFingerprint,
    GenomeBuild,
    GenomicEvent,
    InputKind,
    Locus,
    ModuleRunStatus,
    SampleManifest,
    ToolRecord,
    Verdict,
)
from .reference import sha256_file

_BND_MATE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
_VERSION = re.compile(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)")
_SV_TYPES = {
    "DEL": EventType.DELETION,
    "DUP": EventType.DUPLICATION,
    "INV": EventType.INVERSION,
    "INS": EventType.INSERTION,
    "BND": EventType.TRANSLOCATION,
    "TRA": EventType.TRANSLOCATION,
}


class _RejectedRecord(ValueError):
    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def _open_vcf(path: Path) -> Iterator[str]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        yield from handle


def _info(raw: str) -> dict[str, str | bool]:
    result: dict[str, str | bool] = {}
    if raw == ".":
        return result
    for item in raw.split(";"):
        key, separator, value = item.partition("=")
        if not key:
            raise _RejectedRecord("malformed_info")
        result[key] = value if separator else True
    return result


def _first(value: str | bool | None) -> str | None:
    if not isinstance(value, str):
        return None
    first = value.split(",", maxsplit=1)[0]
    return None if first in {"", "."} else first


def _integer(value: str | bool | None, reason: str) -> int:
    raw = _first(value)
    if raw is None:
        raise _RejectedRecord(reason)
    try:
        return int(raw)
    except ValueError as exc:
        raise _RejectedRecord(reason) from exc


def _quality(raw: str) -> float | None:
    if raw in {"", "."}:
        return None
    try:
        value = float(raw)
    except ValueError as exc:
        raise _RejectedRecord("malformed_quality") from exc
    if value < 0:
        raise _RejectedRecord("malformed_quality")
    return value


def _sample_values(fields: list[str]) -> dict[str, str]:
    if len(fields) < 10 or fields[8] in {"", "."} or fields[9] in {"", "."}:
        return {}
    return dict(zip(fields[8].split(":"), fields[9].split(":"), strict=False))


def _support(info: dict[str, str | bool], sample: dict[str, str]) -> int:
    for key in ("RE", "SUPPORT"):
        if key in info:
            return _integer(info[key], "missing_or_invalid_support")
    for key in ("DV", "RE"):
        try:
            return int(sample[key])
        except (KeyError, ValueError):
            continue
    raise _RejectedRecord("missing_or_invalid_support")


def _vaf(sample: dict[str, str], support: int) -> float | None:
    try:
        reference = int(sample["DR"])
    except (KeyError, ValueError):
        return None
    total = reference + support
    return support / total if total > 0 else None


def _loci(
    chromosome: str,
    position: int,
    alternate: str,
    event_type: EventType,
    info: dict[str, str | bool],
) -> tuple[Locus, Locus | None, int | None]:
    try:
        if event_type == EventType.TRANSLOCATION:
            mate = _BND_MATE.search(alternate)
            chromosome2 = mate.group(1) if mate else (_first(info.get("CHR2")) or "")
            position2 = int(mate.group(2)) if mate else _integer(info.get("END"), "missing_mate")
            return (
                Locus(chromosome=chromosome, start=position - 1, end=position),
                Locus(chromosome=chromosome2, start=position2 - 1, end=position2),
                None,
            )
        raw_length = _first(info.get("SVLEN"))
        if raw_length is not None:
            length = abs(int(raw_length))
        elif event_type == EventType.INSERTION:
            raise _RejectedRecord("missing_svlen")
        else:
            length = abs(_integer(info.get("END"), "missing_end") - position + 1)
        if length < 1:
            raise _RejectedRecord("malformed_svlen")
        start = position - 1
        end = start + 1 if event_type == EventType.INSERTION else start + length
        return Locus(chromosome=chromosome, start=start, end=end), None, length
    except (ValidationError, ValueError) as exc:
        if isinstance(exc, _RejectedRecord):
            raise
        raise _RejectedRecord("invalid_locus") from exc


def _normalize_record(
    fields: list[str], *, record_number: int, policy: CuteSvPolicy, caller_version: str
) -> GenomicEvent:
    if len(fields) < 8:
        raise _RejectedRecord("malformed_record")
    chromosome, raw_position, _record_id, _ref, alt, raw_quality, raw_filter, raw_info = fields[:8]
    try:
        position = int(raw_position)
    except ValueError as exc:
        raise _RejectedRecord("malformed_position") from exc
    if position < 1:
        raise _RejectedRecord("malformed_position")
    filters = raw_filter.split(";") if raw_filter not in {"", "."} else []
    if policy.pass_only and filters != ["PASS"]:
        raise _RejectedRecord("filter_not_pass")
    info = _info(raw_info)
    raw_type = _first(info.get("SVTYPE"))
    event_type = _SV_TYPES.get((raw_type or "").upper().split(":", maxsplit=1)[0])
    if event_type is None:
        raise _RejectedRecord("unsupported_svtype")
    if event_type not in policy.allowed_sv_types:
        raise _RejectedRecord("svtype_not_allowed")
    sample = _sample_values(fields)
    support = _support(info, sample)
    if support < policy.min_support:
        raise _RejectedRecord("support_below_policy")
    quality = _quality(raw_quality)
    if policy.minimum_quality is not None and (quality is None or quality < policy.minimum_quality):
        raise _RejectedRecord("quality_below_policy")
    primary, secondary, length = _loci(chromosome, position, alt, event_type, info)
    if length is not None and length < policy.min_sv_length:
        raise _RejectedRecord("sv_length_below_policy")
    strand = _first(info.get("STRAND"))
    if strand not in {None, "+", "-", "+-"}:
        strand = None
    precise = True if "PRECISE" in info else False if "IMPRECISE" in info else None
    return GenomicEvent(
        event_id=f"CUTESV-{record_number:06d}",
        event_type=event_type,
        primary=primary,
        secondary=secondary,
        length_bp=length,
        evidence=[
            Evidence(
                caller="cuteSV",
                caller_version=caller_version,
                support_reads=support,
                variant_allele_fraction=_vaf(sample, support),
                quality=quality,
                filters=filters,
                supporting_read_strands=strand,
                precise=precise,
            )
        ],
        reportable=False,
        notes=[
            "Normalized from a cuteSV VCF into 0-based, half-open coordinates.",
            "Candidate only; assay-specific analytical validation and expert review are pending.",
        ],
    )


def normalize_cutesv_vcf(
    path: Path,
    *,
    sample_id: str,
    genome_build: GenomeBuild,
    policy: CuteSvPolicy,
    tool: ToolRecord,
) -> CuteSvCallReport:
    if not path.is_file():
        raise ValueError("cuteSV VCF is missing or unreadable")
    if tool.version != policy.expected_version:
        raise ValueError(
            f"cuteSV version {tool.version!r} does not match policy lock "
            f"{policy.expected_version!r}"
        )
    raw_count = 0
    events: list[GenomicEvent] = []
    rejections: Counter[str] = Counter()
    saw_fileformat = False
    saw_columns = False
    sample_column_count: int | None = None
    for raw_line in _open_vcf(path):
        line = raw_line.rstrip("\r\n")
        if line.startswith("##fileformat=VCF"):
            saw_fileformat = True
            continue
        if line.startswith("#CHROM"):
            saw_columns = True
            sample_column_count = max(0, len(line.split("\t")) - 9)
            continue
        if not line or line.startswith("#"):
            continue
        raw_count += 1
        try:
            events.append(
                _normalize_record(
                    line.split("\t"),
                    record_number=raw_count,
                    policy=policy,
                    caller_version=tool.version,
                )
            )
        except _RejectedRecord as exc:
            rejections[exc.reason] += 1
    if not saw_fileformat or not saw_columns:
        raise ValueError("cuteSV output is not a complete VCF document")
    if sample_column_count not in {0, 1}:
        raise ValueError("Single-sample analysis requires at most one VCF sample column")
    warnings = [policy.note] if policy.status != "validated" else []
    if rejections:
        warnings.append("One or more cuteSV records were rejected by normalization policy.")
    if not events:
        warnings.append(
            "NO_CALL means no cuteSV record passed this technical policy; it is not a "
            "clinical negative result."
        )
    return CuteSvCallReport(
        sample_id=sample_id,
        genome_build=genome_build,
        status=ModuleRunStatus.COMPLETED if events else ModuleRunStatus.NO_CALL,
        policy=policy,
        events=events,
        raw_record_count=raw_count,
        accepted_record_count=len(events),
        rejected_record_count=sum(rejections.values()),
        rejection_counts=dict(sorted(rejections.items())),
        tool=tool,
        vcf_fingerprint=FileFingerprint(size_bytes=path.stat().st_size, sha256=sha256_file(path)),
        warnings=warnings,
        limitations=[
            "Only canonical chromosomes supported by the normalized locus contract are retained.",
            "cuteSV candidates are not validated clinical findings or fusion assertions.",
        ],
    )


def cutesv_version(text: str) -> str:
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    return text.splitlines()[0].strip()[:80] if text.strip() else "unknown"


def run_cutesv(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    policy: CuteSvPolicy,
    *,
    reference_fasta: Path,
    output_vcf: Path,
    runner: CommandRunner | None = None,
    cutesv: str = "cuteSV",
    threads: int = 4,
) -> CuteSvCallReport:
    if manifest.input.kind != InputKind.ALIGNED_BAM:
        raise ValueError("cuteSV requires input.kind=aligned_bam")
    if manifest.sample_id != intake.sample_id or manifest.assay.genome_build != intake.genome_build:
        raise ValueError("Manifest and intake artifact do not describe the same sample/build")
    if manifest.assay.reference_id != intake.reference_id:
        raise ValueError("Manifest and intake artifact use different reference IDs")
    if intake.verdict == Verdict.FAIL:
        raise ValueError("cuteSV cannot run after a failed aligned-BAM intake gate")
    if not reference_fasta.is_file():
        raise ValueError("cuteSV requires the locked reference FASTA")
    if output_vcf.exists():
        raise ValueError("Refusing to overwrite an existing cuteSV VCF")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    command_runner = runner or SubprocessRunner()
    probe = command_runner.run([cutesv, "--version"], timeout_seconds=30)
    if probe.returncode != 0:
        raise ValueError("cuteSV version probe returned a non-zero exit code")
    version = cutesv_version(f"{probe.stdout}\n{probe.stderr}")
    if version != policy.expected_version:
        raise ValueError(
            f"cuteSV version {version!r} does not match policy lock {policy.expected_version!r}"
        )
    staged_dir = Path(tempfile.mkdtemp(prefix=".cutesv-", dir=output_vcf.parent))
    staged_vcf = staged_dir / output_vcf.name
    (staged_dir / "work").mkdir(parents=True, exist_ok=True)
    parameters: dict[str, str | int | float | bool] = {
        "threads": threads,
        "min_support": policy.min_support,
        "min_size": policy.min_sv_length,
        "max_cluster_bias_INS": policy.max_cluster_bias_ins,
        "diff_ratio_merging_INS": policy.diff_ratio_merging_ins,
        "max_cluster_bias_DEL": policy.max_cluster_bias_del,
        "diff_ratio_merging_DEL": policy.diff_ratio_merging_del,
        "expected_version": policy.expected_version,
        "normalizer_pass_only": True,
    }
    argv = [
        cutesv,
        manifest.input.path,
        str(reference_fasta),
        str(staged_dir / "work"),
        str(staged_vcf),
        "--threads",
        str(threads),
        "--min_support",
        str(policy.min_support),
        "--min_size",
        str(policy.min_sv_length),
        "--max_cluster_bias_INS",
        str(policy.max_cluster_bias_ins),
        "--diff_ratio_merging_INS",
        str(policy.diff_ratio_merging_ins),
        "--max_cluster_bias_DEL",
        str(policy.max_cluster_bias_del),
        "--diff_ratio_merging_DEL",
        str(policy.diff_ratio_merging_del),
    ]
    try:
        result = command_runner.run(argv, timeout_seconds=7200)
        if result.returncode != 0:
            detail = result.stderr.strip()[-2000:]
            suffix = f": {detail}" if detail else ""
            raise ValueError(f"cuteSV failed with exit code {result.returncode}{suffix}")
        if not staged_vcf.is_file():
            raise ValueError("cuteSV returned success but produced no VCF")
        report = normalize_cutesv_vcf(
            staged_vcf,
            sample_id=manifest.sample_id,
            genome_build=manifest.assay.genome_build,
            policy=policy,
            tool=ToolRecord(name="cuteSV", version=version, parameters=parameters),
        )
        os.replace(staged_vcf, output_vcf)
        return report
    finally:
        shutil.rmtree(staged_dir, ignore_errors=True)

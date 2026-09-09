from __future__ import annotations

import gzip
import os
import re
import tempfile
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

from pydantic import ValidationError

from .breakends import BreakendParseError, resolve_breakend
from .execution import CommandRunner, SubprocessRunner
from .models import (
    AlignedBamIntakeReport,
    EventType,
    Evidence,
    FileFingerprint,
    GenomeBuild,
    GenomicEvent,
    InputKind,
    Locus,
    ModuleRunStatus,
    SampleManifest,
    SnifflesCallReport,
    SnifflesMode,
    SnifflesPolicy,
    ToolRecord,
    Verdict,
)
from .reference import sha256_file

_VERSION = re.compile(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)")
_SV_TYPE_MAP = {
    "DEL": EventType.DELETION,
    "DUP": EventType.DUPLICATION,
    "INV": EventType.INVERSION,
    "INS": EventType.INSERTION,
    "BND": EventType.TRANSLOCATION,
    "TRA": EventType.TRANSLOCATION,
}
_KNOWN_SNIFFLES_FILTERS = {
    "ALN_NM",
    "COV_CHANGE_DEL",
    "COV_CHANGE_DUP",
    "COV_CHANGE_FRAC_CE",
    "COV_CHANGE_FRAC_ED",
    "COV_CHANGE_FRAC_SC",
    "COV_CHANGE_FRAC_US",
    "COV_CHANGE_INS",
    "COV_MIN",
    "COV_MIN_GT",
    "COV_VAR",
    "GT",
    "GT_FAILED",
    "INLINE_SA",
    "MOSAIC_SV_CLOSE_EDGE",
    "MOSAIC_VAF",
    "NOT_MOSAIC_VAF",
    "SINGLE_BREAK",
    "STDEV_LEN",
    "STDEV_POS",
    "STRAND",
    "STRAND_BND",
    "STRAND_MOSAIC",
    "SUPPORT_MIN",
    "SVLEN_MAX_MOSAIC",
    "SVLEN_MIN",
    "SVLEN_MIN_MOSAIC",
}


class _RejectedRecord(ValueError):
    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def _open_vcf(path: Path) -> Iterator[str]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield from handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield from handle


def _parse_info(raw: str) -> dict[str, str | bool]:
    parsed: dict[str, str | bool] = {}
    if raw == ".":
        return parsed
    for item in raw.split(";"):
        key, separator, value = item.partition("=")
        if not key:
            raise _RejectedRecord("malformed_info")
        parsed[key] = value if separator else True
    return parsed


def _first_value(value: str | bool | None) -> str | None:
    if not isinstance(value, str):
        return None
    first = value.split(",", maxsplit=1)[0]
    return None if first in {"", "."} else first


def _integer(value: str | bool | None, *, reason: str) -> int:
    first = _first_value(value)
    if first is None:
        raise _RejectedRecord(reason)
    try:
        return int(first)
    except ValueError as exc:
        raise _RejectedRecord(reason) from exc


def _optional_float(value: str | bool | None, *, reason: str) -> float | None:
    first = _first_value(value)
    if first is None:
        return None
    try:
        return float(first)
    except ValueError as exc:
        raise _RejectedRecord(reason) from exc


def _optional_nonnegative_float(value: str | bool | None, *, reason: str) -> float | None:
    number = _optional_float(value, reason=reason)
    if number is not None and number < 0:
        raise _RejectedRecord(reason)
    return number


def _coverage_context(value: str | bool | None) -> list[float]:
    if not isinstance(value, str):
        return []
    result: list[float] = []
    for item in value.split(","):
        if item in {"", "."}:
            continue
        try:
            number = float(item)
        except ValueError as exc:
            raise _RejectedRecord("malformed_coverage") from exc
        if number < 0:
            raise _RejectedRecord("malformed_coverage")
        result.append(number)
    return result


def _sv_length(
    event_type: EventType,
    primary: Locus,
    info: dict[str, str | bool],
) -> int | None:
    if event_type == EventType.TRANSLOCATION:
        return None
    raw_length = _first_value(info.get("SVLEN"))
    if raw_length is not None:
        try:
            length = abs(int(raw_length))
        except ValueError as exc:
            raise _RejectedRecord("malformed_svlen") from exc
    elif event_type == EventType.INSERTION:
        raise _RejectedRecord("missing_svlen")
    else:
        length = primary.end - primary.start
    if length < 1:
        raise _RejectedRecord("malformed_svlen")
    return length


def _format_vaf(fields: list[str]) -> float | None:
    if len(fields) < 10 or fields[8] in {"", "."} or fields[9] in {"", "."}:
        return None
    names = fields[8].split(":")
    values = fields[9].split(":")
    sample = dict(zip(names, values, strict=False))
    try:
        reference_reads = int(sample["DR"])
        variant_reads = int(sample["DV"])
    except (KeyError, ValueError):
        return None
    total = reference_reads + variant_reads
    return variant_reads / total if total > 0 else None


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


def _caller_filter_reason(filters: list[str]) -> str:
    recognized = sorted(set(filters).intersection(_KNOWN_SNIFFLES_FILTERS))
    return "filter_not_pass:" + (",".join(recognized) if recognized else "OTHER")


def _simple_locus(
    chromosome: str,
    position: int,
    event_type: EventType,
    info: dict[str, str | bool],
) -> Locus:
    start = position - 1
    if event_type == EventType.INSERTION:
        end = start + 1
    else:
        raw_end = _first_value(info.get("END"))
        if raw_end is not None:
            try:
                end = int(raw_end)
            except ValueError as exc:
                raise _RejectedRecord("malformed_end") from exc
        else:
            length = abs(_integer(info.get("SVLEN"), reason="missing_end_and_svlen"))
            end = start + max(length, 1)
    try:
        return Locus(chromosome=chromosome, start=start, end=end)
    except ValidationError as exc:
        raise _RejectedRecord("invalid_primary_locus") from exc


def _breakend_loci(
    chromosome: str,
    position: int,
    alternate: str,
    info: dict[str, str | bool],
) -> tuple[Locus, Locus]:
    try:
        resolved = resolve_breakend(
            alternate,
            declared_chromosome=info.get("CHR2"),
            declared_position=info.get("END"),
        )
    except BreakendParseError as exc:
        raise _RejectedRecord(exc.reason) from exc
    try:
        primary = Locus(chromosome=chromosome, start=position - 1, end=position)
        secondary = Locus(
            chromosome=resolved.mate_chromosome,
            start=resolved.mate_position_0based,
            end=resolved.mate_position_0based + 1,
        )
    except ValidationError as exc:
        raise _RejectedRecord("invalid_breakend_locus") from exc
    return primary, secondary


def _normalize_record(
    fields: list[str],
    *,
    record_number: int,
    policy: SnifflesPolicy,
    caller_version: str,
) -> GenomicEvent:
    if len(fields) < 8:
        raise _RejectedRecord("malformed_record")
    (
        chromosome,
        raw_position,
        _record_id,
        _reference,
        alternate,
        raw_quality,
        raw_filter,
        raw_info,
    ) = fields[:8]
    try:
        position = int(raw_position)
    except ValueError as exc:
        raise _RejectedRecord("malformed_position") from exc
    if position < 1:
        raise _RejectedRecord("malformed_position")

    filters = raw_filter.split(";") if raw_filter not in {"", "."} else []
    if policy.pass_only and filters != ["PASS"]:
        raise _RejectedRecord(_caller_filter_reason(filters))

    info = _parse_info(raw_info)
    raw_sv_type = _first_value(info.get("SVTYPE"))
    if raw_sv_type is None:
        raise _RejectedRecord("missing_svtype")
    event_type = _SV_TYPE_MAP.get(raw_sv_type.upper().split(":", maxsplit=1)[0])
    if event_type is None:
        raise _RejectedRecord("unsupported_svtype")
    if event_type not in policy.allowed_sv_types:
        raise _RejectedRecord("svtype_not_allowed")

    support = _integer(info.get("SUPPORT"), reason="missing_or_invalid_support")
    if support < policy.min_support:
        raise _RejectedRecord("support_below_policy")
    quality = _quality(raw_quality)
    if policy.minimum_quality is not None:
        if quality is None:
            raise _RejectedRecord("quality_missing")
        if quality < policy.minimum_quality:
            raise _RejectedRecord("quality_below_policy")
    if policy.require_precise and "PRECISE" not in info:
        raise _RejectedRecord("imprecise_record")

    if event_type == EventType.TRANSLOCATION:
        primary, secondary = _breakend_loci(chromosome, position, alternate, info)
    else:
        primary = _simple_locus(chromosome, position, event_type, info)
        secondary = None
    length_bp = _sv_length(event_type, primary, info)
    if length_bp is not None and length_bp < policy.min_sv_length:
        raise _RejectedRecord("sv_length_below_policy")

    coverage = _coverage_context(info.get("COVERAGE"))
    vaf = _optional_float(info.get("VAF"), reason="malformed_vaf")
    if vaf is None:
        vaf = _format_vaf(fields)
    if vaf is not None and not 0 <= vaf <= 1:
        raise _RejectedRecord("invalid_vaf")
    strand = _first_value(info.get("STRAND"))
    if strand is not None and strand not in {"+", "-", "+-"}:
        raise _RejectedRecord("malformed_supporting_read_strands")
    mean_nm = _optional_nonnegative_float(info.get("NM"), reason="malformed_nm")
    position_stdev = _optional_nonnegative_float(
        info.get("STDEV_POS"), reason="malformed_position_standard_deviation"
    )
    length_stdev = _optional_nonnegative_float(
        info.get("STDEV_LEN"), reason="malformed_length_standard_deviation"
    )
    if "PRECISE" in info and "IMPRECISE" in info:
        raise _RejectedRecord("conflicting_precision_flags")
    if "PRECISE" in info:
        precise = True
    elif "IMPRECISE" in info:
        precise = False
    else:
        precise = None
    notes = [
        "Normalized from a Sniffles2 VCF into 0-based, half-open coordinates.",
        "Candidate only; assay-specific analytical validation and expert review are pending.",
    ]
    if event_type == EventType.TRANSLOCATION:
        notes.append("Breakend evidence is not equivalent to an annotated gene fusion.")
    return GenomicEvent(
        event_id=f"SNIFFLES2-{record_number:06d}",
        event_type=event_type,
        primary=primary,
        secondary=secondary,
        length_bp=length_bp,
        evidence=[
            Evidence(
                caller="Sniffles2",
                caller_version=caller_version,
                support_reads=support,
                local_coverage=(sum(coverage) / len(coverage) if coverage else None),
                variant_allele_fraction=vaf,
                quality=quality,
                filters=filters,
                supporting_read_strands=strand,
                coverage_context=coverage,
                mean_alignment_nm=mean_nm,
                position_standard_deviation=position_stdev,
                length_standard_deviation=length_stdev,
                precise=precise,
            )
        ],
        confidence="unclassified",
        reportable=False,
        notes=notes,
    )


def normalize_sniffles_vcf(
    path: Path,
    *,
    sample_id: str,
    genome_build: GenomeBuild,
    policy: SnifflesPolicy,
    tool: ToolRecord,
) -> SnifflesCallReport:
    if not path.is_file():
        raise ValueError("Sniffles VCF is missing or unreadable")
    if tool.version != policy.expected_version:
        raise ValueError(
            f"Sniffles2 version {tool.version!r} does not match policy lock "
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
            event = _normalize_record(
                line.split("\t"),
                record_number=raw_count,
                policy=policy,
                caller_version=tool.version,
            )
        except _RejectedRecord as exc:
            rejections[exc.reason] += 1
        else:
            events.append(event)
    if not saw_fileformat or not saw_columns:
        raise ValueError("Sniffles output is not a complete VCF document")
    if sample_column_count not in {0, 1}:
        raise ValueError("Single-sample analysis requires at most one VCF sample column")

    warnings: list[str] = []
    if policy.status != "validated":
        warnings.append(policy.note)
    if rejections:
        warnings.append(
            "One or more VCF records were rejected by the explicit normalization policy."
        )
    if not events:
        warnings.append(
            "NO_CALL means no Sniffles2 record passed this technical policy; it is not a "
            "clinical negative result."
        )
    return SnifflesCallReport(
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
        vcf_fingerprint=FileFingerprint(
            size_bytes=path.stat().st_size,
            sha256=sha256_file(path),
        ),
        warnings=warnings,
        limitations=[
            "Only canonical chromosomes supported by the current normalized locus contract are "
            "retained.",
            "Read names, inserted sequences and source paths are deliberately excluded.",
            "Sniffles2 candidate calls are not validated clinical findings or ISCN assertions.",
            "The caller-level mapping-quality threshold is provenance; Sniffles2 does not emit "
            "per-event mapping quality in its standard VCF fields.",
        ],
    )


def sniffles_version(text: str) -> str:
    """Parse a Sniffles2 version from its probe output.

    Public because preflight has to reach the same answer this module will. A preflight
    that parsed versions differently from the run it precedes could clear a run that then
    fails on the version lock, which is worse than not checking at all.
    """
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    first_line = text.splitlines()[0].strip() if text.strip() else "unknown"
    return first_line[:80]


def _staged_vcf_path(output_vcf: Path) -> Path:
    """Reserve a unique same-filesystem path, then leave it absent for Sniffles to create."""
    if output_vcf.name.lower().endswith(".vcf.gz"):
        suffix = ".vcf.gz"
    else:
        suffix = output_vcf.suffix or ".vcf"
    descriptor, staged_name = tempfile.mkstemp(
        dir=output_vcf.parent,
        prefix=f".{output_vcf.name}.",
        suffix=suffix,
    )
    os.close(descriptor)
    staged = Path(staged_name)
    staged.unlink()
    return staged


def run_sniffles(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    policy: SnifflesPolicy,
    *,
    output_vcf: Path,
    runner: CommandRunner | None = None,
    sniffles: str = "sniffles",
    threads: int = 4,
) -> SnifflesCallReport:
    if manifest.input.kind != InputKind.ALIGNED_BAM:
        raise ValueError("Sniffles2 requires input.kind=aligned_bam")
    if manifest.sample_id != intake.sample_id:
        raise ValueError("Manifest and intake artifact must refer to the same sample")
    if manifest.assay.genome_build != intake.genome_build:
        raise ValueError("Manifest and intake artifact use different genome builds")
    if manifest.assay.reference_id != intake.reference_id:
        raise ValueError("Manifest and intake artifact use different reference IDs")
    if intake.verdict == Verdict.FAIL:
        raise ValueError("Sniffles2 cannot run after a failed aligned-BAM intake gate")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if output_vcf.exists():
        raise ValueError("Refusing to overwrite an existing Sniffles VCF")
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    command_runner = runner or SubprocessRunner()
    version_result = command_runner.run([sniffles, "--version"], timeout_seconds=30)
    if version_result.returncode != 0:
        raise ValueError("Sniffles2 version probe returned a non-zero exit code")
    version = sniffles_version(f"{version_result.stdout}\n{version_result.stderr}")
    if version != policy.expected_version:
        raise ValueError(
            f"Sniffles2 version {version!r} does not match policy lock {policy.expected_version!r}"
        )
    parameters: dict[str, str | int | bool] = {
        "threads": threads,
        "minsupport": policy.min_support,
        "minsvlen": policy.min_sv_length,
        "mapq": policy.mapq,
        "caller_pass_only": False,
        "normalizer_pass_only": True,
        "symbolic": True,
        "mosaic": policy.mode == SnifflesMode.MOSAIC,
        "output_read_names": False,
        "expected_version": policy.expected_version,
    }
    staged_vcf = _staged_vcf_path(output_vcf)
    argv = [
        sniffles,
        "--input",
        manifest.input.path,
        "--vcf",
        str(staged_vcf),
        "--threads",
        str(threads),
        "--minsupport",
        str(policy.min_support),
        "--minsvlen",
        str(policy.min_sv_length),
        "--mapq",
        str(policy.mapq),
        "--symbolic",
        "--no-progress",
    ]
    if policy.mode == SnifflesMode.MOSAIC:
        argv.append("--mosaic")
    try:
        result = command_runner.run(argv, timeout_seconds=7200)
        if result.returncode != 0:
            diagnostic = result.stderr.strip()[-2000:]
            detail = f": {diagnostic}" if diagnostic else ""
            raise ValueError(f"Sniffles2 failed with exit code {result.returncode}{detail}")
        if not staged_vcf.is_file():
            raise ValueError("Sniffles2 returned success but produced no VCF")
        tool = ToolRecord(name="Sniffles2", version=version, parameters=parameters)
        report = normalize_sniffles_vcf(
            staged_vcf,
            sample_id=manifest.sample_id,
            genome_build=manifest.assay.genome_build,
            policy=policy,
            tool=tool,
        )
        os.replace(staged_vcf, output_vcf)
        return report
    except BaseException:
        staged_vcf.unlink(missing_ok=True)
        raise

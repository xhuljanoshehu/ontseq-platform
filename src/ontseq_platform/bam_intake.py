from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .execution import CommandRunner, SubprocessRunner, ToolExecutionError
from .models import (
    AlignedBamIntakeReport,
    BamHeaderSummary,
    CheckStatus,
    FileFingerprint,
    InputKind,
    ReferenceLock,
    SampleManifest,
    ToolRecord,
    ValidationCheck,
    Verdict,
)
from .reference import contig_signature, sha256_file


@dataclass(frozen=True)
class ParsedBamHeader:
    sort_order: str | None
    contigs: tuple[tuple[str, int], ...]
    read_group_count: int
    sample_tag_count: int
    program_count: int


def _tags(fields: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for field in fields:
        key, separator, value = field.partition(":")
        if separator:
            parsed[key] = value
    return parsed


def parse_sam_header(text: str) -> ParsedBamHeader:
    sort_order: str | None = None
    contigs: list[tuple[str, int]] = []
    read_group_count = 0
    sample_tag_count = 0
    program_count = 0
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        if not raw_line.startswith("@"):
            continue
        fields = raw_line.split("\t")
        record_type = fields[0]
        tags = _tags(fields[1:])
        if record_type == "@HD":
            sort_order = tags.get("SO")
        elif record_type == "@SQ":
            name = tags.get("SN")
            raw_length = tags.get("LN")
            if not name or not raw_length:
                raise ValueError(f"Malformed @SQ record on header line {line_number}")
            try:
                length = int(raw_length)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid @SQ length on header line {line_number}: {raw_length!r}"
                ) from exc
            contigs.append((name, length))
        elif record_type == "@RG":
            read_group_count += 1
            if tags.get("SM"):
                sample_tag_count += 1
        elif record_type == "@PG":
            program_count += 1
    names = [name for name, _ in contigs]
    if len(names) != len(set(names)):
        raise ValueError("BAM header contains duplicate @SQ sequence names")
    return ParsedBamHeader(
        sort_order=sort_order,
        contigs=tuple(contigs),
        read_group_count=read_group_count,
        sample_tag_count=sample_tag_count,
        program_count=program_count,
    )


def _fingerprint(path: Path, *, checksum: bool) -> FileFingerprint:
    return FileFingerprint(
        size_bytes=path.stat().st_size, sha256=sha256_file(path) if checksum else None
    )


def _verdict(checks: list[ValidationCheck]) -> Verdict:
    if any(check.status == CheckStatus.FAIL for check in checks):
        return Verdict.FAIL
    if any(check.status == CheckStatus.WARN for check in checks):
        return Verdict.WARN
    return Verdict.PASS


def _check(
    checks: list[ValidationCheck],
    name: str,
    status: CheckStatus,
    message: str,
    **details: int | float | str | bool | None,
) -> None:
    checks.append(ValidationCheck(name=name, status=status, message=message, details=details))


def _reference_check(
    header: ParsedBamHeader, reference_lock: ReferenceLock
) -> tuple[CheckStatus, str, dict[str, int | str | bool | None]]:
    observed = dict(header.contigs)
    expected = {item.name: item.length for item in reference_lock.contigs}
    missing = [name for name in expected if name not in observed]
    length_mismatches = [
        name for name, length in expected.items() if name in observed and observed[name] != length
    ]
    extras = [name for name in observed if name not in expected]
    expected_present_order = [name for name in expected if name in observed]
    observed_locked_order = [name for name in observed if name in expected]
    order_mismatch = expected_present_order != observed_locked_order
    details: dict[str, int | str | bool | None] = {
        "expected_contigs": len(expected),
        "observed_contigs": len(observed),
        "expected_reference_bases": sum(expected.values()),
        "observed_reference_bases": sum(observed.values()),
        "expected_dictionary_sha256": contig_signature(expected.items()),
        "observed_dictionary_sha256": contig_signature(observed.items()),
        "missing_contigs": len(missing),
        "length_mismatches": len(length_mismatches),
        "extra_contigs": len(extras),
        "contig_order_mismatch": order_mismatch,
    }
    failed = (
        missing
        or length_mismatches
        or order_mismatch
        or (extras and not reference_lock.allow_extra_contigs)
    )
    if failed:
        order_detail = "; contig order differs" if order_mismatch else ""
        return (
            CheckStatus.FAIL,
            "Reference dictionary mismatch: "
            f"expected {len(expected)}, observed {len(observed)}; "
            f"{len(missing)} missing, {len(length_mismatches)} length mismatches, "
            f"{len(extras)} extra contigs{order_detail}",
            details,
        )
    if extras:
        return (
            CheckStatus.WARN,
            "Reference dictionary matched under the explicit extra-contig policy: "
            f"{len(extras)} extra contigs permitted",
            details,
        )
    return CheckStatus.PASS, "BAM sequence dictionary matches the reference lock", details


class AlignedBamInspector:
    def __init__(self, *, runner: CommandRunner | None = None, samtools: str = "samtools") -> None:
        self.runner = runner or SubprocessRunner()
        self.samtools = samtools

    def inspect(
        self,
        manifest: SampleManifest,
        reference_lock: ReferenceLock,
        *,
        include_checksums: bool = False,
    ) -> AlignedBamIntakeReport:
        if manifest.input.kind != InputKind.ALIGNED_BAM:
            raise ValueError("AlignedBamInspector requires input.kind=aligned_bam")
        if (
            manifest.input.index_path is None
        ):  # defensive; the manifest validator already enforces it
            raise ValueError("aligned_bam requires index_path")

        checks: list[ValidationCheck] = []
        bam_path = Path(manifest.input.path)
        index_path = Path(manifest.input.index_path)
        input_fingerprint: FileFingerprint | None = None
        index_fingerprint: FileFingerprint | None = None
        header_summary: BamHeaderSummary | None = None
        tool: ToolRecord | None = None

        if manifest.assay.reference_id == reference_lock.reference_id:
            _check(checks, "reference_id", CheckStatus.PASS, "Manifest reference ID matches lock")
        else:
            _check(
                checks,
                "reference_id",
                CheckStatus.FAIL,
                "Manifest reference ID does not match lock",
            )
        if manifest.assay.genome_build == reference_lock.genome_build:
            _check(checks, "genome_build", CheckStatus.PASS, "Manifest genome build matches lock")
        else:
            _check(
                checks,
                "genome_build",
                CheckStatus.FAIL,
                "Manifest genome build does not match lock",
            )

        if bam_path.is_file():
            checksum_required = include_checksums or manifest.input.sha256 is not None
            input_fingerprint = _fingerprint(bam_path, checksum=checksum_required)
            _check(
                checks,
                "bam_readable",
                CheckStatus.PASS,
                "Aligned BAM exists and is a regular file",
                size_bytes=input_fingerprint.size_bytes,
            )
            if bam_path.suffix.lower() != ".bam":
                _check(
                    checks,
                    "bam_extension",
                    CheckStatus.WARN,
                    "Input is declared as BAM but does not use the .bam extension",
                )
        else:
            _check(checks, "bam_readable", CheckStatus.FAIL, "Aligned BAM is missing or unreadable")

        if index_path.is_file():
            index_fingerprint = _fingerprint(index_path, checksum=include_checksums)
            _check(
                checks,
                "index_readable",
                CheckStatus.PASS,
                "BAM index exists and is a regular file",
                size_bytes=index_fingerprint.size_bytes,
            )
        else:
            _check(checks, "index_readable", CheckStatus.FAIL, "BAM index is missing or unreadable")

        if manifest.input.sha256 and input_fingerprint:
            if input_fingerprint.sha256 == manifest.input.sha256:
                _check(checks, "bam_sha256", CheckStatus.PASS, "BAM checksum matches manifest")
            else:
                _check(
                    checks, "bam_sha256", CheckStatus.FAIL, "BAM checksum does not match manifest"
                )

        if not bam_path.is_file() or not index_path.is_file():
            return self._report(
                manifest,
                reference_lock,
                checks,
                input_fingerprint,
                index_fingerprint,
                header_summary,
                tool,
            )

        try:
            version_result = self.runner.run([self.samtools, "--version"], timeout_seconds=30)
        except ToolExecutionError as exc:
            _check(checks, "samtools_available", CheckStatus.FAIL, str(exc))
            return self._report(
                manifest,
                reference_lock,
                checks,
                input_fingerprint,
                index_fingerprint,
                header_summary,
                tool,
            )
        if version_result.returncode != 0:
            _check(
                checks,
                "samtools_available",
                CheckStatus.FAIL,
                "samtools version probe returned a non-zero exit code",
                returncode=version_result.returncode,
            )
            return self._report(
                manifest,
                reference_lock,
                checks,
                input_fingerprint,
                index_fingerprint,
                header_summary,
                tool,
            )
        version_line = (
            version_result.stdout.splitlines()[0] if version_result.stdout else "samtools unknown"
        )
        version = version_line.removeprefix("samtools ").strip()
        tool = ToolRecord(
            name="samtools",
            version=version,
            parameters={"checks": ["quickcheck", "view -H", "idxstats -X"]},
        )
        _check(checks, "samtools_available", CheckStatus.PASS, "samtools is available")

        quickcheck = self.runner.run(
            [self.samtools, "quickcheck", "-v", str(bam_path)], timeout_seconds=120
        )
        if quickcheck.returncode == 0:
            _check(
                checks,
                "bam_quickcheck",
                CheckStatus.PASS,
                "BAM header and EOF block passed samtools quickcheck",
            )
        else:
            _check(
                checks,
                "bam_quickcheck",
                CheckStatus.FAIL,
                "samtools quickcheck detected a missing header or EOF problem",
                returncode=quickcheck.returncode,
            )

        header_result = self.runner.run(
            [self.samtools, "view", "-H", str(bam_path)], timeout_seconds=120
        )
        if header_result.returncode != 0:
            _check(
                checks,
                "bam_header",
                CheckStatus.FAIL,
                "samtools could not read the BAM header",
                returncode=header_result.returncode,
            )
        else:
            try:
                parsed_header = parse_sam_header(header_result.stdout)
            except ValueError as exc:
                _check(checks, "bam_header", CheckStatus.FAIL, str(exc))
            else:
                header_summary = BamHeaderSummary(
                    sort_order=parsed_header.sort_order,
                    sequence_count=len(parsed_header.contigs),
                    total_reference_bases=sum(length for _, length in parsed_header.contigs),
                    contig_signature_sha256=contig_signature(parsed_header.contigs),
                    read_group_count=parsed_header.read_group_count,
                    sample_tag_count=parsed_header.sample_tag_count,
                    program_count=parsed_header.program_count,
                )
                if parsed_header.contigs:
                    _check(
                        checks,
                        "bam_header",
                        CheckStatus.PASS,
                        "BAM header contains a sequence dictionary",
                        sequence_count=len(parsed_header.contigs),
                    )
                else:
                    _check(
                        checks,
                        "bam_header",
                        CheckStatus.FAIL,
                        "BAM header contains no @SQ sequence records",
                    )
                if parsed_header.sort_order == "coordinate":
                    _check(
                        checks,
                        "coordinate_sort",
                        CheckStatus.PASS,
                        "BAM header declares coordinate sorting",
                    )
                else:
                    _check(
                        checks,
                        "coordinate_sort",
                        CheckStatus.FAIL,
                        "BAM must declare SO:coordinate for indexed analysis",
                    )
                if parsed_header.read_group_count:
                    status = (
                        CheckStatus.PASS
                        if parsed_header.sample_tag_count == parsed_header.read_group_count
                        else CheckStatus.WARN
                    )
                    _check(
                        checks,
                        "read_groups",
                        status,
                        "BAM read-group metadata inspected",
                        read_group_count=parsed_header.read_group_count,
                        sample_tag_count=parsed_header.sample_tag_count,
                    )
                else:
                    _check(
                        checks,
                        "read_groups",
                        CheckStatus.WARN,
                        "BAM header contains no read groups; provenance may be incomplete",
                    )
                reference_status, message, details = _reference_check(parsed_header, reference_lock)
                _check(checks, "sequence_dictionary", reference_status, message, **details)

        idxstats = self.runner.run(
            [self.samtools, "idxstats", "-X", str(bam_path), str(index_path)],
            timeout_seconds=300,
        )
        if idxstats.returncode == 0 and idxstats.stdout.strip():
            malformed_rows = 0
            rows = 0
            for line in idxstats.stdout.splitlines():
                fields = line.split("\t")
                rows += 1
                if len(fields) != 4:
                    malformed_rows += 1
                    continue
                try:
                    int(fields[1])
                    int(fields[2])
                    int(fields[3])
                except ValueError:
                    malformed_rows += 1
            if malformed_rows:
                _check(
                    checks,
                    "bam_index",
                    CheckStatus.FAIL,
                    "samtools idxstats returned malformed rows",
                    rows=rows,
                    malformed_rows=malformed_rows,
                )
            else:
                _check(
                    checks,
                    "bam_index",
                    CheckStatus.PASS,
                    "The declared BAM index is readable and compatible with the BAM",
                    rows=rows,
                )
        else:
            _check(
                checks,
                "bam_index",
                CheckStatus.FAIL,
                "samtools idxstats could not use the declared BAM index",
                returncode=idxstats.returncode,
            )

        return self._report(
            manifest,
            reference_lock,
            checks,
            input_fingerprint,
            index_fingerprint,
            header_summary,
            tool,
        )

    @staticmethod
    def _report(
        manifest: SampleManifest,
        reference_lock: ReferenceLock,
        checks: list[ValidationCheck],
        input_fingerprint: FileFingerprint | None,
        index_fingerprint: FileFingerprint | None,
        header: BamHeaderSummary | None,
        tool: ToolRecord | None,
    ) -> AlignedBamIntakeReport:
        return AlignedBamIntakeReport(
            sample_id=manifest.sample_id,
            reference_id=reference_lock.reference_id,
            genome_build=reference_lock.genome_build,
            input_fingerprint=input_fingerprint,
            index_fingerprint=index_fingerprint,
            header=header,
            checks=checks,
            verdict=_verdict(checks),
            tool=tool,
            limitations=[
                "samtools quickcheck verifies the header and EOF block but not "
                "internal BAM corruption.",
                "A technical PASS is not evidence of adequate biological coverage "
                "or clinical validity.",
                "No patient or read-level content is copied into the intake report.",
            ],
        )

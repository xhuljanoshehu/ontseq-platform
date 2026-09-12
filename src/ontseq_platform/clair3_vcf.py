"""Read a Clair3 VCF into policy inputs, and refuse to run against the wrong binary.

Two jobs, both of them about not guessing.

*Reading.* A record only reaches :mod:`ontseq_platform.small_variants` if it carries the read
counts that policy needs. Clair3 reports depth as ``DP`` and allele support as ``AD``, and a
record missing either cannot be depth-gated at all. Such a record is counted as malformed and
dropped -- never accepted with an assumed depth, and never silently skipped, because a file
that lost most of its records to malformation looks identical to one that had few records and
the two mean opposite things. Multi-allelic records are refused rather than split: choosing
which alternate allele a single ``AD`` pair refers to is a decision this module is not in a
position to make correctly.

*Probing.* The version check exists because a caller that is not the pinned one produces
results the pinned parameters do not describe. :func:`clair3_version` is public for the same
reason its Sniffles2 counterpart is: preflight has to reach the same answer the run will, and
a preflight that parsed versions differently could clear a run that then fails the lock.

The model gate is the sharper of the two. Clair3 requires a model matched to the basecaller
chemistry, and a mismatched model does not error -- it quietly produces worse calls. So a
missing model pin stops the run rather than falling back to whatever the installation ships,
which is the one behaviour that would make the failure invisible.

Nothing here has been executed against the real Clair3 binary, so no verification status is
claimed for it (ADR-015). Research use only.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

from ontseq_platform.execution import CommandRunner
from ontseq_platform.small_variants import Clair3Policy, SmallVariant, SmallVariantError

_VERSION = re.compile(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)")

MALFORMED_NO_FORMAT = "record_has_no_format_or_sample_column"
MALFORMED_MISSING_DEPTH = "record_carries_no_dp_field"
MALFORMED_MISSING_ALLELE_DEPTH = "record_carries_no_ad_field"
MALFORMED_UNPARSABLE_COUNTS = "dp_or_ad_could_not_be_read_as_numbers"
MALFORMED_MULTI_ALLELIC = "record_is_multi_allelic_and_was_not_split"
MALFORMED_ZERO_DEPTH = "record_reports_zero_depth"
MALFORMED_INCONSISTENT_COUNTS = "variant_reads_exceed_reported_depth"


class Clair3ReadError(ValueError):
    """Raised when a VCF cannot be read as a Clair3 call set at all."""


@dataclass(frozen=True)
class Clair3VcfContents:
    """Everything a VCF yielded, including what it could not yield."""

    variants: tuple[SmallVariant, ...]
    malformed: tuple[tuple[int, str], ...]

    @property
    def malformed_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for _line, reason in self.malformed:
            counts[reason] = counts.get(reason, 0) + 1
        return counts

    @property
    def total_records(self) -> int:
        return len(self.variants) + len(self.malformed)


def _open_vcf(path: Path) -> Iterator[str]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield from handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield from handle


def _sample_fields(format_column: str, sample_column: str) -> dict[str, str]:
    keys = format_column.split(":")
    values = sample_column.split(":")
    return dict(zip(keys, values, strict=False))


def _variant_reads(allele_depth: str) -> int:
    """Alternate-allele support from a Clair3 ``AD`` field.

    ``AD`` is ``ref,alt``. The alternate count is taken rather than derived from depth minus
    reference, because at a site with a third allele those differ and the derived figure would
    silently overstate support.
    """
    parts = allele_depth.split(",")
    if len(parts) < 2:
        raise ValueError("AD does not carry a reference and alternate count")
    return int(parts[1])


def read_clair3_vcf(path: Path) -> Clair3VcfContents:
    """Parse a Clair3 VCF, keeping every record it could not use and why."""
    try:
        lines = list(_open_vcf(path))
    except OSError as error:
        raise Clair3ReadError(f"clair3 VCF is unreadable: {path}") from error

    variants: list[SmallVariant] = []
    malformed: list[tuple[int, str]] = []

    for number, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        columns = line.rstrip("\n").split("\t")
        if len(columns) < 10:
            malformed.append((number, MALFORMED_NO_FORMAT))
            continue

        chrom, position, _identifier, reference, alternate = columns[:5]
        quality_text, filter_status = columns[5], columns[6]
        fields = _sample_fields(columns[8], columns[9])

        if "," in alternate:
            malformed.append((number, MALFORMED_MULTI_ALLELIC))
            continue
        if "DP" not in fields:
            malformed.append((number, MALFORMED_MISSING_DEPTH))
            continue
        if "AD" not in fields:
            malformed.append((number, MALFORMED_MISSING_ALLELE_DEPTH))
            continue

        try:
            depth = int(fields["DP"])
            variant_reads = _variant_reads(fields["AD"])
            quality = 0.0 if quality_text == "." else float(quality_text)
            start = int(position)
        except ValueError:
            malformed.append((number, MALFORMED_UNPARSABLE_COUNTS))
            continue

        if depth <= 0:
            malformed.append((number, MALFORMED_ZERO_DEPTH))
            continue
        if variant_reads > depth:
            malformed.append((number, MALFORMED_INCONSISTENT_COUNTS))
            continue

        try:
            variants.append(
                SmallVariant(
                    chrom=chrom,
                    position=start,
                    reference=reference,
                    alternate=alternate,
                    depth=depth,
                    variant_reads=max(0, variant_reads),
                    quality=quality,
                    filter_status=filter_status,
                )
            )
        except SmallVariantError:
            malformed.append((number, MALFORMED_UNPARSABLE_COUNTS))

    return Clair3VcfContents(variants=tuple(variants), malformed=tuple(malformed))


def clair3_version(text: str) -> str:
    """Parse a Clair3 version from its probe output.

    Public because preflight has to reach the same answer the run will. A preflight that
    parsed versions differently from the run it precedes could clear a run that then fails the
    version lock, which is worse than not checking at all.
    """
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    first_line = text.splitlines()[0].strip() if text.strip() else "unknown"
    return first_line[:80]


@dataclass(frozen=True)
class PreconditionFailure:
    reason: str
    detail: str


def check_preconditions(
    command_runner: CommandRunner,
    *,
    binary: str,
    policy: Clair3Policy,
    argv: Sequence[str] = ("--version",),
    timeout_seconds: int = 30,
) -> PreconditionFailure | None:
    """Return the reason a run must not start, or ``None`` when it may.

    Returns rather than raises so a caller can report every blocked stage in one pass instead
    of stopping at the first, which is what preflight wants.
    """
    if policy.required_model_id is None:
        return PreconditionFailure(
            reason="model_not_pinned",
            detail=(
                "No Clair3 model is pinned. Clair3 requires a model matched to the basecaller "
                "chemistry, and a mismatched model does not error -- it quietly produces "
                "worse calls. Falling back to whatever the installation ships is the one "
                "behaviour that would make that failure invisible, so the run stops here."
            ),
        )

    result = command_runner.run([binary, *argv], timeout_seconds=timeout_seconds)
    if result.returncode != 0:
        return PreconditionFailure(
            reason="version_probe_failed",
            detail=f"clair3 version probe exited {result.returncode}",
        )

    observed = clair3_version(f"{result.stdout}\n{result.stderr}")
    if observed != policy.expected_version:
        return PreconditionFailure(
            reason="version_mismatch",
            detail=(
                f"clair3 version {observed!r} does not match policy lock "
                f"{policy.expected_version!r}"
            ),
        )
    return None

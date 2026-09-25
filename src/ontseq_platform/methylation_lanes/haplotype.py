"""Haplotype-resolved methylation from a haplotagged BAM (research only, issue #97).

The lane reports per region, modification code and haplotype (HP=1, HP=2, unphased) the
same count classes and fractions as the aggregate lane, so allele-specific methylation
becomes *observable*. It does not classify it: no imprinting, clonality or promoter-silencing
call is made, and a haplotype difference is shown only where both haplotypes are assessable.

What the pinned tool does, established with modkit 0.6.4 on synthetic haplotagged MM/ML
BAMs (``tests/test_haplotype_methylation_real_tool.py``) rather than inferred from its
documentation:

* ``modkit pileup --phased <out>`` writes exactly ``combined.bedmethyl``,
  ``hp1.bedmethyl`` and ``hp2.bedmethyl`` into the output directory; ``combined`` counts
  every record, ``hp1``/``hp2`` only records tagged ``HP=1``/``HP=2``;
* a haplotype without reads at a site has no row (absent, not zero), so unphased counts are
  ``combined - hp1 - hp2`` per site and count class, and must never be negative;
* a record with any other ``HP`` value (for example ``HP=3``) makes the pinned binary panic,
  and ``HP=0`` is silently counted as unphased — both are refused before the tool runs;
* phase sets are ignored: ``HP=1`` reads from two phase blocks are pooled into one
  haplotype, although nothing makes them the same parental allele.

The last point is load-bearing. Haplotype labels are only comparable inside one phase block,
so the lane reads the ``PS`` values of the haplotagged reads overlapping each region and
declares a region with more than one phase block ``NOT_ASSESSABLE`` instead of reporting a
pooled "haplotype" that may mix both alleles. Reads carrying ``HP`` without ``PS`` cannot be
placed in a block and are refused for the whole BAM.

The lane never phases on its own. It requires a haplotagging step recorded in the BAM
header by an accepted program; the command line is kept only as a SHA-256 because it
usually contains file paths.
"""

from __future__ import annotations

import hashlib
import importlib
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Any, Final, Literal

from pydantic import Field, model_validator

from ..execution import CommandRunner, SubprocessRunner
from ..methylation import (
    MODIFICATION_NAMES,
    MODKIT_MODIFICATION_NAMES,
    MethylationRegionSummary,
    ModificationCode,
    _bed_regions,
    _modkit_failed_processing_count,
    _modkit_include_bed,
    _Region,
    _Site,
    _summarize_region,
    count_reads_with_independent_cytosine_mod_groups,
    count_reads_with_modified_base_tags,
    modkit_version,
    parse_bedmethyl,
)
from ..methylation import _assign_sites as assign_sites
from ..models import (
    AlignedBamIntakeReport,
    FileFingerprint,
    GenomeBuild,
    InputKind,
    ModuleRunStatus,
    SampleManifest,
    StrictModel,
    ToolRecord,
    Verdict,
)
from ..modkit_build import PR709_BUILD_ID, identify_modkit_binary
from ..reference import sha256_file

#: The qualified ``--phased`` output layout. Anything else in the output directory, or any
#: of these missing, means the tool no longer behaves as it did when the lane was qualified.
PHASED_OUTPUT_FILES: tuple[str, str, str] = (
    "combined.bedmethyl",
    "hp1.bedmethyl",
    "hp2.bedmethyl",
)
PHASED_OUTPUT_LAYOUT = "modkit-0.6.4-phased-combined-hp1-hp2-v1"
InformativeReadDefinition = Literal[
    "mapped primary alignment overlapping the region, not secondary, supplementary, "
    "QC-failed or duplicate-flagged"
]
#: What counts as an informative read for the assessability rule. Recorded in the report.
INFORMATIVE_READ_DEFINITION: Final[InformativeReadDefinition] = (
    "mapped primary alignment overlapping the region, not secondary, supplementary, "
    "QC-failed or duplicate-flagged"
)
_EXCLUDED_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
_UNSUPPORTED_HP_EXPR = "exists([HP]) && [HP] != 1 && [HP] != 2"
_HP_WITHOUT_PS_EXPR = "exists([HP]) && !exists([PS])"
_HAS_HP_EXPR = "exists([HP])"


class HaplotypeMethylationPolicy(StrictModel):
    """Versioned technical settings of the haplotype lane. None is a validated threshold."""

    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    #: Pinned modkit version. The ``--phased`` layout above was qualified on 0.6.4 only.
    expected_version: str = Field(default="0.6.4", pattern=r"^\d+\.\d+\.\d+$")
    modification_codes: list[ModificationCode] = Field(
        default_factory=lambda: [ModificationCode.FIVE_MC], min_length=1
    )
    cpg_only: bool = True
    combine_strands: bool = True
    filter_threshold: float = Field(default=0.8, ge=0, le=1)
    #: Valid calls a site needs *within one haplotype* to enter that haplotype's aggregate.
    minimum_valid_coverage: int = Field(default=5, ge=1)
    #: Assessability rule: informative reads per haplotype in the region ...
    minimum_informative_reads_per_haplotype: int = Field(default=5, ge=1)
    #: ... and CpG sites at the coverage floor per haplotype. Below either: NOT_ASSESSABLE.
    minimum_sites_at_floor_per_haplotype: int = Field(default=3, ge=1)
    #: Programs whose haplotagging step, recorded in the BAM header, the lane accepts.
    accepted_phasing_programs: list[str] = Field(
        default_factory=lambda: ["whatshap", "longphase"], min_length=1
    )
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def policy_is_consistent(self) -> HaplotypeMethylationPolicy:
        if len(set(self.modification_codes)) != len(self.modification_codes):
            raise ValueError("Haplotype policy contains duplicate modification codes")
        if self.combine_strands and not self.cpg_only:
            raise ValueError("combine_strands requires cpg_only; strands are folded per CpG")
        if any(
            not name.strip() or name != name.strip().lower()
            for name in self.accepted_phasing_programs
        ):
            raise ValueError("Accepted phasing programs are lower-case program names")
        return self


class Haplotype(StrEnum):
    HP1 = "HP1"
    HP2 = "HP2"


class HaplotypeAssessability(StrEnum):
    ASSESSABLE = "ASSESSABLE"
    NOT_ASSESSABLE = "NOT_ASSESSABLE"


class NotAssessableReason(StrEnum):
    #: No haplotagged read overlaps the region.
    NO_PHASED_READS = "NO_PHASED_READS"
    #: Haplotagged reads come from more than one phase block; HP labels are not comparable.
    MULTIPLE_PHASE_BLOCKS = "MULTIPLE_PHASE_BLOCKS"
    #: Fewer informative reads for this haplotype than the policy requires.
    INSUFFICIENT_READS = "INSUFFICIENT_READS"
    #: Fewer CpG sites at the coverage floor for this haplotype than the policy requires.
    INSUFFICIENT_SITES = "INSUFFICIENT_SITES"


class PhasingProvenance(StrictModel):
    """A haplotagging step recorded in the BAM header (``@PG``)."""

    program: str = Field(min_length=1)
    version: str | None = None
    #: The command line usually names local files, so only its digest is retained.
    command_line_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class HaplotypeCall(StrictModel):
    """One haplotype of one region and modification code."""

    haplotype: Haplotype
    informative_reads: int = Field(ge=0)
    summary: MethylationRegionSummary
    assessability: HaplotypeAssessability
    reasons: list[NotAssessableReason] = Field(default_factory=list)

    @model_validator(mode="after")
    def assessability_is_explained(self) -> HaplotypeCall:
        assessable = self.assessability is HaplotypeAssessability.ASSESSABLE
        if assessable == bool(self.reasons):
            raise ValueError("A haplotype is assessable exactly when no reason blocks it")
        if assessable and self.summary.mean_modified_fraction is None:
            raise ValueError("An assessable haplotype must carry a measured fraction")
        return self


class HaplotypeRegionRow(StrictModel):
    """One region and modification code, split by haplotype."""

    region_id: str = Field(min_length=1)
    chromosome: str
    start: int = Field(ge=0)
    end: int = Field(gt=0)
    modification_code: ModificationCode
    modification_name: str
    #: Distinct phase sets (``PS``) of the haplotagged reads overlapping the region.
    phase_sets: list[int]
    hp1: HaplotypeCall
    hp2: HaplotypeCall
    #: Calls from reads without an HP tag. Always reported; never compared.
    unphased: MethylationRegionSummary
    unphased_informative_reads: int = Field(ge=0)
    #: HP1 minus HP2 call-weighted fraction, only when both haplotypes are assessable.
    haplotype_difference: float | None = Field(default=None, ge=-1, le=1)

    @model_validator(mode="after")
    def row_is_consistent(self) -> HaplotypeRegionRow:
        if self.end <= self.start:
            raise ValueError("Haplotype region end must be greater than start")
        if (self.hp1.haplotype, self.hp2.haplotype) != (Haplotype.HP1, Haplotype.HP2):
            raise ValueError("Haplotype calls are in HP1, HP2 order")
        if self.phase_sets != sorted(set(self.phase_sets)):
            raise ValueError("Phase sets are sorted and unique")
        both = all(
            call.assessability is HaplotypeAssessability.ASSESSABLE for call in (self.hp1, self.hp2)
        )
        if both != (self.haplotype_difference is not None):
            raise ValueError("A haplotype difference exists exactly when both are assessable")
        if self.haplotype_difference is not None:
            first = self.hp1.summary.mean_modified_fraction
            second = self.hp2.summary.mean_modified_fraction
            assert first is not None and second is not None
            if not math.isclose(self.haplotype_difference, first - second, abs_tol=1e-12):
                raise ValueError("Haplotype difference disagrees with the haplotype fractions")
        if both and len(self.phase_sets) != 1:
            raise ValueError("Haplotypes are only comparable inside a single phase block")
        return self


class HaplotypeMethylationReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: HaplotypeMethylationPolicy
    phasing: list[PhasingProvenance] = Field(min_length=1)
    reads_with_modified_base_tags: int | None = Field(default=None, ge=0)
    reads_with_haplotype_tags: int = Field(ge=1)
    informative_read_definition: InformativeReadDefinition = INFORMATIVE_READ_DEFINITION
    rows: list[HaplotypeRegionRow] = Field(default_factory=list)
    summary_metrics: dict[str, int] = Field(default_factory=dict)
    bedmethyl_fingerprints: dict[str, FileFingerprint]
    region_bed_fingerprint: FileFingerprint
    tool: ToolRecord
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def report_is_consistent(self) -> HaplotypeMethylationReport:
        if self.status not in {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL}:
            raise ValueError("A normalized haplotype report is either COMPLETED or NO_CALL")
        if self.tool.version != self.policy.expected_version:
            raise ValueError("modkit tool version does not match the haplotype policy lock")
        if set(self.bedmethyl_fingerprints) != {name.split(".")[0] for name in PHASED_OUTPUT_FILES}:
            raise ValueError("A haplotype report fingerprints exactly combined, hp1 and hp2")
        allowed = set(self.policy.modification_codes)
        if any(row.modification_code not in allowed for row in self.rows):
            raise ValueError("A row reports a modification code the policy does not allow")
        compared = any(row.haplotype_difference is not None for row in self.rows)
        expected = ModuleRunStatus.COMPLETED if compared else ModuleRunStatus.NO_CALL
        if self.status != expected:
            raise ValueError("Haplotype status must follow whether any region was comparable")
        return self


@dataclass(frozen=True)
class RegionPhasing:
    """What the haplotagged reads overlapping one region look like; no read identity."""

    phase_sets: frozenset[int] = frozenset()
    hp1_reads: int = 0
    hp2_reads: int = 0
    unphased_reads: int = 0


#: ``scan(bam, regions) -> {region: RegionPhasing}``; replaceable in unit tests.
PhasingScanner = Callable[[Path, Sequence[_Region]], Mapping[_Region, RegionPhasing]]


@dataclass
class _Tally:
    phase_sets: set[int] = field(default_factory=set)
    hp1: int = 0
    hp2: int = 0
    unphased: int = 0


def _pysam() -> Any:
    try:
        return importlib.import_module("pysam")
    except ImportError as exc:
        raise ValueError(
            "The haplotype lane reads HP/PS tags and the BAM header with pysam, which is not "
            "installed in this runtime"
        ) from exc


def scan_region_phasing(bam: Path, regions: Sequence[_Region]) -> dict[_Region, RegionPhasing]:
    """Count informative reads per haplotype and collect phase sets, region by region.

    Only tag values are read; read names and sequences are never retained.
    """
    pysam = _pysam()
    result: dict[_Region, RegionPhasing] = {}
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        for region in regions:
            tally = _Tally()
            for read in handle.fetch(region.chromosome, region.start, region.end):
                if read.flag & _EXCLUDED_FLAGS:
                    continue
                if not read.has_tag("HP"):
                    tally.unphased += 1
                    continue
                haplotype = read.get_tag("HP")
                if haplotype == 1:
                    tally.hp1 += 1
                elif haplotype == 2:
                    tally.hp2 += 1
                else:
                    raise ValueError("A read carries an HP value other than 1 or 2")
                if not read.has_tag("PS"):
                    raise ValueError("A haplotagged read carries no PS phase set")
                tally.phase_sets.add(int(read.get_tag("PS")))
            result[region] = RegionPhasing(
                phase_sets=frozenset(tally.phase_sets),
                hp1_reads=tally.hp1,
                hp2_reads=tally.hp2,
                unphased_reads=tally.unphased,
            )
    return result


def read_phasing_provenance(bam: Path, accepted_programs: Sequence[str]) -> list[PhasingProvenance]:
    """Return the accepted haplotagging steps recorded in the BAM header, or refuse."""
    pysam = _pysam()
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        records = handle.header.to_dict().get("PG", [])
    accepted = set(accepted_programs)
    found: list[PhasingProvenance] = []
    for record in records:
        program = str(record.get("PN") or record.get("ID") or "").strip()
        command = str(record.get("CL") or "")
        if program.lower() not in accepted or "haplotag" not in command.split():
            continue
        version = record.get("VN")
        found.append(
            PhasingProvenance(
                program=program.lower(),
                version=str(version) if version is not None else None,
                command_line_sha256=hashlib.sha256(command.encode("utf-8")).hexdigest(),
            )
        )
    if not found:
        raise ValueError(
            "The BAM header records no haplotagging step by an accepted phasing program "
            f"({', '.join(sorted(accepted))}). The lane never phases on its own and does not "
            "trust HP tags of unknown provenance"
        )
    return found


def _count(
    bam: Path, expression: str, *, runner: CommandRunner, samtools: str, threads: int
) -> int:
    result = runner.run(
        [samtools, "view", "-c", "-@", str(threads), "-e", expression, str(bam)],
        timeout_seconds=7200,
    )
    text = result.stdout.strip()
    if result.returncode != 0 or not text.isdigit():
        raise ValueError(
            "samtools could not evaluate a haplotype tag expression; refusing an unverified "
            "haplotype pileup"
        )
    return int(text)


_COUNT_FIELDS = (
    "valid_coverage",
    "modified_calls",
    "canonical_calls",
    "other_mod_calls",
    "fail_calls",
    "nocall_calls",
    "delete_calls",
    "diff_calls",
)


def split_unphased(
    combined: Sequence[_Site], hp1: Sequence[_Site], hp2: Sequence[_Site]
) -> list[_Site]:
    """Unphased calls per site: ``combined - hp1 - hp2``, refusing any negative remainder."""
    by_key = {(site.chromosome, site.start, site.code): site for site in combined}
    phased: dict[tuple[str, int, ModificationCode], list[_Site]] = {}
    for site in (*hp1, *hp2):
        key = (site.chromosome, site.start, site.code)
        if key not in by_key:
            raise ValueError("A haplotype pileup site is missing from the combined pileup")
        phased.setdefault(key, []).append(site)
    unphased: list[_Site] = []
    for key, total in by_key.items():
        parts = phased.get(key, [])
        values = {
            name: getattr(total, name) - sum(getattr(part, name) for part in parts)
            for name in _COUNT_FIELDS
        }
        if any(value < 0 for value in values.values()):
            raise ValueError(
                "Haplotype pileups exceed the combined pileup at one site; the phased "
                "output is inconsistent"
            )
        if not any(values.values()):
            continue
        unphased.append(
            _Site(
                chromosome=total.chromosome,
                start=total.start,
                end=total.end,
                code=total.code,
                **values,
            )
        )
    return unphased


def _assess(
    summary: MethylationRegionSummary,
    reads: int,
    phasing: RegionPhasing,
    policy: HaplotypeMethylationPolicy,
) -> tuple[HaplotypeAssessability, list[NotAssessableReason]]:
    reasons: list[NotAssessableReason] = []
    if not phasing.phase_sets:
        reasons.append(NotAssessableReason.NO_PHASED_READS)
    elif len(phasing.phase_sets) > 1:
        reasons.append(NotAssessableReason.MULTIPLE_PHASE_BLOCKS)
    if reads < policy.minimum_informative_reads_per_haplotype:
        reasons.append(NotAssessableReason.INSUFFICIENT_READS)
    if summary.sites_at_minimum_coverage < policy.minimum_sites_at_floor_per_haplotype:
        reasons.append(NotAssessableReason.INSUFFICIENT_SITES)
    status = HaplotypeAssessability.NOT_ASSESSABLE if reasons else HaplotypeAssessability.ASSESSABLE
    return status, reasons


_LIMITATIONS: tuple[str, ...] = (
    "Research use only. Haplotype-resolved fractions are descriptive; no imprinting, "
    "allele-specific-silencing or clonality statement is made, and parent of origin is not "
    "assigned.",
    "Haplotype labels are only as good as the recorded haplotagging step. A switch error "
    "inside a phase block mixes both alleles within one label and is not detectable here.",
    "The pinned modkit --phased mode ignores phase sets. Regions whose haplotagged reads span "
    "more than one phase block are therefore NOT_ASSESSABLE rather than pooled.",
    "Assessability thresholds (informative reads and CpG sites per haplotype) are technical "
    "defaults, not validated limits; below them a haplotype is NOT_ASSESSABLE, which is not "
    "a statement about its methylation.",
    "The --phased output layout was qualified against the pinned binary on synthetic "
    "haplotagged MM/ML fixtures (both strands). That is tool interoperability, not "
    "analytical recovery on biological data.",
    "No read names, per-read probabilities, source BAM path or phasing command line are "
    "copied into this report.",
)


def normalize_haplotype_methylation(
    *,
    sample_id: str,
    genome_build: GenomeBuild,
    phased_dir: Path,
    region_bed: Path,
    regions: Sequence[_Region],
    phasing_by_region: Mapping[_Region, RegionPhasing],
    policy: HaplotypeMethylationPolicy,
    tool: ToolRecord,
    phasing: Sequence[PhasingProvenance],
    reads_with_modified_base_tags: int | None,
    reads_with_haplotype_tags: int,
    warnings: Sequence[str] = (),
) -> HaplotypeMethylationReport:
    """Turn the qualified ``--phased`` output directory into the validated report."""
    present = sorted(path.name for path in phased_dir.iterdir()) if phased_dir.is_dir() else []
    if present != sorted(PHASED_OUTPUT_FILES):
        raise ValueError(
            "modkit --phased did not produce exactly the qualified combined/hp1/hp2 layout; "
            "refusing an unqualified output"
        )
    parsed = {
        name.split(".")[0]: parse_bedmethyl(
            phased_dir / name, allowed_codes=policy.modification_codes
        )
        for name in PHASED_OUTPUT_FILES
    }
    combined, _skipped = parsed["combined"]
    hp1_sites, _ = parsed["hp1"]
    hp2_sites, _ = parsed["hp2"]
    unphased_sites = split_unphased(combined, hp1_sites, hp2_sites)
    assigned = {
        label: assign_sites(sites, regions)
        for label, sites in (
            ("hp1", hp1_sites),
            ("hp2", hp2_sites),
            ("unphased", unphased_sites),
        )
    }
    rows: list[HaplotypeRegionRow] = []
    for region in regions:
        region_phasing = phasing_by_region.get(region)
        if region_phasing is None:
            raise ValueError(f"No phasing scan was recorded for region {region.region_id!r}")
        for code in policy.modification_codes:
            summaries = {
                label: _summarize_region(
                    region,
                    code,
                    buckets.get((region, code), []),
                    minimum_valid_coverage=policy.minimum_valid_coverage,
                )
                for label, buckets in assigned.items()
            }
            calls: dict[Haplotype, HaplotypeCall] = {}
            for haplotype, label, reads in (
                (Haplotype.HP1, "hp1", region_phasing.hp1_reads),
                (Haplotype.HP2, "hp2", region_phasing.hp2_reads),
            ):
                status, reasons = _assess(summaries[label], reads, region_phasing, policy)
                calls[haplotype] = HaplotypeCall(
                    haplotype=haplotype,
                    informative_reads=reads,
                    summary=summaries[label],
                    assessability=status,
                    reasons=reasons,
                )
            both = all(
                call.assessability is HaplotypeAssessability.ASSESSABLE for call in calls.values()
            )
            first = calls[Haplotype.HP1].summary.mean_modified_fraction
            second = calls[Haplotype.HP2].summary.mean_modified_fraction
            difference = (
                first - second if both and first is not None and second is not None else None
            )
            assert region.start is not None and region.end is not None
            rows.append(
                HaplotypeRegionRow(
                    region_id=region.region_id,
                    chromosome=region.chromosome,
                    start=region.start,
                    end=region.end,
                    modification_code=code,
                    modification_name=MODIFICATION_NAMES[code],
                    phase_sets=sorted(region_phasing.phase_sets),
                    hp1=calls[Haplotype.HP1],
                    hp2=calls[Haplotype.HP2],
                    unphased=summaries["unphased"],
                    unphased_informative_reads=region_phasing.unphased_reads,
                    haplotype_difference=difference,
                )
            )
    compared = sum(row.haplotype_difference is not None for row in rows)
    collected = [policy.note, *warnings]
    if not compared:
        collected.append(
            "No region met the assessability rule on both haplotypes inside one phase block. "
            "This NO_CALL says allele-specific methylation was not assessable at this depth "
            "and phasing; it does not say there is none."
        )
    blocked = {
        reason.value: sum(reason in call.reasons for row in rows for call in (row.hp1, row.hp2))
        for reason in NotAssessableReason
    }
    return HaplotypeMethylationReport(
        sample_id=sample_id,
        genome_build=genome_build,
        status=ModuleRunStatus.COMPLETED if compared else ModuleRunStatus.NO_CALL,
        policy=policy,
        phasing=list(phasing),
        reads_with_modified_base_tags=reads_with_modified_base_tags,
        reads_with_haplotype_tags=reads_with_haplotype_tags,
        rows=rows,
        summary_metrics={
            "region_count": len(regions),
            "row_count": len(rows),
            "comparable_row_count": compared,
            "combined_site_count": len(combined),
            **{f"haplotype_calls_{key.lower()}": value for key, value in blocked.items()},
        },
        bedmethyl_fingerprints={
            name.split(".")[0]: FileFingerprint(
                size_bytes=(phased_dir / name).stat().st_size,
                sha256=sha256_file(phased_dir / name),
            )
            for name in PHASED_OUTPUT_FILES
        },
        region_bed_fingerprint=FileFingerprint(
            size_bytes=region_bed.stat().st_size, sha256=sha256_file(region_bed)
        ),
        tool=tool,
        warnings=collected,
        limitations=list(_LIMITATIONS),
    )


def _discard_phased_output(phased_dir: Path) -> None:
    """Remove a partial ``--phased`` result this call created; nothing else is touched."""
    for name in PHASED_OUTPUT_FILES:
        (phased_dir / name).unlink(missing_ok=True)
    if phased_dir.is_dir() and not any(phased_dir.iterdir()):
        phased_dir.rmdir()


def run_haplotype_methylation(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    policy: HaplotypeMethylationPolicy,
    *,
    region_bed: Path,
    output_dir: Path,
    reference_fasta: Path | None,
    runner: CommandRunner | None = None,
    modkit: str = "modkit",
    samtools: str = "samtools",
    threads: int = 4,
    scan: PhasingScanner = scan_region_phasing,
    provenance: Callable[[Path, Sequence[str]], list[PhasingProvenance]] = read_phasing_provenance,
) -> HaplotypeMethylationReport:
    """Refuse unsupported inputs, run ``modkit pileup --phased`` and normalize the result."""
    if manifest.input.kind != InputKind.ALIGNED_BAM:
        raise ValueError("Haplotype methylation requires input.kind=aligned_bam")
    if (manifest.sample_id, manifest.assay.genome_build, manifest.assay.reference_id) != (
        intake.sample_id,
        intake.genome_build,
        intake.reference_id,
    ):
        raise ValueError("Manifest and intake artifact describe different samples or references")
    if intake.verdict == Verdict.FAIL:
        raise ValueError("Haplotype methylation cannot run after a failed aligned-BAM intake gate")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if reference_fasta is None or not reference_fasta.is_file():
        raise ValueError("modkit --modified-bases requires the locked reference FASTA")
    regions = _bed_regions(region_bed)
    bam = Path(manifest.input.path)

    output_dir.mkdir(parents=True, exist_ok=True)
    phased_dir = output_dir / f"{manifest.sample_id}.modkit-phased"
    log_path = output_dir / f"{manifest.sample_id}.modkit-phased.log"
    include_path = output_dir / f"{manifest.sample_id}.modkit-phased.include.bed"
    if phased_dir.exists() or log_path.exists() or include_path.exists():
        raise ValueError("Refusing to overwrite existing haplotype methylation outputs")

    command_runner = runner or SubprocessRunner()
    binary = identify_modkit_binary(modkit)
    modkit = binary.executable
    version_result = command_runner.run([modkit, "--version"], timeout_seconds=30)
    if version_result.returncode != 0:
        raise ValueError("modkit version probe returned a non-zero exit code")
    version = modkit_version(f"{version_result.stdout}\n{version_result.stderr}")
    if version != policy.expected_version:
        raise ValueError(
            f"modkit version {version!r} does not match policy lock {policy.expected_version!r}"
        )

    phasing = provenance(bam, policy.accepted_phasing_programs)
    tagged = count_reads_with_modified_base_tags(
        bam, runner=command_runner, samtools=samtools, threads=threads
    )
    if tagged == 0:
        raise ValueError(
            "The aligned BAM carries no MM modified-base tags; refusing a pileup that would "
            "read as unmethylated on both haplotypes"
        )
    haplotagged = _count(
        bam, _HAS_HP_EXPR, runner=command_runner, samtools=samtools, threads=threads
    )
    if haplotagged == 0:
        raise ValueError("The aligned BAM carries no HP haplotype tags; nothing to resolve")
    unsupported = _count(
        bam, _UNSUPPORTED_HP_EXPR, runner=command_runner, samtools=samtools, threads=threads
    )
    if unsupported:
        raise ValueError(
            f"{unsupported} read(s) carry an HP value other than 1 or 2. The pinned modkit "
            "--phased mode supports diploid labels only (other values abort it or are "
            "silently treated as unphased)"
        )
    without_phase_set = _count(
        bam, _HP_WITHOUT_PS_EXPR, runner=command_runner, samtools=samtools, threads=threads
    )
    if without_phase_set:
        raise ValueError(
            f"{without_phase_set} haplotagged read(s) carry no PS phase set, so their "
            "haplotype cannot be placed in a phase block"
        )
    if version == "0.6.4" and binary.build_id != PR709_BUILD_ID:
        independent = count_reads_with_independent_cytosine_mod_groups(
            bam, runner=command_runner, samtools=samtools, threads=threads
        )
        if independent is None or independent > 0:
            raise ValueError(
                "modkit 0.6.4 can silently lose or misassign calls for independent cytosine "
                "MM groups, and this BAM contains them or could not be checked"
            )
    phasing_by_region = scan(bam, regions)

    include_bed = _modkit_include_bed(region_bed, include_path)
    argv = [
        modkit,
        "pileup",
        str(bam),
        str(phased_dir),
        "--phased",
        "--threads",
        str(threads),
        "--filter-threshold",
        f"{policy.filter_threshold:g}",
        "--suppress-progress",
        "--log-filepath",
        str(log_path),
        "--modified-bases",
        *(MODKIT_MODIFICATION_NAMES[code] for code in policy.modification_codes),
        "--ref",
        str(reference_fasta),
        "--include-bed",
        str(include_bed),
    ]
    if policy.cpg_only:
        argv.append("--cpg")
        if policy.combine_strands:
            argv.append("--combine-strands")
    result = command_runner.run(argv, timeout_seconds=14400)
    if identify_modkit_binary(modkit).sha256 != binary.sha256:
        _discard_phased_output(phased_dir)
        raise ValueError("modkit executable changed during the haplotype pileup")
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        tail = detail[-1] if detail else "no diagnostic output"
        raise ValueError(
            f"modkit pileup --phased failed with exit code {result.returncode}: {tail}"
        )
    failed = _modkit_failed_processing_count(log_path, result.stdout, result.stderr)
    if failed:
        _discard_phased_output(phased_dir)
        raise ValueError(
            f"modkit pileup --phased reported {failed} failed record(s); refusing partial output"
        )
    parameters: dict[str, object] = {
        **binary.parameters(),
        "subcommand": "pileup",
        "phased": True,
        "phased_output_layout": PHASED_OUTPUT_LAYOUT,
        "threads": threads,
        "filter_threshold": policy.filter_threshold,
        "cpg_only": policy.cpg_only,
        "combine_strands": policy.combine_strands,
        "modification_codes": [code.value for code in policy.modification_codes],
        "minimum_valid_coverage": policy.minimum_valid_coverage,
        "minimum_informative_reads_per_haplotype": policy.minimum_informative_reads_per_haplotype,
        "minimum_sites_at_floor_per_haplotype": policy.minimum_sites_at_floor_per_haplotype,
        "include_bed_format": "bed3-or-bed6-v1",
        "expected_version": policy.expected_version,
    }
    return normalize_haplotype_methylation(
        sample_id=manifest.sample_id,
        genome_build=manifest.assay.genome_build,
        phased_dir=phased_dir,
        region_bed=region_bed,
        regions=regions,
        phasing_by_region=phasing_by_region,
        policy=policy,
        tool=ToolRecord(name="modkit", version=version, parameters=parameters),
        phasing=phasing,
        reads_with_modified_base_tags=tagged,
        reads_with_haplotype_tags=haplotagged,
        warnings=(
            []
            if tagged is not None
            else [
                "The installed samtools could not evaluate a tag filter expression, so the "
                "presence of MM/ML tags was not verified before the pileup ran."
            ]
        ),
    )

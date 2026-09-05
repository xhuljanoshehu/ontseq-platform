"""Modified-base pileup: from MM/ML tags to a normalized, region-aggregated methylation report.

The lane exists because the information is already in the data and is thrown away by
default. Dorado writes ``MM``/``ML`` tags, alignment carries them through unchanged
(:mod:`ontseq_platform.align`), and nothing downstream has ever read them. This adapter
reads them, via modkit's ``pileup``, and normalizes the result the same way every other
lane in this repository is normalized: a versioned policy, a locked tool version, an
explicit status, and numbers that state what they are.

Three decisions are load-bearing and deliberately not configurable away:

**An empty pileup is never "unmethylated".** A BAM basecalled without a modified-base
model carries no ``MM``/``ML`` tags at all, and modkit answers that with an empty file —
which looks exactly like a sample with no methylation. The two are separated before the
tool runs (a tag probe) and again afterwards (``NO_CALL`` with the reason named), because
a plausible-looking zero is the most dangerous output this lane could produce.

**The confidence threshold is pinned, not estimated.** modkit can pick a filter threshold
from the data it is looking at. That makes the parameter a function of the sample, so two
runs of the same pipeline are no longer running the same pipeline. The policy carries an
explicit threshold instead.

**A region with no site that met the coverage floor reports ``None``, not ``0.0``.** A
fraction of zero is a measurement; the absence of a measurement is not one.
"""

from __future__ import annotations

import gzip
import re
import statistics
from bisect import bisect_right
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
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

#: Columns in modkit's bedMethyl output (BED9+9). Fixed by the format, not by policy.
_BEDMETHYL_COLUMNS = 18


class ModificationCode(StrEnum):
    """The modified-base codes this adapter is prepared to normalize.

    Deliberately a closed set. modkit reports whatever the basecalling model emitted, and
    a code nobody planned for must surface as a refusal rather than be dropped on the
    floor: silently discarding ``h`` rows from a 5mC+5hmC model would understate total
    modification without anything in the report saying so.
    """

    FIVE_MC = "m"
    FIVE_HMC = "h"
    SIX_MA = "a"


#: What each code means, carried into the report so a reader never has to decode a letter.
MODIFICATION_NAMES: dict[ModificationCode, str] = {
    ModificationCode.FIVE_MC: "5-methylcytosine",
    ModificationCode.FIVE_HMC: "5-hydroxymethylcytosine",
    ModificationCode.SIX_MA: "N6-methyladenine",
}


class MethylationRegionSource(StrEnum):
    """What the rows of the report are aggregated over."""

    #: The locked Adaptive Sampling target design; one row per target and code.
    TARGET_BED = "target_bed"
    #: Canonical chromosomes; one row per chromosome and code.
    CHROMOSOME = "chromosome"


class MethylationPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    #: modkit is version-locked like every other caller. This default is an engineering
    #: pin, not a recommendation: re-pin it deliberately against the installed binary and
    #: record the change, because a pileup is not comparable across modkit majors.
    expected_version: str = Field(default="0.4.1", pattern=r"^\d+\.\d+\.\d+$")
    #: Codes the report is allowed to contain. A code in the pileup and not here fails the
    #: run rather than being discarded.
    modification_codes: list[ModificationCode] = Field(
        default_factory=lambda: [ModificationCode.FIVE_MC], min_length=1
    )
    #: Restrict the pileup to CpG dinucleotides. Requires the locked reference FASTA,
    #: because the motif is a property of the reference rather than of the reads.
    cpg_only: bool = True
    #: Fold the reverse-strand call of a CpG onto its forward-strand partner. Only
    #: meaningful with ``cpg_only``.
    combine_strands: bool = True
    #: Codes handed to ``--ignore``: their calls are folded into the canonical count
    #: rather than reported. Must not overlap ``modification_codes``.
    ignored_codes: list[ModificationCode] = Field(default_factory=list)
    #: modkit's per-call confidence threshold. Explicit on purpose — see the module
    #: docstring. 0.8 is a technical starting point and not a validated cut-off.
    filter_threshold: float = Field(default=0.8, ge=0, le=1)
    #: Sites below this valid coverage are excluded from every aggregate. They are still
    #: counted, so the report can say how much of the design was not measurable.
    minimum_valid_coverage: int = Field(default=5, ge=1)
    #: Aggregate over the target design or over chromosomes.
    region_source: MethylationRegionSource = MethylationRegionSource.CHROMOSOME
    #: Probe the BAM for MM/ML tags before running modkit. Costs one pass over the file
    #: and buys the difference between "no methylation" and "no modified-base data".
    verify_modified_base_tags: bool = True
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def policy_is_internally_consistent(self) -> MethylationPolicy:
        if len(set(self.modification_codes)) != len(self.modification_codes):
            raise ValueError("Methylation policy contains duplicate modification codes")
        if len(set(self.ignored_codes)) != len(self.ignored_codes):
            raise ValueError("Methylation policy contains duplicate ignored codes")
        overlap = set(self.ignored_codes) & set(self.modification_codes)
        if overlap:
            raise ValueError(
                "A modification code cannot be both reported and ignored: "
                + ", ".join(sorted(item.value for item in overlap))
            )
        if self.combine_strands and not self.cpg_only:
            raise ValueError("combine_strands requires cpg_only; strands are folded per CpG")
        return self


class MethylationRegionSummary(StrictModel):
    """One region, one modification code.

    ``start``/``end`` are present for BED-derived rows and absent for chromosome rows: the
    adapter aggregates the pileup it was given and does not read contig lengths, so it
    cannot honestly claim a chromosome span it never saw.
    """

    region_id: str = Field(min_length=1)
    chromosome: str = Field(pattern=r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
    start: int | None = Field(default=None, ge=0)
    end: int | None = Field(default=None, gt=0)
    modification_code: ModificationCode
    modification_name: str = Field(min_length=1)
    #: Pileup rows falling in this region for this code, before the coverage floor.
    sites_total: int = Field(ge=0)
    #: Rows that met ``minimum_valid_coverage`` and therefore entered the aggregates.
    sites_at_minimum_coverage: int = Field(ge=0)
    #: Summed valid calls and modified calls over the sites that entered the aggregates.
    valid_call_count: int = Field(ge=0)
    modified_call_count: int = Field(ge=0)
    #: Call-weighted modified fraction. ``None`` when no site met the coverage floor —
    #: never 0.0, which would read as a measured absence of methylation.
    mean_modified_fraction: float | None = Field(default=None, ge=0, le=1)
    #: Unweighted median across qualifying sites, which a single deep site cannot skew.
    median_site_modified_fraction: float | None = Field(default=None, ge=0, le=1)
    mean_valid_coverage: float | None = Field(default=None, ge=0)

    @model_validator(mode="after")
    def summary_is_consistent(self) -> MethylationRegionSummary:
        if (self.start is None) != (self.end is None):
            raise ValueError("A methylation region declares both coordinates or neither")
        if self.start is not None and self.end is not None and self.end <= self.start:
            raise ValueError("Methylation region end must be greater than start")
        if self.sites_at_minimum_coverage > self.sites_total:
            raise ValueError("More qualifying sites than sites were counted")
        if self.modified_call_count > self.valid_call_count:
            raise ValueError("Modified calls cannot exceed valid calls")
        measured = self.sites_at_minimum_coverage > 0
        if measured != (self.mean_modified_fraction is not None):
            raise ValueError("mean_modified_fraction must be present exactly when a site qualified")
        if measured != (self.median_site_modified_fraction is not None):
            raise ValueError(
                "median_site_modified_fraction must be present exactly when a site qualified"
            )
        if measured != (self.mean_valid_coverage is not None):
            raise ValueError("mean_valid_coverage must be present exactly when a site qualified")
        if not measured and self.valid_call_count:
            raise ValueError("A region with no qualifying site cannot carry valid calls")
        return self


class MethylationReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    sample_id: str
    genome_build: GenomeBuild
    status: ModuleRunStatus
    policy: MethylationPolicy
    region_source: MethylationRegionSource
    #: How many reads carried MM/ML tags, when the probe ran. ``None`` means the probe was
    #: disabled or the installed samtools could not answer — not that the answer was zero.
    reads_with_modified_base_tags: int | None = Field(default=None, ge=0)
    summary_metrics: dict[str, float | int] = Field(default_factory=dict)
    regions: list[MethylationRegionSummary] = Field(default_factory=list)
    bedmethyl_fingerprint: FileFingerprint
    target_bed_fingerprint: FileFingerprint | None = None
    tool: ToolRecord
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def report_is_consistent(self) -> MethylationReport:
        if self.status not in {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL}:
            raise ValueError("A normalized methylation report is either COMPLETED or NO_CALL")
        if self.tool.version != self.policy.expected_version:
            raise ValueError("modkit tool version does not match the methylation policy lock")
        allowed = set(self.policy.modification_codes)
        if any(region.modification_code not in allowed for region in self.regions):
            raise ValueError("A region reports a modification code the policy does not allow")
        measured = any(region.sites_at_minimum_coverage > 0 for region in self.regions)
        expected_status = ModuleRunStatus.COMPLETED if measured else ModuleRunStatus.NO_CALL
        if self.status != expected_status:
            raise ValueError("Methylation status is inconsistent with the aggregated sites")
        if self.region_source != self.policy.region_source:
            raise ValueError("Report region source does not match the policy")
        if self.summary_metrics.get("region_row_count") != len(self.regions):
            raise ValueError("Methylation region_row_count is inconsistent")
        return self


@dataclass(frozen=True)
class _Region:
    """An aggregation bucket: a named interval, or a whole chromosome when unbounded."""

    region_id: str
    chromosome: str
    start: int | None
    end: int | None


@dataclass(frozen=True)
class _Site:
    chromosome: str
    start: int
    end: int
    code: ModificationCode
    valid_coverage: int
    modified_calls: int

    @property
    def modified_fraction(self) -> float:
        return self.modified_calls / self.valid_coverage


def modkit_version(text: str) -> str:
    """Parse modkit's ``--version`` banner (``mod_kit 0.4.1``).

    Public for the same reason :func:`ontseq_platform.target_coverage.mosdepth_version` is:
    a preflight that reads the version differently from the adapter can clear a run the
    adapter then refuses.
    """
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    first_line = text.splitlines()[0].strip() if text.strip() else "unknown"
    return first_line[:80]


def _open_text(path: Path) -> list[str]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return handle.read().splitlines()
    return path.read_text(encoding="utf-8").splitlines()


def _split_bedmethyl(line: str) -> list[str]:
    """Split one bedMethyl row.

    The adapter passes ``--only-tabs`` so every column is tab-separated. Older modkit
    builds separate the nine trailing count columns with spaces instead, and a file
    produced that way is still perfectly readable — so a row that does not split into the
    expected column count on tabs is retried on arbitrary whitespace. No bedMethyl column
    may contain a space, which is what makes the retry safe rather than a guess.
    """
    fields = line.split("\t")
    if len(fields) == _BEDMETHYL_COLUMNS:
        return fields
    return line.split()


def parse_bedmethyl(
    path: Path,
    *,
    allowed_codes: Sequence[ModificationCode],
) -> tuple[list[_Site], int]:
    """Read a modkit bedMethyl file into sites, returning the non-canonical rows skipped.

    A modification code outside ``allowed_codes`` is refused rather than ignored: the
    policy declares what the run is reporting, and a model that also emitted 5hmC changes
    what a 5mC fraction means.
    """
    if not path.is_file():
        raise ValueError("modkit bedMethyl output is missing")
    permitted = set(allowed_codes)
    sites: list[_Site] = []
    skipped_non_canonical = 0
    seen: set[tuple[str, int, str]] = set()
    for line_number, line in enumerate(_open_text(path), start=1):
        if not line or line.startswith(("#", "track ", "browser ")):
            continue
        fields = _split_bedmethyl(line)
        if len(fields) != _BEDMETHYL_COLUMNS:
            raise ValueError(
                f"bedMethyl line {line_number} has {len(fields)} columns, expected "
                f"{_BEDMETHYL_COLUMNS}"
            )
        chromosome = fields[0]
        raw_code = fields[3]
        try:
            code = ModificationCode(raw_code)
        except ValueError as exc:
            raise ValueError(
                f"bedMethyl line {line_number} reports modification code {raw_code!r}, which "
                "this adapter cannot normalize"
            ) from exc
        if code not in permitted:
            raise ValueError(
                f"bedMethyl contains modification code {raw_code!r}, which the policy does not "
                "list. Dropping it would silently change what the reported fractions mean"
            )
        if _CANONICAL_CHROMOSOME.fullmatch(chromosome) is None:
            skipped_non_canonical += 1
            continue
        start = _parse_int(fields[1], field=f"bedMethyl line {line_number} start")
        end = _parse_int(fields[2], field=f"bedMethyl line {line_number} end")
        if end <= start:
            raise ValueError(f"bedMethyl line {line_number} has invalid coordinates")
        valid_coverage = _parse_int(fields[9], field=f"bedMethyl line {line_number} valid coverage")
        modified_calls = _parse_int(fields[11], field=f"bedMethyl line {line_number} Nmod")
        if modified_calls > valid_coverage:
            raise ValueError(
                f"bedMethyl line {line_number} reports more modified calls than valid calls"
            )
        key = (chromosome, start, raw_code)
        if key in seen:
            raise ValueError("bedMethyl contains a duplicate site for one modification code")
        seen.add(key)
        sites.append(
            _Site(
                chromosome=chromosome,
                start=start,
                end=end,
                code=code,
                valid_coverage=valid_coverage,
                modified_calls=modified_calls,
            )
        )
    return sites, skipped_non_canonical


def _parse_int(raw: str, *, field: str) -> int:
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"Invalid integer in {field}") from exc
    if value < 0:
        raise ValueError(f"Negative value in {field}")
    return value


def _chromosome_regions(sites: Sequence[_Site]) -> list[_Region]:
    ordered = sorted({site.chromosome for site in sites})
    return [
        _Region(region_id=chromosome, chromosome=chromosome, start=None, end=None)
        for chromosome in ordered
    ]


def _bed_regions(target_bed: Path) -> list[_Region]:
    """Reuse the target-coverage BED reader so both lanes agree on what a target is."""
    from .target_coverage import load_target_bed

    return [
        _Region(
            region_id=item.region_id,
            chromosome=item.chromosome,
            start=item.start,
            end=item.end,
        )
        for item in load_target_bed(target_bed)
    ]


def _assign_sites(
    sites: Sequence[_Site], regions: Sequence[_Region]
) -> dict[tuple[_Region, ModificationCode], list[_Site]]:
    """Bucket sites by region and code.

    Regions may overlap — a buffered panel routinely does — so a site is counted in every
    region containing it. That inflates nothing: each row is an independent statement
    about one design interval, exactly as the target-coverage lane treats overlapping
    intervals.
    """
    # Keyed by the region itself, not by ``region_id``: a BED names intervals and a panel
    # routinely gives several the same name (two exons of one gene). ``load_target_bed``
    # rejects duplicate *coordinates*, not duplicate names, so keying on the name merged
    # those intervals into one bucket and then reported it once per interval — each row
    # carrying sites from outside its own declared range, with the counts double-counted
    # across the report.
    buckets: dict[tuple[_Region, ModificationCode], list[_Site]] = {}
    bounded = [region for region in regions if region.start is not None]
    unbounded = {region.chromosome: region for region in regions if region.start is None}

    by_chromosome: dict[str, list[_Region]] = {}
    for region in bounded:
        by_chromosome.setdefault(region.chromosome, []).append(region)
    starts: dict[str, list[int]] = {}
    for chromosome, items in by_chromosome.items():
        items.sort(key=lambda item: (item.start or 0, item.end or 0))
        starts[chromosome] = [item.start or 0 for item in items]

    for site in sites:
        whole = unbounded.get(site.chromosome)
        if whole is not None:
            buckets.setdefault((whole, site.code), []).append(site)
        candidates = by_chromosome.get(site.chromosome)
        if not candidates:
            continue
        # Every interval starting at or before the site is a candidate; the end is what
        # decides. Intervals are few (a panel, not a genome), so the linear scan back from
        # the insertion point stays cheap and cannot miss a long enclosing interval.
        limit = bisect_right(starts[site.chromosome], site.start)
        for region in candidates[:limit]:
            if region.end is not None and site.start < region.end:
                buckets.setdefault((region, site.code), []).append(site)
    return buckets


def _summarize_region(
    region: _Region,
    code: ModificationCode,
    sites: Sequence[_Site],
    *,
    minimum_valid_coverage: int,
) -> MethylationRegionSummary:
    # The policy floor is at least 1, so a qualifying site always carries at least one
    # valid call. That is what lets "a site qualified" and "a fraction exists" be the same
    # condition, which the report model then enforces.
    qualifying = [site for site in sites if site.valid_coverage >= minimum_valid_coverage]
    valid_calls = sum(site.valid_coverage for site in qualifying)
    modified_calls = sum(site.modified_calls for site in qualifying)
    mean_fraction: float | None = None
    median_fraction: float | None = None
    mean_coverage: float | None = None
    if qualifying:
        mean_fraction = modified_calls / valid_calls
        median_fraction = statistics.median(site.modified_fraction for site in qualifying)
        mean_coverage = valid_calls / len(qualifying)
    return MethylationRegionSummary(
        region_id=region.region_id,
        chromosome=region.chromosome,
        start=region.start,
        end=region.end,
        modification_code=code,
        modification_name=MODIFICATION_NAMES[code],
        sites_total=len(sites),
        sites_at_minimum_coverage=len(qualifying),
        valid_call_count=valid_calls,
        modified_call_count=modified_calls,
        mean_modified_fraction=mean_fraction,
        median_site_modified_fraction=median_fraction,
        mean_valid_coverage=mean_coverage,
    )


_LIMITATIONS: tuple[str, ...] = (
    "Modified-base fractions are descriptive technical measurements. No methylation "
    "threshold, region set or classifier in this repository is analytically validated.",
    "The modkit adapter has not been executed against the real binary in this repository's "
    "continuous integration; its behaviour on real modified-base data is an assumption.",
    "Aggregated fractions depend on the basecalling model that produced the MM/ML tags. "
    "Runs basecalled with different models are not comparable.",
    "Strand-folded CpG values combine both strands of one dinucleotide; they are not "
    "per-strand measurements.",
    "No read names, per-read modification probabilities or source BAM path are copied into "
    "this report.",
)


def normalize_methylation(
    *,
    sample_id: str,
    genome_build: GenomeBuild,
    bedmethyl_path: Path,
    policy: MethylationPolicy,
    tool: ToolRecord,
    regions: Sequence[_Region] | None = None,
    target_bed: Path | None = None,
    reads_with_modified_base_tags: int | None = None,
    warnings: Sequence[str] = (),
) -> MethylationReport:
    """Turn a bedMethyl file into the validated report contract."""
    if tool.version != policy.expected_version:
        raise ValueError(
            f"modkit version {tool.version!r} does not match policy lock "
            f"{policy.expected_version!r}"
        )
    sites, skipped_non_canonical = parse_bedmethyl(
        bedmethyl_path, allowed_codes=policy.modification_codes
    )
    buckets_regions = list(regions) if regions is not None else _chromosome_regions(sites)
    assigned = _assign_sites(sites, buckets_regions)

    summaries: list[MethylationRegionSummary] = []
    for region in buckets_regions:
        for code in policy.modification_codes:
            summaries.append(
                _summarize_region(
                    region,
                    code,
                    assigned.get((region, code), []),
                    minimum_valid_coverage=policy.minimum_valid_coverage,
                )
            )

    qualifying_sites = sum(item.sites_at_minimum_coverage for item in summaries)
    collected: list[str] = [policy.note, *warnings]
    summary_metrics: dict[str, float | int] = {
        "region_row_count": len(summaries),
        "site_count": len(sites),
        "skipped_non_canonical_site_count": skipped_non_canonical,
        "minimum_valid_coverage": policy.minimum_valid_coverage,
    }
    for code in policy.modification_codes:
        code_sites = [site for site in sites if site.code == code]
        qualifying = [
            site for site in code_sites if site.valid_coverage >= policy.minimum_valid_coverage
        ]
        summary_metrics[f"site_count_{code.value}"] = len(code_sites)
        summary_metrics[f"site_count_at_minimum_coverage_{code.value}"] = len(qualifying)
        valid_calls = sum(site.valid_coverage for site in qualifying)
        if qualifying:
            modified = sum(site.modified_calls for site in qualifying)
            summary_metrics[f"mean_modified_fraction_{code.value}"] = modified / valid_calls
            summary_metrics[f"mean_valid_coverage_{code.value}"] = valid_calls / len(qualifying)
    if skipped_non_canonical:
        collected.append(
            f"{skipped_non_canonical} pileup row(s) on non-canonical contigs were excluded; "
            "this report describes chromosomes 1-22, X and Y only."
        )
    if not qualifying_sites:
        collected.append(
            "No site reached the configured valid-coverage floor. This NO_CALL says the "
            "sample was not measurable at this depth; it does not say the DNA is "
            "unmethylated."
        )

    return MethylationReport(
        sample_id=sample_id,
        genome_build=genome_build,
        status=ModuleRunStatus.COMPLETED if qualifying_sites else ModuleRunStatus.NO_CALL,
        policy=policy,
        region_source=policy.region_source,
        reads_with_modified_base_tags=reads_with_modified_base_tags,
        summary_metrics=summary_metrics,
        regions=summaries,
        bedmethyl_fingerprint=FileFingerprint(
            size_bytes=bedmethyl_path.stat().st_size,
            sha256=sha256_file(bedmethyl_path),
        ),
        target_bed_fingerprint=(
            FileFingerprint(size_bytes=target_bed.stat().st_size, sha256=sha256_file(target_bed))
            if target_bed is not None
            else None
        ),
        tool=tool,
        warnings=collected,
        limitations=list(_LIMITATIONS),
    )


def count_reads_with_modified_base_tags(
    bam: Path,
    *,
    runner: CommandRunner,
    samtools: str = "samtools",
    threads: int = 4,
) -> int | None:
    """Count reads carrying an ``MM`` tag, or return ``None`` when samtools cannot say.

    A samtools without filter expressions cannot answer this, and that is not a reason to
    refuse the run — it is a reason to say the question went unanswered. ``None`` and
    ``0`` therefore mean different things here, and the caller treats only ``0`` as fatal.
    """
    result = runner.run(
        [samtools, "view", "-c", "-@", str(threads), "-e", "[MM]", str(bam)],
        timeout_seconds=7200,
    )
    if result.returncode != 0:
        return None
    text = result.stdout.strip()
    if not text.isdigit():
        return None
    return int(text)


def _build_argv(
    *,
    modkit: str,
    bam: Path,
    output_path: Path,
    log_path: Path,
    policy: MethylationPolicy,
    reference_fasta: Path | None,
    include_bed: Path | None,
    threads: int,
) -> list[str]:
    argv = [
        modkit,
        "pileup",
        str(bam),
        str(output_path),
        "--threads",
        str(threads),
        "--filter-threshold",
        f"{policy.filter_threshold:g}",
        "--only-tabs",
        "--suppress-progress",
        "--log-filepath",
        str(log_path),
    ]
    if policy.cpg_only:
        if reference_fasta is None:
            raise ValueError("cpg_only requires the locked reference FASTA")
        argv.extend(["--cpg", "--ref", str(reference_fasta)])
        if policy.combine_strands:
            argv.append("--combine-strands")
    for code in policy.ignored_codes:
        argv.extend(["--ignore", code.value])
    if include_bed is not None:
        argv.extend(["--include-bed", str(include_bed)])
    return argv


def run_methylation(
    manifest: SampleManifest,
    intake: AlignedBamIntakeReport,
    policy: MethylationPolicy,
    *,
    output_dir: Path,
    reference_fasta: Path | None = None,
    runner: CommandRunner | None = None,
    modkit: str = "modkit",
    samtools: str = "samtools",
    threads: int = 4,
) -> MethylationReport:
    """Run ``modkit pileup`` over an aligned BAM and normalize the result."""
    if manifest.input.kind != InputKind.ALIGNED_BAM:
        raise ValueError("Methylation requires input.kind=aligned_bam")
    if manifest.sample_id != intake.sample_id:
        raise ValueError("Manifest and intake artifact must refer to the same sample")
    if manifest.assay.genome_build != intake.genome_build:
        raise ValueError("Manifest and intake artifact use different genome builds")
    if manifest.assay.reference_id != intake.reference_id:
        raise ValueError("Manifest and intake artifact use different reference IDs")
    if intake.verdict == Verdict.FAIL:
        raise ValueError("Methylation cannot run after a failed aligned-BAM intake gate")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if policy.cpg_only and reference_fasta is None:
        raise ValueError(
            "cpg_only restricts the pileup to a reference motif and therefore requires the "
            "locked reference FASTA"
        )
    if reference_fasta is not None and not reference_fasta.is_file():
        raise ValueError("Reference FASTA is missing or unreadable")

    warnings: list[str] = []
    target_bed: Path | None = None
    regions: list[_Region] | None = None
    if policy.region_source == MethylationRegionSource.TARGET_BED:
        if not manifest.assay.target_bed:
            raise ValueError(
                "region_source=target_bed was selected but the manifest declares no target "
                "BED. Falling back to chromosomes would report a different question than "
                "the one that was asked"
            )
        target_bed = Path(manifest.assay.target_bed)
        regions = _bed_regions(target_bed)
    elif manifest.assay.mode == AssayMode.ADAPTIVE_SAMPLING and manifest.assay.target_bed:
        # Enrichment leaves the off-target genome at a depth where a chromosome-wide
        # fraction mixes measured targets with barely-observed background. The design is
        # recorded because it is real context for the run, but it did not constrain this
        # pileup — and a bare checksum in the report reads as though it had, so say so.
        target_bed = Path(manifest.assay.target_bed)
        warnings.append(
            "This run is enriched and the report records its target design, but the policy "
            "aggregates over chromosomes, so the pileup was not restricted to that design. "
            "Every fraction here mixes enriched targets with off-target background at a "
            "depth the design never intended. Set region_source=target_bed to aggregate "
            "over the design instead."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    bedmethyl_path = output_dir / f"{manifest.sample_id}.modkit.bedmethyl"
    log_path = output_dir / f"{manifest.sample_id}.modkit.log"
    if bedmethyl_path.exists() or log_path.exists():
        raise ValueError("Refusing to overwrite existing modkit methylation outputs")

    command_runner = runner or SubprocessRunner()
    version_result = command_runner.run([modkit, "--version"], timeout_seconds=30)
    if version_result.returncode != 0:
        raise ValueError("modkit version probe returned a non-zero exit code")
    version = modkit_version(f"{version_result.stdout}\n{version_result.stderr}")
    if version != policy.expected_version:
        raise ValueError(
            f"modkit version {version!r} does not match policy lock {policy.expected_version!r}"
        )

    tagged_reads: int | None = None
    if policy.verify_modified_base_tags:
        tagged_reads = count_reads_with_modified_base_tags(
            Path(manifest.input.path),
            runner=command_runner,
            samtools=samtools,
            threads=threads,
        )
        if tagged_reads == 0:
            raise ValueError(
                "The aligned BAM carries no MM modified-base tags, so no methylation can be "
                "inferred from it. Refusing to produce an empty pileup, which would read as "
                "an unmethylated sample. Re-basecall with a modified-base model"
            )
        if tagged_reads is None:
            warnings.append(
                "The installed samtools could not evaluate a tag filter expression, so the "
                "presence of MM/ML tags was not verified before the pileup ran."
            )

    include_bed = target_bed if policy.region_source == MethylationRegionSource.TARGET_BED else None
    argv = _build_argv(
        modkit=modkit,
        bam=Path(manifest.input.path),
        output_path=bedmethyl_path,
        log_path=log_path,
        policy=policy,
        reference_fasta=reference_fasta,
        include_bed=include_bed,
        threads=threads,
    )
    result = command_runner.run(argv, timeout_seconds=14400)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        tail = detail[-1] if detail else "no diagnostic output"
        raise ValueError(f"modkit pileup failed with exit code {result.returncode}: {tail}")
    if not bedmethyl_path.is_file():
        raise ValueError("modkit pileup reported success but produced no bedMethyl output")

    parameters: dict[str, object] = {
        "subcommand": "pileup",
        "threads": threads,
        "filter_threshold": policy.filter_threshold,
        "cpg_only": policy.cpg_only,
        "combine_strands": policy.combine_strands,
        "modification_codes": [code.value for code in policy.modification_codes],
        "ignored_codes": [code.value for code in policy.ignored_codes],
        "minimum_valid_coverage": policy.minimum_valid_coverage,
        "region_source": policy.region_source.value,
        "include_bed": include_bed.name if include_bed is not None else None,
        "expected_version": policy.expected_version,
    }
    return normalize_methylation(
        sample_id=manifest.sample_id,
        genome_build=manifest.assay.genome_build,
        bedmethyl_path=bedmethyl_path,
        policy=policy,
        tool=ToolRecord(name="modkit", version=version, parameters=parameters),
        regions=regions,
        target_bed=target_bed,
        reads_with_modified_base_tags=tagged_reads,
        warnings=warnings,
    )

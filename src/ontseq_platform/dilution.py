"""In-silico tumour dilution series and the limit-of-detection estimate derived from one.

``docs/BENCHMARKING.md`` puts coverage and tumour/blast-fraction dilution on the third rung
of the validation ladder, and ``docs/ROADMAP.md`` has carried "benchmark CNV candidates
across dilution series" as an open item. This module is the machinery for that rung: mix a
characterised tumour BAM into a normal BAM at declared fractions, run the pipeline on each
level, and turn the resulting benchmark reports into a limit-of-detection statement that
says exactly what it does and does not establish.

Three things are deliberate.

**The plan is pure and separate from the mixing.** Which reads go into which level is
arithmetic over two read counts, and arithmetic can be reviewed, diffed and unit tested
without a BAM. :func:`plan_dilution_series` therefore produces the whole series as data —
seeds, subsample fractions, per-level read budgets — and :func:`execute_dilution_series`
only carries it out. A reviewer can check the design of a titration without running one.

**Nominal is not observed.** Subsampling is random, so a level asked for 5 % tumour reads
does not contain exactly 5 %. Both numbers are recorded per level, and a level that drifts
past the policy's tolerance fails the series rather than being quietly reported at its
nominal fraction — a limit of detection computed from mislabelled levels is worse than no
limit of detection.

**An unbracketed limit is not a limit.** If the lowest fraction tested still detects, the
series has not found where detection stops; it has only shown that the stopping point is
somewhere below the lowest level. :class:`LodReport` says so in a field rather than
reporting the lowest tested fraction as though it were the answer.

Nothing here is a clinical limit of detection. It is a technical characterisation of one
pipeline configuration on one pair of inputs, and it inherits every limitation of the
material it was computed from.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .execution import CommandRunner, StreamingCommandRunner, SubprocessRunner
from .models import (
    BenchmarkKind,
    BenchmarkReport,
    FileFingerprint,
    GenomeBuild,
    ModuleRunStatus,
    StrictModel,
    ToolRecord,
)
from .reference import sha256_file

_VERSION = re.compile(r"(?<!\d)(\d+\.\d+(?:\.\d+)?)(?!\d)")
_IDENTIFIER = r"^[A-Za-z0-9][A-Za-z0-9._-]{2,79}$"

#: samtools expresses "seed and fraction" as one float, so the fraction is emitted with a
#: fixed number of digits. Six is enough to name a 1-in-a-million subsample exactly and
#: keeps the argument reproducible across platforms' float formatting.
_FRACTION_DIGITS = 6

#: Strata keys a benchmark report must carry to take part in a dilution evaluation. They
#: match the vocabulary ``examples/benchmarks/synthetic_cnv.yaml`` already uses, so a case
#: does not need a second, parallel labelling scheme.
TUMOR_FRACTION_KEY = "tumor_fraction"
REPLICATE_KEY = "replicate"
SERIES_KEY = "dilution_series_id"


class DilutionDetection(StrEnum):
    """What one replicate at one level said.

    ``NO_CALL`` is separated from ``NOT_DETECTED`` for the reason it always is here: a
    level whose truth set was empty asked no question, and folding that into "not
    detected" would depress the detection rate with evidence that does not exist.
    """

    DETECTED = "detected"
    NOT_DETECTED = "not_detected"
    NO_CALL = "no_call"


class DilutionPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    expected_samtools_version: str = Field(default="1.24", pattern=r"^\d+\.\d+(?:\.\d+)?$")
    #: Nominal tumour read fractions, highest first. 1.0 is the undiluted control.
    tumor_fractions: list[float] = Field(
        default_factory=lambda: [1.0, 0.5, 0.2, 0.1, 0.05, 0.025], min_length=1
    )
    #: A pure-normal level. Kept separate from ``tumor_fractions`` because it is a negative
    #: control rather than a dilution step, and it must not be scored as one.
    include_normal_only_control: bool = True
    replicates: int = Field(default=1, ge=1)
    #: Base seed for samtools subsampling. Every level derives a distinct seed from it, so
    #: a whole series is reproducible from one integer.
    seed: int = Field(default=42, ge=0, le=999_999)
    #: Reads per mixed level. ``None`` derives the largest budget every level can afford,
    #: which keeps depth constant across the series — the point of a titration is to vary
    #: tumour fraction and nothing else.
    total_read_target: int | None = Field(default=None, ge=1)
    #: How far an observed fraction may drift from its nominal label before the series is
    #: refused. Subsampling is random; mislabelled levels are not acceptable.
    observed_fraction_tolerance: float = Field(default=0.02, gt=0, le=1)
    note: str = Field(min_length=1)

    @model_validator(mode="after")
    def fractions_are_ordered_and_usable(self) -> DilutionPolicy:
        if any(not 0 < value <= 1 for value in self.tumor_fractions):
            raise ValueError("Dilution fractions must be greater than 0 and at most 1")
        if self.tumor_fractions != sorted(set(self.tumor_fractions), reverse=True):
            raise ValueError("Dilution fractions must be unique and ordered highest first")
        for value in self.tumor_fractions:
            if round(value, _FRACTION_DIGITS) != value:
                raise ValueError(
                    f"Dilution fraction {value} cannot be expressed exactly in "
                    f"{_FRACTION_DIGITS} decimal digits, which is what samtools is given"
                )
        return self


class DilutionLevelPlan(StrictModel):
    """One mixed BAM, fully specified before anything runs."""

    level_id: str = Field(pattern=_IDENTIFIER)
    nominal_tumor_fraction: float = Field(ge=0, le=1)
    replicate: int = Field(ge=1)
    is_negative_control: bool = False
    total_reads_target: int = Field(ge=0)
    tumor_reads_target: int = Field(ge=0)
    normal_reads_target: int = Field(ge=0)
    #: Share of each source BAM this level keeps. 0.0 means the source is not used at all
    #: and 1.0 means it is used whole; only strictly-between values are subsampled.
    tumor_subsample_fraction: float = Field(ge=0, le=1)
    normal_subsample_fraction: float = Field(ge=0, le=1)
    #: The exact ``samtools view -s`` argument, or ``None`` when the source is taken whole
    #: or skipped. Precomputed so the plan is auditable without reimplementing the format.
    tumor_subsample_argument: str | None = None
    normal_subsample_argument: str | None = None

    @model_validator(mode="after")
    def level_is_consistent(self) -> DilutionLevelPlan:
        if self.tumor_reads_target + self.normal_reads_target != self.total_reads_target:
            raise ValueError("Dilution level read budget does not sum to its total")
        for fraction, argument, label in (
            (self.tumor_subsample_fraction, self.tumor_subsample_argument, "tumor"),
            (self.normal_subsample_fraction, self.normal_subsample_argument, "normal"),
        ):
            expected = 0 < fraction < 1
            if expected != (argument is not None):
                raise ValueError(
                    f"{label} subsample argument must be present exactly when the source is "
                    "partially sampled"
                )
        return self


class DilutionSeriesPlan(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    series_id: str = Field(pattern=_IDENTIFIER)
    tumor_sample_id: str = Field(min_length=1)
    normal_sample_id: str = Field(min_length=1)
    genome_build: GenomeBuild
    policy: DilutionPolicy
    tumor_read_count: int = Field(ge=1)
    normal_read_count: int = Field(ge=1)
    levels: list[DilutionLevelPlan] = Field(min_length=1)
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True

    @model_validator(mode="after")
    def level_ids_are_unique(self) -> DilutionSeriesPlan:
        ids = [level.level_id for level in self.levels]
        if len(ids) != len(set(ids)):
            raise ValueError("Dilution plan contains duplicate level identifiers")
        return self


class DilutionLevelResult(StrictModel):
    level_id: str = Field(pattern=_IDENTIFIER)
    nominal_tumor_fraction: float = Field(ge=0, le=1)
    replicate: int = Field(ge=1)
    is_negative_control: bool = False
    tumor_reads_observed: int = Field(ge=0)
    normal_reads_observed: int = Field(ge=0)
    total_reads_observed: int = Field(ge=0)
    #: ``None`` when the level contains no reads at all; a fraction of an empty mixture is
    #: not zero, it is undefined.
    observed_tumor_fraction: float | None = Field(default=None, ge=0, le=1)
    mixed_bam_relative_path: str = Field(min_length=1)
    mixed_bam_fingerprint: FileFingerprint

    @model_validator(mode="after")
    def counts_are_consistent(self) -> DilutionLevelResult:
        if self.tumor_reads_observed + self.normal_reads_observed != self.total_reads_observed:
            raise ValueError("Dilution level read counts do not sum to the observed total")
        measured = self.total_reads_observed > 0
        if measured != (self.observed_tumor_fraction is not None):
            raise ValueError(
                "observed_tumor_fraction must be present exactly when the level has reads"
            )
        return self


class DilutionSeriesReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    series_id: str = Field(pattern=_IDENTIFIER)
    status: ModuleRunStatus
    plan: DilutionSeriesPlan
    levels: list[DilutionLevelResult] = Field(min_length=1)
    tool: ToolRecord
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def report_matches_its_plan(self) -> DilutionSeriesReport:
        if self.status != ModuleRunStatus.COMPLETED:
            raise ValueError("A normalized dilution series report must have COMPLETED status")
        if self.tool.version != self.plan.policy.expected_samtools_version:
            raise ValueError("samtools version does not match the dilution policy lock")
        planned = [level.level_id for level in self.plan.levels]
        produced = [level.level_id for level in self.levels]
        if planned != produced:
            raise ValueError("Dilution report levels do not match the plan they came from")
        return self


class LodPolicy(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    profile_id: str = Field(min_length=1)
    status: Literal["technical_defaults_only", "validated"]
    #: Recall a replicate must reach to count as a detection. 1.0 means every truth event
    #: was recovered, which is the strict reading and the default.
    minimum_recall: float = Field(default=1.0, ge=0, le=1)
    #: Share of evaluable replicates that must detect for a level to pass.
    minimum_detection_rate: float = Field(default=0.95, gt=0, le=1)
    #: Levels evaluated on fewer replicates than this cannot pass. One replicate is enough
    #: to run the machinery and nowhere near enough to characterise a limit.
    minimum_replicates: int = Field(default=1, ge=1)
    #: Require every level above the candidate to pass as well, so a single lucky low
    #: level cannot be reported as the limit.
    require_monotonic: bool = True
    note: str = Field(min_length=1)


class LodLevelOutcome(StrictModel):
    nominal_tumor_fraction: float = Field(ge=0, le=1)
    replicates_total: int = Field(ge=1)
    #: Replicates that asked a question at all, i.e. everything but ``NO_CALL``.
    replicates_evaluated: int = Field(ge=0)
    replicates_detected: int = Field(ge=0)
    replicates_no_call: int = Field(ge=0)
    #: ``None`` when nothing was evaluable at this level.
    detection_rate: float | None = Field(default=None, ge=0, le=1)
    mean_recall: float | None = Field(default=None, ge=0, le=1)
    mean_precision: float | None = Field(default=None, ge=0, le=1)
    meets_criterion: bool
    case_ids: list[str] = Field(min_length=1)

    @model_validator(mode="after")
    def outcome_is_consistent(self) -> LodLevelOutcome:
        if self.replicates_evaluated + self.replicates_no_call != self.replicates_total:
            raise ValueError("Dilution level replicate counts do not sum to the total")
        if self.replicates_detected > self.replicates_evaluated:
            raise ValueError("More detections than evaluable replicates")
        if (self.replicates_evaluated > 0) != (self.detection_rate is not None):
            raise ValueError(
                "detection_rate must be present exactly when a replicate was evaluable"
            )
        if self.meets_criterion and self.detection_rate is None:
            raise ValueError("A level with nothing evaluable cannot meet a detection criterion")
        return self


class LodReport(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    series_id: str = Field(pattern=_IDENTIFIER)
    kind: BenchmarkKind
    genome_build: GenomeBuild
    policy: LodPolicy
    levels: list[LodLevelOutcome] = Field(min_length=1)
    #: Lowest tumour fraction meeting the criterion, or ``None`` when no level did.
    detection_limit_fraction: float | None = Field(default=None, ge=0, le=1)
    #: ``True`` only when a *failing* level was observed below the limit. Without one, the
    #: series has bounded the limit from above and not located it.
    bracketed: bool = False
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    research_only: Literal[True] = True
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))

    @model_validator(mode="after")
    def bracketing_requires_a_limit(self) -> LodReport:
        if self.bracketed and self.detection_limit_fraction is None:
            raise ValueError("A series with no detection limit cannot be bracketed")
        fractions = [level.nominal_tumor_fraction for level in self.levels]
        if fractions != sorted(set(fractions), reverse=True):
            raise ValueError("LoD levels must be unique and ordered highest fraction first")
        return self


_SERIES_LIMITATIONS: tuple[str, ...] = (
    "An in-silico mixture is not a wet-lab dilution. It reproduces read-fraction effects "
    "and reproduces nothing about library preparation, input mass or capture behaviour at "
    "low tumour content.",
    "Levels share reads with one another and with the undiluted control, so replicates "
    "within a series are not statistically independent samples.",
    "Subsampling is random. Nominal and observed fractions are both recorded because they "
    "are not the same number.",
)

_LOD_LIMITATIONS: tuple[str, ...] = (
    "This is a technical detection limit for one pipeline configuration on one pair of "
    "inputs. It is not a validated analytical or clinical limit of detection.",
    "The estimate can be no finer than the fractions that were tested; the true limit lies "
    "between the lowest passing level and the highest failing one.",
    "Replicates drawn from one in-silico series are not independent specimens, so the "
    "detection rate is not a confidence statement about future samples.",
)


def _subsample_argument(seed: int, fraction: float) -> str | None:
    """Render samtools' combined seed-and-fraction argument, or ``None`` for whole/none.

    ``samtools view -s 42.050000`` means "seed 42, keep 5 %". A fraction of exactly 1.0
    has no representation in that format — ``.999999`` is not "all reads" — so taking the
    source whole is expressed by not subsampling it at all, and the caller branches on
    ``None`` rather than emitting an argument that would silently drop reads.
    """
    if fraction <= 0 or fraction >= 1:
        return None
    digits = round(fraction * 10**_FRACTION_DIGITS)
    if digits >= 10**_FRACTION_DIGITS:
        # A fraction just below 1.0 rounds up to a seventh digit, and the format has no
        # place for it: "42.1000000" is read back as seed 42 keeping 10 %, so the level
        # would silently end up at a tenth of its intended depth. At this precision the
        # fraction is indistinguishable from taking the source whole, which is what it means.
        return None
    if digits == 0:
        # Rounding down to ".000000" is samtools for "keep nothing", which is not what a
        # positive fraction asked for. Refuse rather than emit an argument that says the
        # opposite of the plan.
        raise ValueError(
            f"Subsample fraction {fraction!r} is smaller than this format can express "
            f"({_FRACTION_DIGITS} decimal places); the level would keep no reads at all"
        )
    return f"{seed}.{digits:0{_FRACTION_DIGITS}d}"


def _level_id(series_id: str, fraction: float, replicate: int) -> str:
    permille = round(fraction * 1000)
    return f"{series_id}_tf{permille:04d}_r{replicate}"


def _feasible_total_reads(
    fractions: Sequence[float], tumor_read_count: int, normal_read_count: int
) -> int:
    """The largest per-level read budget every level in the series can actually fund."""
    budgets: list[float] = []
    for fraction in fractions:
        if fraction > 0:
            budgets.append(tumor_read_count / fraction)
        if fraction < 1:
            budgets.append(normal_read_count / (1 - fraction))
    return int(min(budgets)) if budgets else 0


def plan_dilution_series(
    policy: DilutionPolicy,
    *,
    series_id: str,
    tumor_sample_id: str,
    normal_sample_id: str,
    genome_build: GenomeBuild,
    tumor_read_count: int,
    normal_read_count: int,
) -> DilutionSeriesPlan:
    """Lay out the whole titration as data, without touching a BAM."""
    if tumor_read_count < 1 or normal_read_count < 1:
        raise ValueError("Both source BAMs must contain at least one read")
    if tumor_sample_id == normal_sample_id:
        raise ValueError("Tumour and normal sources must be distinct samples")

    fractions = list(policy.tumor_fractions)
    if policy.include_normal_only_control:
        fractions.append(0.0)

    feasible = _feasible_total_reads(fractions, tumor_read_count, normal_read_count)
    if feasible < 1:
        raise ValueError(
            "The source BAMs cannot fund a single read at every requested fraction; "
            f"tumour has {tumor_read_count} read(s) and normal has {normal_read_count}"
        )
    if policy.total_read_target is None:
        total = feasible
    elif policy.total_read_target > feasible:
        raise ValueError(
            f"total_read_target {policy.total_read_target} exceeds what the sources can fund "
            f"at every level ({feasible}). Lower the target, drop the most extreme fraction, "
            "or sequence deeper"
        )
    else:
        total = policy.total_read_target

    levels: list[DilutionLevelPlan] = []
    for index, fraction in enumerate(fractions):
        for replicate in range(1, policy.replicates + 1):
            # Each source of each level gets its own seed, so no two subsamples of the same
            # BAM anywhere in the series are the same draw. Derived rather than random:
            # the series is reproducible from policy.seed alone.
            tumor_seed = policy.seed + index * 2 * policy.replicates + (replicate - 1) * 2
            normal_seed = tumor_seed + 1
            tumor_reads = round(total * fraction)
            normal_reads = total - tumor_reads
            if tumor_reads > tumor_read_count or normal_reads > normal_read_count:
                # Rounding a budget up can ask for one read more than a source holds. Say
                # so with both numbers rather than letting a subsample fraction above 1.0
                # surface as an opaque schema error.
                raise ValueError(
                    f"level at fraction {fraction} needs {tumor_reads} tumour and "
                    f"{normal_reads} normal read(s), but the sources hold "
                    f"{tumor_read_count} and {normal_read_count}"
                )
            levels.append(
                DilutionLevelPlan(
                    level_id=_level_id(series_id, fraction, replicate),
                    nominal_tumor_fraction=fraction,
                    replicate=replicate,
                    is_negative_control=fraction == 0.0,
                    total_reads_target=total,
                    tumor_reads_target=tumor_reads,
                    normal_reads_target=normal_reads,
                    tumor_subsample_fraction=tumor_reads / tumor_read_count,
                    normal_subsample_fraction=normal_reads / normal_read_count,
                    tumor_subsample_argument=_subsample_argument(
                        tumor_seed, tumor_reads / tumor_read_count
                    ),
                    normal_subsample_argument=_subsample_argument(
                        normal_seed, normal_reads / normal_read_count
                    ),
                )
            )

    warnings = [policy.note]
    if policy.replicates < 3:
        warnings.append(
            f"{policy.replicates} replicate(s) per level cannot characterise a detection "
            "rate; this series demonstrates behaviour rather than measuring a limit."
        )
    if policy.total_read_target is None:
        warnings.append(
            f"Per-level depth was derived from the sources rather than declared: {total} "
            "read(s) per level. Two series built from different inputs are not comparable "
            "unless total_read_target pins this."
        )
    return DilutionSeriesPlan(
        series_id=series_id,
        tumor_sample_id=tumor_sample_id,
        normal_sample_id=normal_sample_id,
        genome_build=genome_build,
        policy=policy,
        tumor_read_count=tumor_read_count,
        normal_read_count=normal_read_count,
        levels=levels,
        warnings=warnings,
        limitations=list(_SERIES_LIMITATIONS),
    )


def samtools_version(text: str) -> str:
    match = _VERSION.search(text)
    if match:
        return match.group(1)
    first_line = text.splitlines()[0].strip() if text.strip() else "unknown"
    return first_line[:80]


def count_reads(
    bam: Path,
    *,
    runner: CommandRunner,
    samtools: str = "samtools",
    threads: int = 4,
) -> int:
    """Count the alignment records in a BAM, refusing to guess when samtools will not say."""
    result = runner.run(
        [samtools, "view", "-c", "-@", str(threads), str(bam)], timeout_seconds=7200
    )
    if result.returncode != 0:
        raise ValueError(f"samtools view -c failed for {bam.name} with {result.returncode}")
    text = result.stdout.strip()
    if not text.isdigit():
        raise ValueError(f"samtools view -c returned no read count for {bam.name}")
    return int(text)


def _run_checked(
    runner: CommandRunner, argv: list[str], *, label: str, timeout_seconds: int = 7200
) -> None:
    result = runner.run(argv, timeout_seconds=timeout_seconds)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        tail = detail[-1] if detail else "no diagnostic output"
        raise ValueError(f"{label} failed with exit code {result.returncode}: {tail}")


def _subsample(
    runner: StreamingCommandRunner,
    *,
    source: Path,
    destination: Path,
    argument: str | None,
    fraction: float,
    samtools: str,
    threads: int,
) -> bool:
    """Write one source's contribution to a level. Returns whether a file was produced."""
    if fraction <= 0:
        return False
    argv = [samtools, "view", "-b", "-@", str(threads)]
    if argument is not None:
        argv.extend(["-s", argument])
    argv.append(str(source))
    result = runner.run_to_file(argv, destination, timeout_seconds=14400)
    if result.returncode != 0:
        detail = (result.stderr or "").strip().splitlines()
        tail = detail[-1] if detail else "no diagnostic output"
        raise ValueError(f"samtools view subsample failed with {result.returncode}: {tail}")
    return True


def execute_dilution_series(
    plan: DilutionSeriesPlan,
    *,
    tumor_bam: Path,
    normal_bam: Path,
    output_dir: Path,
    runner: StreamingCommandRunner | None = None,
    samtools: str = "samtools",
    threads: int = 4,
) -> DilutionSeriesReport:
    """Carry out a planned series, verifying each level against the fraction it claims."""
    if not tumor_bam.is_file():
        raise ValueError("Tumour BAM is missing or unreadable")
    if not normal_bam.is_file():
        raise ValueError("Normal BAM is missing or unreadable")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    output_dir.mkdir(parents=True, exist_ok=True)

    command_runner: StreamingCommandRunner = runner or SubprocessRunner()
    version_result = command_runner.run([samtools, "--version"], timeout_seconds=30)
    if version_result.returncode != 0:
        raise ValueError("samtools version probe returned a non-zero exit code")
    version = samtools_version(f"{version_result.stdout}\n{version_result.stderr}")
    if version != plan.policy.expected_samtools_version:
        raise ValueError(
            f"samtools version {version!r} does not match policy lock "
            f"{plan.policy.expected_samtools_version!r}"
        )

    results: list[DilutionLevelResult] = []
    warnings: list[str] = list(plan.warnings)
    for level in plan.levels:
        mixed = output_dir / f"{level.level_id}.bam"
        if mixed.exists():
            raise ValueError(f"Refusing to overwrite an existing dilution level: {mixed.name}")
        tumor_part = output_dir / f"{level.level_id}.tumor.part.bam"
        normal_part = output_dir / f"{level.level_id}.normal.part.bam"
        parts: list[Path] = []
        if _subsample(
            command_runner,
            source=tumor_bam,
            destination=tumor_part,
            argument=level.tumor_subsample_argument,
            fraction=level.tumor_subsample_fraction,
            samtools=samtools,
            threads=threads,
        ):
            parts.append(tumor_part)
        if _subsample(
            command_runner,
            source=normal_bam,
            destination=normal_part,
            argument=level.normal_subsample_argument,
            fraction=level.normal_subsample_fraction,
            samtools=samtools,
            threads=threads,
        ):
            parts.append(normal_part)
        if not parts:
            raise ValueError(f"Dilution level {level.level_id} would contain no reads at all")

        # Counted before the merge, because the merge is where a source stops being
        # distinguishable: afterwards nothing in the mixture says which read came from
        # which BAM, and the observed fraction is the number this whole lane rests on.
        tumor_reads = (
            count_reads(tumor_part, runner=command_runner, samtools=samtools, threads=threads)
            if tumor_part.is_file()
            else 0
        )
        normal_reads = (
            count_reads(normal_part, runner=command_runner, samtools=samtools, threads=threads)
            if normal_part.is_file()
            else 0
        )
        total_reads = tumor_reads + normal_reads
        observed = tumor_reads / total_reads if total_reads else None

        if len(parts) == 1:
            parts[0].replace(mixed)
        else:
            _run_checked(
                command_runner,
                [samtools, "merge", "-@", str(threads), str(mixed), *(str(p) for p in parts)],
                label="samtools merge",
                timeout_seconds=14400,
            )
        _run_checked(
            command_runner,
            [samtools, "index", "-@", str(threads), str(mixed)],
            label="samtools index",
        )
        if observed is not None:
            drift = abs(observed - level.nominal_tumor_fraction)
            if drift > plan.policy.observed_fraction_tolerance:
                raise ValueError(
                    f"Dilution level {level.level_id} mixed to {observed:.4f} tumour reads "
                    f"against a nominal {level.nominal_tumor_fraction:.4f}; the drift of "
                    f"{drift:.4f} exceeds the policy tolerance of "
                    f"{plan.policy.observed_fraction_tolerance}. Refusing to label a level "
                    "with a fraction it does not contain"
                )
        for part in (tumor_part, normal_part):
            part.unlink(missing_ok=True)
        results.append(
            DilutionLevelResult(
                level_id=level.level_id,
                nominal_tumor_fraction=level.nominal_tumor_fraction,
                replicate=level.replicate,
                is_negative_control=level.is_negative_control,
                tumor_reads_observed=tumor_reads,
                normal_reads_observed=normal_reads,
                total_reads_observed=total_reads,
                observed_tumor_fraction=observed,
                mixed_bam_relative_path=mixed.name,
                mixed_bam_fingerprint=FileFingerprint(
                    size_bytes=mixed.stat().st_size, sha256=sha256_file(mixed)
                ),
            )
        )

    return DilutionSeriesReport(
        series_id=plan.series_id,
        status=ModuleRunStatus.COMPLETED,
        plan=plan,
        levels=results,
        tool=ToolRecord(
            name="samtools",
            version=version,
            parameters={
                "subcommands": ["view", "merge", "index"],
                "threads": threads,
                "seed": plan.policy.seed,
                "level_count": len(results),
                "expected_version": plan.policy.expected_samtools_version,
            },
        ),
        warnings=warnings,
        limitations=list(_SERIES_LIMITATIONS),
    )


def _strata_float(report: BenchmarkReport, key: str) -> float:
    value = report.strata.get(key)
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError(
            f"benchmark case {report.case_id!r} carries no numeric {key!r} stratum, so it "
            "cannot be placed in a dilution series"
        )
    return float(value)


def _strata_replicate(report: BenchmarkReport) -> int:
    value = report.strata.get(REPLICATE_KEY, 1)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(
            f"benchmark case {report.case_id!r} carries a non-integer {REPLICATE_KEY!r} stratum"
        )
    if value < 1:
        raise ValueError(f"benchmark case {report.case_id!r} declares a replicate below 1")
    return value


def classify_detection(report: BenchmarkReport, policy: LodPolicy) -> DilutionDetection:
    """Decide what one benchmark report says about detection at its level."""
    if report.metrics.recall is None:
        return DilutionDetection.NO_CALL
    if report.metrics.recall >= policy.minimum_recall:
        return DilutionDetection.DETECTED
    return DilutionDetection.NOT_DETECTED


def evaluate_lod(
    reports: Sequence[BenchmarkReport],
    policy: LodPolicy,
    *,
    series_id: str,
) -> LodReport:
    """Turn per-level benchmark reports into a detection-limit statement.

    Every report must say which level and replicate it belongs to. Guessing would let two
    runs of the same level be counted as two replicates, which is precisely the error that
    makes a detection rate look better than the evidence supports.
    """
    if not reports:
        raise ValueError("A limit-of-detection evaluation needs at least one benchmark report")
    kinds = {report.kind for report in reports}
    builds = {report.genome_build for report in reports}
    if len(kinds) != 1:
        raise ValueError("A dilution series must compare one benchmark kind")
    if len(builds) != 1:
        raise ValueError("A dilution series must compare one genome build")

    grouped: dict[float, list[BenchmarkReport]] = {}
    seen: set[tuple[float, int]] = set()
    for report in reports:
        declared = report.strata.get(SERIES_KEY)
        if declared is not None and declared != series_id:
            raise ValueError(
                f"benchmark case {report.case_id!r} belongs to series {declared!r}, not "
                f"{series_id!r}"
            )
        fraction = _strata_float(report, TUMOR_FRACTION_KEY)
        if not 0 <= fraction <= 1:
            raise ValueError(f"benchmark case {report.case_id!r} declares an impossible fraction")
        replicate = _strata_replicate(report)
        key = (fraction, replicate)
        if key in seen:
            raise ValueError(
                f"two benchmark reports describe tumour fraction {fraction} replicate "
                f"{replicate}; a repeated level is not an extra replicate"
            )
        seen.add(key)
        grouped.setdefault(fraction, []).append(report)

    outcomes: list[LodLevelOutcome] = []
    for fraction in sorted(grouped, reverse=True):
        level_reports = grouped[fraction]
        detections = [classify_detection(item, policy) for item in level_reports]
        evaluated = [
            item
            for item, verdict in zip(level_reports, detections, strict=True)
            if verdict is not DilutionDetection.NO_CALL
        ]
        detected = sum(1 for verdict in detections if verdict is DilutionDetection.DETECTED)
        no_call = sum(1 for verdict in detections if verdict is DilutionDetection.NO_CALL)
        rate = detected / len(evaluated) if evaluated else None
        recalls = [item.metrics.recall for item in evaluated if item.metrics.recall is not None]
        precisions = [
            item.metrics.precision for item in evaluated if item.metrics.precision is not None
        ]
        meets = (
            rate is not None
            and len(evaluated) >= policy.minimum_replicates
            and rate >= policy.minimum_detection_rate
        )
        outcomes.append(
            LodLevelOutcome(
                nominal_tumor_fraction=fraction,
                replicates_total=len(level_reports),
                replicates_evaluated=len(evaluated),
                replicates_detected=detected,
                replicates_no_call=no_call,
                detection_rate=rate,
                mean_recall=sum(recalls) / len(recalls) if recalls else None,
                mean_precision=sum(precisions) / len(precisions) if precisions else None,
                meets_criterion=meets,
                case_ids=sorted(item.case_id for item in level_reports),
            )
        )

    limit: float | None = None
    if policy.require_monotonic:
        # Walk down from the highest fraction and stop at the first level that fails, so a
        # level that passes only because a lower one happened to work is never reported.
        for outcome in outcomes:
            if not outcome.meets_criterion:
                break
            limit = outcome.nominal_tumor_fraction
    else:
        passing = [item.nominal_tumor_fraction for item in outcomes if item.meets_criterion]
        limit = min(passing) if passing else None

    # A fraction of 0.0 is the pure-normal negative control, not a dilution step: the policy
    # keeps it out of ``tumor_fractions`` for exactly that reason. Measuring "bracketed"
    # against it made every series look bracketed — the control cannot meet a detection
    # criterion it has no truth set for, so it always sat below the limit and always
    # suppressed the "bounded from above, not located" warning that is the whole point.
    dilution_steps = [item for item in outcomes if item.nominal_tumor_fraction > 0]
    lowest = dilution_steps[-1].nominal_tumor_fraction if dilution_steps else None
    bracketed = limit is not None and lowest is not None and limit > lowest
    warnings: list[str] = [policy.note]
    if limit is None:
        warnings.append(
            "No tested fraction met the detection criterion. This series did not establish "
            "a detection limit; it did not establish that detection is impossible either."
        )
    elif not bracketed:
        warnings.append(
            f"The lowest dilution tested ({lowest}) still met the criterion, so the limit is "
            "bounded from above and not located. Extend the series downwards before quoting "
            "this number as a limit."
        )
    if any(item.replicates_no_call for item in outcomes):
        warnings.append(
            "Some replicates produced no evaluable comparison and were excluded from the "
            "detection rate rather than counted as failures."
        )
    if all(item.replicates_total < 3 for item in outcomes):
        warnings.append(
            "Every level was evaluated on fewer than three replicates; the detection rates "
            "here describe what happened, not what a detection rate is."
        )

    return LodReport(
        series_id=series_id,
        kind=next(iter(kinds)),
        genome_build=next(iter(builds)),
        policy=policy,
        levels=outcomes,
        detection_limit_fraction=limit,
        bracketed=bracketed,
        warnings=warnings,
        limitations=list(_LOD_LIMITATIONS),
    )

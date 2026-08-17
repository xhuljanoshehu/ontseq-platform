"""Turning a drop folder into runs, without inventing anything about the samples.

:mod:`ontseq_platform.pipeline.watch` decides *what* to process and *whether it is ready*.
This module does the rest: work out which file in a ready directory is the input, build a
manifest for it, hand it to the pipeline, and record what happened.

The constraint that shapes all of it: **a watch folder must not guess clinical metadata.**
Sample identity, reference, genome build and assay mode are facts about a patient sample,
and a filename is not evidence of any of them. They come from a manifest template the lab
writes once. The only per-sample things derived here are the sample identifier — from the
directory name, and only when that name already satisfies the manifest contract — and the
path of the input file.

Input kind is declared, not sniffed
-----------------------------------

Which kind of input a drop folder holds is *configuration*, not discovery. A lab's drop
folder receives POD5, or unaligned BAM, or aligned BAM; it does not alternate per sample.

Sniffing would mean rules like "a BAM without an index must be unaligned" — wrong the first
time an index has not finished copying, and wrong *silently*: the run would strip and
re-align an already-aligned BAM and produce a plausible result nobody asked for. Declaring
the kind and failing closed when a directory does not match is the version that cannot be
quietly wrong.

Configuration errors fail the watcher, not the samples
------------------------------------------------------

Policies and the manifest template are resolved once per sweep, before any sample is
touched. A template missing ``assay.reference_id`` is one mistake, and recording it as a
failure against every sample in the folder would both misattribute it and require
``--retry-failed`` for all of them after the one-line fix.

Nothing here decides anything biological. A watcher that reaches a `PASS` verdict has
produced evidence for a human to read, exactly as `ontseq run` does when a person types it.
"""

from __future__ import annotations

import signal
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path
from types import FrameType

from . import __version__
from .align import AlignmentPolicy
from .io import load_mapping, load_model
from .models import InputKind, QCPolicy, ReferenceLock, SampleManifest, SnifflesPolicy
from .pipeline.lock import RunAlreadyRunning
from .pipeline.runner import RunConfiguration, run_pipeline
from .pipeline.watch import (
    LEDGER_FILENAME,
    Attempt,
    Candidate,
    Ledger,
    Outcome,
    discover,
    sample_id_from_directory,
    should_attempt,
)

#: What each declared input kind looks like on disk, matched case-insensitively.
INPUT_SUFFIXES: dict[InputKind, tuple[str, ...]] = {
    InputKind.POD5: (".pod5",),
    InputKind.UNALIGNED_BAM: (".bam",),
    InputKind.ALIGNED_BAM: (".bam",),
}

#: Used only to prove a manifest template validates before any real sample is built from it.
PROBE_SAMPLE_ID = "TEMPLATE_PROBE_000"


class DropRejected(ValueError):
    """Raised when a ready directory cannot be turned into a run, and why."""


class WatchConfigurationError(ValueError):
    """Raised when the watcher's own configuration is unusable.

    Distinct from :class:`DropRejected` because it is nobody's sample's fault and must stop
    the watcher rather than be recorded against whichever directory happened to be first.
    """


@dataclass(frozen=True)
class WatchSettings:
    """Everything the watcher needs that does not change between samples."""

    watch_dir: Path
    output_dir: Path
    manifest_template: Path
    reference_lock: Path
    qc_policy: Path
    input_kind: InputKind
    sniffles_policy: Path | None = None
    alignment_policy: Path | None = None
    reference_fasta: Path | None = None
    run_id_prefix: str = ""
    ready_marker: str | None = None
    quiet_seconds: float = 300.0
    threads: int = 4
    git_commit: str = "UNKNOWN"
    retry_failed: bool = False
    executables: dict[str, str] | None = None

    @property
    def ledger_path(self) -> Path:
        return self.output_dir / LEDGER_FILENAME


@dataclass(frozen=True)
class ResolvedConfiguration:
    """Policies loaded and validated once, before any sample is touched."""

    reference_lock: ReferenceLock
    qc_policy: QCPolicy
    sniffles_policy: SnifflesPolicy | None
    alignment_policy: AlignmentPolicy | None
    manifest_template: dict[str, object]


@dataclass(frozen=True)
class PassResult:
    """What one sweep of the watch folder did."""

    attempted: tuple[Attempt, ...]
    skipped: tuple[tuple[str, str], ...]

    @property
    def failures(self) -> tuple[Attempt, ...]:
        return tuple(item for item in self.attempted if item.outcome is Outcome.FAILED)

    @property
    def completed(self) -> tuple[Attempt, ...]:
        return tuple(item for item in self.attempted if item.outcome is Outcome.COMPLETED)


def resolve(settings: WatchSettings) -> ResolvedConfiguration:
    """Load and validate everything shared between samples, or say precisely what is wrong."""
    if not settings.watch_dir.is_dir():
        raise WatchConfigurationError(f"watch folder does not exist: {settings.watch_dir}")
    try:
        reference_lock = load_model(settings.reference_lock, ReferenceLock)
        qc_policy = load_model(settings.qc_policy, QCPolicy)
        sniffles = (
            load_model(settings.sniffles_policy, SnifflesPolicy)
            if settings.sniffles_policy is not None and settings.sniffles_policy.is_file()
            else None
        )
        alignment = (
            load_model(settings.alignment_policy, AlignmentPolicy)
            if settings.alignment_policy is not None and settings.alignment_policy.is_file()
            else None
        )
        template = load_mapping(settings.manifest_template)
    except (OSError, ValueError) as error:
        raise WatchConfigurationError(f"watcher configuration is unusable: {error}") from error

    if settings.input_kind is not InputKind.ALIGNED_BAM and settings.reference_fasta is None:
        raise WatchConfigurationError(
            f"input kind {settings.input_kind.value} has to be aligned, which needs "
            "--reference-fasta"
        )

    # Prove the template produces a valid manifest now, rather than once per sample later.
    try:
        _manifest_from(
            template,
            sample_id=PROBE_SAMPLE_ID,
            kind=settings.input_kind,
            input_path=Path("/probe/input.bam"),
            index_path=Path("/probe/input.bam.bai"),
        )
    except ValueError as error:
        raise WatchConfigurationError(
            f"manifest template {settings.manifest_template} does not produce a valid "
            f"manifest: {error}"
        ) from error

    return ResolvedConfiguration(
        reference_lock=reference_lock,
        qc_policy=qc_policy,
        sniffles_policy=sniffles,
        alignment_policy=alignment,
        manifest_template=template,
    )


def _matching_files(directory: Path, suffixes: Sequence[str]) -> list[Path]:
    wanted = {suffix.lower() for suffix in suffixes}
    return sorted(
        item
        for item in directory.rglob("*")
        if item.is_file() and item.suffix.lower() in wanted and not item.name.startswith(".")
    )


def find_input(directory: Path, kind: InputKind) -> tuple[Path, Path | None]:
    """Locate the declared input inside a drop directory.

    Returns the input path and, for an aligned BAM, its index. Ambiguity is an error rather
    than a choice: picking one of two BAMs would silently analyse an arbitrary half of what
    somebody delivered.
    """
    matches = _matching_files(directory, INPUT_SUFFIXES[kind])
    if not matches:
        wanted = "/".join(INPUT_SUFFIXES[kind])
        raise DropRejected(
            f"declared input kind {kind.value} but no {wanted} file was found in "
            f"{directory.name}"
        )

    if kind is InputKind.POD5:
        # Dorado reads a directory of POD5 files, so many are expected and the common parent
        # is the input. Spread across several directories that parent is ambiguous.
        parents = {item.parent for item in matches}
        if len(parents) != 1:
            raise DropRejected(
                f"{directory.name} holds POD5 files in {len(parents)} directories; "
                "the basecaller takes exactly one"
            )
        return parents.pop(), None

    if len(matches) != 1:
        names = ", ".join(item.name for item in matches)
        raise DropRejected(
            f"{directory.name} holds {len(matches)} BAM files ({names}); exactly one is "
            "required so that no arbitrary choice is made"
        )
    bam = matches[0]

    if kind is InputKind.UNALIGNED_BAM:
        return bam, None

    index = next(
        (
            candidate
            for candidate in (Path(f"{bam}.bai"), bam.with_suffix(".bai"))
            if candidate.is_file()
        ),
        None,
    )
    if index is None:
        raise DropRejected(
            f"{bam.name} was declared as an aligned BAM but has no .bai beside it. "
            "A missing index is not treated as evidence that the BAM is unaligned."
        )
    return bam, index


def _manifest_from(
    template: dict[str, object],
    *,
    sample_id: str,
    kind: InputKind,
    input_path: Path,
    index_path: Path | None,
) -> SampleManifest:
    """Fill a template in with the two things a directory can legitimately supply.

    Everything else — reference, genome build, assay mode, analysis profile — comes from the
    template unchanged, because none of it is derivable from a directory listing.
    """
    payload = dict(template)
    payload["sample_id"] = sample_id
    payload["input"] = {
        "kind": kind.value,
        "path": str(input_path),
        **({"index_path": str(index_path)} if index_path is not None else {}),
    }
    payload.setdefault("run_id", sample_id)
    return SampleManifest.model_validate(payload)


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _attempt_one(
    candidate: Candidate, settings: WatchSettings, resolved: ResolvedConfiguration
) -> Attempt:
    """Run one drop directory, converting every outcome into a recorded attempt."""
    sample_id = sample_id_from_directory(candidate.name)
    if sample_id is None:
        return Attempt(
            name=candidate.name,
            outcome=Outcome.REJECTED,
            detail=(
                f"{candidate.name!r} is not a usable sample identifier (3-64 characters, "
                "starting alphanumeric, then letters, digits, dot, dash or underscore). "
                "Rename the directory; the watcher will not invent an identifier that "
                "would then appear on a reviewer artifact."
            ),
            attempted_at=_now(),
        )

    run_id = f"{settings.run_id_prefix}{sample_id}"
    base = {"name": candidate.name, "attempted_at": _now(), "run_id": run_id}
    try:
        input_path, index_path = find_input(candidate.path, settings.input_kind)
        manifest = _manifest_from(
            resolved.manifest_template,
            sample_id=sample_id,
            kind=settings.input_kind,
            input_path=input_path.resolve(),
            index_path=index_path.resolve() if index_path is not None else None,
        )
    except DropRejected as error:
        return Attempt(**base, outcome=Outcome.REJECTED, detail=str(error), sample_id=sample_id)
    except (OSError, ValueError) as error:
        return Attempt(**base, outcome=Outcome.FAILED, detail=str(error), sample_id=sample_id)

    configuration = RunConfiguration(
        manifest=manifest,
        reference_lock=resolved.reference_lock,
        output_base=settings.output_dir,
        run_id=run_id,
        pipeline_version=__version__,
        git_commit=settings.git_commit,
        qc_policy=resolved.qc_policy,
        sniffles_policy=resolved.sniffles_policy,
        alignment_policy=resolved.alignment_policy,
        reference_fasta=settings.reference_fasta,
        threads=settings.threads,
    )
    if settings.executables:
        configuration = replace(configuration, executables=settings.executables)

    try:
        report, _ = run_pipeline(configuration)
    except RunAlreadyRunning as error:
        return Attempt(**base, outcome=Outcome.LOCKED, detail=str(error), sample_id=sample_id)
    except (OSError, ValueError) as error:
        return Attempt(
            **base,
            outcome=Outcome.FAILED,
            detail=f"the run could not be executed: {error}",
            sample_id=sample_id,
        )

    return Attempt(
        **base,
        outcome=Outcome.COMPLETED if report.passed else Outcome.FAILED,
        detail=report.verdict_reason,
        sample_id=sample_id,
    )


def sweep(
    settings: WatchSettings,
    resolved: ResolvedConfiguration | None = None,
    *,
    now: float | None = None,
) -> PassResult:
    """Make one pass over the watch folder.

    Every candidate is either attempted and recorded, or skipped with a stated reason. The
    ledger is written after each attempt rather than at the end, so a watcher killed between
    two samples does not repeat the one it just finished.
    """
    configuration = resolved if resolved is not None else resolve(settings)
    ledger = Ledger.load(settings.ledger_path)
    candidates = discover(
        settings.watch_dir,
        ready_marker=settings.ready_marker,
        quiet_seconds=settings.quiet_seconds,
        now=now,
    )
    attempted: list[Attempt] = []
    skipped: list[tuple[str, str]] = []
    for candidate in candidates:
        run, why = should_attempt(
            candidate, ledger.get(candidate.name), retry_failed=settings.retry_failed
        )
        if not run:
            skipped.append((candidate.name, why))
            continue
        attempt = _attempt_one(candidate, settings, configuration)
        ledger.record(attempt)
        attempted.append(attempt)
    return PassResult(attempted=tuple(attempted), skipped=tuple(skipped))


class StopRequested:
    """Turns SIGINT/SIGTERM into "stop after the current sample".

    Killing a watcher mid-run leaves a lock for the next start to reclaim and an envelope
    half-built. Finishing the sample in hand costs one sample's runtime and avoids both.
    """

    def __init__(self) -> None:
        self.requested = False
        self._previous: dict[int, object] = {}

    def __enter__(self) -> StopRequested:
        for number in (signal.SIGINT, signal.SIGTERM):
            self._previous[number] = signal.getsignal(number)
            signal.signal(number, self._handle)
        return self

    def __exit__(self, *_: object) -> None:
        for number, handler in self._previous.items():
            signal.signal(number, handler)  # type: ignore[arg-type]

    def _handle(self, signum: int, frame: FrameType | None) -> None:
        self.requested = True


def watch(
    settings: WatchSettings,
    *,
    once: bool = False,
    poll_seconds: float = 60.0,
    report: Callable[[PassResult], None] | None = None,
) -> list[PassResult]:
    """Sweep repeatedly until asked to stop, or exactly once when ``once`` is set.

    Configuration is resolved once up front so a bad policy stops the watcher immediately
    instead of being discovered one sample at a time. ``report`` is called with each
    :class:`PassResult`; this module prints nothing itself.
    """
    resolved = resolve(settings)
    passes: list[PassResult] = []
    with StopRequested() as stop:
        while True:
            result = sweep(settings, resolved)
            passes.append(result)
            if report is not None:
                report(result)
            if once or stop.requested:
                break
            # Slept in short slices so a stop request is noticed promptly rather than after
            # a whole poll interval.
            slept = 0.0
            while slept < poll_seconds and not stop.requested:
                time.sleep(min(1.0, poll_seconds - slept))
                slept += 1.0
            if stop.requested:
                break
    return passes

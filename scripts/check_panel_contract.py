"""Check every committed panel against its lock, before a run is ever attempted.

The pipeline fingerprints the target BED it was handed, so editing a BED re-runs the stage
that used it. Nothing yet reads the lock beside it, which means the facts that decide how
per-target depth may be read -- assembly, buffered or not, confirmed or not -- are recorded
but never checked. This script checks them.

What it proves:

* the lock parses and carries the fields a reader needs to judge it;
* the BED on disk is byte-for-byte the design the lock describes;
* the interval and label counts agree with the lock.

What it then prints, without deciding anything: what the panel does not establish. An
unconfirmed design is a legitimate research input; a report that does not say so is not.

A BED with no lock beside it is reported rather than failed, because synthetic fixtures
legitimately have no provenance to check. ``--strict`` turns every note into a failure,
which is the right setting for an automated gate.

Exit codes:

* ``0`` -- every lock verified against its BED
* ``1`` -- a lock is unusable, or a BED is not the design its lock describes

Research use only. Nothing here validates a panel for clinical use, and no threshold,
coordinate or gene identity is checked or corrected.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from ontseq_platform.pipeline.panel_lock import (
    PanelLockError,
    load_panel_lock,
    panel_usage_warnings,
    target_labels,
    verify_panel_bed,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
PANEL_DIRECTORY = REPO_ROOT / "configs" / "panels"
LOCK_SUFFIX = ".lock.yaml"


def _resolve_bed(lock_path: Path, recorded: str) -> Path:
    """Locate the BED a lock points at, whether it records a repo path or a bare name."""
    candidate = REPO_ROOT / recorded
    if candidate.is_file():
        return candidate
    return lock_path.parent / Path(recorded).name


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify committed panel BEDs against their locks.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--panel-dir",
        type=Path,
        default=PANEL_DIRECTORY,
        help="directory holding the panel BEDs and locks (default: configs/panels)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="fail when a panel is unconfirmed or a BED carries no lock",
    )
    args = parser.parse_args(argv)

    panel_dir: Path = args.panel_dir
    if not panel_dir.is_dir():
        print(f"FAILED  no such panel directory: {panel_dir}")
        return 1

    locks = sorted(panel_dir.glob(f"*{LOCK_SUFFIX}"))
    if not locks:
        print(f"FAILED  no panel lock found in {panel_dir}")
        return 1

    failures: list[str] = []
    notes = 0
    locked_beds: set[Path] = set()

    for lock_path in locks:
        print(f"panel   {lock_path.name}")
        try:
            lock = load_panel_lock(lock_path)
            bed_path = _resolve_bed(lock_path, lock.bed_path)
            verify_panel_bed(lock, bed_path)
            labels = target_labels(bed_path)
        except PanelLockError as error:
            failures.append(f"{lock_path.name}: {error}")
            print(f"  FAILED  {error}")
            continue
        locked_beds.add(bed_path.resolve())
        print(
            f"  ok      {bed_path.name} matches the lock: {lock.target_count} interval(s), "
            f"{len(set(labels))} distinct label(s), build {lock.genome_build}"
        )
        for warning in panel_usage_warnings(lock, labels=labels):
            notes += 1
            print(f"  note    {warning}")

    for bed_path in sorted(panel_dir.glob("*.bed")):
        if bed_path.resolve() in locked_beds:
            continue
        notes += 1
        print(f"panel   {bed_path.name}")
        print(
            "  note    no lock beside this BED. A manifest pointing at it carries no "
            "provenance that can be checked; acceptable for a synthetic fixture, not for a "
            "design a run is measured against."
        )

    print()
    if failures:
        print(f"FAILED  {len(failures)} panel contract violation(s)")
        return 1
    if notes and args.strict:
        print(f"FAILED  {notes} note(s) and --strict was requested")
        return 1
    print(f"ok      {len(locks)} lock(s) verified, {notes} note(s) recorded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

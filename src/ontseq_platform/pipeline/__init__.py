"""End-to-end pipeline execution.

``stages``
    The declarative stage graph: applicability, dependencies, failure propagation and the
    verification status of each adapter. Dependency-free and independently testable.
``envelope``
    The run directory, atomic writes, artifact fingerprints and resume signatures.
    Dependency-free.
``state``
    Versioned contracts for the run report and the release bundle.
``runner``
    Sequences the existing adapters, records every outcome and resumes safely.
"""

from __future__ import annotations

__all__ = ["envelope", "runner", "stages", "state"]

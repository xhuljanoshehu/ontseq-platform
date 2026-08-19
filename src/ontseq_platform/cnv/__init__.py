"""Copy-number variation benchmarking, truth representation and baseline calling.

The subsystem is layered so that the scientific core can be tested and reasoned about
without the serialization contract:

``intervals``, ``states``, ``stats``, ``core``
    Dependency-free numerical and semantic core. No pydantic, no I/O.
``mask``
    Construction of the evaluable genome from reference, exclusions and coverage.
``models``
    Versioned pydantic contracts for truth sets, call sets and evaluation reports.
``truth``
    Multi-source truth ingestion, including cytogenetic band uncertainty.
``evaluate``
    Orchestration between the contract layer and :mod:`core`.
``simulate``
    Deterministic synthetic bin counts for dilution and coverage series.
``segment``
    A transparent read-depth baseline caller used as an experimental control.
``adapters``
    Parsers that normalize third-party caller output into the shared contract.
``lea_compat`` / ``lea_truth_tables``
    Strict research-only compatibility boundary for the Lea Evers historical comparator.
``strata``
    Cross-run aggregation and limit-of-detection analysis.
"""

from __future__ import annotations

__all__ = [
    "adapters",
    "core",
    "evaluate",
    "intervals",
    "lea_compat",
    "lea_truth_tables",
    "mask",
    "models",
    "segment",
    "simulate",
    "states",
    "stats",
    "strata",
    "truth",
]

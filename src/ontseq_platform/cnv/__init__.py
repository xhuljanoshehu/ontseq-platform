"""Copy-number analysis: the runtime lane, and the benchmark subsystem that judges it.

Two things live here, and keeping them apart is deliberate.

The **runtime lane** is what a patient run executes. ``qdnaseq`` drives QDNAseq +
ACE, and ``extension`` installs it into the execution graph. It produces the
copy-number evidence a physician eventually reads.

The **benchmark subsystem** never runs on a patient sample. It exists to answer
"how good is a caller, on what, and where can it not answer at all" before any
lane is allowed near a report. It is layered so the scientific core can be tested
and reasoned about without the serialization contract:

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
``strata``
    Cross-run aggregation and limit-of-detection analysis.

A caller being measured well here is not a promotion. Nothing in this package
makes a lane clinically reportable; that decision needs real cohort validation
and a human who signs for it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover - import-time typing only
    from .qdnaseq import (
        CnvChromosomeConsensus,
        CnvFit,
        QDNAseqCallReport,
        QDNAseqPolicy,
        run_qdnaseq_ace,
    )

# The runtime names stay reachable from the package root, but they are resolved on first
# access rather than at import. Binding them eagerly would make `cnv.core`, `cnv.intervals`,
# `cnv.states` and `cnv.stats` import pydantic through this file — and those four modules
# are dependency-free on purpose, which is a property that only holds if nothing above them
# quietly breaks it.
_RUNTIME_EXPORTS = frozenset(
    {
        "CnvChromosomeConsensus",
        "CnvFit",
        "QDNAseqCallReport",
        "QDNAseqPolicy",
        "run_qdnaseq_ace",
    }
)


def __getattr__(name: str) -> Any:
    if name in _RUNTIME_EXPORTS:
        from . import qdnaseq

        return getattr(qdnaseq, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(__all__)


__all__ = [
    "CnvChromosomeConsensus",
    "CnvFit",
    "QDNAseqCallReport",
    "QDNAseqPolicy",
    "adapters",
    "core",
    "evaluate",
    "intervals",
    "mask",
    "models",
    "run_qdnaseq_ace",
    "segment",
    "simulate",
    "states",
    "stats",
    "strata",
    "truth",
]

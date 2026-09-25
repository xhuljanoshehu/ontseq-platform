"""Research-only methylation lanes that go beyond region-aggregated fractions.

Each lane is a separate, fail-closed contract with its own versioned policy and report
(issue #96). None of them changes the aggregate methylation lane
(:mod:`ontseq_platform.methylation`), its thresholds or the reporting boundary.

``haplotype``
    Haplotype-resolved aggregates from a haplotagged BAM via the qualified
    ``modkit pileup --phased`` layout, bounded by phase-block identity (issue #97).
"""

from __future__ import annotations

__all__ = ["haplotype"]

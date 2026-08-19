from __future__ import annotations

from pathlib import Path

from .breakends import breakend_descriptors_from_sniffles_vcf
from .fusion import (
    FusionInterpretationReport,
    GeneAnnotationIndex,
    ObservabilityRegion,
    interpret_sniffles_fusions,
)
from .models import SnifflesCallReport


def interpret_sniffles_vcf_fusions(
    report: SnifflesCallReport,
    sniffles_vcf: Path,
    annotation: GeneAnnotationIndex,
    *,
    observability: list[ObservabilityRegion] | None = None,
    known_pairs: set[tuple[str, str]] | None = None,
    flank_bp: int = 0,
) -> FusionInterpretationReport:
    """Join a normalized Sniffles report with privacy-safe BND orientation evidence.

    The VCF is read locally only to derive `BreakendDescriptor` objects. Raw VCF IDs,
    read names, REF/ALT strings, inserted sequence and the source path are not propagated
    into the returned fusion report. Event reconciliation remains exact and fail-closed in
    `FusionCandidate`.

    This helper does not normalize the VCF itself. The caller must supply the
    `SnifflesCallReport` produced from the same VCF and caller run.
    """

    descriptors = breakend_descriptors_from_sniffles_vcf(sniffles_vcf)
    return interpret_sniffles_fusions(
        report,
        annotation,
        observability=observability,
        known_pairs=known_pairs,
        breakend_descriptors=descriptors,
        flank_bp=flank_bp,
    )

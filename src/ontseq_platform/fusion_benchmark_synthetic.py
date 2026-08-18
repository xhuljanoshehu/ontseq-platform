from __future__ import annotations

from pathlib import Path

from .fusion import GeneAnnotationIndex, ObservabilityRegion, ObservabilityStatus
from .fusion_benchmark import (
    FusionBenchmarkCaseResult,
    PartnerExpectation,
    SyntheticFusionFixture,
    evaluate_fusion_benchmark_case,
)
from .fusion_workflow import interpret_sniffles_vcf_fusions
from .models import EventType, GenomeBuild, SnifflesPolicy, StrictModel, ToolRecord
from .sniffles import normalize_sniffles_vcf

_PRIMARY_CHROMOSOME = "chr1"
_SECONDARY_CHROMOSOME = "chr2"
_DUMMY_CHROMOSOME = "chr3"
_PRIMARY_POSITION_0BASED = 100_000
_SECONDARY_POSITION_0BASED = 200_000
_SYNTHETIC_PARTNER = "SYNTHPARTNER"
_SYNTHETIC_GENE = "SYNTHGENE"
_DUMMY_GENE = "SYNTHDUMMY"


class SyntheticFusionExecutionResult(StrictModel):
    fixture_id: str
    benchmark: FusionBenchmarkCaseResult
    normalized_event_count: int
    fusion_candidate_count: int
    privacy_profile: str = "synthetic_nonbiological_local_files_only"


def _genes_for_fixture(fixture: SyntheticFusionFixture) -> tuple[str | None, str | None]:
    if fixture.partner_expectation == PartnerExpectation.EXACT_PAIR:
        assert fixture.expected_gene_pair is not None
        return fixture.expected_gene_pair
    if fixture.partner_expectation == PartnerExpectation.ANY_PARTNER:
        assert fixture.anchor_gene is not None
        return fixture.anchor_gene, _SYNTHETIC_PARTNER
    if fixture.expected_classification.value == "gene_intergenic":
        return _SYNTHETIC_GENE, None
    return None, None


def _write_synthetic_vcf(path: Path, fixture: SyntheticFusionFixture) -> None:
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    if fixture.expected_candidate_count == 0:
        path.write_text(header, encoding="utf-8")
        return
    if fixture.expected_candidate_count != 1:
        raise ValueError("synthetic fusion harness currently supports zero or one candidate")
    primary_position_1based = _PRIMARY_POSITION_0BASED + 1
    secondary_position_1based = _SECONDARY_POSITION_0BASED + 1
    record = (
        f"{_PRIMARY_CHROMOSOME}\t{primary_position_1based}\tSYNTHETIC_BND\tN\t"
        f"N]{_SECONDARY_CHROMOSOME}:{secondary_position_1based}]\t60\tPASS\t"
        "PRECISE;SVTYPE=BND;SUPPORT=12;STRAND=+-\n"
    )
    path.write_text(header + record, encoding="utf-8")


def _bed_feature(chromosome: str, position: int, gene: str, strand: str) -> str:
    return (
        f"{chromosome}\t{position - 50}\t{position + 51}\t{gene}\t0\t{strand}\tENST_SYNTH_{gene}\n"
    )


def _write_synthetic_annotation(
    path: Path,
    fixture: SyntheticFusionFixture,
) -> set[tuple[str, str]]:
    primary_gene, secondary_gene = _genes_for_fixture(fixture)
    rows: list[str] = []
    if primary_gene is not None:
        rows.append(
            _bed_feature(
                _PRIMARY_CHROMOSOME,
                _PRIMARY_POSITION_0BASED,
                primary_gene,
                "+",
            )
        )
    if secondary_gene is not None:
        rows.append(
            _bed_feature(
                _SECONDARY_CHROMOSOME,
                _SECONDARY_POSITION_0BASED,
                secondary_gene,
                "-",
            )
        )
    if not rows:
        rows.append(_bed_feature(_DUMMY_CHROMOSOME, 300_000, _DUMMY_GENE, "+"))
    path.write_text("".join(rows), encoding="utf-8")

    if fixture.partner_expectation == PartnerExpectation.EXACT_PAIR:
        assert fixture.expected_gene_pair is not None
        return {fixture.expected_gene_pair}
    return set()


def _observability_region(
    chromosome: str,
    position: int,
    status: ObservabilityStatus | None,
) -> ObservabilityRegion | None:
    if status is None or status == ObservabilityStatus.UNKNOWN:
        return None
    return ObservabilityRegion(
        chromosome=chromosome,
        start=position - 100,
        end=position + 101,
        status=status,
        reason="deterministic synthetic benchmark observability",
    )


def _observability_for_fixture(fixture: SyntheticFusionFixture) -> list[ObservabilityRegion]:
    regions = [
        _observability_region(
            _PRIMARY_CHROMOSOME,
            _PRIMARY_POSITION_0BASED,
            fixture.expected_primary_observability,
        ),
        _observability_region(
            _SECONDARY_CHROMOSOME,
            _SECONDARY_POSITION_0BASED,
            fixture.expected_secondary_observability,
        ),
    ]
    return [region for region in regions if region is not None]


def execute_synthetic_fusion_fixture(
    fixture: SyntheticFusionFixture,
    workdir: Path,
) -> SyntheticFusionExecutionResult:
    """Execute one fixture through the real local VCF-to-fusion software path.

    All generated coordinates, genes without a family label, VCF IDs, transcripts and file
    contents are deterministic synthetic test material. The function never asserts clinical
    validity and does not invoke Sniffles2 itself; it exercises the production VCF normalizer,
    privacy-safe BND parser, annotation, observability, fusion interpretation and benchmark
    evaluator using a synthetic Sniffles-compatible VCF.
    """

    workdir.mkdir(parents=True, exist_ok=True)
    vcf_path = workdir / f"{fixture.fixture_id}.synthetic.vcf"
    annotation_path = workdir / f"{fixture.fixture_id}.synthetic.genes.bed"
    _write_synthetic_vcf(vcf_path, fixture)
    known_pairs = _write_synthetic_annotation(annotation_path, fixture)

    policy = SnifflesPolicy(
        profile_id="synthetic_fusion_benchmark_v1",
        status="technical_defaults_only",
        expected_version="2.8.0",
        min_support=1,
        allowed_sv_types=[EventType.TRANSLOCATION],
        note="Synthetic software benchmark only; not an analytical or clinical policy.",
    )
    tool = ToolRecord(name="Sniffles2", version="2.8.0")
    normalized = normalize_sniffles_vcf(
        vcf_path,
        sample_id=f"SYNTHETIC_{fixture.fixture_id.upper()}",
        genome_build=GenomeBuild.GRCH38,
        policy=policy,
        tool=tool,
    )
    annotation = GeneAnnotationIndex.from_bed(
        annotation_path,
        resource_id="synthetic_fusion_benchmark_genes",
        resource_version="v1",
        genome_build=GenomeBuild.GRCH38,
    )
    fusion_report = interpret_sniffles_vcf_fusions(
        normalized,
        vcf_path,
        annotation,
        observability=_observability_for_fixture(fixture),
        known_pairs=known_pairs,
    )
    benchmark = evaluate_fusion_benchmark_case(fixture, fusion_report)
    return SyntheticFusionExecutionResult(
        fixture_id=fixture.fixture_id,
        benchmark=benchmark,
        normalized_event_count=len(normalized.events),
        fusion_candidate_count=len(fusion_report.candidates),
    )

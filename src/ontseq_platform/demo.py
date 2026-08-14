from __future__ import annotations

from .iscn import build_iscn_proposal
from .models import (
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    EventType,
    Evidence,
    GenomeBuild,
    GenomicEvent,
    InputKind,
    InputSpec,
    Locus,
    ModuleOutcome,
    ModuleRunStatus,
    PipelineResult,
    Provenance,
    QCMetrics,
    SampleManifest,
    ToolRecord,
    Verdict,
)


def build_demo_result() -> PipelineResult:
    manifest = SampleManifest(
        sample_id="SYNTHETIC_AML_001",
        run_id="DEMO_RUN_001",
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path="/secure-input/SYNTHETIC_AML_001.bam",
            index_path="/secure-input/SYNTHETIC_AML_001.bam.bai",
        ),
        assay=AssaySpec(
            mode=AssayMode.ADAPTIVE_SAMPLING,
            genome_build=GenomeBuild.GRCH38,
            reference_id="GRCh38-demo-not-for-analysis",
            target_bed="configs/panels/aml_fusions.synthetic.bed",
            target_bed_version="synthetic-v1",
        ),
        analysis=AnalysisSpec(
            profile="adaptive-sampling-demo",
            modules=[
                AnalysisModule.QC,
                AnalysisModule.CNV,
                AnalysisModule.SV,
                AnalysisModule.FUSION,
                AnalysisModule.ISCN,
                AnalysisModule.REPORT,
            ],
        ),
    )
    events = [
        GenomicEvent(
            event_id="CNV-001",
            event_type=EventType.CHROMOSOME_GAIN,
            primary=Locus(chromosome="chr8", start=0, end=145_138_636),
            copy_number=3,
            confidence="moderate",
            reportable=True,
            evidence=[
                Evidence(
                    caller="QDNAseq+ACE",
                    caller_version="demo",
                    local_coverage=3.2,
                    filters=["SYNTHETIC"],
                )
            ],
        ),
        GenomicEvent(
            event_id="CNV-002",
            event_type=EventType.DELETION,
            primary=Locus(
                chromosome="chr5",
                start=70_000_000,
                end=160_000_000,
                cytoband_start="q13",
                cytoband_end="q34",
            ),
            copy_number=1,
            confidence="moderate",
            reportable=True,
            evidence=[
                Evidence(
                    caller="QDNAseq+ACE",
                    caller_version="demo",
                    local_coverage=3.0,
                    filters=["SYNTHETIC"],
                )
            ],
        ),
        GenomicEvent(
            event_id="FUS-001",
            event_type=EventType.FUSION,
            primary=Locus(
                chromosome="chr8",
                start=92_000_000,
                end=92_000_100,
                cytoband_start="q22",
                gene="RUNX1T1",
            ),
            secondary=Locus(
                chromosome="chr21",
                start=36_000_000,
                end=36_000_100,
                cytoband_start="q22.1",
                gene="RUNX1",
            ),
            genes=["RUNX1", "RUNX1T1"],
            confidence="high",
            reportable=False,
            evidence=[
                Evidence(
                    caller="Sniffles2",
                    caller_version="demo",
                    support_reads=12,
                    local_coverage=24,
                    variant_allele_fraction=0.50,
                    filters=["SYNTHETIC"],
                ),
                Evidence(
                    caller="cuteSV",
                    caller_version="demo",
                    support_reads=11,
                    local_coverage=24,
                    variant_allele_fraction=0.46,
                    filters=["SYNTHETIC"],
                ),
            ],
            notes=["Illustrative fusion only; coordinates are synthetic."],
        ),
    ]
    proposal = build_iscn_proposal(events, chromosome_count=46, sex_chromosomes="XX")
    return PipelineResult(
        manifest=manifest,
        qc=QCMetrics(
            verdict=Verdict.PASS,
            metrics={
                "number_of_reads": 626_650,
                "total_yield_gb": 9.41,
                "n50_bp": 22_401,
                "median_identity_percent": 93.3,
                "aligned_percent": 73.7,
                "mean_coverage_x": 3.04,
                "target_coverage_x": 24.0,
            },
            warnings=["Synthetic QC metrics; no biological interpretation is permitted."],
        ),
        events=events,
        iscn=proposal,
        provenance=Provenance(
            pipeline_version="0.3.0",
            git_commit="UNCOMMITTED-DEMO",
            tools=[
                ToolRecord(name="Cramino", version="demo"),
                ToolRecord(
                    name="QDNAseq",
                    version="demo",
                    parameters={"bins_kb": [100, 500, 1000]},
                ),
                ToolRecord(name="ACE", version="demo", parameters={"penalty": 0.6}),
                ToolRecord(name="Sniffles2", version="demo"),
                ToolRecord(name="cuteSV", version="demo"),
            ],
            reference_checksums={"reference": "synthetic-not-a-real-checksum"},
        ),
        modules=[
            ModuleOutcome(
                module=module,
                status=ModuleRunStatus.COMPLETED,
                reason="Synthetic demonstration only; no scientific tool was executed",
            )
            for module in manifest.analysis.modules
        ],
        warnings=[
            "All data in this report are synthetic.",
            "The ISCN renderer implements only an unvalidated subset of ISCN 2024.",
            "No result may be used for diagnosis or treatment decisions.",
        ],
    )

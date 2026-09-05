"""One assemble stage, which extensions may add to but never replace.

The CNV lane used to arrive by registration and *replace* the assemble stage, leaving two
implementations of one contract. The failure mode was silent: whichever report the copy
forgot simply did not reach the result, while the artifact it came from still sat in the
envelope looking as though it was used. The SV consensus went missing exactly that way —
it landed after the extension was written, the runner's assemble was updated to pass it and
the replacement was not, so every run with CNV registered assembled reviewer artifacts from
raw Sniffles calls.

Extensions now contribute instead: they receive the finished result and return an enriched
one. There is only one assemble body, so there is nothing left to forget.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cnv import extension as cnv_extension
from ontseq_platform.cnv.qdnaseq import QDNAseqPolicy
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CraminoQCReport,
    EventType,
    GenomeBuild,
    GenomicEvent,
    InputKind,
    InputSpec,
    Locus,
    ModuleRunStatus,
    PipelineResult,
    QCMetrics,
    QCPolicy,
    ReferenceContig,
    ReferenceLock,
    SampleManifest,
    SnifflesCallReport,
    SnifflesPolicy,
    SvConsensusPolicy,
    SvConsensusReport,
    ToolRecord,
    Verdict,
)
from ontseq_platform.pipeline import runner as pipeline_runner
from ontseq_platform.pipeline.envelope import RunEnvelope
from ontseq_platform.pipeline.runner import (
    ASSEMBLE_SOURCE_ARTIFACTS,
    RunConfiguration,
    RunContext,
    StagePlan,
    load_assemble_inputs,
)
from ontseq_platform.pipeline.stages import (
    SPEC_BY_STAGE,
    StageId,
    VerificationStatus,
    verification_of,
)

SAMPLE = "ASSEMBLE_001"

#: Ids chosen so a result says which report it was built from.
SNIFFLES_EVENT_ID = "sniffles-raw-call"
CONSENSUS_EVENT_ID = "consensus-annotated-call"


class _NullRunner:
    def run(self, argv, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        raise AssertionError("assemble is pure Python and must run no tool")


def _event(event_id: str) -> GenomicEvent:
    return GenomicEvent(
        event_id=event_id,
        event_type=EventType.DELETION,
        primary=Locus(chromosome="chr1", start=1000, end=5000),
    )


def _manifest() -> SampleManifest:
    return SampleManifest(
        sample_id=SAMPLE,
        run_id="ASSEMBLE_RUN_001",
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path="/nowhere/input.bam",
            index_path="/nowhere/input.bam.bai",
        ),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id="ASSEMBLE_REFERENCE_V1",
        ),
        analysis=AnalysisSpec(profile="assemble", modules=[AnalysisModule.QC, AnalysisModule.SV]),
    )


def _intake() -> AlignedBamIntakeReport:
    return AlignedBamIntakeReport(
        sample_id=SAMPLE,
        reference_id="ASSEMBLE_REFERENCE_V1",
        genome_build=GenomeBuild.GRCH38,
        checks=[],
        verdict=Verdict.PASS,
    )


def _qc() -> CraminoQCReport:
    return CraminoQCReport(
        sample_id=SAMPLE,
        qc=QCMetrics(verdict=Verdict.PASS, metrics={"reads": 100}),
        tool=ToolRecord(name="cramino", version="1.3.0"),
    )


def _sniffles() -> SnifflesCallReport:
    from ontseq_platform.models import FileFingerprint

    return SnifflesCallReport(
        sample_id=SAMPLE,
        genome_build=GenomeBuild.GRCH38,
        status=ModuleRunStatus.COMPLETED,
        policy=SnifflesPolicy(profile_id="assemble", status="technical_defaults_only", note="test"),
        events=[_event(SNIFFLES_EVENT_ID)],
        raw_record_count=1,
        accepted_record_count=1,
        rejected_record_count=0,
        tool=ToolRecord(name="sniffles", version="2.8.0"),
        vcf_fingerprint=FileFingerprint(size_bytes=1),
    )


def _consensus() -> SvConsensusReport:
    return SvConsensusReport(
        sample_id=SAMPLE,
        genome_build=GenomeBuild.GRCH38,
        status=ModuleRunStatus.COMPLETED,
        policy=SvConsensusPolicy(
            profile_id="assemble", status="technical_defaults_only", note="test"
        ),
        events=[_event(CONSENSUS_EVENT_ID)],
        input_event_count=2,
        consolidated_event_count=1,
        caller_names=["sniffles", "cutesv"],
    )


class AssembleCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.base = Path(self._tmp.name)
        self.manifest = _manifest()
        config = RunConfiguration(
            manifest=self.manifest,
            reference_lock=ReferenceLock(
                reference_id="ASSEMBLE_REFERENCE_V1",
                genome_build=GenomeBuild.GRCH38,
                contigs=[ReferenceContig(name="chr1", length=100000)],
                source_fai_sha256="a" * 64,
            ),
            output_base=self.base / "runs",
            run_id="ASSEMBLE_RUN_001",
            pipeline_version="0.0.0-test",
            git_commit="0" * 40,
            qc_policy=QCPolicy(status="technical_defaults_only", note="test"),
        )
        envelope = RunEnvelope.create(config.output_base, run_id=config.run_id, sample_id=SAMPLE)
        self.ctx = RunContext(
            config=config, envelope=envelope, runner=_NullRunner(), manifest=self.manifest
        )
        self._write(pipeline_runner.INTAKE_REPORT, _intake().model_dump_json(indent=2))
        self._write(pipeline_runner.QC_REPORT, _qc().model_dump_json(indent=2))
        self._write(self.ctx.path(pipeline_runner.SV_REPORT), _sniffles().model_dump_json(indent=2))
        self._write(
            self.ctx.path(pipeline_runner.SV_CONSENSUS_REPORT),
            _consensus().model_dump_json(indent=2),
        )

    def _write(self, relative: str, content: str) -> None:
        path = self.ctx.envelope.path(relative)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content + "\n", encoding="utf-8")

    def _result_event_ids(self) -> list[str]:
        payload = self.ctx.envelope.path(self.ctx.path(pipeline_runner.RESULT_JSON))
        result = PipelineResult.model_validate_json(payload.read_text(encoding="utf-8"))
        return [item.event_id for item in result.events]

    def _plan(self) -> StagePlan:
        return StagePlan(parameters={}, tool_versions={})


class _CnvRegistration:
    """Register the CNV lane for one test and undo it afterwards.

    Registration is process-global by design — an extension is installed for the life of a
    process — so a test that installs one has to put the runner back as it found it.
    """

    def __enter__(self) -> _CnvRegistration:
        # CNV has no entry until something registers one: the graph declares the stage,
        # nothing implements it by default. Restoring means removing, not overwriting.
        self.implementation = pipeline_runner.IMPLEMENTATIONS.get(StageId.CNV)
        self.result_contributions = list(pipeline_runner.RESULT_CONTRIBUTIONS)
        self.report_contributions = list(pipeline_runner.REPORT_CONTRIBUTIONS)
        cnv_extension.register_qdnaseq_extension(
            cnv_extension.QDNAseqExtensionSettings(
                policy=QDNAseqPolicy(profile_id="assemble-test", note="test")
            )
        )
        return self

    def __exit__(self, *exc: object) -> None:
        if self.implementation is None:
            pipeline_runner.IMPLEMENTATIONS.pop(StageId.CNV, None)
        else:
            pipeline_runner.IMPLEMENTATIONS[StageId.CNV] = self.implementation
        pipeline_runner.RESULT_CONTRIBUTIONS[:] = self.result_contributions
        pipeline_runner.REPORT_CONTRIBUTIONS[:] = self.report_contributions


class ConsensusReachesTheResultTests(AssembleCase):
    def test_the_result_is_assembled_from_the_consensus(self) -> None:
        pipeline_runner._assemble_execute(self.ctx, self._plan())
        self.assertEqual(self._result_event_ids(), [CONSENSUS_EVENT_ID])

    def test_the_consensus_survives_a_registered_extension(self) -> None:
        """The regression: with CNV registered, the result used to lose the consensus.

        Assembling from `sniffles.events` instead dropped cuteSV entirely along with
        consensus merging, gene and cytoband annotation, artifact context flags, Adaptive
        Sampling observability and AML prioritization — invisibly, because the consensus
        artifact was still written to the envelope.
        """
        with _CnvRegistration():
            pipeline_runner._assemble_execute(self.ctx, self._plan())
        self.assertEqual(self._result_event_ids(), [CONSENSUS_EVENT_ID])

    def test_without_a_consensus_the_result_falls_back_to_sniffles(self) -> None:
        self.ctx.envelope.path(self.ctx.path(pipeline_runner.SV_CONSENSUS_REPORT)).unlink()
        pipeline_runner._assemble_execute(self.ctx, self._plan())
        self.assertEqual(self._result_event_ids(), [SNIFFLES_EVENT_ID])


class ExtensionsMayNotReplaceTests(AssembleCase):
    """An extension owns its own stage and contributes to everyone else's."""

    def test_registration_leaves_the_assemble_and_report_stages_alone(self) -> None:
        before = (
            pipeline_runner.IMPLEMENTATIONS[StageId.ASSEMBLE],
            pipeline_runner.IMPLEMENTATIONS[StageId.REPORT],
        )
        with _CnvRegistration():
            after = (
                pipeline_runner.IMPLEMENTATIONS[StageId.ASSEMBLE],
                pipeline_runner.IMPLEMENTATIONS[StageId.REPORT],
            )
            self.assertEqual(before, after)
            self.assertEqual(
                [item.extension_id for item in pipeline_runner.RESULT_CONTRIBUTIONS],
                [cnv_extension.EXTENSION_ID],
            )

    def test_registration_does_not_mutate_the_stage_graph(self) -> None:
        """The graph is data. A registration that edits it makes honesty order-dependent."""
        before = SPEC_BY_STAGE[StageId.CNV]
        with _CnvRegistration():
            self.assertIs(SPEC_BY_STAGE[StageId.CNV], before)
            self.assertIs(
                SPEC_BY_STAGE[StageId.CNV].verification,
                VerificationStatus.NOT_IMPLEMENTED,
            )

    def test_the_implementation_carries_its_own_verification(self) -> None:
        with _CnvRegistration():
            self.assertIs(
                pipeline_runner.implementation_verification()[StageId.CNV],
                VerificationStatus.VERIFIED_WITH_REAL_TOOL,
            )
            self.assertIs(
                verification_of(StageId.CNV, pipeline_runner.implementation_verification()),
                VerificationStatus.VERIFIED_WITH_REAL_TOOL,
            )
        self.assertNotIn(StageId.CNV, pipeline_runner.implementation_verification())

    def test_registering_twice_does_not_duplicate_a_contribution(self) -> None:
        with _CnvRegistration(), _CnvRegistration():
            ids = [item.extension_id for item in pipeline_runner.RESULT_CONTRIBUTIONS]
            self.assertEqual(ids, [cnv_extension.EXTENSION_ID])


class LoaderTests(AssembleCase):
    def test_the_loader_returns_every_evidence_report_present(self) -> None:
        inputs = load_assemble_inputs(self.ctx)
        self.assertIsNotNone(inputs.sniffles)
        self.assertIsNotNone(inputs.sv_consensus)
        self.assertIsNone(inputs.methylation)
        assert inputs.sv_consensus is not None
        self.assertEqual(
            [item.event_id for item in inputs.sv_consensus.events], [CONSENSUS_EVENT_ID]
        )


class ResumeSignatureTests(AssembleCase):
    """Assemble declares only QC as a dependency, so upstream artifacts do not cover it.

    `stage_signature` hashes the artifacts of a stage's *declared* dependencies. Assemble
    reads the SV, consensus and methylation reports without depending on those stages, so
    each one has to be fingerprinted explicitly or a changed report resumes a stale result.
    """

    def test_every_source_artifact_is_fingerprinted(self) -> None:
        present = {
            Path(self.ctx.path(relative)).name
            for relative in ASSEMBLE_SOURCE_ARTIFACTS
            if self.ctx.envelope.path(self.ctx.path(relative)).is_file()
        }
        self.assertTrue(present, "fixture writes at least one evidence report")
        fingerprinted = {
            name for name, _ in pipeline_runner._assemble_plan(self.ctx).external_inputs
        }
        self.assertTrue(
            present <= fingerprinted,
            f"{sorted(present - fingerprinted)} not covered by the resume signature",
        )

    def test_a_contribution_puts_its_own_sources_in_the_signature(self) -> None:
        """A contribution reads artifacts the stage does not depend on; they must count."""
        cnv_report = self.ctx.path(cnv_extension.CNV_REPORT)
        self._write(cnv_report, "{}")
        with _CnvRegistration():
            plan = pipeline_runner._assemble_plan(self.ctx)
        fingerprinted = {name for name, _ in plan.external_inputs}
        self.assertIn(Path(cnv_report).name, fingerprinted)
        self.assertIn(cnv_extension.EXTENSION_ID, plan.parameters["result_contributions"])

    def test_a_changed_consensus_changes_the_signature(self) -> None:
        before = pipeline_runner._assemble_plan(self.ctx).external_inputs
        self._write(
            self.ctx.path(pipeline_runner.SV_CONSENSUS_REPORT),
            _consensus()
            .model_copy(update={"caller_names": ["sniffles"]})
            .model_dump_json(indent=2),
        )
        self.assertNotEqual(before, pipeline_runner._assemble_plan(self.ctx).external_inputs)


if __name__ == "__main__":
    unittest.main()

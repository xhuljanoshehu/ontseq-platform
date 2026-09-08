"""A stage runs what the manifest asked for, and says so when it was not asked.

The manifest is a run's scope contract. Two failure modes this guards against:

1. A stage that runs anyway produces evidence nobody requested, and can kill a run over a
   tool the operator had no reason to install — which is how a CNV-only run came to die on
   a missing cuteSV reference.
2. A stage that skips silently looks the same as one that looked and found nothing. Every
   skip here carries a reason saying it is a scope statement, not a negative result.
"""

from __future__ import annotations

import unittest
from dataclasses import dataclass, field
from pathlib import Path

from ontseq_platform.methylation import MethylationPolicy
from ontseq_platform.models import (
    AnalysisModule,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    GenomeBuild,
    InputKind,
    InputSpec,
    ModuleRunStatus,
    SampleManifest,
    SnifflesPolicy,
)
from ontseq_platform.pipeline.runner import (
    StageFailure,
    StagePlan,
    _methylation_execute,
    _methylation_plan,
    _sv_execute,
    _sv_plan,
)


class _ExplodingRunner:
    """Any probe at all is a failure: a skipped stage must not touch its tools."""

    def run(self, argv, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        raise AssertionError(f"a stage that was not requested probed a tool: {list(argv)}")


@dataclass
class _Config:
    sniffles_policy: SnifflesPolicy | None = None
    cutesv_policy: None = None
    sv_consensus_policy: None = None
    sv_evidence_policy: None = None
    methylation_policy: MethylationPolicy | None = None
    reference_fasta: Path | None = None
    gene_annotation: None = None
    cytoband_annotation: None = None
    sv_context_resources: tuple[()] = ()
    aml_knowledge: None = None
    sv_minimum_mean_depth: float = 10.0
    threads: int = 2
    executables: dict[str, str] = field(default_factory=dict)

    def executable(self, name: str) -> str:
        return self.executables.get(name, name)


@dataclass
class _Ctx:
    manifest: SampleManifest
    config: _Config
    runner: _ExplodingRunner = field(default_factory=_ExplodingRunner)


def _manifest(*modules: AnalysisModule) -> SampleManifest:
    return SampleManifest(
        sample_id="GATING_001",
        run_id="GATING_RUN_001",
        input=InputSpec(
            kind=InputKind.ALIGNED_BAM,
            path="/nowhere/input.bam",
            index_path="/nowhere/input.bam.bai",
        ),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id="GATING_REFERENCE_V1",
        ),
        analysis=AnalysisSpec(profile="gating", modules=list(modules)),
    )


def _sniffles_policy() -> SnifflesPolicy:
    return SnifflesPolicy(
        profile_id="gating",
        status="technical_defaults_only",
        note="Synthetic technical policy",
    )


class StructuralVariantGatingTests(unittest.TestCase):
    def test_a_run_without_the_module_never_probes_a_caller(self) -> None:
        ctx = _Ctx(
            manifest=_manifest(AnalysisModule.QC, AnalysisModule.CNV),
            config=_Config(sniffles_policy=_sniffles_policy()),
        )
        plan = _sv_plan(ctx)  # type: ignore[arg-type]
        self.assertEqual(plan.parameters, {"requested": False})
        self.assertEqual(plan.tool_versions, {})

    def test_a_run_without_the_module_records_a_scope_statement(self) -> None:
        ctx = _Ctx(manifest=_manifest(AnalysisModule.QC), config=_Config())
        result = _sv_execute(ctx, StagePlan(parameters={"requested": False}, tool_versions={}))  # type: ignore[arg-type]
        self.assertEqual(result.status, ModuleRunStatus.NOT_RUN)
        self.assertIn("does not request", result.reason)
        self.assertIn("scope statement", result.reason)
        self.assertEqual(result.outputs, [])

    def test_a_cnv_only_run_is_unaffected_by_a_missing_cutesv_reference(self) -> None:
        """The regression this gate exists for: no reference FASTA, and no failure."""
        ctx = _Ctx(
            manifest=_manifest(AnalysisModule.QC, AnalysisModule.CNV, AnalysisModule.REPORT),
            config=_Config(sniffles_policy=_sniffles_policy(), reference_fasta=None),
        )
        self.assertEqual(_sv_plan(ctx).parameters, {"requested": False})  # type: ignore[arg-type]

    def test_requesting_the_module_without_a_policy_still_fails_closed(self) -> None:
        """Asking for SV evidence that cannot be produced is a configuration error."""
        ctx = _Ctx(manifest=_manifest(AnalysisModule.QC, AnalysisModule.SV), config=_Config())
        with self.assertRaises(StageFailure) as raised:
            _sv_plan(ctx)  # type: ignore[arg-type]
        self.assertIn("requires Sniffles2 and/or cuteSV policy", str(raised.exception))


class MethylationGatingTests(unittest.TestCase):
    def test_a_run_without_the_module_never_probes_modkit(self) -> None:
        ctx = _Ctx(manifest=_manifest(AnalysisModule.QC), config=_Config())
        plan = _methylation_plan(ctx)  # type: ignore[arg-type]
        self.assertEqual(plan.parameters, {"requested": False})
        self.assertEqual(plan.tool_versions, {})

    def test_a_run_without_the_module_records_a_scope_statement(self) -> None:
        ctx = _Ctx(manifest=_manifest(AnalysisModule.QC), config=_Config())
        result = _methylation_execute(  # type: ignore[arg-type]
            ctx, StagePlan(parameters={"requested": False}, tool_versions={})
        )
        self.assertEqual(result.status, ModuleRunStatus.NOT_RUN)
        self.assertIn("does not request", result.reason)
        self.assertIn("scope statement", result.reason)

    def test_requesting_the_module_without_a_policy_still_fails_closed(self) -> None:
        ctx = _Ctx(manifest=_manifest(AnalysisModule.METHYLATION), config=_Config())
        with self.assertRaises(StageFailure) as raised:
            _methylation_plan(ctx)  # type: ignore[arg-type]
        self.assertIn("no methylation policy was supplied", str(raised.exception))


class VocabularyTests(unittest.TestCase):
    def test_both_module_gated_stages_use_the_same_key(self) -> None:
        """`requested` means gated on the manifest; `applicable` means gated on the assay.

        Keeping the two words apart is what lets a reader tell "this run did not ask for
        it" from "this assay has nothing for it to measure".
        """
        ctx = _Ctx(manifest=_manifest(AnalysisModule.QC), config=_Config())
        self.assertEqual(
            set(_sv_plan(ctx).parameters),  # type: ignore[arg-type]
            set(_methylation_plan(ctx).parameters),  # type: ignore[arg-type]
        )


if __name__ == "__main__":
    unittest.main()

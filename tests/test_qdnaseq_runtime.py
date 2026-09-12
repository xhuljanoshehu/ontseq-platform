from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from ontseq_platform.cnv.extension import CNV_REPORT, _assemble_execute, _assemble_plan
from ontseq_platform.cnv.qdnaseq import (
    CnvFit,
    QDNAseqCallReport,
    QDNAseqPolicy,
    _events_from_primary_segments,
    run_qdnaseq_ace,
)
from ontseq_platform.demo import build_demo_result
from ontseq_platform.execution import CommandResult
from ontseq_platform.iscn import ISCN_RULE_PROFILE, build_iscn_proposal
from ontseq_platform.models import (
    EventType,
    GenomeBuild,
    ISCNResourceProvenance,
    ISCNSelectionPolicy,
    ModuleRunStatus,
    ReferenceContig,
    ReferenceDictionaryContract,
    ReferenceLock,
    ToolRecord,
)
from ontseq_platform.pipeline.envelope import sha256_file, stage_signature
from ontseq_platform.pipeline.runner import (
    INTAKE_REPORT,
    QC_REPORT,
    SV_CONSENSUS_REPORT,
    SV_REPORT,
    StagePlan,
)
from ontseq_platform.pipeline.runner import _assemble_plan as _core_assemble_plan


class FakeQDNAseqRunner:
    def __init__(self, *, fail: bool = False, event_call: float = -1.0) -> None:
        self.fail = fail
        self.event_call = event_call

    def run(self, argv, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        del timeout_seconds
        args = [str(item) for item in argv]
        if self.fail:
            return CommandResult(tuple(args), 9, "", "synthetic R failure")
        out = Path(args[args.index("--output-dir") + 1])
        sample = args[args.index("--sample-id") + 1]
        out.mkdir(parents=True, exist_ok=True)
        runs = []
        for bin_size, cellularity, ploidy in (
            (100, 0.55, 2.0),
            (500, 0.57, 2.0),
            (1000, 0.56, 2.0),
        ):
            segment_name = f"{sample}.{bin_size}kbp.segments.tsv"
            chromosome_name = f"{sample}.{bin_size}kbp.chromosomes.tsv"
            fit_name = f"{sample}.{bin_size}kbp.ace-fit.png"
            copy_name = f"{sample}.{bin_size}kbp.copy-number.png"
            rds_name = f"{sample}.{bin_size}kbp.segmented.rds"
            bins_name = f"{sample}.{bin_size}kbp.bins.tsv"
            model_name = f"{sample}.{bin_size}kbp.ace-models.tsv"
            (out / segment_name).write_text(
                "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\tqnorm_log10\tcoordinate_system\n"
                f"chr7\t0\t950\t2\t1.0\t{self.event_call}\t-5.0\tzero_based_half_open\n",
                encoding="utf-8",
            )
            (out / chromosome_name).write_text(
                "chromosome\tcopy_number\nchr7\t1.0\n",
                encoding="utf-8",
            )
            (out / fit_name).write_bytes(b"PNG-fit")
            (out / copy_name).write_bytes(b"PNG-copy")
            (out / rds_name).write_bytes(b"RDS")
            (out / bins_name).write_text(
                "chromosome\tstart\tend\treads\tcoordinate_system\n"
                "chr7\t0\t500\t12\tzero_based_half_open\n",
                encoding="utf-8",
            )
            (out / model_name).write_text(
                "cellularity\tploidy\terror\n0.57\t2.0\t0.1\n",
                encoding="utf-8",
            )
            runs.append(
                {
                    "bin_size_kbp": bin_size,
                    "cellularity": cellularity,
                    "ploidy": ploidy,
                    "fit_error": 0.1 + bin_size / 10000,
                    "candidate_count": 3,
                    "alternatives": [
                        {
                            "cellularity": cellularity,
                            "ploidy": ploidy,
                            "fit_error": 0.1,
                        }
                    ],
                    "segment_file": segment_name,
                    "chromosome_file": chromosome_name,
                    "bins_file": bins_name,
                    "model_file": model_name,
                    "fit_plot": fit_name,
                    "copy_number_plot": copy_name,
                    "rds_file": rds_name,
                    "segment_count": 1,
                }
            )
        consensus_name = f"{sample}.consensus.chromosomes.tsv"
        (out / consensus_name).write_text(
            "chromosome\tmedian_copy_number\trounded_copy_number\tagreeing_bins\t"
            "contributing_bins\tmin_copy_number\tmax_copy_number\n"
            "chr7\t1.0\t1\t3\t3\t1.0\t1.0\n",
            encoding="utf-8",
        )
        summary = {
            "schema_version": "0.1.0",
            "sample_id": sample,
            "genome_build": "GRCh37",
            "genome_annotation": "hg19",
            "coordinate_system": "zero_based_half_open",
            "source_coordinate_system": "one_based_inclusive",
            "coordinate_normalization": "qdnaseq_start_minus_one_end_unchanged",
            "primary_bin_size_kbp": 500,
            "ace_penalty": 0.6,
            "ploidy_search": {"min": 1.5, "max": 4.5, "step": 0.05},
            "consensus_strategy": "median_rounded_across_bins",
            "package_versions": {
                "R": "4.4.0",
                "QDNAseq": "1.42.0",
                "ACE": "1.24.0",
                "DNAcopy": "1.80.0",
                "QDNAseq_hg19": "1.36.0",
                "QDNAseq_hg38": None,
            },
            "runs": runs,
            "primary": runs[1],
            "consensus_file": consensus_name,
        }
        (out / f"{sample}.qdnaseq-ace.summary.json").write_text(
            json.dumps(summary),
            encoding="utf-8",
        )
        return CommandResult(tuple(args), 0, "", "")

    def run_to_file(self, argv, output_path, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        del argv, output_path, timeout_seconds
        raise AssertionError("QDNAseq adapter should not use run_to_file")


def _lock() -> ReferenceLock:
    return ReferenceLock(
        reference_id="hg19-test",
        genome_build=GenomeBuild.GRCH37,
        contigs=[ReferenceContig(name="chr7", length=1000)],
        source_fai_sha256="a" * 64,
    )


def _policy() -> QDNAseqPolicy:
    return QDNAseqPolicy(
        profile_id="test",
        cytoband_affected_fraction=0.66,
        expected_qdnaseq_version="1.42.0",
        expected_ace_version="1.24.0",
        note="test policy",
    )


def _empty_cnv_report(*, sample_id: str, genome_build: GenomeBuild) -> QDNAseqCallReport:
    fit = CnvFit(
        bin_size_kbp=500,
        cellularity=1.0,
        ploidy=2.0,
        fit_error=0.01,
        candidate_count=1,
        segment_count=0,
        alternatives=[],
        segment_file="synthetic.segments.tsv",
        chromosome_file="synthetic.chromosomes.tsv",
        bins_file="synthetic.bins.tsv",
        model_file="synthetic.models.tsv",
        fit_plot="synthetic.fit.png",
        copy_number_plot="synthetic.copy-number.png",
        rds_file="synthetic.rds",
    )
    return QDNAseqCallReport(
        sample_id=sample_id,
        genome_build=genome_build,
        status=ModuleRunStatus.NO_CALL,
        primary_fit=fit,
        fits=[fit],
        chromosome_consensus=[],
        events=[],
        tools=[],
        output_files=[],
    )


class QDNAseqRuntimeTests(unittest.TestCase):
    def test_assemble_rejects_cnv_artifact_from_another_sample_or_build(self) -> None:
        baseline = build_demo_result()
        manifest = baseline.manifest
        mismatches = (
            (
                _empty_cnv_report(
                    sample_id="SYNTHETIC_OTHER_SAMPLE",
                    genome_build=manifest.assay.genome_build,
                ),
                "same sample",
            ),
            (
                _empty_cnv_report(
                    sample_id=manifest.sample_id,
                    genome_build=(
                        GenomeBuild.GRCH37
                        if manifest.assay.genome_build == GenomeBuild.GRCH38
                        else GenomeBuild.GRCH38
                    ),
                ),
                "different genome builds",
            ),
        )

        for cnv_report, expected_error in mismatches:
            with (
                self.subTest(expected_error=expected_error),
                tempfile.TemporaryDirectory() as directory,
            ):
                root = Path(directory)
                for relative in (INTAKE_REPORT, QC_REPORT):
                    path = root / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text("{}", encoding="utf-8")
                context = SimpleNamespace(
                    manifest=manifest,
                    config=SimpleNamespace(
                        pipeline_version=baseline.provenance.pipeline_version,
                        git_commit=baseline.provenance.git_commit,
                        resource_context=baseline.reference_context,
                    ),
                    envelope=SimpleNamespace(
                        root=root,
                        path=lambda relative, root=root: root / relative,
                    ),
                    path=lambda template: template.format(sample=manifest.sample_id),
                )

                with (
                    patch(
                        "ontseq_platform.cnv.extension.AlignedBamIntakeReport.model_validate_json",
                        return_value=object(),
                    ),
                    patch(
                        "ontseq_platform.cnv.extension.CraminoQCReport.model_validate_json",
                        return_value=object(),
                    ),
                    patch(
                        "ontseq_platform.cnv.extension.assemble_aligned_bam_mvp",
                        return_value=baseline,
                    ),
                    patch("ontseq_platform.cnv.extension._load_cnv", return_value=cnv_report),
                    self.assertRaisesRegex(ValueError, expected_error),
                ):
                    _assemble_execute(context, StagePlan(parameters={}, tool_versions={}))

    def test_core_assemble_plan_fingerprints_iscn_resources_and_actual_cache(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "annotation.sqlite"
            reference_lock = root / "reference.lock.json"
            cytobands = root / "cytobands.tsv"
            cache.write_bytes(b"synthetic-cache-v1")
            reference_lock.write_text("synthetic-reference-lock", encoding="utf-8")
            cytobands.write_text("synthetic-cytobands", encoding="utf-8")
            checksums = {
                "reference.reference_lock": "a" * 64,
                "reference.cytobands": "b" * 64,
                "reference.annotation_cache": "c" * 64,
            }
            context = SimpleNamespace(
                config=SimpleNamespace(
                    pipeline_version="test",
                    git_commit="0" * 40,
                    annotation_cache=cache,
                    resource_context=SimpleNamespace(
                        resource_checksums=checksums,
                        resource_paths={
                            "reference.reference_lock": str(reference_lock),
                            "reference.cytobands": str(cytobands),
                        },
                    ),
                )
            )

            first = _core_assemble_plan(context)
            first_signature = stage_signature(
                stage="assemble",
                upstream=[],
                parameters=first.parameters,
                tool_versions=first.tool_versions,
                external_inputs=first.external_inputs,
            )

            self.assertEqual(first.parameters["iscn_rule_profile"], ISCN_RULE_PROFILE)
            self.assertEqual(
                first.parameters["iscn_selection_policy"],
                ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1.value,
            )
            self.assertEqual(first.parameters["iscn_reference_lock_sha256"], "a" * 64)
            self.assertEqual(first.parameters["iscn_cytoband_sha256"], "b" * 64)
            self.assertEqual(first.parameters["iscn_annotation_cache_sha256"], "c" * 64)
            self.assertEqual(
                dict(first.external_inputs),
                {
                    "iscn_annotation_cache": sha256_file(cache),
                    "iscn_reference_lock": sha256_file(reference_lock),
                    "iscn_cytobands": sha256_file(cytobands),
                },
            )

            checksums["reference.cytobands"] = "d" * 64
            changed_expected_resource = _core_assemble_plan(context)
            expected_resource_signature = stage_signature(
                stage="assemble",
                upstream=[],
                parameters=changed_expected_resource.parameters,
                tool_versions=changed_expected_resource.tool_versions,
                external_inputs=changed_expected_resource.external_inputs,
            )
            self.assertNotEqual(first_signature, expected_resource_signature)

            checksums["reference.cytobands"] = "b" * 64
            cache.write_bytes(b"synthetic-cache-v2")
            changed = _core_assemble_plan(context)
            changed_signature = stage_signature(
                stage="assemble",
                upstream=[],
                parameters=changed.parameters,
                tool_versions=changed.tool_versions,
                external_inputs=changed.external_inputs,
            )
            self.assertNotEqual(first_signature, changed_signature)

    def test_assemble_resume_fingerprints_raw_and_consensus_sv_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample_id = "SAMPLE_000"
            paths = (
                CNV_REPORT.format(sample=sample_id),
                SV_REPORT.format(sample=sample_id),
                SV_CONSENSUS_REPORT.format(sample=sample_id),
            )
            for relative in paths:
                target = root / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(relative, encoding="utf-8")
            context = SimpleNamespace(
                config=SimpleNamespace(
                    pipeline_version="test",
                    git_commit="0" * 40,
                    annotation_cache=None,
                    resource_context=None,
                ),
                envelope=SimpleNamespace(path=lambda relative: root / relative),
                path=lambda template: template.format(sample=sample_id),
            )

            with patch(
                "ontseq_platform.cnv.extension._settings",
                return_value=SimpleNamespace(policy=_policy()),
            ):
                plan = _assemble_plan(context)

        self.assertEqual(
            {name for name, _checksum in plan.external_inputs},
            {
                Path(CNV_REPORT.format(sample=sample_id)).name,
                Path(SV_REPORT.format(sample=sample_id)).name,
                Path(SV_CONSENSUS_REPORT.format(sample=sample_id)).name,
            },
        )
        self.assertEqual(plan.parameters["iscn_rule_profile"], ISCN_RULE_PROFILE)
        self.assertEqual(
            plan.parameters["iscn_selection_policy"],
            ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1.value,
        )
        self.assertEqual(plan.parameters["iscn_reference_lock_sha256"], "UNAVAILABLE")
        self.assertEqual(plan.parameters["iscn_cytoband_sha256"], "UNAVAILABLE")
        self.assertEqual(plan.parameters["iscn_annotation_cache_sha256"], "UNAVAILABLE")

    def test_assemble_resume_signature_tracks_iscn_resources_and_actual_cache(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "annotation.sqlite"
            reference_lock = root / "reference.lock.json"
            cytobands = root / "cytobands.tsv"
            cache.write_bytes(b"synthetic-cache-v1")
            reference_lock.write_text("synthetic-reference-lock", encoding="utf-8")
            cytobands.write_text("synthetic-cytobands", encoding="utf-8")
            checksums = {
                "reference.reference_lock": "a" * 64,
                "reference.cytobands": "b" * 64,
                "reference.annotation_cache": "c" * 64,
            }
            context = SimpleNamespace(
                config=SimpleNamespace(
                    pipeline_version="test",
                    git_commit="0" * 40,
                    annotation_cache=cache,
                    resource_context=SimpleNamespace(
                        resource_checksums=checksums,
                        resource_paths={
                            "reference.reference_lock": str(reference_lock),
                            "reference.cytobands": str(cytobands),
                        },
                    ),
                ),
                envelope=SimpleNamespace(path=lambda relative: root / relative),
                path=lambda template: template.format(sample="SAMPLE_000"),
            )

            with patch(
                "ontseq_platform.cnv.extension._settings",
                return_value=SimpleNamespace(policy=_policy()),
            ):
                first = _assemble_plan(context)

                def signature(plan):  # noqa: ANN001, ANN202
                    return stage_signature(
                        stage="assemble",
                        upstream=[],
                        parameters=plan.parameters,
                        tool_versions=plan.tool_versions,
                        external_inputs=plan.external_inputs,
                    )

                first_signature = signature(first)
                self.assertEqual(first.parameters["iscn_reference_lock_sha256"], "a" * 64)
                self.assertEqual(first.parameters["iscn_cytoband_sha256"], "b" * 64)
                self.assertEqual(first.parameters["iscn_annotation_cache_sha256"], "c" * 64)
                self.assertEqual(first.parameters["iscn_whole_chromosome_fraction"], 0.90)
                self.assertEqual(first.parameters["iscn_cytoband_affected_fraction"], 0.66)
                self.assertEqual(
                    dict(first.external_inputs),
                    {
                        "iscn_annotation_cache": sha256_file(cache),
                        "iscn_reference_lock": sha256_file(reference_lock),
                        "iscn_cytobands": sha256_file(cytobands),
                    },
                )

                checksums["reference.cytobands"] = "d" * 64
                changed_expected_resource = _assemble_plan(context)
                self.assertNotEqual(first_signature, signature(changed_expected_resource))

                checksums["reference.cytobands"] = "b" * 64
                cache.write_bytes(b"synthetic-cache-v2")
                changed_actual_cache = _assemble_plan(context)
                self.assertNotEqual(first_signature, signature(changed_actual_cache))

            cache.write_bytes(b"synthetic-cache-v1")
            changed_policy = _policy().model_copy(update={"cytoband_affected_fraction": 0.75})
            with patch(
                "ontseq_platform.cnv.extension._settings",
                return_value=SimpleNamespace(policy=changed_policy),
            ):
                changed_policy_plan = _assemble_plan(context)
            self.assertNotEqual(first_signature, signature(changed_policy_plan))

    def test_promotes_complete_result_and_normalizes_event(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "sample.bam"
            bam.write_bytes(b"bam")
            script = root / "runner.R"
            script.write_text("# synthetic", encoding="utf-8")
            output = root / "cnv"

            report = run_qdnaseq_ace(
                bam=bam,
                sample_id="SAMPLE_001",
                genome_build=GenomeBuild.GRCH37,
                reference_lock=_lock(),
                policy=_policy(),
                output_dir=output,
                script=script,
                runner=FakeQDNAseqRunner(),
                threads=2,
            )

            self.assertTrue(output.is_dir())
            self.assertEqual(report.primary_fit.bin_size_kbp, 500)
            self.assertAlmostEqual(report.primary_fit.cellularity, 0.57)
            self.assertAlmostEqual(report.primary_fit.ploidy, 2.0)
            self.assertEqual(len(report.events), 1)
            self.assertEqual(report.events[0].event_type, EventType.CHROMOSOME_LOSS)
            self.assertFalse(report.events[0].whole_chromosome_span_confirmed)
            self.assertEqual(
                (report.events[0].primary.start, report.events[0].primary.end), (0, 950)
            )
            self.assertEqual(report.events[0].length_bp, 950)
            self.assertAlmostEqual(report.events[0].copy_number or -1, 1.0)
            self.assertEqual(report.chromosome_consensus[0].agreeing_bins, 3)
            self.assertTrue(any(name.endswith(".rds") for name in report.output_files))
            self.assertIn(report.primary_fit.bins_file, report.output_files)
            self.assertIn(report.primary_fit.model_file, report.output_files)

    def test_reference_lock_exact_spans_gate_whole_chromosome_iscn(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            segments = Path(directory) / "segments.tsv"
            segments.write_text(
                "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\t"
                "qnorm_log10\tcoordinate_system\n"
                "chr7\t0\t1000\t10\t1.0\t-1.0\t-5.0\tzero_based_half_open\n"
                "chr8\t0\t2000\t20\t3.0\t1.0\t5.0\tzero_based_half_open\n"
                "chr8\t0\t1900\t19\t3.0\t1.0\t5.0\tzero_based_half_open\n",
                encoding="utf-8",
            )
            reference_lock = ReferenceLock(
                reference_id="synthetic-exact-span",
                genome_build=GenomeBuild.GRCH38,
                contigs=[
                    ReferenceContig(name="chr7", length=1000),
                    ReferenceContig(name="chr8", length=2000),
                ],
                source_fai_sha256="a" * 64,
            )
            fit = CnvFit(
                bin_size_kbp=500,
                cellularity=0.5,
                ploidy=2.0,
                fit_error=0.1,
                candidate_count=1,
                segment_count=3,
                alternatives=[],
                segment_file=segments.name,
                chromosome_file="synthetic.chromosomes.tsv",
                bins_file="synthetic.bins.tsv",
                model_file="synthetic.models.tsv",
                fit_plot="synthetic.fit.png",
                copy_number_plot="synthetic.copy-number.png",
                rds_file="synthetic.rds",
            )
            events, warnings = _events_from_primary_segments(
                segments,
                sample_id="SYNTHETIC",
                fit=fit,
                tools=[
                    ToolRecord(name="QDNAseq", version="1.42.0"),
                    ToolRecord(name="ACE", version="1.24.0"),
                ],
                reference_lock=reference_lock,
                minimum_segment_bins=1,
                whole_chromosome_fraction=0.9,
                consensus={},
            )

            self.assertEqual(warnings, [])
            self.assertEqual(
                [event.whole_chromosome_span_confirmed for event in events],
                [True, True, False],
            )
            proposal = build_iscn_proposal(
                events,
                genome_build=GenomeBuild.GRCH38,
                selection_policy=ISCNSelectionPolicy.TECHNICAL_CANDIDATES_V1,
                resource_provenance=ISCNResourceProvenance(
                    genome_build=GenomeBuild.GRCH38,
                    reference_dictionary_contract=ReferenceDictionaryContract.EXACT_FULL,
                    reference_bundle_id="SYNTHETIC_GRCH38",
                    reference_bundle_version="test-v1",
                    reference_lock_sha256="a" * 64,
                    annotation_cache_sha256="c" * 64,
                    cytoband_release="synthetic-test",
                    cytoband_sha256="b" * 64,
                ),
                cnv_source_event_ids={event.event_id for event in events},
            )

            self.assertEqual(proposal.notation, "-7,+8")
            self.assertIn("ISCN rendering is therefore suppressed", events[2].notes[-1])

    def test_returns_no_call_when_ace_call_is_neutral(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "sample.bam"
            bam.write_bytes(b"bam")
            script = root / "runner.R"
            script.write_text("# synthetic", encoding="utf-8")

            report = run_qdnaseq_ace(
                bam=bam,
                sample_id="SAMPLE_002",
                genome_build=GenomeBuild.GRCH37,
                reference_lock=_lock(),
                policy=_policy(),
                output_dir=root / "cnv",
                script=script,
                runner=FakeQDNAseqRunner(event_call=0.0),
            )

            self.assertEqual(report.events, [])
            self.assertEqual(report.status.value, "NO_CALL")

    def test_failure_never_promotes_partial_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "sample.bam"
            bam.write_bytes(b"bam")
            script = root / "runner.R"
            script.write_text("# synthetic", encoding="utf-8")
            output = root / "cnv"

            with self.assertRaisesRegex(ValueError, "exited with code 9"):
                run_qdnaseq_ace(
                    bam=bam,
                    sample_id="SAMPLE_003",
                    genome_build=GenomeBuild.GRCH37,
                    reference_lock=_lock(),
                    policy=_policy(),
                    output_dir=output,
                    script=script,
                    runner=FakeQDNAseqRunner(fail=True),
                )

            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".cnv.*")), [])


if __name__ == "__main__":
    unittest.main()

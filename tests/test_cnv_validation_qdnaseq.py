from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cnv.qdnaseq import (
    CnvChromosomeConsensus,
    CnvFit,
    QDNAseqCallReport,
    QDNAseqPolicy,
)
from ontseq_platform.cnv_validation_contracts import (
    CnvAssessabilityMask,
    CnvCallerLock,
    CnvDataBasis,
    CnvRepeatKind,
    CnvTruthSource,
    CnvValidationSpecimen,
)
from ontseq_platform.cnv_validation_evidence import (
    CnvContributionStatus,
    CnvEvidenceRecordKind,
)
from ontseq_platform.models import (
    EventType,
    GenomeBuild,
    GenomicEvent,
    Locus,
    ModuleRunStatus,
    ToolRecord,
)


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _policy() -> QDNAseqPolicy:
    return QDNAseqPolicy(
        profile_id="qdnaseq-ace-multibin-v1",
        bin_sizes_kbp=[100, 500, 1000],
        primary_bin_size_kbp=500,
        cytoband_affected_fraction=0.66,
        expected_qdnaseq_version="1.42.0",
        expected_ace_version="1.24.0",
        note="Synthetic Block-4 adapter contract.",
    )


def _fit(bin_size_kbp: int) -> CnvFit:
    selected_cellularity = {100: 0.41, 500: 0.52, 1000: 0.63}[bin_size_kbp]
    selected_ploidy = {100: 2.0, 500: 2.2, 1000: 2.4}[bin_size_kbp]
    return CnvFit(
        bin_size_kbp=bin_size_kbp,
        cellularity=selected_cellularity,
        ploidy=selected_ploidy,
        fit_error=0.1 + bin_size_kbp / 10000,
        candidate_count=2,
        segment_count=1,
        alternatives=[
            {
                "cellularity": selected_cellularity - 0.05,
                "ploidy": selected_ploidy + 0.1,
                "fit_error": 0.2 + bin_size_kbp / 10000,
            }
        ],
        segment_file=f"SYNTHETIC_CNV_001.{bin_size_kbp}kbp.segments.tsv",
        chromosome_file=f"SYNTHETIC_CNV_001.{bin_size_kbp}kbp.chromosomes.tsv",
        bins_file=f"SYNTHETIC_CNV_001.{bin_size_kbp}kbp.bins.tsv",
        model_file=f"SYNTHETIC_CNV_001.{bin_size_kbp}kbp.ace-models.tsv",
        fit_plot=f"SYNTHETIC_CNV_001.{bin_size_kbp}kbp.ace-fit.png",
        copy_number_plot=f"SYNTHETIC_CNV_001.{bin_size_kbp}kbp.copy-number.png",
        rds_file=f"SYNTHETIC_CNV_001.{bin_size_kbp}kbp.segmented.rds",
    )


def _report() -> QDNAseqCallReport:
    fits = [_fit(100), _fit(500), _fit(1000)]
    output_files = [
        "SYNTHETIC_CNV_001.qdnaseq-ace.summary.json",
        "SYNTHETIC_CNV_001.consensus.chromosomes.tsv",
    ]
    for fit in fits:
        output_files.extend(
            [
                fit.segment_file,
                fit.chromosome_file,
                fit.bins_file,
                fit.model_file,
                fit.fit_plot,
                fit.copy_number_plot,
                fit.rds_file,
            ]
        )
    return QDNAseqCallReport(
        sample_id="SYNTHETIC_CNV_001",
        genome_build=GenomeBuild.GRCH38,
        status=ModuleRunStatus.COMPLETED,
        primary_fit=fits[1],
        fits=fits,
        chromosome_consensus=[
            CnvChromosomeConsensus(
                chromosome="chr7",
                median_copy_number=1.0,
                rounded_copy_number=1,
                agreeing_bins=3,
                contributing_bins=3,
                min_copy_number=1.0,
                max_copy_number=1.0,
            )
        ],
        events=[
            GenomicEvent(
                event_id="CNV_SYNTHETIC_0001",
                event_type=EventType.DELETION,
                primary=Locus(chromosome="chr7", start=0, end=950),
                copy_number=1.0,
            )
        ],
        tools=[
            ToolRecord(name="QDNAseq", version="1.42.0"),
            ToolRecord(name="ACE", version="1.24.0"),
            ToolRecord(name="R", version="4.4.0"),
        ],
        output_files=sorted(output_files),
    )


def _specimen() -> CnvValidationSpecimen:
    return CnvValidationSpecimen(
        specimen_id="SYNTHETIC_CNV_001",
        biological_specimen_id="SYNTHETIC_BIO_001",
        repeat_kind=CnvRepeatKind.INDEPENDENT,
        access_basis="public",
        material_type="synthetic-DNA",
        genome_build=GenomeBuild.GRCH38,
        data_basis=CnvDataBasis.LCWGS_GENOME_WIDE,
        reference_id="synthetic-grch38",
        reference_sha256=_sha("reference"),
        input_sha256=_sha("input"),
        coverage_x=5.0,
        coverage_definition="synthetic mean autosomal depth",
        tumor_fraction=0.25,
        tumor_fraction_method="synthetic orthogonal fraction",
        tumor_fraction_timepoint="synthetic same aliquot",
        truth_sources=[
            CnvTruthSource(
                method_name="synthetic-karyotype",
                method_version="v1",
                resource_id="synthetic-truth-v1",
                resource_sha256=_sha("truth"),
                provenance_reference="Synthetic fixture only.",
            )
        ],
        assessability_mask=CnvAssessabilityMask(
            resource_id="synthetic-assessability-v1",
            resource_sha256=_sha("assessability"),
            unit="regions",
            definition="Synthetic assessability mask for Block-4 adapter contract tests.",
        ),
    )


def _caller_lock() -> CnvCallerLock:
    return CnvCallerLock(
        caller_id="qdnaseq_ace",
        caller_version="QDNAseq=1.42.0;ACE=1.24.0",
        adapter_policy_sha256=_sha("adapter-policy"),
        execution_identity_sha256=_sha("runtime"),
    )


def _write_native_files(output_dir: Path, report: QDNAseqCallReport) -> None:
    for fit in report.fits:
        (output_dir / fit.segment_file).write_text(
            "chromosome\tstart\tend\tcoordinate_system\tbin_count\t"
            "absolute_copy_number\tcall\tqnorm_log10\n"
            "chr7\t0\t950\tzero_based_half_open\t2\t1.0\t-1.0\t-5.0\n",
            encoding="utf-8",
        )
        (output_dir / fit.chromosome_file).write_text(
            "chromosome\tcopy_number\nchr7\t1.0\n",
            encoding="utf-8",
        )
        assert fit.bins_file is not None
        (output_dir / fit.bins_file).write_text(
            "chromosome\tstart\tend\treads\tgc\tuse\tcoordinate_system\n"
            "chr7\t0\t500\t12\t0.41\tTRUE\tzero_based_half_open\n",
            encoding="utf-8",
        )
        assert fit.model_file is not None
        (output_dir / fit.model_file).write_text(
            "cellularity\tploidy\terror\n0.52\t2.2\t0.1\n",
            encoding="utf-8",
        )
        (output_dir / fit.fit_plot).write_bytes(b"PNG-fit")
        (output_dir / fit.copy_number_plot).write_bytes(b"PNG-copy")
        (output_dir / fit.rds_file).write_bytes(b"RDS")

    for relative_path in report.output_files:
        path = output_dir / relative_path
        if path.exists():
            continue
        if relative_path.endswith(".qdnaseq-ace.summary.json"):
            path.write_text("{}\n", encoding="utf-8")
        elif relative_path.endswith(".consensus.chromosomes.tsv"):
            path.write_text(
                "chromosome\tmedian_copy_number\trounded_copy_number\tagreeing_bins\t"
                "contributing_bins\tmin_copy_number\tmax_copy_number\n"
                "chr7\t1.0\t1\t3\t3\t1.0\t1.0\n",
                encoding="utf-8",
            )
        else:
            raise AssertionError(f"Unhandled synthetic QDNAseq artifact: {relative_path}")


class QDNAseqValidationAdapterTests(unittest.TestCase):
    def test_maps_all_bin_sizes_and_retained_ace_alternatives_without_dropping_native_artifacts(
        self,
    ) -> None:
        from ontseq_platform.cnv_validation_qdnaseq import map_qdnaseq_ace_full_evidence

        report = _report()
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)
            _write_native_files(output_dir, report)

            artifacts, full_evidence = map_qdnaseq_ace_full_evidence(
                report=report,
                policy=_policy(),
                output_dir=output_dir,
                registration_sha256=_sha("registration"),
                specimen=_specimen(),
                caller_lock=_caller_lock(),
                replicate_id="replicate-1",
            )

        self.assertEqual(
            {item.relative_path for item in artifacts},
            set(report.output_files),
        )
        artifact_ids = {item.artifact_id for item in artifacts}
        referenced_artifact_ids = {
            artifact_id for record in full_evidence for artifact_id in record.native_artifact_ids
        }
        self.assertEqual(referenced_artifact_ids, artifact_ids)

        fit_records = [
            record
            for record in full_evidence
            if record.record_kind == CnvEvidenceRecordKind.CALLER_FIT
        ]
        self.assertEqual(len(fit_records), 6)
        self.assertEqual({record.bin_size_kbp for record in fit_records}, {100, 500, 1000})

        for bin_size_kbp in (100, 500, 1000):
            records = [record for record in fit_records if record.bin_size_kbp == bin_size_kbp]
            self.assertEqual(len(records), 2)
            selected = [record for record in records if record.selected_fit]
            alternatives = [record for record in records if record.selected_fit is False]
            self.assertEqual(len(selected), 1)
            self.assertEqual(len(alternatives), 1)
            self.assertEqual(
                selected[0].fit_group_id,
                alternatives[0].fit_group_id,
            )

        primary = next(
            record for record in fit_records if record.bin_size_kbp == 500 and record.selected_fit
        )
        self.assertEqual(
            primary.contribution_status,
            CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
        )
        self.assertTrue(all(record.research_only for record in full_evidence))

    def test_maps_bins_segments_and_chromosome_summaries_for_every_resolution(self) -> None:
        from ontseq_platform.cnv_validation_qdnaseq import map_qdnaseq_ace_full_evidence

        report = _report()
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)
            _write_native_files(output_dir, report)

            _artifacts, full_evidence = map_qdnaseq_ace_full_evidence(
                report=report,
                policy=_policy(),
                output_dir=output_dir,
                registration_sha256=_sha("registration"),
                specimen=_specimen(),
                caller_lock=_caller_lock(),
                replicate_id="replicate-1",
            )

        bin_records = [
            record
            for record in full_evidence
            if record.record_kind == CnvEvidenceRecordKind.CALLER_BIN
        ]
        segment_records = [
            record
            for record in full_evidence
            if record.record_kind == CnvEvidenceRecordKind.CALLER_SEGMENT
        ]
        chromosome_records = [
            record
            for record in full_evidence
            if record.record_kind == CnvEvidenceRecordKind.CHROMOSOME_SUMMARY
        ]

        self.assertEqual({record.bin_size_kbp for record in bin_records}, {100, 500, 1000})
        self.assertEqual({record.bin_size_kbp for record in segment_records}, {100, 500, 1000})
        self.assertEqual(len(bin_records), 3)
        self.assertEqual(len(segment_records), 3)
        self.assertEqual(len(chromosome_records), 4)

        primary_segment = next(record for record in segment_records if record.bin_size_kbp == 500)
        self.assertEqual(
            primary_segment.contribution_status,
            CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
        )
        self.assertIsNotNone(primary_segment.primary)
        assert primary_segment.primary is not None
        self.assertEqual(primary_segment.primary.chromosome, "chr7")
        self.assertEqual(primary_segment.primary.start, 0)
        self.assertEqual(primary_segment.primary.end, 950)
        self.assertEqual(primary_segment.copy_number, 1.0)
        segment_measurements = {
            item.name: item.value for item in primary_segment.numeric_measurements
        }
        self.assertEqual(segment_measurements["bin_count"], 2.0)
        self.assertEqual(segment_measurements["call"], -1.0)
        self.assertEqual(segment_measurements["qnorm_log10"], -5.0)

        primary_bin = next(record for record in bin_records if record.bin_size_kbp == 500)
        bin_measurements = {item.name: item.value for item in primary_bin.numeric_measurements}
        self.assertEqual(bin_measurements["reads"], 12.0)
        self.assertEqual(bin_measurements["gc"], 0.41)

        consensus = next(record for record in chromosome_records if record.bin_size_kbp is None)
        consensus_measurements = {item.name: item.value for item in consensus.numeric_measurements}
        self.assertEqual(consensus_measurements["agreeing_bins"], 3.0)
        self.assertEqual(consensus_measurements["contributing_bins"], 3.0)


    def test_seals_traceable_manifest_without_discarding_multiresolution_evidence(self) -> None:
        from ontseq_platform.cnv_validation_qdnaseq import (
            build_qdnaseq_ace_validation_manifest,
        )

        report = _report()
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)
            _write_native_files(output_dir, report)

            manifest = build_qdnaseq_ace_validation_manifest(
                manifest_id="SYNTHETIC_QDNASEQ_EVIDENCE_001",
                report=report,
                policy=_policy(),
                output_dir=output_dir,
                registration_sha256=_sha("registration"),
                specimen=_specimen(),
                caller_lock=_caller_lock(),
                replicate_id="replicate-1",
            )

        self.assertTrue(manifest.retain_all_evidence)
        self.assertEqual(len(manifest.native_artifacts), len(report.output_files))
        self.assertEqual(len(manifest.normalized_events), 1)
        normalized = manifest.normalized_events[0]
        self.assertEqual(normalized.bin_size_kbp, 500)
        self.assertEqual(
            normalized.contribution_status,
            CnvContributionStatus.USED_FOR_PRIMARY_ANALYSIS,
        )
        self.assertEqual(normalized.primary, report.events[0].primary)
        self.assertEqual(normalized.normalized_copy_number, report.events[0].copy_number)
        self.assertEqual(normalized.normalization_policy_sha256, _caller_lock().adapter_policy_sha256)

        source_ids = set(normalized.source_full_evidence_ids)
        source_records = [
            record for record in manifest.full_evidence if record.record_id in source_ids
        ]
        self.assertEqual(len(source_records), 1)
        self.assertEqual(source_records[0].record_kind, CnvEvidenceRecordKind.CALLER_SEGMENT)
        self.assertEqual(source_records[0].bin_size_kbp, 500)

        referenced_artifacts = {
            artifact_id
            for record in manifest.full_evidence
            for artifact_id in record.native_artifact_ids
        }
        self.assertEqual(
            referenced_artifacts,
            {artifact.artifact_id for artifact in manifest.native_artifacts},
        )

    def test_manifest_build_fails_closed_if_a_declared_native_artifact_is_missing(self) -> None:
        from ontseq_platform.cnv_validation_qdnaseq import (
            build_qdnaseq_ace_validation_manifest,
        )

        report = _report()
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)
            _write_native_files(output_dir, report)
            missing = output_dir / report.primary_fit.segment_file
            missing.unlink()

            with self.assertRaisesRegex(ValueError, "artifact is missing"):
                build_qdnaseq_ace_validation_manifest(
                    manifest_id="SYNTHETIC_QDNASEQ_EVIDENCE_001",
                    report=report,
                    policy=_policy(),
                    output_dir=output_dir,
                    registration_sha256=_sha("registration"),
                    specimen=_specimen(),
                    caller_lock=_caller_lock(),
                    replicate_id="replicate-1",
                )


if __name__ == "__main__":
    unittest.main()

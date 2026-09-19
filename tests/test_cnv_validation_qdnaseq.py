from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from ontseq_platform.cnv.qdnaseq import CnvFit, QDNAseqCallReport, QDNAseqPolicy
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
        chromosome_consensus=[],
        events=[
            GenomicEvent(
                event_id="CNV_SYNTHETIC_0001",
                event_type=EventType.DELETION,
                primary=Locus(chromosome="chr7", start=1_000_000, end=6_000_000),
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


class QDNAseqValidationAdapterTests(unittest.TestCase):
    def test_maps_all_bin_sizes_and_retained_ace_alternatives_without_dropping_native_artifacts(
        self,
    ) -> None:
        from ontseq_platform.cnv_validation_qdnaseq import map_qdnaseq_ace_full_evidence

        report = _report()
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)
            for relative_path in report.output_files:
                (output_dir / relative_path).write_bytes(
                    f"synthetic-native-artifact:{relative_path}".encode()
                )

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


if __name__ == "__main__":
    unittest.main()

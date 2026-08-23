from __future__ import annotations

import unittest

from ontseq_platform.cnv.cytobands import Cytoband, CytobandTable
from ontseq_platform.cnv.lea_compat import (
    LEA_ACE_2026_HG19,
    LeaCompatibilityError,
    lea_ace_call_set_from_outputs,
    parse_lea_cn_csv,
    parse_lea_dels_dups_csv,
)
from ontseq_platform.cnv.lea_truth_tables import (
    LeaTruthTableError,
    parse_lea_gt_full_csv,
    parse_lea_gt_tsv,
)
from ontseq_platform.cnv.models import CnvDataBasis
from ontseq_platform.cnv.states import CopyNumberState
from ontseq_platform.models import GenomeBuild, ModuleRunStatus


def _cn_lines(*, overrides: dict[str, tuple[float, float]] | None = None) -> list[str]:
    """Return a complete synthetic CN.csv for chromosomes 1-22,X,Y.

    overrides maps contig -> (Copies, Ploidy).
    """
    overrides = overrides or {}
    lines = ["Chromosome,Copies,Ploidy,CNA"]
    for index in range(1, 23):
        contig = str(index)
        copies, ploidy = overrides.get(contig, (2.0, 2.0))
        lines.append(f"{contig},{copies:g},{ploidy:g},{copies - ploidy:g}")
    for contig, defaults in (("X", (1.0, 1.0)), ("Y", (1.0, 1.0))):
        copies, ploidy = overrides.get(contig, defaults)
        lines.append(f"{contig},{copies:g},{ploidy:g},{copies - ploidy:g}")
    return lines


def _contig_lengths() -> dict[str, int]:
    return {str(index): 1000 for index in range(1, 23)} | {"X": 1000, "Y": 1000}


def _cytobands(build: GenomeBuild = GenomeBuild.GRCH37) -> CytobandTable:
    return CytobandTable(
        resource_id="synthetic-hg19-cytobands",
        genome_build=build,
        source="synthetic",
        bands=[
            Cytoband(contig="chr5", start=0, end=400, name="p15", stain="gneg"),
            Cytoband(contig="chr5", start=400, end=600, name="q13.1", stain="gpos25"),
            Cytoband(contig="chr5", start=600, end=800, name="q13.2", stain="gneg"),
            Cytoband(contig="chr5", start=800, end=1000, name="q13.3", stain="gpos50"),
        ],
    )


class LeaCnCsvTests(unittest.TestCase):
    def test_cna_identity_is_checked(self) -> None:
        with self.assertRaises(LeaCompatibilityError) as context:
            parse_lea_cn_csv(["Chromosome,Copies,Ploidy,CNA", "5,1,2,0"])
        self.assertIn("Copies - Ploidy", str(context.exception))

    def test_duplicate_chromosome_is_rejected(self) -> None:
        with self.assertRaises(LeaCompatibilityError):
            parse_lea_cn_csv(["Chromosome,Copies,Ploidy,CNA", "5,2,2,0", "chr5,2,2,0"])

    def test_empty_cn_table_is_not_a_negative_result(self) -> None:
        with self.assertRaises(LeaCompatibilityError):
            parse_lea_cn_csv(["Chromosome,Copies,Ploidy,CNA"])


class LeaBandCsvTests(unittest.TestCase):
    def test_historical_fraction_threshold_is_enforced(self) -> None:
        with self.assertRaises(LeaCompatibilityError) as context:
            parse_lea_dels_dups_csv(["chromosome,name,event,frac_abr", "chr5,q13.1,del,0.65"])
        self.assertIn("threshold", str(context.exception))

    def test_conflicting_event_for_same_band_is_rejected(self) -> None:
        with self.assertRaises(LeaCompatibilityError):
            parse_lea_dels_dups_csv(
                [
                    "chromosome,name,event,frac_abr",
                    "chr5,q13.1,del,0.70",
                    "chr5,q13.1,dup,0.80",
                ]
            )


class LeaAceCallSetTests(unittest.TestCase):
    def test_partial_band_call_maps_through_locked_cytobands(self) -> None:
        result = lea_ace_call_set_from_outputs(
            cn_lines=_cn_lines(),
            dels_dups_lines=[
                "chromosome,name,event,frac_abr",
                "chr5,q13.2,del,0.75",
            ],
            cytobands=_cytobands(),
            contig_lengths=_contig_lengths(),
            call_set_id="LEA_SYNTHETIC_001",
            sample_id="SYNTHETIC_AML_001",
            genome_build=GenomeBuild.GRCH37,
            data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
        )
        self.assertEqual(result.status, ModuleRunStatus.COMPLETED)
        self.assertFalse(result.reportable)
        self.assertTrue(result.research_only)
        self.assertEqual(len(result.segments), 1)
        event = result.segments[0]
        self.assertEqual((event.contig, event.start, event.end), ("5", 600, 800))
        self.assertEqual(event.state, CopyNumberState.LOSS)
        self.assertIsNone(event.copy_number)

    def test_whole_chromosome_call_preserves_absolute_copies(self) -> None:
        result = lea_ace_call_set_from_outputs(
            cn_lines=_cn_lines(overrides={"7": (3.0, 2.0)}),
            dels_dups_lines=["chromosome,name,event,frac_abr"],
            cytobands=_cytobands(),
            contig_lengths=_contig_lengths(),
            call_set_id="LEA_SYNTHETIC_002",
            sample_id="SYNTHETIC_AML_002",
            genome_build=GenomeBuild.GRCH37,
            data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
        )
        self.assertEqual(len(result.segments), 1)
        event = result.segments[0]
        self.assertEqual((event.contig, event.start, event.end), ("7", 0, 1000))
        self.assertEqual(event.state, CopyNumberState.GAIN)
        self.assertEqual(event.copy_number, 3.0)

    def test_whole_chromosome_and_partial_call_conflict_fails_closed(self) -> None:
        with self.assertRaises(LeaCompatibilityError) as context:
            lea_ace_call_set_from_outputs(
                cn_lines=_cn_lines(overrides={"5": (1.0, 2.0)}),
                dels_dups_lines=[
                    "chromosome,name,event,frac_abr",
                    "chr5,q13.2,del,0.75",
                ],
                cytobands=_cytobands(),
                contig_lengths=_contig_lengths(),
                call_set_id="LEA_SYNTHETIC_003",
                sample_id="SYNTHETIC_AML_003",
                genome_build=GenomeBuild.GRCH37,
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
            )
        self.assertIn("artifacts are inconsistent", str(context.exception))

    def test_grch38_is_never_silently_lifted(self) -> None:
        with self.assertRaises(LeaCompatibilityError) as context:
            lea_ace_call_set_from_outputs(
                cn_lines=_cn_lines(),
                dels_dups_lines=["chromosome,name,event,frac_abr"],
                cytobands=_cytobands(),
                contig_lengths=_contig_lengths(),
                call_set_id="LEA_SYNTHETIC_004",
                sample_id="SYNTHETIC_AML_004",
                genome_build=GenomeBuild.GRCH38,
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
            )
        self.assertIn("never performs lift-over", str(context.exception))

    def test_cytoband_build_must_match(self) -> None:
        with self.assertRaises(LeaCompatibilityError):
            lea_ace_call_set_from_outputs(
                cn_lines=_cn_lines(),
                dels_dups_lines=["chromosome,name,event,frac_abr"],
                cytobands=_cytobands(GenomeBuild.GRCH38),
                contig_lengths=_contig_lengths(),
                call_set_id="LEA_SYNTHETIC_005",
                sample_id="SYNTHETIC_AML_005",
                genome_build=GenomeBuild.GRCH37,
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
            )

    def test_missing_autosome_refuses_the_frozen_whole_genome_profile(self) -> None:
        lines = [line for line in _cn_lines() if not line.startswith("12,")]
        with self.assertRaises(LeaCompatibilityError) as context:
            lea_ace_call_set_from_outputs(
                cn_lines=lines,
                dels_dups_lines=["chromosome,name,event,frac_abr"],
                cytobands=_cytobands(),
                contig_lengths=_contig_lengths(),
                call_set_id="LEA_SYNTHETIC_006",
                sample_id="SYNTHETIC_AML_006",
                genome_build=GenomeBuild.GRCH37,
                data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
            )
        self.assertIn("lacks autosome", str(context.exception))

    def test_valid_empty_alteration_set_is_completed_but_bounded(self) -> None:
        result = lea_ace_call_set_from_outputs(
            cn_lines=_cn_lines(),
            dels_dups_lines=["chromosome,name,event,frac_abr"],
            cytobands=_cytobands(),
            contig_lengths=_contig_lengths(),
            call_set_id="LEA_SYNTHETIC_007",
            sample_id="SYNTHETIC_AML_007",
            genome_build=GenomeBuild.GRCH37,
            data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
        )
        self.assertEqual(result.status, ModuleRunStatus.COMPLETED)
        self.assertTrue(result.reports_biological_negative)
        self.assertEqual(result.segments, [])
        self.assertTrue(
            any("not an assay-wide or clinical negative" in x for x in result.limitations)
        )

    def test_runtime_versions_are_not_claimed_as_verified(self) -> None:
        result = lea_ace_call_set_from_outputs(
            cn_lines=_cn_lines(),
            dels_dups_lines=["chromosome,name,event,frac_abr"],
            cytobands=_cytobands(),
            contig_lengths=_contig_lengths(),
            call_set_id="LEA_SYNTHETIC_008",
            sample_id="SYNTHETIC_AML_008",
            genome_build=GenomeBuild.GRCH37,
            data_basis=CnvDataBasis.LOW_COVERAGE_WGS,
        )
        self.assertFalse(LEA_ACE_2026_HG19.runtime_versions_verified)
        self.assertTrue(any("remain unverified" in x for x in result.warnings))


class HistoricalTruthTableTests(unittest.TestCase):
    def test_gt_tsv_preserves_uncertain_iscn_verbatim(self) -> None:
        rows = parse_lea_gt_tsv(["sample\tiscn", "SYNTHETIC_001\t45,XY,-7,del(7)(q?31)"])
        self.assertEqual(rows[0].karyotype, "45,XY,-7,del(7)(q?31)")

    def test_gt_full_preserves_labels_without_recomputing_them(self) -> None:
        rows = parse_lea_gt_full_csv(
            [
                "sample_name,karyotype_cg,karyotype_ont,cellularity,complex,monosomal,mrc,mrca,mra",
                'SYNTHETIC_001,"46,XX,del(5)(q13q33)","46,XX,del(5)(q13q33)",0.60,1,0,1,1,1',
            ]
        )
        row = rows[0]
        self.assertEqual(row.cellularity, 0.60)
        self.assertTrue(row.complex_karyotype)
        self.assertFalse(row.monosomal_karyotype)

    def test_gt_full_rejects_nonbinary_historical_flag(self) -> None:
        with self.assertRaises(LeaTruthTableError):
            parse_lea_gt_full_csv(
                [
                    "sample_name,karyotype_cg,karyotype_ont,cellularity,complex,"
                    "monosomal,mrc,mrca,mra",
                    'SYNTHETIC_001,"46,XX","46,XX",0.60,2,0,1,1,1',
                ]
            )


if __name__ == "__main__":
    unittest.main()

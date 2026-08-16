from __future__ import annotations

import unittest
from pathlib import Path

from ontseq_platform.cnv.cytobands import load_cytoband_file, parse_ucsc_cytoband
from ontseq_platform.cnv.models import CnvSegment, CnvTruthSource
from ontseq_platform.cnv.states import CopyNumberState
from ontseq_platform.cnv.truth import (
    convert_karyotype,
    truth_from_fish,
    truth_from_karyotype,
    truth_from_segments,
)
from ontseq_platform.models import GenomeBuild

CYTOBANDS = Path("examples/references/synthetic.cytobands.txt")


def _table():
    return load_cytoband_file(
        CYTOBANDS, genome_build=GenomeBuild.GRCH38, resource_id="SYNTHETIC_CYTOBANDS_V1"
    )


class CytobandTableTests(unittest.TestCase):
    def test_fixture_parses_and_covers_expected_contigs(self) -> None:
        table = _table()
        contigs = {band.contig for band in table.bands}
        self.assertEqual(contigs, {"chr2", "chr5", "chr7", "chr8", "chr17", "chr20"})
        self.assertEqual(table.contig_length("chr5"), 181_538_259)
        self.assertEqual(table.contig_length("5"), 181_538_259)

    def test_arm_intervals(self) -> None:
        table = _table()
        self.assertEqual(table.arm_interval("chr17", "p"), (0, 25_100_000))
        self.assertEqual(table.arm_interval("chr17", "q"), (25_100_000, 83_257_441))

    def test_band_interval_and_uncertainty(self) -> None:
        table = _table()
        self.assertEqual(table.band_interval("chr5", "q13"), (60_000_000, 80_000_000))
        self.assertEqual(table.band_uncertainty("chr5", "q13"), 20_000_000)

    def test_unknown_band_raises(self) -> None:
        with self.assertRaises(KeyError):
            _table().band_interval("chr5", "q99")

    def test_malformed_band_designation_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _table().band_interval("chr5", "banana")

    def test_overlapping_bands_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            parse_ucsc_cytoband(
                ["chr1\t0\t100\tp11\tgneg", "chr1\t50\t200\tp12\tgpos"],
                genome_build=GenomeBuild.GRCH38,
                resource_id="BAD",
                source="test",
            )

    def test_non_canonical_contigs_are_skipped_not_fatal(self) -> None:
        table = parse_ucsc_cytoband(
            [
                "chr1\t0\t100\tp11\tgneg",
                "chr1_KI270706v1_random\t0\t100\tp11\tgneg",
                "chrUn_KI270302v1\t0\t100\tp11\tgneg",
            ],
            genome_build=GenomeBuild.GRCH38,
            resource_id="OK",
            source="test",
        )
        self.assertEqual(len(table.bands), 1)


class KaryotypeConversionTests(unittest.TestCase):
    def test_whole_chromosome_gain(self) -> None:
        conversion = convert_karyotype("47,XY,+8[20]", _table())
        self.assertTrue(conversion.is_complete)
        self.assertEqual(len(conversion.segments), 1)
        segment = conversion.segments[0]
        self.assertEqual((segment.contig, segment.start, segment.end), ("8", 0, 145_138_636))
        self.assertEqual(segment.state, CopyNumberState.GAIN)
        self.assertEqual(segment.copy_number, 3.0)
        self.assertEqual(conversion.clone_cell_counts, [20])

    def test_whole_chromosome_loss(self) -> None:
        conversion = convert_karyotype("45,XX,-7", _table())
        self.assertEqual(conversion.segments[0].state, CopyNumberState.LOSS)
        self.assertEqual(conversion.segments[0].copy_number, 1.0)

    def test_interstitial_deletion_carries_band_uncertainty(self) -> None:
        conversion = convert_karyotype("46,XY,del(5)(q13q33)", _table())
        segment = conversion.segments[0]
        self.assertEqual((segment.start, segment.end), (60_000_000, 155_000_000))
        self.assertEqual(segment.state, CopyNumberState.LOSS)
        # Uncertainty equals the width of the named bands, not zero.
        self.assertEqual(segment.start_uncertainty_bp, 20_000_000)
        self.assertEqual(segment.end_uncertainty_bp, 15_000_000)

    def test_terminal_deletion_runs_to_the_arm_end(self) -> None:
        conversion = convert_karyotype("46,XX,del(7)(q22)", _table())
        segment = conversion.segments[0]
        self.assertEqual((segment.start, segment.end), (98_000_000, 159_345_973))

    def test_isochromosome_asserts_both_arm_changes(self) -> None:
        conversion = convert_karyotype("46,XY,i(17)(q10)", _table())
        states = {(s.state, s.start, s.end) for s in conversion.segments}
        self.assertIn((CopyNumberState.LOSS, 0, 25_100_000), states)
        self.assertIn((CopyNumberState.GAIN, 25_100_000, 83_257_441), states)

    def test_balanced_translocation_asserts_no_copy_number_change(self) -> None:
        conversion = convert_karyotype("46,XX,t(9;22)(q34;q11)", _table())
        self.assertEqual(conversion.segments, [])
        self.assertEqual(conversion.balanced_constructs, ["t(9;22)(q34;q11)"])
        self.assertTrue(conversion.is_complete)

    def test_unsupported_construct_is_recorded_never_dropped(self) -> None:
        conversion = convert_karyotype("46,XY,der(1;7)(q10;p10),+mar", _table())
        self.assertFalse(conversion.is_complete)
        tokens = {item.token for item in conversion.unsupported}
        self.assertIn("der(1;7)(q10;p10)", tokens)
        self.assertIn("+mar", tokens)
        for item in conversion.unsupported:
            self.assertTrue(item.reason)

    def test_uncertainty_marker_is_refused(self) -> None:
        conversion = convert_karyotype("46,XY,del(5)(q13q3?3)", _table())
        self.assertFalse(conversion.is_complete)

    def test_unusual_sex_complement_is_flagged(self) -> None:
        conversion = convert_karyotype("45,X,-7", _table())
        tokens = {item.token for item in conversion.unsupported}
        self.assertIn("X", tokens)

    def test_multiple_clones_are_merged_with_a_warning_and_deduplicated(self) -> None:
        conversion = convert_karyotype("47,XY,+8[15]/47,XY,+8[5]", _table())
        # The same alteration in two clones is one genomic claim, not two segments.
        self.assertEqual(len(conversion.segments), 1)
        self.assertTrue(any("clone" in warning.lower() for warning in conversion.warnings))
        self.assertEqual(conversion.clone_cell_counts, [15, 5])

    def test_isochromosome_in_two_clones_does_not_duplicate_segments(self) -> None:
        conversion = convert_karyotype("46,XY,i(17)(q10)[10]/46,XY,i(17)(q10)[10]", _table())
        self.assertEqual(len(conversion.segments), 2)

    def test_complex_aml_karyotype(self) -> None:
        conversion = convert_karyotype("46,XY,del(5)(q13q33),-7,+8,t(9;22)(q34;q11)", _table())
        self.assertTrue(conversion.is_complete)
        by_contig = {segment.contig: segment.state for segment in conversion.segments}
        self.assertEqual(by_contig["5"], CopyNumberState.LOSS)
        self.assertEqual(by_contig["7"], CopyNumberState.LOSS)
        self.assertEqual(by_contig["8"], CopyNumberState.GAIN)
        self.assertEqual(conversion.balanced_constructs, ["t(9;22)(q34;q11)"])

    def test_segments_never_overlap(self) -> None:
        conversion = convert_karyotype("46,XY,i(17)(q10),del(5)(q13q33),+8", _table())
        by_contig: dict[str, list[CnvSegment]] = {}
        for segment in conversion.segments:
            by_contig.setdefault(segment.contig, []).append(segment)
        for items in by_contig.values():
            items.sort(key=lambda item: item.start)
            for previous, current in zip(items, items[1:], strict=False):
                self.assertLessEqual(previous.end, current.start)


class TruthSetConstructionTests(unittest.TestCase):
    def test_karyotype_truth_declares_resolution_and_limitations(self) -> None:
        truth, conversion = truth_from_karyotype(
            truth_id="KARYO_001",
            sample_id="SYNTHETIC_AML_001",
            karyotype="46,XY,del(5)(q13q33),+8",
            cytobands=_table(),
            source_version="synthetic-cytobands-v1",
        )
        self.assertEqual(truth.source, CnvTruthSource.ISCN_KARYOTYPE)
        self.assertEqual(truth.background_state, CopyNumberState.NEUTRAL)
        self.assertEqual(truth.resolution_bp, 10_000_000)
        self.assertTrue(any("band" in item.lower() for item in truth.limitations))
        self.assertTrue(any("loss of heterozygosity" in item.lower() for item in truth.limitations))
        self.assertTrue(conversion.is_complete)

    def test_karyotype_truth_surfaces_unconverted_constructs_in_limitations(self) -> None:
        truth, _ = truth_from_karyotype(
            truth_id="KARYO_002",
            sample_id="SYNTHETIC_AML_002",
            karyotype="46,XY,der(1;7)(q10;p10)",
            cytobands=_table(),
            source_version="synthetic-cytobands-v1",
        )
        self.assertTrue(any("refused to interpret" in item for item in truth.limitations))

    def test_fish_truth_is_open_world_and_probe_scoped(self) -> None:
        truth = truth_from_fish(
            truth_id="FISH_001",
            sample_id="SYNTHETIC_AML_001",
            genome_build=GenomeBuild.GRCH38,
            source_version="probe-set-v1",
            probes=[("chr8", 90_000_000, 90_100_000, "MYC", 3.0)],
        )
        self.assertEqual(truth.background_state, CopyNumberState.NO_CALL)
        self.assertEqual(len(truth.informative_regions), 1)
        self.assertTrue(any("unexamined" in item for item in truth.limitations))

    def test_closed_world_truth_requires_a_resolution(self) -> None:
        with self.assertRaises(ValueError):
            truth_from_segments(
                truth_id="ARRAY_BAD",
                sample_id="SYNTHETIC_AML_001",
                genome_build=GenomeBuild.GRCH38,
                source=CnvTruthSource.SNP_ARRAY,
                source_version="v1",
                segments=[],
                resolution_bp=0,
            )

    def test_array_truth_can_carry_copy_neutral_loh(self) -> None:
        truth = truth_from_segments(
            truth_id="ARRAY_001",
            sample_id="SYNTHETIC_AML_001",
            genome_build=GenomeBuild.GRCH38,
            source=CnvTruthSource.SNP_ARRAY,
            source_version="v1",
            resolution_bp=50_000,
            segments=[
                CnvSegment(
                    contig="chr17",
                    start=0,
                    end=25_100_000,
                    state=CopyNumberState.COPY_NEUTRAL_LOH,
                    copy_number=2.0,
                )
            ],
        )
        self.assertEqual(truth.segments[0].state, CopyNumberState.COPY_NEUTRAL_LOH)
        self.assertEqual(truth.resolution_bp, 50_000)

    def test_fish_truth_cannot_claim_a_genome_wide_neutral_background(self) -> None:
        with self.assertRaises(ValueError):
            truth_from_segments(
                truth_id="FISH_BAD",
                sample_id="SYNTHETIC_AML_001",
                genome_build=GenomeBuild.GRCH38,
                source=CnvTruthSource.FISH,
                source_version="v1",
                segments=[],
                resolution_bp=1000,
                closed_world=True,
            )


if __name__ == "__main__":
    unittest.main()

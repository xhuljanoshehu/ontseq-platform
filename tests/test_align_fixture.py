from __future__ import annotations

import unittest

from ontseq_platform.align_fixture import (
    CONTIG_LENGTHS,
    DELETION_BP,
    DELETION_READS,
    READ_ARM_BP,
    READ_GROUP_ID,
    READ_START,
    REFERENCE_READS,
    REVERSE_READS,
    SAMPLE_ID,
    SECOND_READ_GROUP_ID,
    deterministic_sequence,
    format_fasta,
    reverse_complement,
    synthetic_reference,
    unaligned_sam_text,
)


def _header_field(line: str, tag: str) -> str:
    """Read one ``TAG:value`` field out of a SAM header line, or fail loudly."""
    for field in line.split("\t")[1:]:
        key, _, value = field.partition(":")
        if key == tag:
            return value
    raise AssertionError(f"header line carries no {tag} field: {line}")


def _read_group_of(record: str) -> str:
    """The record's ``RG`` tag as a whole value, never as a substring match."""
    for field in record.split("\t")[11:]:
        if field.startswith("RG:Z:"):
            return field[len("RG:Z:") :]
    raise AssertionError(f"record carries no RG tag: {record[:60]}")


class SequenceTests(unittest.TestCase):
    def test_sequence_is_reproducible_for_a_seed(self) -> None:
        self.assertEqual(
            deterministic_sequence(500, seed="x"), deterministic_sequence(500, seed="x")
        )

    def test_different_seeds_give_different_sequences(self) -> None:
        self.assertNotEqual(
            deterministic_sequence(500, seed="x"), deterministic_sequence(500, seed="y")
        )

    def test_a_prefix_is_stable_as_the_length_grows(self) -> None:
        """The keystream must extend, not reshuffle, or fixtures change size to size."""
        self.assertTrue(
            deterministic_sequence(900, seed="x").startswith(deterministic_sequence(300, seed="x"))
        )

    def test_only_dna_bases_are_emitted(self) -> None:
        self.assertEqual(set(deterministic_sequence(2000, seed="x")) - set("ACGT"), set())

    def test_length_is_exact(self) -> None:
        for length in (1, 31, 32, 33, 1000):
            self.assertEqual(len(deterministic_sequence(length, seed="x")), length)

    def test_zero_length_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            deterministic_sequence(0, seed="x")

    def test_reverse_complement_is_an_involution(self) -> None:
        sequence = deterministic_sequence(200, seed="x")
        self.assertEqual(reverse_complement(reverse_complement(sequence)), sequence)


class ReferenceTests(unittest.TestCase):
    def test_contig_lengths_match_the_declared_layout(self) -> None:
        reference = synthetic_reference()
        for name, length in CONTIG_LENGTHS:
            self.assertEqual(len(reference[name]), length)

    def test_contigs_are_not_copies_of_each_other(self) -> None:
        reference = synthetic_reference()
        self.assertNotEqual(reference["chr1"][:1000], reference["chr2"][:1000])

    def test_fasta_wraps_at_the_requested_width(self) -> None:
        rendered = format_fasta({"c": "ACGT" * 40}, line_length=60)
        lines = rendered.splitlines()
        self.assertEqual(lines[0], ">c")
        self.assertEqual(len(lines[1]), 60)
        self.assertEqual("".join(lines[1:]), "ACGT" * 40)

    def test_fasta_rejects_a_nonsense_width(self) -> None:
        with self.assertRaises(ValueError):
            format_fasta({"c": "ACGT"}, line_length=0)


class UnalignedSamTests(unittest.TestCase):
    def setUp(self) -> None:
        self.reference = synthetic_reference()
        self.text = unaligned_sam_text(self.reference)
        self.lines = self.text.splitlines()
        self.records = [line for line in self.lines if not line.startswith("@")]

    def test_no_sequence_dictionary_is_declared(self) -> None:
        """An unaligned BAM must not carry @SQ records; otherwise it is not unaligned."""
        self.assertFalse([line for line in self.lines if line.startswith("@SQ")])

    def test_two_read_groups_are_declared_each_with_a_sample(self) -> None:
        """Two, because one cannot distinguish "preserved" from "collapsed onto the first".

        An aligner that dropped every read onto whichever group it saw first would look
        perfect against a single-group fixture. The second group is what makes the claim
        that read groups survive alignment testable at all.
        """
        read_groups = [line for line in self.lines if line.startswith("@RG")]
        self.assertEqual(len(read_groups), 2)
        self.assertEqual(
            [_header_field(line, "ID") for line in read_groups],
            [READ_GROUP_ID, SECOND_READ_GROUP_ID],
        )
        for line in read_groups:
            self.assertEqual(_header_field(line, "SM"), SAMPLE_ID)

    def test_every_record_is_unmapped(self) -> None:
        for record in self.records:
            fields = record.split("\t")
            self.assertEqual(fields[1], "4", record[:40])
            self.assertEqual(fields[2], "*")
            self.assertEqual(fields[5], "*")

    def test_expected_number_of_reads(self) -> None:
        self.assertEqual(len(self.records), DELETION_READS + REFERENCE_READS + REVERSE_READS)

    def test_every_record_carries_the_tags_alignment_must_preserve(self) -> None:
        for record in self.records:
            # Compared as a whole field, not as a substring: ``SYNTHETIC_ALIGN_RG`` is a
            # prefix of ``SYNTHETIC_ALIGN_RG2``, so a containment check would pass on a
            # record belonging to the other group and prove nothing.
            self.assertIn(_read_group_of(record), {READ_GROUP_ID, SECOND_READ_GROUP_ID})
            self.assertIn("MM:Z:C+m?", record)
            self.assertIn("ML:B:C,", record)

    def test_both_read_groups_actually_carry_reads(self) -> None:
        """A declared-but-empty second group would restore the blind spot silently."""
        counts: dict[str, int] = {}
        for record in self.records:
            group = _read_group_of(record)
            counts[group] = counts.get(group, 0) + 1
        self.assertEqual(counts[SECOND_READ_GROUP_ID], REVERSE_READS)
        self.assertEqual(counts[READ_GROUP_ID], DELETION_READS + REFERENCE_READS)

    def test_the_reverse_reads_are_the_ones_on_the_second_group(self) -> None:
        """Which reads sit on which group has to be stated, or the counts prove nothing."""
        for record in self.records:
            reverse = record.startswith("SYNTH_ALIGN_REV_")
            expected = SECOND_READ_GROUP_ID if reverse else READ_GROUP_ID
            self.assertEqual(_read_group_of(record), expected, record[:40])

    def test_sequence_and_quality_lengths_agree(self) -> None:
        for record in self.records:
            fields = record.split("\t")
            self.assertEqual(len(fields[9]), len(fields[10]))

    def test_reference_reads_are_exact_substrings_of_the_reference(self) -> None:
        """If they were not, a mapped-read assertion would prove nothing."""
        chromosome = self.reference["chr1"]
        for record in self.records:
            if not record.startswith("SYNTH_ALIGN_REF_"):
                continue
            self.assertIn(record.split("\t")[9], chromosome)

    def test_deletion_reads_span_the_deleted_interval(self) -> None:
        chromosome = self.reference["chr1"]
        deletion_records = [r for r in self.records if r.startswith("SYNTH_ALIGN_DEL_")]
        self.assertEqual(len(deletion_records), DELETION_READS)
        for record in deletion_records:
            sequence = record.split("\t")[9]
            self.assertEqual(len(sequence), READ_ARM_BP * 2)
            # The two arms exist in the reference; the joined read does not.
            self.assertIn(sequence[:READ_ARM_BP], chromosome)
            self.assertIn(sequence[READ_ARM_BP:], chromosome)
            self.assertNotIn(sequence, chromosome)

    def test_the_deleted_bases_are_the_ones_that_were_skipped(self) -> None:
        chromosome = self.reference["chr1"]
        record = next(r for r in self.records if r.startswith("SYNTH_ALIGN_DEL_001"))
        sequence = record.split("\t")[9]
        left = chromosome[READ_START : READ_START + READ_ARM_BP]
        right_start = READ_START + READ_ARM_BP + DELETION_BP
        right = chromosome[right_start : right_start + READ_ARM_BP]
        self.assertEqual(sequence, left + right)

    def test_reverse_reads_are_reverse_complements_of_reference_sequence(self) -> None:
        chromosome = self.reference["chr1"]
        reverse_records = [r for r in self.records if r.startswith("SYNTH_ALIGN_REV_")]
        self.assertEqual(len(reverse_records), REVERSE_READS)
        for record in reverse_records:
            self.assertIn(reverse_complement(record.split("\t")[9]), chromosome)

    def test_read_names_are_unique(self) -> None:
        names = [record.split("\t")[0] for record in self.records]
        self.assertEqual(len(names), len(set(names)))

    def test_the_fixture_is_byte_identical_across_calls(self) -> None:
        self.assertEqual(self.text, unaligned_sam_text())


if __name__ == "__main__":
    unittest.main()

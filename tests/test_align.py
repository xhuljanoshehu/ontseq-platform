from __future__ import annotations

import unittest
from pathlib import Path

from ontseq_platform.align import (
    AlignmentPolicy,
    build_minimap2_argv,
    header_with_read_groups,
    read_group_lines,
)


def _policy(**overrides: object) -> AlignmentPolicy:
    values: dict[str, object] = {
        "profile_id": "test_profile",
        "status": "technical_defaults_only",
        "note": "test",
    }
    values.update(overrides)
    return AlignmentPolicy(**values)  # type: ignore[arg-type]


class PolicyTests(unittest.TestCase):
    def test_hard_clipping_is_refused(self) -> None:
        """SV callers need the clipped sequence, so the policy cannot opt out."""
        with self.assertRaises(ValueError):
            _policy(soft_clip_supplementary=False)

    def test_the_version_lock_must_look_like_a_version(self) -> None:
        """A free-text lock would make the fail-closed version gate unenforceable."""
        with self.assertRaises(ValueError):
            _policy(expected_minimap2_version="latest")


class ArgvTests(unittest.TestCase):
    def _argv(self, policy: AlignmentPolicy) -> list[str]:
        return build_minimap2_argv(
            minimap2="minimap2",
            policy=policy,
            reference_fasta=Path("/ref.fa"),
            reads_fastq=Path("/reads.fastq"),
            output_sam=Path("/out.sam"),
            threads=3,
        )

    def test_the_ont_preset_and_thread_count_are_explicit(self) -> None:
        argv = self._argv(_policy())
        self.assertEqual(argv[:6], ["minimap2", "-a", "-x", "map-ont", "-t", "3"])

    def test_tag_carrying_flags_are_present_by_default(self) -> None:
        argv = self._argv(_policy())
        for flag in ("--MD", "-Y", "-y"):
            self.assertIn(flag, argv)

    def test_dropping_modified_bases_drops_only_that_flag(self) -> None:
        argv = self._argv(_policy(preserve_modified_base_tags=False))
        self.assertNotIn("-y", argv)
        self.assertIn("-Y", argv)

    def test_reference_precedes_reads(self) -> None:
        """minimap2 is positional; swapping these silently aligns the reference."""
        argv = self._argv(_policy())
        self.assertEqual(argv[-2:], ["/ref.fa", "/reads.fastq"])


class ReadGroupHeaderTests(unittest.TestCase):
    HEADER = (
        "@HD\tVN:1.6\tSO:unknown\n"
        "@RG\tID:RG1\tSM:S1\tPL:ONT\n"
        "@RG\tID:RG2\tSM:S1\tPL:ONT\n"
        "@PG\tID:dorado\tPN:dorado\n"
    )

    def test_every_read_group_is_extracted_in_order(self) -> None:
        self.assertEqual(
            read_group_lines(self.HEADER),
            ["@RG\tID:RG1\tSM:S1\tPL:ONT", "@RG\tID:RG2\tSM:S1\tPL:ONT"],
        )

    def test_a_header_without_read_groups_yields_nothing(self) -> None:
        self.assertEqual(read_group_lines("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n"), [])

    def test_a_read_name_starting_with_rg_is_not_mistaken_for_a_header(self) -> None:
        self.assertEqual(read_group_lines("@RGX\tID:no\n"), [])

    def test_read_groups_land_after_the_sequence_dictionary(self) -> None:
        aligned = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n@PG\tID:minimap2\n"
        merged = header_with_read_groups(aligned, read_group_lines(self.HEADER)).splitlines()
        self.assertEqual(
            merged,
            [
                "@HD\tVN:1.6\tSO:coordinate",
                "@SQ\tSN:chr1\tLN:10",
                "@RG\tID:RG1\tSM:S1\tPL:ONT",
                "@RG\tID:RG2\tSM:S1\tPL:ONT",
                "@PG\tID:minimap2",
            ],
        )

    def test_all_sequence_records_stay_above_the_read_groups(self) -> None:
        aligned = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n@SQ\tSN:chr2\tLN:20\n@PG\tID:minimap2\n"
        merged = header_with_read_groups(aligned, ["@RG\tID:RG1"]).splitlines()
        self.assertEqual(merged.index("@RG\tID:RG1"), 3)

    def test_reapplying_the_same_read_groups_is_idempotent(self) -> None:
        aligned = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n@PG\tID:minimap2\n"
        groups = read_group_lines(self.HEADER)
        once = header_with_read_groups(aligned, groups)
        self.assertEqual(header_with_read_groups(once, groups), once)

    def test_a_header_with_no_sequence_records_still_places_read_groups(self) -> None:
        merged = header_with_read_groups("@PG\tID:minimap2\n", ["@RG\tID:RG1"]).splitlines()
        self.assertEqual(merged, ["@RG\tID:RG1", "@PG\tID:minimap2"])


if __name__ == "__main__":
    unittest.main()

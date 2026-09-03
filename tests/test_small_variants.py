from __future__ import annotations

import unittest
from pathlib import Path

from ontseq_platform.small_variants import (
    COMPLETED,
    DELETION,
    INSERTION,
    MNV,
    NO_CALL,
    REJECT_BELOW_DETECTION_FLOOR,
    REJECT_FEW_VARIANT_READS,
    REJECT_LOW_DEPTH,
    REJECT_LOW_QUALITY,
    REJECT_NOT_PASS,
    REJECT_UNRESOLVABLE_DEPTH,
    SNV,
    AcceptedCall,
    Clair3Policy,
    RejectedCall,
    SmallVariant,
    SmallVariantError,
    SomaticStatus,
    apply_policy,
    evaluate,
)

POLICY = Clair3Policy(profile_id="test-v1", expected_version="2.0.2")


def _variant(
    *,
    reference: str = "A",
    alternate: str = "G",
    depth: int = 80,
    variant_reads: int = 30,
    quality: float = 30.0,
    filter_status: str = "PASS",
) -> SmallVariant:
    return SmallVariant(
        chrom="chr5",
        position=171_390_000,
        reference=reference,
        alternate=alternate,
        depth=depth,
        variant_reads=variant_reads,
        quality=quality,
        filter_status=filter_status,
    )


class SomaticBoundaryTests(unittest.TestCase):
    def test_there_is_exactly_one_somatic_status(self) -> None:
        """The boundary is structural, not conventional.

        Adding a member is a decision about the assay -- a matched normal, a
        population-frequency filter, or a tumour-only caller -- not a refactor.
        """
        self.assertEqual(len(list(SomaticStatus)), 1)
        self.assertIs(SomaticStatus.NOT_DETERMINED, next(iter(SomaticStatus)))

    def test_an_accepted_call_is_never_somatic(self) -> None:
        outcome = evaluate(_variant(), POLICY)
        assert isinstance(outcome, AcceptedCall)
        self.assertIs(outcome.somatic_status, SomaticStatus.NOT_DETERMINED)

    def test_an_accepted_call_is_never_reportable(self) -> None:
        outcome = evaluate(_variant(), POLICY)
        assert isinstance(outcome, AcceptedCall)
        self.assertFalse(outcome.reportable)

    def test_the_caveat_says_why_somatic_status_is_undetermined(self) -> None:
        outcome = evaluate(_variant(), POLICY)
        assert isinstance(outcome, AcceptedCall)
        caveat = outcome.caveat()
        self.assertIn("matched normal", caveat)
        self.assertIn("germline", caveat)

    def test_a_profile_claiming_a_somatic_vocabulary_is_refused(self) -> None:
        with self.assertRaises(SmallVariantError):
            Clair3Policy.from_mapping(
                {
                    "profile_id": "wrong",
                    "expected_version": "2.0.2",
                    "caller_vocabulary": "somatic",
                }
            )


class VariantClassTests(unittest.TestCase):
    def test_substitution_is_an_snv(self) -> None:
        self.assertEqual(_variant(reference="A", alternate="G").variant_class, SNV)

    def test_equal_length_multi_base_is_an_mnv(self) -> None:
        self.assertEqual(_variant(reference="AC", alternate="GT").variant_class, MNV)

    def test_npm1_style_four_base_insertion(self) -> None:
        """The canonical NPM1 alteration is a 4 bp insertion, so it must classify as one."""
        variant = _variant(reference="C", alternate="CTCTG")
        self.assertEqual(variant.variant_class, INSERTION)
        self.assertEqual(variant.length_change, 4)
        self.assertTrue(variant.is_indel)

    def test_deletion(self) -> None:
        variant = _variant(reference="CTCTG", alternate="C")
        self.assertEqual(variant.variant_class, DELETION)
        self.assertEqual(variant.length_change, -4)
        self.assertTrue(variant.is_indel)

    def test_an_snv_is_not_an_indel(self) -> None:
        self.assertFalse(_variant().is_indel)


class IndelConfirmationTests(unittest.TestCase):
    def test_an_indel_carries_the_confirmation_flag(self) -> None:
        outcome = evaluate(_variant(reference="C", alternate="CTCTG"), POLICY)
        assert isinstance(outcome, AcceptedCall)
        self.assertTrue(outcome.requires_orthogonal_confirmation)
        self.assertIn("orthogonal confirmation", outcome.caveat())

    def test_an_snv_does_not(self) -> None:
        outcome = evaluate(_variant(), POLICY)
        assert isinstance(outcome, AcceptedCall)
        self.assertFalse(outcome.requires_orthogonal_confirmation)

    def test_the_flag_can_be_switched_off_by_policy(self) -> None:
        policy = Clair3Policy(
            profile_id="t", expected_version="2.0.2", indels_require_orthogonal_confirmation=False
        )
        outcome = evaluate(_variant(reference="C", alternate="CTCTG"), policy)
        assert isinstance(outcome, AcceptedCall)
        self.assertFalse(outcome.requires_orthogonal_confirmation)


class DepthGatingTests(unittest.TestCase):
    def test_the_floor_is_computed_per_variant_not_per_run(self) -> None:
        """An on-target and an off-target variant in one BAM get different floors."""
        on_target = evaluate(_variant(depth=80, variant_reads=16), POLICY)
        assert isinstance(on_target, AcceptedCall)
        deep = evaluate(_variant(depth=400, variant_reads=80), POLICY)
        assert isinstance(deep, AcceptedCall)
        self.assertLess(deep.detection_floor, on_target.detection_floor)

    def test_an_allele_fraction_below_its_own_floor_is_rejected(self) -> None:
        outcome = evaluate(_variant(depth=80, variant_reads=5), POLICY)
        assert isinstance(outcome, RejectedCall)
        self.assertEqual(outcome.reason, REJECT_BELOW_DETECTION_FLOOR)

    def test_the_same_allele_fraction_passes_at_greater_depth(self) -> None:
        """6.25% is below the floor at 80x and above it at 800x. Depth is the difference."""
        shallow = evaluate(_variant(depth=80, variant_reads=5), POLICY)
        deep = evaluate(_variant(depth=800, variant_reads=50), POLICY)
        self.assertIsInstance(shallow, RejectedCall)
        self.assertIsInstance(deep, AcceptedCall)

    def test_a_depth_resolving_nothing_is_named_as_such(self) -> None:
        """Reached when the error rate swallows the depth, not at low depth alone.

        At 1x and a 1% error rate an allele fraction of 95% is still resolvable, so shallow
        coverage on its own does not trigger this. It takes an error rate high enough that no
        read count separates signal from noise -- a defensive branch under realistic
        parameters, and the honest answer when a caller is handed an error model that bad.
        """
        policy = Clair3Policy(
            profile_id="t",
            expected_version="2.0.2",
            min_depth=1,
            min_variant_reads=1,
            error_rate=0.5,
        )
        outcome = evaluate(_variant(depth=2, variant_reads=2), policy)
        assert isinstance(outcome, RejectedCall)
        self.assertEqual(outcome.reason, REJECT_UNRESOLVABLE_DEPTH)

    def test_one_read_still_resolves_a_very_high_allele_fraction(self) -> None:
        """The companion fact, pinned so the branch above is not misread as "shallow"."""
        policy = Clair3Policy(
            profile_id="t", expected_version="2.0.2", min_depth=1, min_variant_reads=1
        )
        outcome = evaluate(_variant(depth=1, variant_reads=1), policy)
        assert isinstance(outcome, AcceptedCall)
        self.assertGreater(outcome.detection_floor, 0.9)


class PolicyRejectionTests(unittest.TestCase):
    def test_non_pass_records_are_rejected_when_pass_only(self) -> None:
        outcome = evaluate(_variant(filter_status="LowQual"), POLICY)
        assert isinstance(outcome, RejectedCall)
        self.assertEqual(outcome.reason, REJECT_NOT_PASS)

    def test_low_quality_is_rejected(self) -> None:
        outcome = evaluate(_variant(quality=2.0), POLICY)
        assert isinstance(outcome, RejectedCall)
        self.assertEqual(outcome.reason, REJECT_LOW_QUALITY)

    def test_low_depth_is_rejected(self) -> None:
        outcome = evaluate(_variant(depth=9, variant_reads=4), POLICY)
        assert isinstance(outcome, RejectedCall)
        self.assertEqual(outcome.reason, REJECT_LOW_DEPTH)

    def test_too_few_variant_reads_is_rejected(self) -> None:
        outcome = evaluate(_variant(depth=80, variant_reads=2), POLICY)
        assert isinstance(outcome, RejectedCall)
        self.assertEqual(outcome.reason, REJECT_FEW_VARIANT_READS)


class MalformedRecordTests(unittest.TestCase):
    def test_zero_depth_is_refused(self) -> None:
        with self.assertRaises(SmallVariantError):
            _variant(depth=0)

    def test_more_variant_reads_than_depth_is_refused(self) -> None:
        with self.assertRaises(SmallVariantError):
            _variant(depth=10, variant_reads=11)

    def test_missing_allele_is_refused(self) -> None:
        with self.assertRaises(SmallVariantError):
            _variant(alternate="")

    def test_zero_position_is_refused(self) -> None:
        with self.assertRaises(SmallVariantError):
            SmallVariant(
                chrom="chr5",
                position=0,
                reference="A",
                alternate="G",
                depth=80,
                variant_reads=30,
                quality=30.0,
                filter_status="PASS",
            )


class CallSetTests(unittest.TestCase):
    def test_an_empty_accepted_set_is_no_call_not_completed(self) -> None:
        """The distinction ADR-007 already draws for Sniffles2.

        An empty result is not a negative finding, and COMPLETED with nothing in it reads as
        one.
        """
        call_set = apply_policy([_variant(filter_status="LowQual")], POLICY)
        self.assertEqual(call_set.status, NO_CALL)
        self.assertEqual(call_set.accepted, ())

    def test_a_surviving_call_gives_completed(self) -> None:
        self.assertEqual(apply_policy([_variant()], POLICY).status, COMPLETED)

    def test_rejections_are_counted_by_reason_not_discarded(self) -> None:
        call_set = apply_policy(
            [
                _variant(filter_status="LowQual"),
                _variant(filter_status="LowQual"),
                _variant(quality=1.0),
                _variant(),
            ],
            POLICY,
        )
        self.assertEqual(call_set.status, COMPLETED)
        self.assertEqual(call_set.rejection_counts(), {REJECT_NOT_PASS: 2, REJECT_LOW_QUALITY: 1})

    def test_indels_needing_confirmation_are_listed_separately(self) -> None:
        call_set = apply_policy([_variant(), _variant(reference="C", alternate="CTCTG")], POLICY)
        self.assertEqual(len(call_set.accepted), 2)
        self.assertEqual(len(call_set.indels_needing_confirmation), 1)

    def test_no_accepted_call_in_a_set_is_ever_reportable(self) -> None:
        call_set = apply_policy(
            [_variant(), _variant(reference="C", alternate="CTCTG"), _variant(quality=1.0)],
            POLICY,
        )
        self.assertTrue(all(not call.reportable for call in call_set.accepted))


class PolicyLoadingTests(unittest.TestCase):
    def test_round_trip_from_a_mapping(self) -> None:
        policy = Clair3Policy.from_mapping(
            {
                "profile_id": "clair3-tumour-only-technical-v1",
                "expected_version": "2.0.2",
                "caller_vocabulary": "germline",
                "min_depth": 20,
                "error_rate": 0.02,
            }
        )
        self.assertEqual(policy.expected_version, "2.0.2")
        self.assertEqual(policy.min_depth, 20)
        self.assertAlmostEqual(policy.error_rate, 0.02)

    def test_a_profile_without_a_pinned_version_is_refused(self) -> None:
        with self.assertRaises(SmallVariantError):
            Clair3Policy.from_mapping({"profile_id": "unpinned"})

    def test_the_shipped_profile_loads_and_pins_the_bioconda_version(self) -> None:
        """The drift guard: the config on disk must satisfy the policy that reads it.

        Skipped where yaml is unavailable, which is why this module parses a mapping rather
        than a path. Clair3 2.0.2 was verified present on the live Bioconda linux-64 index on
        2026-09-02; ClairS-TO was verified absent from it on the same day.
        """
        try:
            import yaml
        except ImportError as error:  # pragma: no cover - only where yaml is absent
            self.skipTest(f"yaml unavailable: {error}")
        path = Path("configs/variants/clair3.technical.yaml")
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
        policy = Clair3Policy.from_mapping(raw)
        self.assertEqual(policy.expected_version, "2.0.2")
        self.assertEqual(raw["status"], "technical_defaults_only")
        self.assertEqual(raw["caller_vocabulary"], "germline")
        self.assertTrue(policy.indels_require_orthogonal_confirmation)
        self.assertIsNone(
            policy.required_model_id, "no model is pinned yet; a run must fail closed"
        )


if __name__ == "__main__":
    unittest.main()

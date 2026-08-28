from __future__ import annotations

import unittest
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory

from ontseq_platform.bam_resolution import (
    default_run_id,
    locate_bam_index,
    resolve_bam_header,
    sample_id_from_bam,
)
from ontseq_platform.models import (
    GenomeBuild,
    ReferenceContig,
    ReferenceDictionaryContract,
    ReferenceLock,
)
from ontseq_platform.reference import canonical_contigs, grch38_canonical_25_contigs


def _lock(build: GenomeBuild, contigs: tuple[tuple[str, int], ...]) -> ReferenceLock:
    return ReferenceLock(
        reference_id=f"{build.value}_TEST",
        genome_build=build,
        contigs=[ReferenceContig(name=name, length=length) for name, length in contigs],
        source_fai_sha256="0" * 64,
    )


def _header(contigs: tuple[tuple[str, int], ...]) -> str:
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines.extend(f"@SQ\tSN:{name}\tLN:{length}" for name, length in contigs)
    return "\n".join(lines) + "\n"


class BamResolutionTests(unittest.TestCase):
    def test_bam_dot_bai_precedes_short_bai(self) -> None:
        with TemporaryDirectory() as raw:
            bam = Path(raw) / "sample.bam"
            preferred = Path(f"{bam}.bai")
            alternate = bam.with_suffix(".bai")
            for path in (bam, preferred, alternate):
                path.write_bytes(b"x")
            self.assertEqual(locate_bam_index(bam), preferred.resolve())

    def test_grch38_dictionary_resolves_exactly(self) -> None:
        contigs = canonical_contigs(GenomeBuild.GRCH38)
        with TemporaryDirectory() as raw:
            bam = Path(raw) / "AML 001.bam"
            bam.write_bytes(b"BAM")
            Path(f"{bam}.bai").write_bytes(b"BAI")
            resolved = resolve_bam_header(
                bam_path=bam,
                header_text=_header(contigs),
                reference_lock=_lock(GenomeBuild.GRCH38, contigs),
            )
            self.assertEqual(resolved.genome_build, GenomeBuild.GRCH38)
            self.assertEqual(resolved.naming_style, "chr-prefixed")
            self.assertEqual(resolved.sample_id, "AML_001")

    def test_grch37_bam_is_rejected_before_analysis(self) -> None:
        grch37 = canonical_contigs(GenomeBuild.GRCH37)
        with TemporaryDirectory() as raw:
            bam = Path(raw) / "sample.bam"
            bam.write_bytes(b"BAM")
            Path(f"{bam}.bai").write_bytes(b"BAI")
            with self.assertRaisesRegex(ValueError, "profile requires GRCh38"):
                resolve_bam_header(
                    bam_path=bam,
                    header_text=_header(grch37),
                    reference_lock=_lock(GenomeBuild.GRCH38, canonical_contigs(GenomeBuild.GRCH38)),
                )

    def test_partial_and_mixed_style_dictionaries_are_rejected(self) -> None:
        full = canonical_contigs(GenomeBuild.GRCH38)
        dictionaries = (full[:-1], (full[0], ("2", full[1][1]), *full[2:]))
        for contigs in dictionaries:
            with self.subTest(contigs=len(contigs)), TemporaryDirectory() as raw:
                bam = Path(raw) / "sample.bam"
                bam.write_bytes(b"BAM")
                Path(f"{bam}.bai").write_bytes(b"BAI")
                with self.assertRaisesRegex(ValueError, "partial, mixed-style"):
                    resolve_bam_header(
                        bam_path=bam,
                        header_text=_header(contigs),
                        reference_lock=_lock(GenomeBuild.GRCH38, full),
                    )

    def test_full_dictionary_must_match_lock_order_and_optional_contigs(self) -> None:
        canonical = canonical_contigs(GenomeBuild.GRCH38)
        lock_contigs = (*canonical, ("chrM", 16569))
        with TemporaryDirectory() as raw:
            bam = Path(raw) / "sample.bam"
            bam.write_bytes(b"BAM")
            Path(f"{bam}.bai").write_bytes(b"BAI")
            with self.assertRaisesRegex(ValueError, "exactly match"):
                resolve_bam_header(
                    bam_path=bam,
                    header_text=_header(canonical),
                    reference_lock=_lock(GenomeBuild.GRCH38, lock_contigs),
                )

    def test_explicit_canonical_25_profile_accepts_real_world_dictionary_shape(self) -> None:
        canonical = grch38_canonical_25_contigs()
        pinned_full = (*canonical, ("GL000008.2", 209709))
        with TemporaryDirectory() as raw:
            bam = Path(raw) / "sample.bam"
            bam.write_bytes(b"BAM")
            Path(f"{bam}.bai").write_bytes(b"BAI")

            resolved = resolve_bam_header(
                bam_path=bam,
                header_text=_header(canonical),
                reference_lock=_lock(GenomeBuild.GRCH38, pinned_full),
                dictionary_contract=ReferenceDictionaryContract.GRCH38_CANONICAL_25,
            )

        self.assertEqual(
            resolved.reference_dictionary_contract,
            ReferenceDictionaryContract.GRCH38_CANONICAL_25,
        )
        self.assertEqual(resolved.contigs, canonical)

    def test_dictionary_contracts_never_fallback_between_25_and_full(self) -> None:
        canonical = grch38_canonical_25_contigs()
        pinned_full = (*canonical, ("GL000008.2", 209709))
        cases = (
            (
                ReferenceDictionaryContract.EXACT_FULL,
                canonical,
                "exactly match",
            ),
            (
                ReferenceDictionaryContract.GRCH38_CANONICAL_25,
                pinned_full,
                "exactly match",
            ),
            (
                ReferenceDictionaryContract.GRCH38_CANONICAL_25,
                canonical[:-1],
                "exactly match",
            ),
        )
        for contract, observed, message in cases:
            with (
                self.subTest(contract=contract, observed=len(observed)),
                TemporaryDirectory() as raw,
            ):
                bam = Path(raw) / "sample.bam"
                bam.write_bytes(b"BAM")
                Path(f"{bam}.bai").write_bytes(b"BAI")
                with self.assertRaisesRegex(ValueError, message):
                    resolve_bam_header(
                        bam_path=bam,
                        header_text=_header(observed),
                        reference_lock=_lock(GenomeBuild.GRCH38, pinned_full),
                        dictionary_contract=contract,
                    )

    def test_filename_and_default_run_id_are_manifest_safe(self) -> None:
        self.assertEqual(sample_id_from_bam(Path("a.bam")), "sample_a")
        moment = datetime(2026, 8, 27, 12, 34, 56, tzinfo=UTC)
        self.assertEqual(
            default_run_id("AML_001", now=moment),
            "AML_001-20260827T123456Z",
        )


if __name__ == "__main__":
    unittest.main()

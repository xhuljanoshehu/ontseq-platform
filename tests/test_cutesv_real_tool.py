"""Opt-in synthetic equivalence check for cuteSV's worker count."""

import os
import random
import shutil
from pathlib import Path

import pytest

from ontseq_platform.cutesv import run_cutesv
from ontseq_platform.models import (
    AlignedBamIntakeReport,
    AnalysisSpec,
    AssayMode,
    AssaySpec,
    CuteSvPolicy,
    GenomeBuild,
    InputKind,
    InputSpec,
    SampleManifest,
    Verdict,
)


@pytest.mark.skipif(os.environ.get("ONTSEQ_CUTESV_REAL_TOOL") != "1", reason="Opt-in real cuteSV")
@pytest.mark.parametrize("first_start", [97_000, 90_909])
def test_one_and_four_workers_recover_the_same_synthetic_deletion(
    tmp_path: Path, first_start: int
) -> None:
    import pysam

    executable = shutil.which("cuteSV")
    assert executable is not None, "Requested real-tool test requires cuteSV"
    rng = random.Random(739)
    from ontseq_platform.cutesv_build import BUILD_ID, executable_identity

    original_identity = executable_identity(executable)
    reference = "".join(rng.choices("ACGT", k=200_000))
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + reference + "\n")
    pysam.faidx(str(fasta))
    bam = tmp_path / "synthetic.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(reference)}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i in range(80):
            start = first_start + i * 10
            deleted = i % 2 == 0
            deletion_start = 100_000 if first_start == 97_000 else 93_000
            left = deletion_start - start
            read = pysam.AlignedSegment()
            read.query_name = f"SYNTHETIC_{i}"
            read.reference_id = 0
            read.reference_start = start
            read.flag = 0
            read.mapping_quality = 60
            read.cigartuples = [(0, left), (2, 300), (0, 5000 - left)] if deleted else [(0, 5000)]
            read.query_sequence = (
                (
                    reference[start:deletion_start]
                    + reference[deletion_start + 300 : deletion_start + 300 + 5000 - left]
                )
                if deleted
                else reference[start : start + 5000]
            )
            read.query_qualities = pysam.qualitystring_to_array("I" * 5000)
            out.write(read)
    pysam.index(str(bam))
    manifest = SampleManifest(
        sample_id="SYNTHETIC_CUTESV_THREADS",
        run_id="SYNTHETIC_THREADS_RUN",
        input=InputSpec(kind=InputKind.ALIGNED_BAM, path=str(bam), index_path=str(bam) + ".bai"),
        assay=AssaySpec(
            mode=AssayMode.LOW_COVERAGE_WGS,
            genome_build=GenomeBuild.GRCH38,
            reference_id="SYNTHETIC_REF",
        ),
        analysis=AnalysisSpec(profile="synthetic", modules=["sv"]),
    )
    intake = AlignedBamIntakeReport(
        sample_id=manifest.sample_id,
        reference_id="SYNTHETIC_REF",
        genome_build=GenomeBuild.GRCH38,
        checks=[],
        verdict=Verdict.PASS,
    )
    policy = CuteSvPolicy(
        profile_id="synthetic",
        status="technical_defaults_only",
        note="Synthetic equivalence check only",
    )
    reports = [
        run_cutesv(
            manifest,
            intake,
            policy,
            reference_fasta=fasta,
            output_vcf=tmp_path / f"threads-{threads}.vcf",
            cutesv=executable,
            threads=threads,
        )
        for threads in (1, 4)
    ]
    assert len(reports[0].events) == 1
    event = reports[0].events[0]
    assert event.primary.chromosome == "chr1"
    # Verify the literal synthetic truth in the caller VCF. The existing ONTSeq
    # normalizer uses an anchor-based locus; this concurrency test does not change
    # or qualify that coordinate convention.
    for threads in (1, 4):
        records = [
            line.split("\t")
            for line in (tmp_path / f"threads-{threads}.vcf").read_text().splitlines()
            if not line.startswith("#")
        ]
        assert len(records) == 1
        assert records[0][1] == ("100000" if first_start == 97_000 else "93000")
        info = dict(item.split("=", 1) for item in records[0][7].split(";") if "=" in item)
        assert info["END"] == ("100300" if first_start == 97_000 else "93300")
        assert info["SVLEN"] == "-300"
        assert info["RE"] == "40"
    assert reports[0].events == reports[1].events
    assert [r.tool.parameters["threads"] for r in reports] == [1, 4]

    assert executable_identity(executable) == original_identity
    assert all(r.tool.parameters["cutesv_build_id"] == BUILD_ID for r in reports)

"""Independent literal oracles for the reviewed modkit scanner, using synthetic reads."""

from __future__ import annotations

import array
import hashlib
import importlib
import os
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from ontseq_platform.methylation import (
    MethylationPolicy,
    _modkit_failed_processing_count,
    run_methylation,
)
from ontseq_platform.models import AlignedBamIntakeReport, ModuleRunStatus, SampleManifest

pytestmark = pytest.mark.skipif(
    os.environ.get("ONTSEQ_MODKIT_PR709_REAL_TOOL") != "1",
    reason="explicit modkit PR709 candidate qualification lane",
)

# Each site: valid calls, 5mC, 5hmC, canonical calls, filtered calls.
# AACGCGTT is reverse-complement symmetric; reversing the read exchanges sites 2 and 4.
CASES = [
    (
        "independent-forward",
        "C+m?,0;C+h?,1;",
        [255, 255],
        [0] * 5,
        {2: (5, 5, 0, 0, 0), 4: (5, 0, 5, 0, 0)},
    ),
    (
        "independent-reverse",
        "C+m?,0;C+h?,1;",
        [255, 255],
        [16] * 5,
        {2: (5, 0, 5, 0, 0), 4: (5, 5, 0, 0, 0)},
    ),
    (
        "independent-both",
        "C+m?,0;C+h?,1;",
        [255, 255],
        [0] * 5 + [16] * 5,
        {2: (10, 5, 5, 0, 0), 4: (10, 5, 5, 0, 0)},
    ),
    (
        "independent-implicit",
        "C+m.,0;C+h.,1;",
        [255, 255],
        [0] * 5,
        {2: (5, 5, 0, 0, 0), 4: (5, 0, 5, 0, 0)},
    ),
    (
        "combined-forward",
        "C+mh?,0,0;",
        [255, 0, 0, 255],
        [0] * 5,
        {2: (5, 5, 0, 0, 0), 4: (5, 0, 5, 0, 0)},
    ),
    (
        "combined-reverse",
        "C+mh?,0,0;",
        [255, 0, 0, 255],
        [16] * 5,
        {2: (5, 0, 5, 0, 0), 4: (5, 5, 0, 0, 0)},
    ),
    (
        "combined-both",
        "C+mh?,0,0;",
        [255, 0, 0, 255],
        [0] * 5 + [16] * 5,
        {2: (10, 5, 5, 0, 0), 4: (10, 5, 5, 0, 0)},
    ),
    (
        "measured-zero",
        "C+mh?,0,0;",
        [0, 0, 0, 0],
        [0] * 5 + [16] * 5,
        {2: (10, 0, 0, 10, 0), 4: (10, 0, 0, 10, 0)},
    ),
    (
        "low-confidence",
        "C+m?,0;C+h?,1;",
        [150, 250],
        [0] * 5,
        {2: (0, 0, 0, 0, 5), 4: (5, 0, 5, 0, 0)},
    ),
    (
        "combined-low-confidence",
        "C+mh?,0,0;",
        [150, 0, 0, 250],
        [0] * 5,
        {2: (0, 0, 0, 0, 5), 4: (5, 0, 5, 0, 0)},
    ),
    (
        "only-5mc",
        "C+m?,0,0;",
        [255, 255],
        [0] * 5 + [16] * 5,
        {2: (10, 10, 0, 0, 0), 4: (10, 10, 0, 0, 0)},
    ),
    ("low-depth", "C+m?,0;C+h?,1;", [255, 255], [0] * 2, {2: (2, 2, 0, 0, 0), 4: (2, 0, 2, 0, 0)}),
]


def _inputs(tmp_path, case, missing, clipped):
    pysam = importlib.import_module("pysam")
    _, mm, probabilities, flags, expected = case
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nAACGCGTT\n", encoding="utf-8")
    pysam.faidx(str(reference))
    bam = tmp_path / "synthetic.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 8}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i, flag in enumerate([*flags, *([0] * missing)]):
            read = pysam.AlignedSegment()
            read.query_name = f"SYNTHETIC_{i}"
            read.query_sequence = "TTAACGCGTTAA" if clipped else "AACGCGTT"
            read.query_qualities = array.array("B", [40] * len(read.query_sequence))
            read.flag, read.reference_id, read.reference_start = flag, 0, 0
            read.mapping_quality = 60
            read.cigarstring = "2S8M2S" if clipped else "8M"
            if i < len(flags):
                read.set_tag("MN", len(read.query_sequence))
                read.set_tag("MM", mm)
                read.set_tag("ML", array.array("B", probabilities))
            out.write(read)
    pysam.index(str(bam))
    return bam, reference


@pytest.mark.parametrize("case", CASES, ids=[case[0] for case in CASES])
@pytest.mark.parametrize("missing,clipped", [(0, False), (3, False), (3, True)])
def test_exact_cpg_positions_counts_and_probabilities(tmp_path, case, missing, clipped):
    bam, reference = _inputs(tmp_path, case, missing, clipped)
    expected = case[-1]
    binary = os.environ["ONTSEQ_MODKIT_PR709_BINARY"]
    output, log = tmp_path / "output.bed", tmp_path / "pileup.log"
    command = [
        binary,
        "pileup",
        str(bam),
        str(output),
        "--ref",
        str(reference),
        "--modified-bases",
        "5mC",
        "5hmC",
        "--cpg",
        "--combine-strands",
        "--filter-threshold",
        "0.8",
        "--threads",
        "1",
        "--io-threads",
        "1",
        "--suppress-progress",
        "--log-filepath",
        str(log),
    ]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60, check=False)
    assert result.returncode == 0, result.stderr
    assert output.is_file(), result.stderr
    # Untagged reads are refused by this pinned tool; partial counts are never released.
    assert _modkit_failed_processing_count(log, result.stdout, result.stderr) == missing
    rows = {}
    for line in output.read_text().splitlines():
        fields = line.split()
        assert fields[0] == "chr1" and fields[5] == "."
        position, code = int(fields[1]), fields[3]
        assert int(fields[2]) == position + 1
        assert position in expected and code in ("m", "h")
        assert (position, code) not in rows
        rows[position, code] = [int(fields[9]), *map(int, fields[11:18])]
    wanted = {}
    for position, (valid, mc, hmc, canonical, failed) in expected.items():
        # The upstream writer omits sites with zero valid calls (into_counts).
        # A combined-tag control checks this independently of the PR709 scanner.
        # Their absence is unmeasured, never a reported zero methylation fraction.
        if valid == 0:
            continue
        wanted[position, "m"] = [valid, mc, canonical, hmc, 0, failed, 0, 0]
        wanted[position, "h"] = [valid, hmc, canonical, mc, 0, failed, 0, 0]
    assert rows == wanted


def test_original_two_base_counterexample(tmp_path: Path) -> None:
    pysam = importlib.import_module("pysam")
    sam = tmp_path / "synthetic.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2\n"
        "SYNTHETIC\t0\tchr1\t1\t60\t2M\t*\t0\t0\tCC\t??\t"
        "MM:Z:C+m?,0;C+h?,1;\tML:B:C,200,250\tMN:i:2\tNM:i:0\n",
        encoding="utf-8",
    )
    bam = tmp_path / "synthetic.bam"
    pysam.view("-b", "-o", str(bam), str(sam), catch_stdout=False)
    pysam.index(str(bam))
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nCC\n", encoding="utf-8")
    pysam.faidx(str(reference))
    output = tmp_path / "output.bed"
    result = subprocess.run(
        [
            os.environ["ONTSEQ_MODKIT_PR709_BINARY"],
            "pileup",
            str(bam),
            str(output),
            "--ref",
            str(reference),
            "--modified-bases",
            "C:m",
            "C:h",
            "--no-filtering",
            "--threads",
            "1",
            "--io-threads",
            "1",
            "--suppress-progress",
        ],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    positive = sorted(
        (int(row[1]), row[3], int(row[11]))
        for row in (line.split() for line in output.read_text().splitlines())
        if int(row[11]) > 0
    )
    assert positive == [(0, "m", 1), (1, "h", 1)]


@pytest.mark.parametrize("case", CASES, ids=[case[0] for case in CASES])
@pytest.mark.parametrize("missing,clipped", [(0, False), (3, False), (3, True)])
def test_actual_adapter_acceptance_and_refusal(tmp_path, case, missing, clipped):
    bam, reference = _inputs(tmp_path, case, missing, clipped)
    manifest = SampleManifest.model_validate(
        {
            "sample_id": "SYNTHETIC_PR709",
            "run_id": "SYNTHETIC_RUN",
            "input": {"kind": "aligned_bam", "path": str(bam), "index_path": str(bam) + ".bai"},
            "assay": {
                "mode": "lcwgs",
                "genome_build": "GRCh38",
                "reference_id": "synthetic-pr709",
            },
            "analysis": {"profile": "synthetic", "modules": ["methylation"]},
        }
    )
    intake = AlignedBamIntakeReport(
        sample_id=manifest.sample_id,
        reference_id=manifest.assay.reference_id,
        genome_build=manifest.assay.genome_build,
        checks=[],
        verdict="PASS",
    )
    policy = MethylationPolicy(
        profile_id="SYNTHETIC_PR709",
        expected_version="0.6.4",
        status="technical_defaults_only",
        modification_codes=["m", "h"],
        minimum_valid_coverage=5,
        cpg_only=True,
        combine_strands=True,
        filter_threshold=0.8,
        region_source="chromosome",
        note="Synthetic implementation qualification; no clinical performance claim",
    )
    binary = os.environ["ONTSEQ_MODKIT_PR709_BINARY"]
    digest = hashlib.sha256(Path(binary).read_bytes()).hexdigest()
    # Exercise the complete adapter before this candidate is added to the production list.
    # This test-only override does not make arbitrary runtime binaries trusted.
    with patch("ontseq_platform.modkit_build.PR709_BINARY_SHA256", frozenset({digest})):
        if missing:
            with pytest.raises(ValueError, match=f"reported {missing} failed record"):
                run_methylation(
                    manifest,
                    intake,
                    policy,
                    output_dir=tmp_path / "adapter",
                    reference_fasta=reference,
                    modkit=binary,
                    threads=1,
                )
            assert not (tmp_path / "adapter/SYNTHETIC_PR709.modkit.bedmethyl").exists()
            return
        report = run_methylation(
            manifest,
            intake,
            policy,
            output_dir=tmp_path / "adapter",
            reference_fasta=reference,
            modkit=binary,
            threads=1,
        )
    assert report.tool.parameters["modkit_binary_sha256"] == digest
    assert report.tool.parameters["modkit_build_id"] == "ontseq-modkit-pr709-v1"
    assert report.reads_with_modified_base_tags == len(case[3])
    qualifying = [counts for counts in case[-1].values() if counts[0] >= 5]
    assert report.status == (ModuleRunStatus.COMPLETED if qualifying else ModuleRunStatus.NO_CALL)
    for region in report.regions:
        index = 1 if region.modification_code.value == "m" else 2
        expected_valid = sum(counts[0] for counts in qualifying)
        expected_modified = sum(counts[index] for counts in qualifying)
        assert region.valid_call_count == expected_valid
        assert region.modified_call_count == expected_modified
        assert region.mean_modified_fraction == (
            expected_modified / expected_valid if expected_valid else None
        )

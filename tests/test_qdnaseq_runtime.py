from __future__ import annotations

import json
from pathlib import Path

import pytest

from ontseq_platform.cnv.qdnaseq import QDNAseqPolicy, run_qdnaseq_ace
from ontseq_platform.execution import CommandResult
from ontseq_platform.models import EventType, GenomeBuild, ReferenceContig, ReferenceLock


class FakeQDNAseqRunner:
    def __init__(self, *, fail: bool = False, event_call: float = -1.0) -> None:
        self.fail = fail
        self.event_call = event_call

    def run(self, argv, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        args = [str(item) for item in argv]
        if self.fail:
            return CommandResult(tuple(args), 9, "", "synthetic R failure")
        out = Path(args[args.index("--output-dir") + 1])
        sample = args[args.index("--sample-id") + 1]
        out.mkdir(parents=True, exist_ok=True)
        runs = []
        for bin_size, cellularity, ploidy in (
            (100, 0.55, 2.0),
            (500, 0.57, 2.0),
            (1000, 0.56, 2.0),
        ):
            seg = f"{sample}.{bin_size}kbp.segments.tsv"
            chrom = f"{sample}.{bin_size}kbp.chromosomes.tsv"
            fit = f"{sample}.{bin_size}kbp.ace-fit.png"
            copy = f"{sample}.{bin_size}kbp.copy-number.png"
            rds = f"{sample}.{bin_size}kbp.segmented.rds"
            (out / seg).write_text(
                "chromosome\tstart\tend\tbin_count\tabsolute_copy_number\tcall\tqnorm_log10\n"
                f"chr7\t0\t950\t2\t1.0\t{self.event_call}\t-5.0\n",
                encoding="utf-8",
            )
            (out / chrom).write_text(
                "chromosome\tcopy_number\nchr7\t1.0\n",
                encoding="utf-8",
            )
            (out / fit).write_bytes(b"PNG-fit")
            (out / copy).write_bytes(b"PNG-copy")
            (out / rds).write_bytes(b"RDS")
            runs.append(
                {
                    "bin_size_kbp": bin_size,
                    "cellularity": cellularity,
                    "ploidy": ploidy,
                    "fit_error": 0.1 + bin_size / 10000,
                    "candidate_count": 3,
                    "alternatives": [
                        {"cellularity": cellularity, "ploidy": ploidy, "fit_error": 0.1}
                    ],
                    "segment_file": seg,
                    "chromosome_file": chrom,
                    "fit_plot": fit,
                    "copy_number_plot": copy,
                    "rds_file": rds,
                    "segment_count": 1,
                }
            )
        consensus = f"{sample}.consensus.chromosomes.tsv"
        (out / consensus).write_text(
            "chromosome\tmedian_copy_number\trounded_copy_number\tagreeing_bins\tcontributing_bins\tmin_copy_number\tmax_copy_number\n"
            "chr7\t1.0\t1\t3\t3\t1.0\t1.0\n",
            encoding="utf-8",
        )
        summary = {
            "schema_version": "0.1.0",
            "sample_id": sample,
            "genome_build": "GRCh37",
            "genome_annotation": "hg19",
            "primary_bin_size_kbp": 500,
            "ace_penalty": 0.6,
            "ploidy_search": {"min": 1.5, "max": 4.5, "step": 0.05},
            "consensus_strategy": "median_rounded_across_bins",
            "package_versions": {
                "R": "4.4.0",
                "QDNAseq": "1.42.0",
                "ACE": "1.24.0",
                "DNAcopy": "1.80.0",
                "QDNAseq_hg19": "1.36.0",
                "QDNAseq_hg38": None,
            },
            "runs": runs,
            "primary": runs[1],
            "consensus_file": consensus,
        }
        (out / f"{sample}.qdnaseq-ace.summary.json").write_text(
            json.dumps(summary), encoding="utf-8"
        )
        return CommandResult(tuple(args), 0, "", "")

    def run_to_file(self, argv, output_path, *, timeout_seconds: int = 300):  # noqa: ANN001, ANN201
        raise AssertionError("QDNAseq adapter should not use run_to_file")


def _lock() -> ReferenceLock:
    return ReferenceLock(
        reference_id="hg19-test",
        genome_build=GenomeBuild.GRCH37,
        contigs=[ReferenceContig(name="chr7", length=1000)],
        source_fai_sha256="a" * 64,
    )


def _policy() -> QDNAseqPolicy:
    return QDNAseqPolicy(
        profile_id="test",
        expected_qdnaseq_version="1.42.0",
        expected_ace_version="1.24.0",
        note="test policy",
    )


def test_qdnaseq_adapter_promotes_complete_result_and_normalizes_event(tmp_path: Path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"bam")
    script = tmp_path / "runner.R"
    script.write_text("# synthetic", encoding="utf-8")
    output = tmp_path / "cnv"

    report = run_qdnaseq_ace(
        bam=bam,
        sample_id="SAMPLE_001",
        genome_build=GenomeBuild.GRCH37,
        reference_lock=_lock(),
        policy=_policy(),
        output_dir=output,
        script=script,
        runner=FakeQDNAseqRunner(),
        threads=2,
    )

    assert output.is_dir()
    assert report.primary_fit.bin_size_kbp == 500
    assert report.primary_fit.cellularity == pytest.approx(0.57)
    assert report.primary_fit.ploidy == pytest.approx(2.0)
    assert len(report.events) == 1
    assert report.events[0].event_type == EventType.CHROMOSOME_LOSS
    assert report.events[0].copy_number == pytest.approx(1.0)
    assert report.chromosome_consensus[0].agreeing_bins == 3
    assert not any(name.endswith(".rds") for name in report.output_files if name == "")


def test_qdnaseq_adapter_returns_no_call_when_ace_call_is_neutral(tmp_path: Path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"bam")
    script = tmp_path / "runner.R"
    script.write_text("# synthetic", encoding="utf-8")

    report = run_qdnaseq_ace(
        bam=bam,
        sample_id="SAMPLE_002",
        genome_build=GenomeBuild.GRCH37,
        reference_lock=_lock(),
        policy=_policy(),
        output_dir=tmp_path / "cnv",
        script=script,
        runner=FakeQDNAseqRunner(event_call=0.0),
    )

    assert report.events == []
    assert report.status.value == "NO_CALL"


def test_qdnaseq_failure_never_promotes_partial_directory(tmp_path: Path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"bam")
    script = tmp_path / "runner.R"
    script.write_text("# synthetic", encoding="utf-8")
    output = tmp_path / "cnv"

    with pytest.raises(ValueError, match="exited with code 9"):
        run_qdnaseq_ace(
            bam=bam,
            sample_id="SAMPLE_003",
            genome_build=GenomeBuild.GRCH37,
            reference_lock=_lock(),
            policy=_policy(),
            output_dir=output,
            script=script,
            runner=FakeQDNAseqRunner(fail=True),
        )

    assert not output.exists()
    assert not list(tmp_path.glob(".cnv.*"))

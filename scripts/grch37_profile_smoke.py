"""Run an explicit GRCh37 dictionary profile on deterministic non-patient depth data.

This fixture is deliberately not a model of biological long-read data. It proves
resource binding and real tool execution; it cannot validate sensitivity or a
negative SV result. No clinical threshold or release gate is relaxed here.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path

from ontseq_platform import __version__
from ontseq_platform.cnv.extension import (
    QDNAseqExtensionSettings,
    register_qdnaseq_extension,
)
from ontseq_platform.cnv.qdnaseq import QDNAseqCallReport, QDNAseqPolicy
from ontseq_platform.io import load_model
from ontseq_platform.models import (
    CheckStatus,
    GenomeBuild,
    ModuleRunStatus,
    PipelineResult,
    ReferenceLock,
)
from ontseq_platform.pipeline.envelope import sha256_file
from ontseq_platform.pipeline.runner import run_pipeline
from ontseq_platform.pipeline.stages import StageId
from ontseq_platform.profile_analysis import AnalyzeSettings, build_profile_run_configuration
from ontseq_platform.reference import reference_lock_for_dictionary_contract
from ontseq_platform.resource_bootstrap import GRCH37_PROFILE_IDS
from ontseq_platform.resource_registry import ResourceRegistry
from ontseq_platform.system_smoke import (
    _verify_release_checksums,
    _write_cnv_sam,
    cnv_truth_checks,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROFILE = "AML_LCWGS_GRCh37"
SAMPLE = "SYNTHETIC_GRCH37_NATIVE"


def make_fixture(output: Path, lock: ReferenceLock) -> Path:
    """Use the profile's exact effective dictionary without aliasing or reheadering."""
    output.mkdir(parents=True, exist_ok=False)
    depth_source = output / "synthetic-depth.sam"
    native_sam = output / "synthetic-native.sam"
    bam = output / f"{SAMPLE}.bam"
    _write_cnv_sam(depth_source)
    with native_sam.open("w", encoding="utf-8", newline="\n") as target:
        target.write("@HD\tVN:1.6\tSO:unsorted\n")
        for contig in lock.contigs:
            target.write(f"@SQ\tSN:{contig.name}\tLN:{contig.length}\n")
        target.write(f"@RG\tID:SYSTEM_CNV_SMOKE\tSM:{SAMPLE}\tPL:ONT\n")
        with depth_source.open(encoding="utf-8") as source:
            for line in source:
                if not line.startswith("@"):
                    target.write(line)
    subprocess.run(["samtools", "sort", "-@", "2", "-o", str(bam), str(native_sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)
    return bam


def coordinate_checks(
    envelope: Path, cnv: QDNAseqCallReport, result: PipelineResult, lock: ReferenceLock
) -> tuple[bool, dict[str, object]]:
    """Verify the real R TSV boundary and its normalized event consumers."""
    detail: dict[str, object] = {"coordinate_system": "zero_based_half_open", "passed": False}
    lengths = {item.name: item.length for item in lock.contigs}
    folder = envelope / "evidence/cnv/qdnaseq"
    try:
        raw_summary = json.loads((folder / f"{SAMPLE}.qdnaseq-ace.summary.json").read_text())
        expected = {
            "coordinate_system": "zero_based_half_open",
            "source_coordinate_system": "one_based_inclusive",
            "coordinate_normalization": "qdnaseq_start_minus_one_end_unchanged",
            "native_rds_coordinate_system": "one_based_inclusive",
        }
        if any(raw_summary.get(key) != value for key, value in expected.items()):
            raise ValueError("R summary does not declare the corrected export/native-RDS contract")
        tables: list[dict[str, object]] = []
        for fit in cnv.fits:
            for name in (fit.bins_file, fit.segment_file):
                if not name:
                    raise ValueError("CNV fit does not expose both coordinate-bearing TSV files")
                with (folder / name).open(newline="", encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                if not rows:
                    raise ValueError(f"synthetic CNV coordinate table is empty: {name}")
                for row in rows:
                    chromosome = row["chromosome"]
                    if not chromosome.startswith("chr"):
                        chromosome = "chr" + chromosome
                    start, end = int(row["start"]), int(row["end"])
                    if row.get("coordinate_system") != "zero_based_half_open":
                        raise ValueError(f"missing/wrong TSV coordinate declaration: {name}")
                    if start < 0 or end <= start or end > lengths[chromosome]:
                        raise ValueError(f"out-of-bounds half-open interval: {name}")
                    if start % (fit.bin_size_kbp * 1000):
                        raise ValueError(f"start is not zero-based QDNAseq-bin aligned: {name}")
                tables.append({"file": name, "rows_checked": len(rows)})
        by_id = {event.event_id: event for event in result.events}
        if not cnv.events:
            raise ValueError("synthetic fixture did not produce CNV events")
        for event in cnv.events:
            locus = event.primary
            if event.length_bp != locus.end - locus.start:
                raise ValueError("CNV event length is inconsistent with its half-open locus")
            if locus.start < 0 or locus.end > lengths[locus.chromosome]:
                raise ValueError("normalized CNV event is outside the reference chromosome")
            if locus.start % (cnv.primary_fit.bin_size_kbp * 1000):
                raise ValueError("normalized CNV event start is not zero-based bin aligned")
            normalized = by_id[event.event_id]
            if (
                normalized.primary.chromosome,
                normalized.primary.start,
                normalized.primary.end,
                normalized.length_bp,
            ) != (locus.chromosome, locus.start, locus.end, event.length_bp):
                raise ValueError("CNV coordinates changed between caller report and PipelineResult")
        chr7_starts = [
            event.primary.start for event in cnv.events if event.primary.chromosome == "chr7"
        ]
        chr8_starts = [
            event.primary.start for event in cnv.events if event.primary.chromosome == "chr8"
        ]
        if not chr7_starts or min(chr7_starts) != 0:
            raise ValueError("synthetic chromosome-7 event must begin at zero")
        detail.update(
            passed=True,
            tsv_files_checked=tables,
            events_checked=len(cnv.events),
            first_chr7_event_start=min(chr7_starts),
            first_chr8_event_start=min(chr8_starts) if chr8_starts else None,
        )
    except (OSError, ValueError, KeyError) as error:
        detail["error"] = str(error)
    return detail["passed"] is True, detail


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resource-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--profile", choices=GRCH37_PROFILE_IDS, default=DEFAULT_PROFILE)
    parser.add_argument("--fixture-only", action="store_true")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        parser.error("output directory must be new; existing data is never overwritten")
    registry = ResourceRegistry(args.resource_root, active_build=GenomeBuild.GRCH37)
    context = registry.resolve_profile(args.profile, verify_files=True)
    source_lock = load_model(
        Path(context.resource_paths["reference.reference_lock"]), ReferenceLock
    )
    lock = reference_lock_for_dictionary_contract(
        source_lock, context.reference_dictionary_contract
    )
    bam = make_fixture(output / "inputs", lock)
    if args.fixture_only:
        print(json.dumps({"fixture": str(bam), "patient_data": False}))
        return 0
    policy = load_model(ROOT / "configs/cnv/qdnaseq_ace.technical.yaml", QDNAseqPolicy)
    register_qdnaseq_extension(
        QDNAseqExtensionSettings(policy=policy, script=ROOT / "scripts/run_qdnaseq_ace.R")
    )
    settings = AnalyzeSettings(
        bam=bam,
        profile_id=args.profile,
        resource_root=args.resource_root,
        output_dir=output / "runs",
        sample_id=SAMPLE,
        run_id="NATIVE_GRCH37_PROFILE_SMOKE",
        pipeline_version=__version__,
        git_commit="LOCAL_WORKTREE",
        threads=2,
    )
    config = build_profile_run_configuration(settings)
    # The launcher normally obtains this at intake. Explicitly bind this fixture too.
    config.manifest.input.sha256 = sha256_file(bam)
    report, bundle = run_pipeline(config)
    envelope = output / "runs" / config.run_id / SAMPLE
    summary: dict[str, object] = {
        "pipeline_version": __version__,
        "profile": args.profile,
        "reference_build": lock.genome_build.value,
        "native_dictionary_contigs": len(lock.contigs),
        "patient_data": False,
        "clinical_validation": False,
        "stages": {item.stage.value: item.status.value for item in report.stages},
        "stage_reasons": {item.stage.value: item.reason for item in report.stages},
        "release_bundle_created": bundle is not None,
        "pipeline_passed": report.passed,
        "envelope": str(envelope),
    }
    cnv_path = envelope / f"evidence/cnv/{SAMPLE}.qdnaseq.json"
    result_path = envelope / f"normalized/{SAMPLE}.result.json"
    checks_ok = False
    if cnv_path.is_file() and result_path.is_file():
        cnv = load_model(cnv_path, QDNAseqCallReport)
        truth = cnv_truth_checks(cnv, policy)
        result = load_model(result_path, PipelineResult)
        checks_ok = all(item.status == CheckStatus.PASS for item in truth)
        checks_ok &= result.manifest.assay.genome_build == GenomeBuild.GRCH37
        coordinates_ok, coordinates = coordinate_checks(envelope, cnv, result, lock)
        checks_ok &= coordinates_ok
        summary["cnv_coordinates"] = coordinates
        summary["cnv_truth"] = [item.model_dump(mode="json") for item in truth]
        summary["reference_context"] = result.reference_context.model_dump(mode="json")
    checksums_before = _verify_release_checksums(envelope)
    summary["release_checksums_before_resume"] = checksums_before
    resumed, _ = run_pipeline(config)
    cnv_stage = next(item for item in resumed.stages if item.stage == StageId.CNV)
    sv_stage = next(item for item in report.stages if item.stage == StageId.SV)
    summary["cnv_resumed"] = cnv_stage.resumed
    summary["release_checksums_after_resume"] = _verify_release_checksums(envelope)
    passed = bool(
        checks_ok
        and report.passed
        and resumed.passed
        and sv_stage.status in {ModuleRunStatus.COMPLETED, ModuleRunStatus.NO_CALL}
        and bundle is not None
        and checksums_before[0]
        and cnv_stage.resumed
        and _verify_release_checksums(envelope)[0]
    )
    summary["technical_smoke_passed"] = passed
    (output / "native-grch37-smoke.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

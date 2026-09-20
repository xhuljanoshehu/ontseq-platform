"""Synthetic adapter tests; these do not validate the GATK binary or an ONT assay."""

import gzip
import hashlib
import importlib
import importlib.util
import json
import subprocess
from dataclasses import replace
from pathlib import Path

import pytest


def test_gatk_adapter_exists():
    assert importlib.util.find_spec("ontseq_platform.gatk_mutect2") is not None


@pytest.fixture
def api():
    return importlib.import_module("ontseq_platform.gatk_mutect2")


def lock(api, path, text="synthetic"):
    path.write_text(text)
    return api.LockedFile(str(path), hashlib.sha256(path.read_bytes()).hexdigest())


@pytest.fixture
def config(api, tmp_path):
    ref = lock(api, tmp_path / "ref.fa", ">chr1\nACGTACGTACGTACGTACGT\n")
    fai = lock(api, tmp_path / "ref.fa.fai", "chr1\t20\t6\t20\t21\n")
    dictionary = lock(api, tmp_path / "ref.dict", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20\n")
    bam = lock(api, tmp_path / "tumor.bam")
    bai = lock(api, tmp_path / "tumor.bam.bai")
    bed = lock(api, tmp_path / "targets.bed", "chr1\t0\t20\n")
    return api.GATKConfig(
        profile_id="synthetic-v1", reference=ref, reference_fai=fai,
        reference_dict=dictionary, tumor=api.Sample("TUMOR", bam, bai, ref.sha256),
        intervals=bed,
    )


def vcf_text(sample_names=("TUMOR",), filters="PASS", ad="18,2", af="0.12", alt="T"):
    row = "GT:AD:AF:DP\t0/1:" + ad + ":" + af + ":21"
    if len(sample_names) == 2:
        row = "GT:AD:AF:DP\t0/0:20,0:0:20\t0/1:" + ad + ":" + af + ":21"
    return (
        "##fileformat=VCFv4.2\n##source=Mutect2\n"
        "##contig=<ID=chr1,length=20>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(sample_names) + "\n"
        + "chr1\t2\t.\tC\t" + alt + "\t.\t" + filters + "\tTLOD=12\t" + row + "\n"
    )


class SyntheticRunner:
    def __init__(self, filtered=None, version="4.6.2.0", fail=None, omit_stats=False,
                 normal=False, header_sample="TUMOR"):
        self.calls = []
        self.filtered = filtered or vcf_text()
        self.version = version
        self.fail = fail
        self.omit_stats = omit_stats
        self.normal = normal
        self.header_sample = header_sample

    def __call__(self, argv, cwd, timeout):
        self.calls.append(argv)
        if "--version" in argv:
            return "The Genome Analysis Toolkit (GATK) v" + self.version
        if argv[0] == "samtools":
            if "view" in argv:
                name = "NORMAL" if "normal.bam" in argv[-1] else self.header_sample
                return "@HD\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n@RG\tID:rg\tSM:" + name + "\n"
            return "samtools 1.22.1" if "version" in argv else ""
        tool = next((x for x in argv if x in {
            "Mutect2", "FilterMutectCalls", "GetPileupSummaries", "CalculateContamination"
        }), None)
        if tool == self.fail:
            raise subprocess.CalledProcessError(1, argv, output="synthetic failure")
        out = Path(argv[argv.index("-O") + 1])
        if tool in ("Mutect2", "FilterMutectCalls"):
            with gzip.open(out, "wt") as stream:
                stream.write(self.filtered)
            if tool == "Mutect2" and not self.omit_stats:
                Path(str(out) + ".stats").write_text("synthetic stats")
        else:
            out.write_text("synthetic table")
        return "synthetic completed"


def test_default_is_blocked_without_running_or_creating_output(api, config, tmp_path):
    runner = SyntheticRunner()
    out = tmp_path / "run"
    with pytest.raises(api.GATKError, match="experimental"):
        api.run_gatk(config, out, runner=runner)
    assert not out.exists()
    assert runner.calls == []


def test_plans_both_modes_and_keeps_stats(api, config, tmp_path):
    plan = api.build_commands(config, tmp_path / "run")
    assert len(plan) == 2
    assert "Mutect2" in plan[0]
    assert "FilterMutectCalls" in plan[1]
    assert "--stats" in plan[1]
    assert "-normal" not in plan[0]
    assert not any("BQSR" in x or "MarkDuplicates" in x or "f1r2" in x for p in plan for x in p)
    normal = api.Sample("NORMAL", lock(api, tmp_path / "normal.bam", "synthetic-normal"),
                        lock(api, tmp_path / "normal.bam.bai"), config.reference.sha256)
    paired = api.build_commands(replace(config, normal=normal), tmp_path / "paired")
    assert paired[0][paired[0].index("-normal") + 1] == "NORMAL"
    assert paired[0].count("-I") == 2


def test_success_is_research_only_with_provenance(api, config, tmp_path):
    runner = SyntheticRunner()
    result = api.run_gatk(config, tmp_path / "run", allow_experimental_ont=True, runner=runner)
    assert result["status"] == "COMPLETED"
    assert result["reportable"] is False
    assert result["analytical_validation"] == "NOT_VALIDATED"
    assert result["variants"][0]["somatic_status"] == "NOT_VALIDATED"
    assert result["variants"][0]["quality"] is None
    assert result["variants"][0]["caller_af"] == 0.12
    assert result["variants"][0]["ad_fraction"] == 0.1
    assert result["variants"][0]["depth"] == 21
    assert result["reference_sha256"] == config.reference.sha256
    assert result["input_sha256"]
    assert result["artifact_sha256"]["filtered.vcf.gz"]
    assert json.loads((tmp_path / "run/evidence.json").read_text())["reportable"] is False


@pytest.mark.parametrize("version", ["4.6.1.0", "4.6.2.01", "4.7.0.0", "garbage"])
def test_version_mismatch_blocks_calling(api, config, tmp_path, version):
    runner = SyntheticRunner(version=version)
    with pytest.raises(api.GATKError, match="version"):
        api.run_gatk(config, tmp_path / "run", allow_experimental_ont=True, runner=runner)
    assert not any("Mutect2" in call for call in runner.calls)
    assert json.loads((tmp_path / "run/failure.json").read_text())["status"] == "FAILED"


@pytest.mark.parametrize("change", ["hash", "build", "sample", "index", "bed"])
def test_bad_inputs_fail_closed(api, config, tmp_path, change):
    if change == "hash":
        Path(config.tumor.bam.path).write_text("changed")
    elif change == "build":
        config = replace(config, tumor=replace(config.tumor, reference_sha256="0" * 64))
    elif change == "sample":
        config = replace(config, tumor=replace(config.tumor, name="WRONG"))
    elif change == "index":
        Path(config.tumor.index.path).unlink()
    else:
        config = replace(config, intervals=lock(api, tmp_path / "bad.bed", "chr1\t0\t99\n"))
    runner = SyntheticRunner()
    with pytest.raises(api.GATKError):
        api.run_gatk(config, tmp_path / "run", allow_experimental_ont=True, runner=runner)
    assert not any("Mutect2" in call for call in runner.calls)


@pytest.mark.parametrize("failure", ["Mutect2", "FilterMutectCalls", "stats"])
def test_failure_never_reuses_or_publishes_evidence(api, config, tmp_path, failure):
    runner = SyntheticRunner(fail=failure, omit_stats=failure == "stats")
    out = tmp_path / "run"
    with pytest.raises(api.GATKError):
        api.run_gatk(config, out, allow_experimental_ont=True, runner=runner)
    assert not (out / "evidence.json").exists()
    assert json.loads((out / "failure.json").read_text())["status"] == "FAILED"
    with pytest.raises(api.GATKError, match="exist"):
        api.run_gatk(config, out, allow_experimental_ont=True, runner=runner)


@pytest.mark.parametrize("filters", [".", "germline", "weak_evidence;strand_bias"])
def test_filtered_records_are_preserved_but_not_pass(api, config, tmp_path, filters):
    result = api.run_gatk(config, tmp_path / "run", allow_experimental_ont=True,
                          runner=SyntheticRunner(filtered=vcf_text(filters=filters)))
    assert result["status"] == "NO_CALL"
    assert result["variants"][0]["native_filter"] == filters
    assert result["variants"][0]["technical_pass"] is False


def test_multiallelic_depth_and_paired_sample_selection(api, config, tmp_path):
    parser = importlib.import_module("ontseq_platform.gatk_vcf")
    path = tmp_path / "test.vcf.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(vcf_text(("NORMAL", "TUMOR"), ad="15,3,2", af="0.15,0.10", alt="T,G"))
    result = parser.read_candidates(path, "TUMOR", ("NORMAL", "TUMOR"))
    assert [row["alternate_reads"] for row in result] == [3, 2]
    assert [row["ad_fraction"] for row in result] == [0.15, 0.10]
    assert [row["caller_af"] for row in result] == [0.15, 0.10]


@pytest.mark.parametrize("text", [
    "", vcf_text(sample_names=("WRONG",)), vcf_text(ad="18"),
    vcf_text(af="nan"), vcf_text(af="1.1"), vcf_text(ad="-1,2"),
    vcf_text(alt="<DEL>"), vcf_text().replace("chr1\t2", "chr1\t0"),
])
def test_malformed_vcf_is_failure_not_negative(api, config, tmp_path, text):
    parser = importlib.import_module("ontseq_platform.gatk_vcf")
    path = tmp_path / "bad.vcf"
    path.write_text(text)
    with pytest.raises(ValueError):
        parser.read_candidates(path, "TUMOR", ("TUMOR",))


def test_header_only_vcf_is_no_call(api, config, tmp_path):
    text = vcf_text().rsplit("\n", 2)[0] + "\n"
    result = api.run_gatk(config, tmp_path / "run", allow_experimental_ont=True,
                          runner=SyntheticRunner(filtered=text))
    assert result["status"] == "NO_CALL"
    assert result["variants"] == []


def test_same_sample_cannot_be_tumor_and_normal(api, config, tmp_path):
    with pytest.raises(api.GATKError):
        api.build_commands(replace(config, normal=config.tumor), tmp_path / "run")


def test_subprocess_runner_uses_argv_and_timeout(api, tmp_path, monkeypatch):
    captured = {}
    def fake_run(argv, **kwargs):
        captured.update(kwargs)
        captured["argv"] = argv
        return subprocess.CompletedProcess(argv, 0, stdout="ok")
    monkeypatch.setattr(subprocess, "run", fake_run)
    argv = ["gatk", "-I", "/tmp/data with spaces;echo unsafe.bam"]
    assert api.subprocess_runner(argv, tmp_path, 123) == "ok"
    assert captured["shell"] is False
    assert captured["timeout"] == 123
    assert captured["argv"] == argv


def resource(api, config, path, technology="not_applicable", assay_id=""):
    with gzip.open(path, "wt") as stream:
        stream.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=20>\n"
                     "##INFO=<ID=AF,Number=A,Type=Float,Description=\"frequency\">\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    locked = api.LockedFile(str(path), hashlib.sha256(path.read_bytes()).hexdigest())
    return api.VCFResource(locked, lock(api, Path(str(path) + ".tbi")),
                            config.reference.sha256, technology, assay_id)


@pytest.mark.parametrize("paired", [False, True])
def test_resources_and_contamination_path(api, config, tmp_path, paired):
    germline = resource(api, config, tmp_path / "germline.vcf.gz")
    pon = resource(api, config, tmp_path / "pon.vcf.gz", "ONT", config.profile_id)
    sites = resource(api, config, tmp_path / "sites.vcf.gz")
    config = replace(config, germline=germline, panel_of_normals=pon, population_sites=sites)
    if paired:
        normal = api.Sample("NORMAL", lock(api, tmp_path / "normal.bam", "normal"),
                            lock(api, tmp_path / "normal.bam.bai"), config.reference.sha256)
        config = replace(config, normal=normal)
    runner = SyntheticRunner(filtered=vcf_text(("NORMAL", "TUMOR")) if paired else None)
    result = api.run_gatk(config, tmp_path / "run", allow_experimental_ont=True, runner=runner)
    assert result["status"] == "COMPLETED"
    assert not any("not estimated" in warning for warning in result["warnings"])
    assert sum("GetPileupSummaries" in call for call in runner.calls) == (2 if paired else 1)
    contamination = next(c for c in runner.calls if "CalculateContamination" in c)
    assert ("--matched-normal" in contamination) is paired
    filtering = next(c for c in runner.calls if "FilterMutectCalls" in c)
    assert "--contamination-table" in filtering


@pytest.mark.parametrize("technology,assay", [("ILLUMINA", "synthetic-v1"), ("ONT", "other")])
def test_incompatible_pon_is_blocked(api, config, tmp_path, technology, assay):
    pon = resource(api, config, tmp_path / "pon.vcf.gz", technology, assay)
    with pytest.raises(api.GATKError, match="PoN"):
        api.build_commands(replace(config, panel_of_normals=pon), tmp_path / "run")


def test_allele_specific_filters_are_not_promoted_by_site_pass(api, tmp_path):
    parser = importlib.import_module("ontseq_platform.gatk_vcf")
    path = tmp_path / "as.vcf"
    path.write_text(vcf_text(ad="15,3,2", af="0.15,0.10", alt="T,G").replace(
        "TLOD=12", "TLOD=12,11;AS_FilterStatus=weak_evidence|SITE"))
    rows = parser.read_candidates(path, "TUMOR", ("TUMOR",))
    assert [v["technical_pass"] for v in rows] == [False, True]
    assert rows[0]["allele_filter"] == "weak_evidence"


def test_ref_allele_mismatch_fails_publication(api, config, tmp_path):
    text = vcf_text().replace("\tC\tT\t", "\tG\tT\t")
    with pytest.raises(api.GATKError, match="REF"):
        api.run_gatk(config, tmp_path / "run", allow_experimental_ont=True,
                      runner=SyntheticRunner(filtered=text))
    assert not (tmp_path / "run/evidence.json").exists()


def test_input_changed_during_execution_fails_publication(api, config, tmp_path):
    fake = SyntheticRunner()
    def runner(argv, cwd, timeout):
        result = fake(argv, cwd, timeout)
        if "FilterMutectCalls" in argv:
            Path(config.tumor.bam.path).write_text("changed during run")
        return result
    with pytest.raises(api.GATKError):
        api.run_gatk(config, tmp_path / "run", allow_experimental_ont=True, runner=runner)
    assert not (tmp_path / "run/evidence.json").exists()


def test_missing_quantities_stay_null(api, tmp_path):
    parser = importlib.import_module("ontseq_platform.gatk_vcf")
    path = tmp_path / "missing.vcf"
    path.write_text(vcf_text().replace("0/1:18,2:0.12:21", "0/1:.:.:."))
    row = parser.read_candidates(path, "TUMOR", ("TUMOR",))[0]
    assert all(row[k] is None for k in ["depth", "alternate_reads", "ad_fraction", "caller_af"])


def test_config_round_trip_and_unknown_fields(api, config, tmp_path):
    from dataclasses import asdict
    path = tmp_path / "manifest.json"
    raw = asdict(config)
    path.write_text(json.dumps(raw))
    assert api.load_config(path) == config
    raw["clinical_release"] = True
    path.write_text(json.dumps(raw))
    with pytest.raises(TypeError):
        api.load_config(path)


def test_invalid_json_top_level_is_rejected(api, tmp_path):
    path = tmp_path / "invalid.json"
    path.write_text('[]')
    with pytest.raises((ValueError, TypeError)):
        api.load_config(path)


@pytest.mark.parametrize("key", ["normal", "germline", "panel_of_normals", "population_sites"])
@pytest.mark.parametrize("value", [{}, [], False, ""])
def test_malformed_optional_inputs_are_not_silently_omitted(api, config, tmp_path, key, value):
    from dataclasses import asdict

    raw = asdict(config)
    raw[key] = value
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(raw))
    with pytest.raises((ValueError, TypeError, KeyError)):
        api.load_config(path)


def test_ref_allele_check_handles_wrapped_fasta(api, config, tmp_path):
    config = replace(config,
                     reference=lock(api, tmp_path / "ref.fa", ">chr1\nACGT\nACGT\nACGT\nACGT\nACGT\n"),
                     reference_fai=lock(api, tmp_path / "ref.fa.fai", "chr1\t20\t6\t4\t5\n"))
    api._verify_ref_alleles(config, [{"chromosome": "chr1", "position": 3, "reference": "GTACGT"}])

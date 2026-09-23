from pathlib import Path

import pytest

from ontseq_platform.marlin_native import combined_probe_fractions, encode_native_features
from ontseq_platform.marlin_native_contracts import NativeMarlinReport
from ontseq_platform.models import GenomeBuild, ModuleRunStatus


def row(start, modified, valid, strand="+"):
    return (
        f"chr1\t{start}\t{start + 1}\tC\t{valid}\t{strand}\t{start}\t{start + 1}\t0,0,0"
        f"\t{valid}\t{100 * modified / valid if valid else 0}\t{modified}"
        f"\t{valid - modified}\t0\t0\t0\t0\t0\n"
    )


def test_weighted_pooling_and_zero_distinct_from_missing(tmp_path):
    probes = tmp_path / "probes.bed"
    probes.write_text("chr1\t10\t11\tp1\nchr1\t11\t12\tp1\nchr1\t20\t21\tp2\n")
    calls = tmp_path / "calls.bed"
    calls.write_text(row(10, 1, 1) + row(11, 0, 9, "-") + row(20, 0, 5))
    fractions = combined_probe_fractions(calls, probes)
    assert fractions == {"p1": 0.1, "p2": 0.0}
    assert encode_native_features(("p1", "p2", "missing"), fractions) == (-1, -1, 0)
    assert encode_native_features(("a", "b", "c", "d"), {"a": 0, "b": 0.49, "c": 0.5, "d": 1}) == (
        -1,
        -1,
        1,
        1,
    )


@pytest.mark.parametrize("payload", [row(10, 1, 2) * 2, row(10, 1, 2, ".") + row(10, 1, 2)])
def test_duplicate_or_mixed_combined_strand_rejected(tmp_path, payload):
    probes = tmp_path / "probes.bed"
    probes.write_text("chr1\t10\t11\tp1\n")
    calls = tmp_path / "calls.bed"
    calls.write_text(payload)
    with pytest.raises(ValueError, match="strand|duplicate"):
        combined_probe_fractions(calls, probes)


def test_report_refuses_stale_prediction_on_failure():
    blocked = dict(
        run_id="run",
        sample_id="sample",
        genome_build=GenomeBuild.GRCH38,
        status=ModuleRunStatus.FAILED,
        reason="worker failed",
    )
    assert NativeMarlinReport(**blocked).top_class is None
    with pytest.raises(ValueError):
        NativeMarlinReport(**blocked, top_class="stale")


def test_report_refuses_success_without_evidence():
    with pytest.raises(ValueError):
        NativeMarlinReport(
            run_id="r",
            sample_id="s",
            genome_build=GenomeBuild.GRCH38,
            status=ModuleRunStatus.COMPLETED,
            reason="ok",
        )


def test_unknown_probe_position_and_nonfinite_counts_refused(tmp_path):
    probes = tmp_path / "probes.bed"
    probes.write_text("chr1\t10\t11\tp1\n")
    calls = tmp_path / "calls.bed"
    calls.write_text(row(11, 1, 2))
    with pytest.raises(ValueError, match="outside"):
        combined_probe_fractions(calls, probes)
    calls.write_text(row(10, 1, 2).replace("\t50.0\t", "\tnan\t"))
    with pytest.raises(ValueError, match="percent"):
        combined_probe_fractions(calls, probes)


def test_actual_network_confinement_is_kernel_enforced():
    import subprocess
    import sys

    from ontseq_platform.marlin_native import _WORKER

    result = subprocess.run(
        [
            sys.executable,
            "-I",
            "-B",
            "-c",
            "import runpy; d=runpy.run_path(" + repr(str(_WORKER)) + '); d["confine_network"]()',
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def installation_fixture(tmp_path, monkeypatch):
    import json

    from openpyxl import Workbook

    from ontseq_platform import marlin_native as native
    from ontseq_platform.reference import sha256_file

    runtime = tmp_path / "runtime"
    (runtime / "bin").mkdir(parents=True)
    python = runtime / "bin/python"
    python.write_bytes(b"synthetic executable")
    model = tmp_path / "model.hdf5"
    model.write_bytes(b"synthetic model only for fake runner")
    features = tmp_path / "features.txt"
    features.write_text("".join(f"p{i}\n" for i in range(357340)))
    classes = tmp_path / "classes.xlsx"
    workbook = Workbook()
    workbook.active.append(["model_id", "class_name_current", "mcf", "lineage"])
    for i in range(1, 43):
        workbook.active.append([i, f"class{i}", f"family{i % 3}", f"lineage{i % 2}"])
    workbook.save(classes)
    probes = tmp_path / "probes.bed"
    probes.write_text("chr1\t10\t11\tp0\n")
    runtime_manifest = tmp_path / "runtime.json"
    runtime_manifest.write_text(json.dumps({str(python): sha256_file(python)}))

    def asset(p):
        return {"path": str(p), "sha256": sha256_file(p), "source_uri": "synthetic://test"}

    config = {
        "model": asset(model),
        "features": asset(features),
        "classes": asset(classes),
        "probe_maps": {"GRCh38": {**asset(probes), "genome_build": "GRCh38"}},
        "python_executable": str(python),
        "python_version": "3.10.21",
        "runtime_manifest": asset(runtime_manifest),
        "runtime_archive_sha256": "a" * 64,
    }
    path = tmp_path / "installation.json"
    path.write_text(json.dumps(config))
    monkeypatch.setattr(native, "MODEL_SHA256", sha256_file(model))
    monkeypatch.setattr(
        native, "_APPROVED_RUNTIME_ARCHIVE_SHA256", config["runtime_archive_sha256"]
    )
    monkeypatch.setattr(
        native,
        "_QUALIFIED_RUNTIME_PROFILES",
        {sha256_file(runtime_manifest): ("synthetic-qualified", "3.10.21")},
    )
    monkeypatch.setattr(native, "_FEATURES_SHA256", sha256_file(features))
    monkeypatch.setattr(native, "_CLASSES_SHA256", sha256_file(classes))
    monkeypatch.setattr(
        native,
        "_PROBE_SHA256",
        {GenomeBuild.GRCH38: sha256_file(probes), GenomeBuild.GRCH37: "b" * 64},
    )
    monkeypatch.setattr(native, "_QUALIFIED_DIRECTORY_LINKS", {})
    monkeypatch.setattr(native, "_preflight", lambda installation: None)
    return path, config


def test_installation_mislabeled_build_and_checksum_refused(tmp_path, monkeypatch):
    import json

    from ontseq_platform.marlin_native import native_marlin_signature

    path, config = installation_fixture(tmp_path, monkeypatch)
    assert native_marlin_signature(path, GenomeBuild.GRCH38)["genome_build"] == "GRCh38"
    with pytest.raises(ValueError, match="probe map"):
        native_marlin_signature(path, GenomeBuild.GRCH37)
    config["probe_maps"]["GRCh38"]["genome_build"] = "GRCh37"
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="build"):
        native_marlin_signature(path, GenomeBuild.GRCH38)
    config["probe_maps"]["GRCh38"]["genome_build"] = "GRCh38"
    path.write_text(json.dumps(config))
    Path(config["features"]["path"]).write_text("changed")
    with pytest.raises(ValueError, match="checksum"):
        native_marlin_signature(path, GenomeBuild.GRCH38)


def test_unmanifested_runtime_file_refused(tmp_path, monkeypatch):
    from ontseq_platform.marlin_native import native_marlin_signature

    path, config = installation_fixture(tmp_path, monkeypatch)
    (Path(config["python_executable"]).parent.parent / "injected.py").write_text("injected")
    with pytest.raises(ValueError, match="inventory"):
        native_marlin_signature(path, GenomeBuild.GRCH38)


def test_probe_map_cannot_be_relabelled_to_other_build(tmp_path, monkeypatch):
    import json

    from ontseq_platform.marlin_native import _load_installation

    path, config = installation_fixture(tmp_path, monkeypatch)
    config["probe_maps"]["GRCh37"] = {**config["probe_maps"]["GRCh38"], "genome_build": "GRCh37"}
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="approved.*probe|probe.*identity"):
        _load_installation(path, GenomeBuild.GRCH37)


@pytest.mark.parametrize(
    "mode",
    [
        "empty",
        "worker-failure",
        "low",
        "high",
        "nan",
        "stock-independent",
        "index-missing-proof",
        "index-mismatch",
        "index-before-pileup",
        "index-during-pileup",
        "index-during-worker",
        "index-during-empty",
        "index-missing-sha",
        "reference-fai-missing",
        "reference-fai-before-pileup",
        "reference-fai-during-pileup",
        "reference-fai-during-worker",
        "reference-fai-during-empty",
        "compressed-reference",
    ],
)
def test_adapter_states_and_full_grouping(tmp_path, monkeypatch, mode):
    import hashlib
    import json
    import struct

    from test_methylation import AdapterTests

    from ontseq_platform import marlin_native as native
    from ontseq_platform.execution import CommandResult

    installation_path, config = installation_fixture(tmp_path, monkeypatch)
    manifest, intake = AdapterTests()._fixture(tmp_path)
    index = Path(manifest.input.index_path)
    index.write_bytes(b"original index bytes")
    intake = intake.model_copy(update={"index_fingerprint": native._fingerprint(index)})
    if mode == "index-missing-proof":
        intake = intake.model_copy(update={"index_fingerprint": None})
    elif mode == "index-mismatch":
        intake = intake.model_copy(
            update={
                "index_fingerprint": intake.index_fingerprint.model_copy(
                    update={"sha256": "0" * 64}
                )
            }
        )
    if mode == "index-missing-sha":
        intake = intake.model_copy(
            update={
                "index_fingerprint": intake.index_fingerprint.model_copy(update={"sha256": None})
            }
        )
    fai = tmp_path / "reference.fa.fai"
    fai.write_text("chr1\t6\t6\t6\t7\n")
    if mode == "reference-fai-missing":
        fai.unlink()
    if mode == "compressed-reference":
        (tmp_path / "reference.fa").write_bytes(b"\x1f\x8bnot-qualified")
    output = tmp_path / "output"

    def mutate_index():
        import os

        dependency = fai if mode.startswith("reference-fai-") else index
        before = dependency.stat()
        dependency.write_bytes(b"x" * before.st_size)
        os.utime(dependency, ns=(before.st_atime_ns, before.st_mtime_ns))

    commands = []

    class Runner:
        def run(self, argv, *, timeout_seconds=300):
            commands.append(list(argv))
            if argv[1] == "--version":
                if mode in {"index-before-pileup", "reference-fai-before-pileup"}:
                    mutate_index()
                return CommandResult(tuple(argv), 0, "modkit 0.6.4", "")
            if argv[1] == "view":
                text = "1" if argv[-2] == "[MM]" or mode == "stock-independent" else "0"
                return CommandResult(tuple(argv), 0, text, "")
            if argv[1] == "pileup":
                Path(argv[3]).write_text(
                    ""
                    if mode in {"empty", "index-during-empty", "reference-fai-during-empty"}
                    else row(10, 1, 2)
                )
                if mode in {
                    "index-during-pileup",
                    "index-during-empty",
                    "reference-fai-during-pileup",
                    "reference-fai-during-empty",
                }:
                    mutate_index()
                return CommandResult(tuple(argv), 0, "", "")
            worker_config = json.loads(Path(argv[-1]).read_text())
            values = [1 / 42] * 42 if mode == "low" else [0.8] + [0.2 / 41] * 41
            if mode == "nan":
                values[0] = float("nan")
            payload = {
                "scores": values,
                "tensorflow": "2.13.1",
                "keras": "2.13.1",
                "python": config["python_version"],
                "threads": 2,
                "confinement": "linux-seccomp-network-deny-v1",
                "model_sha256": config["model"]["sha256"],
                "tensor_sha256": worker_config["tensor_sha256"],
                "score_vector_sha256": hashlib.sha256(struct.pack("<42f", *values)).hexdigest(),
            }
            Path(worker_config["output_path"]).write_text(json.dumps(payload))
            if mode in {"index-during-worker", "reference-fai-during-worker"}:
                mutate_index()
            return CommandResult(
                tuple(argv),
                1 if mode == "worker-failure" else 0,
                "",
                "failed" if mode == "worker-failure" else "",
            )

    report = native.run_native_marlin(
        run_id=manifest.run_id,
        manifest=manifest,
        intake=intake,
        reference_fasta=tmp_path / "reference.fa",
        installation_path=installation_path,
        output_dir=output,
        runner=Runner(),
        modkit="absent-synthetic-modkit",
        samtools="synthetic-samtools",
        threads=2,
    )
    assert (output / "native-marlin.json").is_file()
    if mode == "empty":
        assert report.status == ModuleRunStatus.NO_CALL, report.reason
        assert report.feature_summary.observed_model_feature_count == 0
        assert not (output / "worker-config.json").exists()
        assert report.raw_model_scores == []
    elif mode.startswith(("index-", "reference-fai-")) or mode in {
        "worker-failure",
        "nan",
        "stock-independent",
        "compressed-reference",
    }:
        assert report.status == ModuleRunStatus.FAILED, report.reason
        assert report.feature_summary is None
        assert report.top_class is None
        assert report.raw_model_scores == []
        if mode in {
            "stock-independent",
            "index-missing-proof",
            "index-mismatch",
            "index-before-pileup",
            "index-missing-sha",
            "reference-fai-missing",
            "reference-fai-before-pileup",
            "compressed-reference",
        }:
            assert not any(c[1] == "pileup" for c in commands)
        if mode.startswith("index-"):
            assert "index" in report.reason.lower()
    else:
        assert report.status == ModuleRunStatus.COMPLETED, report.reason
        assert report.input_fingerprints["bam_index"] == intake.index_fingerprint
        assert len(report.raw_model_scores) == 42
        assert len(report.class_scores) == 42
        assert len(report.family_scores) == 3
        assert len(report.lineage_scores) == 2
        assert report.decision == ("UNKNOWN" if mode == "low" else "HIGH_CONFIDENCE")
        assert report.feature_summary.observed_model_feature_count == 1
        pileup_command = next(c for c in commands if c[1] == "pileup")
        assert pileup_command[pileup_command.index("--modified-bases") + 1 :][:2] == ["5mC", "5hmC"]
        assert "--combine-mods" in pileup_command
    if mode in {"low", "high", "empty"}:
        assert report.input_fingerprints["reference_fai"] == native._fingerprint(fai)
        assert report.input_fingerprints["bam_index"] == intake.index_fingerprint
        for dependency in ("bam_index", "reference_fai"):
            incomplete = report.model_dump(mode="json")
            del incomplete["input_fingerprints"][dependency]
            with pytest.raises(ValueError, match="input fingerprints"):
                NativeMarlinReport.model_validate(incomplete)
    again = native.run_native_marlin(
        run_id=manifest.run_id,
        manifest=manifest,
        intake=intake,
        reference_fasta=tmp_path / "reference.fa",
        installation_path=installation_path,
        output_dir=output,
        runner=Runner(),
        modkit="absent-synthetic-modkit",
        samtools="synthetic-samtools",
        threads=2,
    )
    assert again.status == ModuleRunStatus.FAILED
    assert again.top_class is None


@pytest.mark.parametrize("change", ["archive", "self-rehashed-runtime", "python-path"])
def test_runtime_identity_cannot_self_qualify(tmp_path, monkeypatch, change):
    import json

    from ontseq_platform.marlin_native import native_marlin_signature
    from ontseq_platform.reference import sha256_file

    path, config = installation_fixture(tmp_path, monkeypatch)
    if change == "archive":
        config["runtime_archive_sha256"] = "c" * 64
    elif change == "python-path":
        selected = Path(config["python_executable"])
        other = selected.with_name("python-other")
        other.write_bytes(selected.read_bytes())
        config["python_executable"] = str(other)
    else:
        python = Path(config["python_executable"])
        python.write_bytes(b"self-consistent replacement runtime with unchanged version metadata")
        manifest = Path(config["runtime_manifest"]["path"])
        manifest.write_text(json.dumps({str(python): sha256_file(python)}))
        config["runtime_manifest"]["sha256"] = sha256_file(manifest)
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="approved runtime archive|qualified runtime"):
        native_marlin_signature(path, GenomeBuild.GRCH38)


def test_selected_index_refuses_competing_or_unrelated_files(tmp_path):
    from test_methylation import AdapterTests

    from ontseq_platform.marlin_native import selected_marlin_bam_index

    manifest, _ = AdapterTests()._fixture(tmp_path)
    expected = Path(manifest.input.index_path)
    assert selected_marlin_bam_index(manifest) == expected.resolve()
    other = Path(manifest.input.path).with_suffix(".csi")
    other.write_bytes(b"different index")
    with pytest.raises(ValueError, match="competing"):
        selected_marlin_bam_index(manifest)
    other.unlink()
    unrelated = tmp_path / "other.bai"
    unrelated.write_bytes(b"other")
    altered = manifest.model_copy(
        update={"input": manifest.input.model_copy(update={"index_path": str(unrelated)})}
    )
    with pytest.raises(ValueError, match="adjacent"):
        selected_marlin_bam_index(altered)


def test_readiness_never_executes_unverified_runtime(tmp_path, monkeypatch):
    from ontseq_platform import marlin_native as native
    from ontseq_platform.execution import SubprocessRunner

    path, config = installation_fixture(tmp_path, monkeypatch)
    Path(config["python_executable"]).write_bytes(b"tampered executable, not yet verified")

    def unexpected_execution(*args, **kwargs):
        raise AssertionError("Readiness must never execute an unverified runtime")

    monkeypatch.setattr(native, "_preflight", unexpected_execution)
    monkeypatch.setattr(SubprocessRunner, "run", unexpected_execution)
    readiness = native.check_native_marlin_readiness(path, GenomeBuild.GRCH38)
    assert readiness.ready
    assert "not yet verified" in readiness.reason
    with pytest.raises(ValueError, match="checksum"):
        native.native_marlin_signature(path, GenomeBuild.GRCH38)


def test_runtime_directory_symlink_injection_refused(tmp_path, monkeypatch):
    from ontseq_platform.marlin_native import native_marlin_signature

    path, config = installation_fixture(tmp_path, monkeypatch)
    target = tmp_path / "unqualified_imports"
    target.mkdir()
    (target / "payload.py").write_text("# synthetic extra importable module")
    (Path(config["python_executable"]).parent.parent / "unexpected").symlink_to(
        target, target_is_directory=True
    )
    with pytest.raises(ValueError, match="directory.*link"):
        native_marlin_signature(path, GenomeBuild.GRCH38)

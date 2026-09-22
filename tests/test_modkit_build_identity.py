"""A reviewed scanner must be identified by its bytes, not its unchanged version banner."""

from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest
import test_methylation as methylation_tests
from test_methylation import _FakeRunner, _policy
from test_methylation_pipeline_integration import _Fixture

from ontseq_platform.methylation import run_methylation
from ontseq_platform.pipeline.runner import RunContext, _methylation_plan


def test_report_records_the_binary_that_actually_produced_the_pileup(tmp_path: Path) -> None:
    binary = tmp_path / "modkit"
    binary.write_bytes(b"synthetic command-runner executable version A")
    binary.chmod(0o755)
    manifest, intake = methylation_tests.AdapterTests()._fixture(tmp_path)
    report = run_methylation(
        manifest,
        intake,
        _policy(cpg_only=False, combine_strands=False),
        output_dir=tmp_path / "out",
        reference_fasta=tmp_path / "reference.fa",
        runner=_FakeRunner(),
        modkit=str(binary),
    )
    assert "modkit_binary_sha256" in report.tool.parameters
    assert (
        report.tool.parameters["modkit_binary_sha256"]
        == hashlib.sha256(binary.read_bytes()).hexdigest()
    )
    assert report.tool.parameters["modkit_build_id"] is None


def test_replacing_a_binary_invalidates_methylation_resume_even_with_same_version(
    tmp_path: Path,
) -> None:
    fixture = _Fixture(tmp_path)
    binary = tmp_path / "modkit"
    binary.write_bytes(b"synthetic command-runner executable version A")
    binary.chmod(0o755)
    config = replace(fixture.config, executables={"modkit": str(binary)})
    context = RunContext(
        config=config,
        envelope=fixture.envelope(),
        runner=_FakeRunner(),
        manifest=config.manifest,
    )
    before = _methylation_plan(context)
    binary.write_bytes(b"synthetic command-runner executable version B")
    after = _methylation_plan(context)
    assert before.tool_versions == after.tool_versions
    assert before.parameters != after.parameters


def test_only_a_reviewed_binary_can_process_independent_groups(tmp_path: Path) -> None:
    binary = tmp_path / "arbitrary-file-name"
    binary.write_bytes(b"synthetic reviewed executable")
    binary.chmod(0o755)
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    manifest, intake = methylation_tests.AdapterTests()._fixture(tmp_path)
    with patch("ontseq_platform.modkit_build.PR709_BINARY_SHA256", frozenset({digest})):
        report = run_methylation(
            manifest,
            intake,
            _policy(cpg_only=False, combine_strands=False),
            output_dir=tmp_path / "reviewed",
            reference_fasta=tmp_path / "reference.fa",
            runner=_FakeRunner(independent_cytosine_groups="2"),
            modkit=str(binary),
        )
    assert report.regions[0].modified_call_count == 15
    assert report.tool.parameters["modkit_build_id"] == "ontseq-modkit-pr709-v1"
    assert report.tool.parameters["modkit_binary_sha256"] == digest


@pytest.mark.parametrize("replace_at", ["view", "pileup"])
def test_a_binary_replaced_during_execution_cannot_produce_results(
    tmp_path: Path, replace_at: str
) -> None:
    binary = tmp_path / "modkit"
    binary.write_bytes(b"synthetic reviewed executable")
    binary.chmod(0o755)
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    manifest, intake = methylation_tests.AdapterTests()._fixture(tmp_path)

    class ReplacingRunner(_FakeRunner):
        def run(self, argv, *, timeout_seconds=300):
            if replace_at in argv:
                binary.write_bytes(b"unreviewed replacement with the same version banner")
            return super().run(argv, timeout_seconds=timeout_seconds)

    with (
        patch("ontseq_platform.modkit_build.PR709_BINARY_SHA256", frozenset({digest})),
        pytest.raises(ValueError, match="modkit executable changed"),
    ):
        run_methylation(
            manifest,
            intake,
            _policy(cpg_only=False, combine_strands=False),
            output_dir=tmp_path / "replaced",
            reference_fasta=tmp_path / "reference.fa",
            runner=ReplacingRunner(),
            modkit=str(binary),
        )
    assert not (tmp_path / "replaced/SYNTHETIC_001.modkit.bedmethyl").exists()


def test_wrong_planned_digest_stops_before_the_pileup(tmp_path: Path) -> None:
    binary = tmp_path / "modkit"
    binary.write_bytes(b"synthetic executable")
    binary.chmod(0o755)
    manifest, intake = methylation_tests.AdapterTests()._fixture(tmp_path)
    with pytest.raises(ValueError, match="modkit executable changed"):
        run_methylation(
            manifest,
            intake,
            _policy(cpg_only=False, combine_strands=False),
            output_dir=tmp_path / "mismatched",
            reference_fasta=tmp_path / "reference.fa",
            runner=_FakeRunner(),
            modkit=str(binary),
            expected_binary_sha256="0" * 64,
        )
    assert not (tmp_path / "mismatched/SYNTHETIC_001.modkit.bedmethyl").exists()

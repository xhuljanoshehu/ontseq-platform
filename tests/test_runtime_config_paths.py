from __future__ import annotations

import os
import shutil
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from ontseq_platform.cli import _parser as _legacy_parser
from ontseq_platform.profile_analysis import configuration_root
from ontseq_platform.runtime_cli import _parser

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = (ROOT / "configs").resolve()


@contextmanager
def _working_directory(path: Path) -> Iterator[None]:
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def test_configuration_root_does_not_trust_an_unrelated_cwd_configs_tree() -> None:
    with TemporaryDirectory() as raw:
        unrelated = Path(raw)
        decoy = unrelated / "configs" / "qc"
        decoy.mkdir(parents=True)
        (decoy / "defaults.yaml").write_text("not: an ontseq release\n", encoding="utf-8")

        with (
            patch("ontseq_platform.profile_analysis.sys.prefix", str(unrelated / "empty-prefix")),
            _working_directory(unrelated),
        ):
            resolved = configuration_root()

    assert resolved == CONFIGS


def test_serve_defaults_resolve_every_release_owned_input_outside_the_cwd() -> None:
    expected_configs = {
        "qc_policy": "qc/defaults.yaml",
        "sniffles_policy": "sv/sniffles2.conservative.technical.yaml",
        "cutesv_policy": "sv/cutesv.conservative.technical.yaml",
        "sv_consensus_policy": "sv/sniffles2_cutesv.consensus.technical.yaml",
        "sv_evidence_policy": "sv/evidence-priority.technical.yaml",
        "target_coverage_policy": "qc/adaptive_target_coverage.technical.yaml",
        "aml_knowledge": "knowledge/aml_rearrangements.v0.1.json",
        "aml_knowledge_lock": "knowledge/aml_rearrangements.v0.1.lock.json",
        "cnv_policy": "cnv/qdnaseq_ace.technical.yaml",
    }
    with TemporaryDirectory() as raw:
        temporary = Path(raw)
        asset_root = temporary / "packed" / "share" / "ontseq"
        shutil.copytree(CONFIGS, asset_root / "configs")
        (asset_root / "scripts").mkdir()
        shutil.copy2(ROOT / "scripts" / "run_qdnaseq_ace.R", asset_root / "scripts")
        with (
            patch("ontseq_platform.profile_analysis.sys.prefix", str(temporary / "packed")),
            _working_directory(temporary / "packed"),
        ):
            args = _parser().parse_args(
                [
                    "serve",
                    "--resource-root",
                    "resources",
                    "--allow-root",
                    "inputs",
                    "--no-browser",
                ]
            )
        for attribute, relative in expected_configs.items():
            path = getattr(args, attribute)
            assert path == asset_root / "configs" / relative
            assert path.is_file()
        assert args.qdnaseq_script == asset_root / "scripts" / "run_qdnaseq_ace.R"
        assert args.qdnaseq_script.is_file()


def test_run_defaults_include_alignment_basecalling_and_release_owned_knowledge() -> None:
    with TemporaryDirectory() as raw:
        unrelated = Path(raw)
        with (
            patch("ontseq_platform.profile_analysis.sys.prefix", str(unrelated / "empty-prefix")),
            _working_directory(unrelated),
        ):
            args = _parser().parse_args(
                [
                    "run",
                    "manifest.json",
                    "--reference-lock",
                    "reference.lock.json",
                    "--run-id",
                    "RUN",
                ]
            )

    assert args.alignment_policy == CONFIGS / "alignment" / "minimap2.ont.technical.yaml"
    assert args.basecall_policy == CONFIGS / "basecalling" / "dorado.technical.yaml"
    assert args.aml_knowledge == CONFIGS / "knowledge" / "aml_rearrangements.v0.1.json"
    assert args.aml_knowledge_lock == (CONFIGS / "knowledge" / "aml_rearrangements.v0.1.lock.json")
    assert args.alignment_policy.is_file()
    assert args.basecall_policy.is_file()


def test_legacy_local_smoke_defaults_are_release_owned_outside_the_cwd() -> None:
    with TemporaryDirectory() as raw:
        unrelated = Path(raw)
        with (
            patch("ontseq_platform.profile_analysis.sys.prefix", str(unrelated / "empty-prefix")),
            _working_directory(unrelated),
        ):
            args = _legacy_parser().parse_args(["local-smoke"])

    assert args.qc_policy == CONFIGS / "qc" / "defaults.yaml"
    assert args.sniffles_policy == CONFIGS / "sv" / "sniffles2.conservative.technical.yaml"
    assert args.qc_policy.is_file()
    assert args.sniffles_policy.is_file()


def test_explicit_relative_config_and_script_paths_remain_operator_paths() -> None:
    args = _parser().parse_args(
        [
            "serve",
            "--resource-root",
            "resources",
            "--allow-root",
            "inputs",
            "--cutesv-policy",
            "operator/cutesv.yaml",
            "--aml-knowledge",
            "operator/knowledge.json",
            "--qdnaseq-script",
            "operator/run.R",
            "--no-browser",
        ]
    )

    assert args.cutesv_policy == Path("operator/cutesv.yaml")
    assert args.aml_knowledge == Path("operator/knowledge.json")
    assert args.qdnaseq_script == Path("operator/run.R")

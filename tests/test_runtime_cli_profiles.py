from __future__ import annotations

import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from ontseq_platform.models import CuteSvPolicy, GenomeBuild
from ontseq_platform.resource_bootstrap import PROFILE_IDS
from ontseq_platform.runtime_cli import _cutesv_policy_for_run, _doctor_checks, main


class OptionalCuteSvPolicyTests(unittest.TestCase):
    def test_legacy_run_without_reference_does_not_implicitly_enable_cutesv(self) -> None:
        with patch("ontseq_platform.runtime_cli._cutesv_policy") as load_policy:
            policy = _cutesv_policy_for_run(Path("default-cutesv.yaml"), None)

        self.assertIsNone(policy)
        load_policy.assert_not_called()

    def test_run_with_reference_keeps_the_dual_caller_policy(self) -> None:
        expected = CuteSvPolicy(
            profile_id="synthetic-cutesv",
            status="technical_defaults_only",
            note="test",
        )
        with patch(
            "ontseq_platform.runtime_cli._cutesv_policy",
            return_value=expected,
        ) as load_policy:
            policy = _cutesv_policy_for_run(
                Path("default-cutesv.yaml"),
                Path("reference.fa"),
            )

        self.assertIs(policy, expected)
        load_policy.assert_called_once_with(Path("default-cutesv.yaml"))


class MethylationCommandIntegrationTests(unittest.TestCase):
    def test_analyze_forwards_only_explicit_methylation_opt_in(self) -> None:
        for include in (False, True):
            argv = [
                "ontseq",
                "analyze",
                "synthetic.bam",
                "--profile",
                PROFILE_IDS[0],
                "--modkit",
                "pinned-modkit",
                *(["--include-methylation"] if include else []),
            ]
            config = SimpleNamespace(
                manifest=SimpleNamespace(assay=SimpleNamespace(genome_build=GenomeBuild.GRCH38))
            )
            report = SimpleNamespace(
                stages=[], passed=True, verdict_reason="synthetic", unverified_stages=[]
            )
            with (
                self.subTest(include=include),
                patch("sys.argv", argv),
                patch("ontseq_platform.runtime_cli._register_cnv"),
                patch(
                    "ontseq_platform.runtime_cli.build_profile_run_configuration",
                    return_value=config,
                ) as build,
                patch("ontseq_platform.runtime_cli.run_pipeline", return_value=(report, None)),
                redirect_stdout(StringIO()),
            ):
                main()
            settings = build.call_args.args[0]
            self.assertIs(settings.include_methylation, include)
            self.assertEqual(settings.executables["modkit"], "pinned-modkit")

    def test_service_keeps_assay_default_and_explicit_tool_paths(self) -> None:
        argv = [
            "ontseq",
            "serve",
            "--resource-root",
            "resources",
            "--allow-root",
            "inputs",
            "--modkit",
            "pinned-modkit",
            "--samtools",
            "pinned-samtools",
            "--no-browser",
        ]
        with (
            patch("sys.argv", argv),
            patch("ontseq_platform.runtime_cli._register_cnv"),
            patch("ontseq_platform.runtime_cli.serve") as serve,
        ):
            main()
        config = serve.call_args.args[0]
        self.assertIsNone(config.methylation_policy)
        self.assertEqual(config.modkit_executable, "pinned-modkit")
        self.assertEqual(config.samtools_executable, "pinned-samtools")
        self.assertFalse(serve.call_args.kwargs["open_browser"])

    def test_doctor_uses_the_installed_configuration_root(self) -> None:
        with patch("ontseq_platform.runtime_cli.shutil.which", return_value="synthetic-tool"):
            result, code = _doctor_checks()
        self.assertEqual(code, 0)
        checks = {item["check"]: item for item in result["checks"]}
        self.assertTrue(Path(checks["project-root"]["detail"]).is_dir())
        self.assertEqual(checks["modkit"]["status"], "pass")


if __name__ == "__main__":
    unittest.main()

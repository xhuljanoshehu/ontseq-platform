"""The command surface: what ``ontseq`` offers, and what it does with a bad invocation.

Both parsers and the dispatcher had no test coverage at all, which is how
``validate-reference`` came to be a working command that ``ontseq`` never listed. The
checks here are the ones that do not need a genome, a tool or a subprocess: that the three
descriptions of the command set agree with each other, that dispatch sends a command to the
parser that owns it, and that a failure leaves through the documented exit path instead of
a traceback.

Nothing here runs a pipeline. What each command *does* is covered by the module tests for
the adapter behind it; this covers the layer between the operator and those.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import tempfile
import unittest
import unittest.mock
from pathlib import Path

from ontseq_platform import cli, entrypoint, runtime_cli
from ontseq_platform.runtime_cli import CNV_REGISTERING_COMMANDS, RUNTIME_COMMANDS


def _subparsers(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    """The subparser object behind each subcommand name."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return dict(action.choices)
    raise AssertionError("parser declares no subcommands")


def _option_strings(parser: argparse.ArgumentParser) -> set[str]:
    return {option for action in parser._actions for option in action.option_strings}


def _cnv_option_strings() -> set[str]:
    """What ``_add_cnv_options`` actually adds, asked of the function rather than listed.

    Derived so the test cannot drift from the code: adding an option there extends what
    every registering command must accept, without anyone having to remember to edit here.
    """
    probe = argparse.ArgumentParser()
    runtime_cli._add_cnv_options(probe)
    return _option_strings(probe) - {"-h", "--help"}


def _subcommands(parser: argparse.ArgumentParser) -> set[str]:
    """The subcommand names a parser accepts."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return set(action.choices)
    raise AssertionError("parser declares no subcommands")


class CommandSetTests(unittest.TestCase):
    """Three places describe the command set. A command missing from one is invisible.

    ``entrypoint`` lists commands for the overview, ``RUNTIME_COMMANDS`` decides which
    parser a command reaches, and the two parsers define what actually exists. Drift
    between them is silent in both directions: a command absent from the overview cannot
    be found, and one absent from the dispatch set falls through to the other parser and
    fails with an unrelated message.
    """

    def test_the_overview_lists_every_runtime_command(self) -> None:
        listed = {name for name, _ in entrypoint._RUNTIME_COMMANDS}
        self.assertEqual(listed, _subcommands(runtime_cli._parser()))

    def test_every_command_that_registers_cnv_accepts_its_options(self) -> None:
        """Registering the CNV lane and accepting its options must move together.

        ``_register_cnv`` reads ``args.cnv_policy`` straight off the namespace, so a
        command listed in ``CNV_REGISTERING_COMMANDS`` whose parser never declared the
        option raises ``AttributeError`` before the command does any work. That is exactly
        how adding ``preflight`` to the set broke every preflight invocation: unit tests
        never build the real namespace, so only CI caught it.
        """
        parsers = _subparsers(runtime_cli._parser())
        expected = _cnv_option_strings()
        self.assertTrue(expected, "the CNV options helper adds at least one option")
        for command in sorted(CNV_REGISTERING_COMMANDS):
            self.assertIn(command, parsers, f"{command} registers CNV but is not a command")
            missing = sorted(expected - _option_strings(parsers[command]))
            self.assertEqual(
                missing, [], f"{command} registers the CNV lane but does not accept {missing}"
            )

    def test_the_overview_lists_every_scientific_command(self) -> None:
        listed = {name for name, _ in entrypoint._SCIENTIFIC_COMMANDS}
        self.assertEqual(listed, _subcommands(cli._parser()))

    def test_dispatch_routes_exactly_the_runtime_parser_s_commands(self) -> None:
        self.assertEqual(RUNTIME_COMMANDS, _subcommands(runtime_cli._parser()))

    def test_the_two_parsers_share_no_command_name(self) -> None:
        """A shared name would be routed by the dispatch set alone, unreadably."""
        self.assertEqual(_subcommands(cli._parser()) & _subcommands(runtime_cli._parser()), set())

    def test_every_listed_command_has_a_summary(self) -> None:
        for name, summary in entrypoint._RUNTIME_COMMANDS + entrypoint._SCIENTIFIC_COMMANDS:
            with self.subTest(command=name):
                self.assertTrue(summary.strip(), f"{name} is listed without a summary")


class OverviewTests(unittest.TestCase):
    def test_a_bare_invocation_names_both_command_groups(self) -> None:
        overview = entrypoint._overview()
        self.assertIn("Execution and operations:", overview)
        self.assertIn("Single-step adapters and reporting:", overview)

    def test_the_overview_says_the_output_is_not_clinical(self) -> None:
        """The research-use statement is the first thing an operator sees."""
        self.assertIn("no output is a clinical result", entrypoint._overview())

    def test_every_command_appears_in_the_rendered_overview(self) -> None:
        overview = entrypoint._overview()
        for name, _ in entrypoint._RUNTIME_COMMANDS + entrypoint._SCIENTIFIC_COMMANDS:
            with self.subTest(command=name):
                self.assertIn(name, overview)

    def test_help_prints_the_overview_rather_than_dispatching(self) -> None:
        for argv in (["ontseq"], ["ontseq", "--help"], ["ontseq", "-h"], ["ontseq", "help"]):
            with self.subTest(argv=argv):
                captured = io.StringIO()
                with unittest.mock.patch("sys.argv", argv), contextlib.redirect_stdout(captured):
                    entrypoint.main()
                self.assertIn("usage: ontseq <command>", captured.getvalue())


class FailureExitTests(unittest.TestCase):
    """A bad invocation must exit through the documented path, not a traceback."""

    def _run(self, argv: list[str]) -> str:
        captured = io.StringIO()
        with (
            unittest.mock.patch("sys.argv", ["ontseq", *argv]),
            contextlib.redirect_stderr(captured),
            self.assertRaises(SystemExit) as raised,
        ):
            cli.main()
        return f"{raised.exception}{captured.getvalue()}"

    def test_a_missing_manifest_is_reported_rather_than_raised(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            message = self._run(["validate-manifest", str(Path(temporary) / "absent.json")])
        self.assertIn("ERROR:", message)

    def test_a_malformed_manifest_is_reported_rather_than_raised(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            manifest = Path(temporary) / "manifest.json"
            manifest.write_text('{"sample_id": "S"}', encoding="utf-8")
            message = self._run(["validate-manifest", str(manifest)])
        self.assertIn("ERROR:", message)

    def test_an_unparseable_document_is_reported_rather_than_raised(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            broken = Path(temporary) / "manifest.json"
            broken.write_text("{ not json", encoding="utf-8")
            message = self._run(["validate-manifest", str(broken)])
        self.assertIn("ERROR:", message)

    def test_an_unknown_command_does_not_reach_a_runtime_import(self) -> None:
        captured = io.StringIO()
        with (
            unittest.mock.patch("sys.argv", ["ontseq", "not-a-command"]),
            contextlib.redirect_stderr(captured),
            self.assertRaises(SystemExit),
        ):
            entrypoint.main()
        self.assertIn("invalid choice", captured.getvalue())


class DemoCommandTests(unittest.TestCase):
    """The one command that runs end to end without a genome, a tool or a real BAM."""

    def test_demo_writes_the_three_reviewer_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "demo"
            captured = io.StringIO()
            with (
                unittest.mock.patch("sys.argv", ["ontseq", "demo", "--output-dir", str(output)]),
                contextlib.redirect_stdout(captured),
            ):
                cli.main()

            written = {path.suffix for path in output.iterdir()}
            self.assertEqual(written, {".json", ".html", ".xlsx"})

    def test_the_demo_report_carries_the_research_use_banner(self) -> None:
        """A rendered report must never leave without saying what it is not."""
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "demo"
            with (
                unittest.mock.patch("sys.argv", ["ontseq", "demo", "--output-dir", str(output)]),
                contextlib.redirect_stdout(io.StringIO()),
            ):
                cli.main()

            report = next(output.glob("*.report.html")).read_text(encoding="utf-8")
            self.assertIn("RESEARCH USE ONLY", report)

            result = json.loads(next(output.glob("*.result.json")).read_text(encoding="utf-8"))
            self.assertIn("manifest", result)


if __name__ == "__main__":
    unittest.main()

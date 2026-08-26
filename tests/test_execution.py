from __future__ import annotations

import os
import subprocess
import sys
import unittest

from ontseq_platform.execution import SubprocessRunner, ToolExecutionError


class SubprocessRunnerTests(unittest.TestCase):
    def test_executes_argument_vector_without_shell(self) -> None:
        result = SubprocessRunner().run(
            [sys.executable, "-c", "import sys; print(sys.argv[1])", "value;not-a-command"]
        )

        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), "value;not-a-command")

    def test_missing_executable_is_a_typed_error(self) -> None:
        with self.assertRaises(ToolExecutionError):
            SubprocessRunner().run(["ontseq-executable-that-does-not-exist"])

    def test_rejects_nul_arguments(self) -> None:
        with self.assertRaises(ValueError):
            SubprocessRunner().run([sys.executable, "bad\x00argument"])

    def test_non_ascii_tool_output_decodes_under_any_locale(self) -> None:
        """A tool banner with a non-ASCII byte must not depend on the operator's locale.

        ``text=True`` would decode with ``locale.getpreferredencoding()``, which is ASCII
        under ``LC_ALL=C`` — a common setting in containers, cron jobs and freshly
        installed WSL distributions. A version probe would then raise mid-run on a
        character in a tool's own banner, and the run would fail on one machine and
        succeed on another from identical inputs.
        """
        emitter = "import sys; sys.stdout.buffer.write(b'tool 1.3.0 \\xe2\\x80\\x94 build\\n')"
        environment = {
            **os.environ,
            "LC_ALL": "C",
            "LANG": "C",
            "PYTHONUTF8": "0",
            "PYTHONCOERCECLOCALE": "0",
        }
        probe = (
            "from ontseq_platform.execution import SubprocessRunner;"
            "import sys;"
            f"print(SubprocessRunner().run([sys.executable, '-c', {emitter!r}])"
            ".stdout.encode('unicode_escape').decode())"
        )
        completed = subprocess.run(
            [sys.executable, "-c", probe],
            check=False,
            capture_output=True,
            encoding="utf-8",
            errors="replace",
            env=environment,
            timeout=60,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("tool 1.3.0", completed.stdout)
        self.assertIn(chr(0x2014).encode("unicode_escape").decode(), completed.stdout)

    def test_undecodable_tool_output_is_replaced_rather_than_raising(self) -> None:
        """Bytes that are not valid UTF-8 at all still have to come back as a string."""
        emitter = "import sys; sys.stdout.buffer.write(b'ok \\xff\\xfe end')"
        result = SubprocessRunner().run([sys.executable, "-c", emitter])

        self.assertEqual(result.returncode, 0)
        self.assertTrue(result.stdout.startswith("ok "))
        self.assertTrue(result.stdout.endswith(" end"))


if __name__ == "__main__":
    unittest.main()

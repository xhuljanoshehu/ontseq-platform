from __future__ import annotations

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


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from pathlib import Path

import post_merge_audit_apply as audit
import post_merge_audit_apply_v2 as audit_v2


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    audit_v2.main()

    path = ROOT / "tests/test_methylation.py"
    audit.replace_once(
        path,
        '''        class _NoFilterExpressions(_FakeRunner):
            def run(self, argv, *, timeout_seconds: int = 300) -> CommandResult:
                argv = [str(item) for item in argv]
                if "view" in argv:
                    self.commands.append(argv)
                    return CommandResult(tuple(argv), 1, "", "unrecognised option -e")
                return super().run(argv, timeout_seconds=timeout_seconds)
''',
        '''        class _NoFilterExpressions(_FakeRunner):
            def run(self, argv, *, timeout_seconds: int = 300) -> CommandResult:
                argv = [str(item) for item in argv]
                if "view" in argv:
                    self.commands.append(argv)
                    if any("C[+-]" in item for item in argv):
                        # The generic MM-tag-presence probe may be unavailable without
                        # making the independent same-base-group corruption-risk check
                        # unknowable. Those are deliberately separate safety questions.
                        return CommandResult(tuple(argv), 0, "0", "")
                    return CommandResult(tuple(argv), 1, "", "unrecognised option -e")
                return super().run(argv, timeout_seconds=timeout_seconds)
''',
        label="separate methylation probes in regression test",
    )


if __name__ == "__main__":
    main()

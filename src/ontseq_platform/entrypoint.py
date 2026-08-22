from __future__ import annotations

import sys

from .runtime_cli import RUNTIME_COMMANDS


_EXECUTION_COMMANDS = frozenset({"run", "serve", "watch"})


def main() -> None:
    """Dispatch execution commands without coupling the legacy/scientific CLI to runtime code."""
    command = sys.argv[1] if len(sys.argv) > 1 else None
    if command in RUNTIME_COMMANDS:
        if command in _EXECUTION_COMMANDS:
            from .runtime_extensions import register_builtin_runtime_extensions

            register_builtin_runtime_extensions()
        from .runtime_cli import main as runtime_main

        runtime_main()
    else:
        from .cli import main as legacy_main

        legacy_main()

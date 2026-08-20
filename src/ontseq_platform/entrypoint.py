from __future__ import annotations

import sys

from .runtime_cli import RUNTIME_COMMANDS


def main() -> None:
    """Dispatch execution commands without coupling the legacy/scientific CLI to runtime code."""
    command = sys.argv[1] if len(sys.argv) > 1 else None
    if command in RUNTIME_COMMANDS:
        from .runtime_cli import main as runtime_main

        runtime_main()
    else:
        from .cli import main as legacy_main

        legacy_main()

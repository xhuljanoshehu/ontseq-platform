"""The local browser interface: an HTTP front door onto ``ontseq run``.

``guard``
    Everything a wrong answer makes dangerous — token comparison, the allowed-root
    boundary, the anti-rebinding host check, and Windows/WSL path translation.
    Dependency-free and unit tested, deliberately: this is the part that must not be
    first executed on a runner.
``app``
    The HTTP transport. Serves the page, starts runs, reports progress. Computes nothing.
"""

from __future__ import annotations

__all__ = ["app", "guard"]

"""A deterministic UTC clock for synthetic experiments, independent of host clock steps."""

from datetime import UTC, datetime, timedelta
from itertools import count
from unittest.mock import patch

SYNTHETIC_REGISTERED_AT = datetime(2026, 9, 5, tzinfo=UTC)


def install_experiment_clock(test_case, *modules):
    ticks = count(1)
    for module in modules:
        replacement = patch.object(module, "datetime", wraps=datetime)
        clock = replacement.start()
        test_case.addCleanup(replacement.stop)
        clock.now.side_effect = lambda tz: SYNTHETIC_REGISTERED_AT + timedelta(seconds=next(ticks))

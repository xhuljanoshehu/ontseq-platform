"""A deterministic UTC clock for synthetic experiments, independent of host clock steps."""

from datetime import UTC, datetime, timedelta
from itertools import count
from unittest.mock import patch

SYNTHETIC_REGISTERED_AT = datetime(2026, 9, 5, tzinfo=UTC)


def install_experiment_clock(test_case, runner_module):
    ticks = count(1)
    replacement = patch.object(runner_module, "datetime", wraps=datetime)
    clock = replacement.start()
    test_case.addCleanup(replacement.stop)
    clock.now.side_effect = lambda tz: SYNTHETIC_REGISTERED_AT + timedelta(seconds=next(ticks))

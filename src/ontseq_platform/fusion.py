from __future__ import annotations

import json
from collections import defaultdict
from enum import StrEnum
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from .models import (
    EventType,
    Evidence,
    GenomeBuild,
    GenomicEvent,
    ModuleRunStatus,
    SnifflesCallReport,
    StrictModel,
)


class ObservabilityStatus(StrEnum):
    OBSERVABLE = "observable"
    LIMITED = "limited"
    NOT
"""Desktop application contracts and backend orchestration.

The Qt user interface is intentionally not imported here so core/server installations do not
require PySide6.
"""

from .backend import (
    DesktopAnalysisRequest,
    DesktopBackend,
    DesktopBackendError,
    DesktopRunResult,
    DesktopStage,
    DesktopStageStatus,
    DiagnosticCheck,
    ProgressEvent,
    locate_bam_index,
    sanitize_sample_id,
)
from .config import (
    DesktopConfig,
    DesktopReferenceProfile,
    default_config_path,
    load_desktop_config,
    save_desktop_config,
)

__all__ = [
    "DesktopAnalysisRequest",
    "DesktopBackend",
    "DesktopBackendError",
    "DesktopConfig",
    "DesktopReferenceProfile",
    "DesktopRunResult",
    "DesktopStage",
    "DesktopStageStatus",
    "DiagnosticCheck",
    "ProgressEvent",
    "default_config_path",
    "load_desktop_config",
    "locate_bam_index",
    "sanitize_sample_id",
    "save_desktop_config",
]

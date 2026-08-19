from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from ..models import GenomeBuild, StrictModel


class DesktopReferenceProfile(StrictModel):
    """Deployment-local reference configuration.

    These paths point to local/on-premises resources only. They must never contain
    patient data or direct identifiers.
    """

    genome_build: GenomeBuild
    reference_id: str = Field(min_length=1)
    reference_lock_path: str = Field(min_length=1)
    adaptive_sampling_target_bed_path: str | None = None
    adaptive_sampling_target_bed_version: str | None = None

    @model_validator(mode="after")
    def target_bed_metadata_is_complete(self) -> DesktopReferenceProfile:
        if bool(self.adaptive_sampling_target_bed_path) != bool(
            self.adaptive_sampling_target_bed_version
        ):
            raise ValueError(
                "adaptive sampling target BED path and version must be configured together"
            )
        return self


class DesktopConfig(StrictModel):
    schema_version: Literal["0.1.0"] = "0.1.0"
    backend_mode: Literal["wsl", "local"] = "wsl"
    wsl_distribution: str | None = None
    wsl_project_root: str = "~/ontseq-platform"
    output_root: str = Field(default_factory=lambda: str(Path.home() / "ONTSeq Results"))
    qc_policy_path: str = "configs/qc/defaults.yaml"
    sniffles_policy_path: str = "configs/sv/sniffles2.conservative.technical.yaml"
    reference_profiles: list[DesktopReferenceProfile] = Field(default_factory=list)

    @model_validator(mode="after")
    def one_profile_per_genome_build(self) -> DesktopConfig:
        builds = [profile.genome_build for profile in self.reference_profiles]
        if len(builds) != len(set(builds)):
            raise ValueError("only one desktop reference profile per genome build is allowed")
        return self

    def reference_for(self, genome_build: GenomeBuild) -> DesktopReferenceProfile:
        for profile in self.reference_profiles:
            if profile.genome_build == genome_build:
                return profile
        raise ValueError(f"No desktop reference profile configured for {genome_build.value}")


def default_config_path() -> Path:
    if os.name == "nt":
        base = Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming"))
        return base / "ONTSeq Desktop" / "config.json"
    xdg = os.environ.get("XDG_CONFIG_HOME")
    base = Path(xdg) if xdg else Path.home() / ".config"
    return base / "ontseq-desktop" / "config.json"


def load_desktop_config(path: Path | None = None) -> DesktopConfig:
    config_path = path or default_config_path()
    if not config_path.exists():
        return DesktopConfig()
    payload = json.loads(config_path.read_text(encoding="utf-8"))
    return DesktopConfig.model_validate(payload)


def save_desktop_config(config: DesktopConfig, path: Path | None = None) -> Path:
    config_path = path or default_config_path()
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        json.dumps(config.model_dump(mode="json"), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return config_path

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, TypeVar

import yaml
from pydantic import BaseModel

ModelT = TypeVar("ModelT", bound=BaseModel)


def load_mapping(path: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    with path.open("r", encoding="utf-8") as handle:
        if suffix in {".yaml", ".yml"}:
            value = yaml.safe_load(handle)
        elif suffix == ".json":
            value = json.load(handle)
        else:
            raise ValueError(f"Unsupported input format: {path}")
    if not isinstance(value, dict):
        raise ValueError(f"Expected a mapping at document root: {path}")
    return value


def load_model(path: Path, model: type[ModelT]) -> ModelT:
    return model.model_validate(load_mapping(path))


def write_json(model: BaseModel, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(model.model_dump_json(indent=2), encoding="utf-8")
    return path

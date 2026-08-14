from __future__ import annotations

import argparse
import json
from pathlib import Path

from ontseq_platform.models import PipelineResult, SampleManifest


def _render() -> dict[Path, str]:
    return {
        Path("schemas/sample-manifest.schema.json"): json.dumps(
            SampleManifest.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
        Path("schemas/pipeline-result.schema.json"): json.dumps(
            PipelineResult.model_json_schema(), indent=2, sort_keys=True
        )
        + "\n",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    mismatches: list[str] = []
    for path, content in _render().items():
        if args.check:
            if not path.exists() or path.read_text(encoding="utf-8") != content:
                mismatches.append(str(path))
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
            print(path)
    if mismatches:
        raise SystemExit(f"Schemas are stale: {', '.join(mismatches)}")


if __name__ == "__main__":
    main()

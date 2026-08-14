from __future__ import annotations

import argparse
from pathlib import Path

from pydantic import ValidationError

from .demo import build_demo_result
from .io import load_model, write_json
from .models import PipelineResult, SampleManifest
from .report import render_html
from .workbook import render_workbook


def _render(result: PipelineResult, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = result.manifest.sample_id
    outputs = [
        write_json(result, output_dir / f"{stem}.result.json"),
        render_html(result, output_dir / f"{stem}.report.html"),
        render_workbook(result, output_dir / f"{stem}.results.xlsx"),
    ]
    return outputs


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ontseq", description="ONTSeq Platform core CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    demo = subparsers.add_parser("demo", help="Generate a synthetic HTML/Excel/JSON report")
    demo.add_argument("--output-dir", type=Path, default=Path("results/demo"))

    validate_manifest = subparsers.add_parser("validate-manifest")
    validate_manifest.add_argument("path", type=Path)

    validate_result = subparsers.add_parser("validate-result")
    validate_result.add_argument("path", type=Path)

    render = subparsers.add_parser("render", help="Render HTML and Excel from a result JSON")
    render.add_argument("result", type=Path)
    render.add_argument("--output-dir", type=Path, required=True)
    return parser


def main() -> None:
    args = _parser().parse_args()
    try:
        if args.command == "demo":
            for path in _render(build_demo_result(), args.output_dir):
                print(path)
        elif args.command == "validate-manifest":
            manifest = load_model(args.path, SampleManifest)
            print(f"VALID manifest: {manifest.sample_id}")
        elif args.command == "validate-result":
            result = load_model(args.path, PipelineResult)
            print(f"VALID result: {result.manifest.sample_id}")
        elif args.command == "render":
            result = load_model(args.result, PipelineResult)
            for path in _render(result, args.output_dir):
                print(path)
    except (OSError, ValueError, ValidationError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

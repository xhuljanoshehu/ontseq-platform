"""Run-bound, offline presentation contract for the integrated Befund interface.

No input paths, example patients, or classification decisions enter this adapter.
The static audit report remains available if JavaScript cannot run.
"""

from __future__ import annotations

import base64
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from .cnv.qdnaseq import QDNAseqCallReport
from .marlin_native_contracts import NativeMarlinReport
from .models import PipelineResult
from .report_formatting import redact_paths
from .report_marlin import validate_marlin_identity
from .report_view import build_report_view

REPORT_PRESENTATION_VERSION = "befund-interactive-v1"


def read_plot(root: Path, name: str) -> str | None:
    """Only embed PNG evidence from inside the explicit CNV evidence directory."""
    path = (root / name).resolve()
    if Path(name).is_absolute() or not path.is_relative_to(root.resolve()):
        raise ValueError("CNV plot must remain inside the evidence directory")
    if not path.is_file():
        return None
    content = path.read_bytes()
    if not content.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError("CNV plot must be a PNG")
    return "data:image/png;base64," + base64.b64encode(content).decode("ascii")


def _redact(value: Any) -> Any:
    if isinstance(value, str):
        return redact_paths(value)
    if isinstance(value, dict):
        return {str(key): _redact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_redact(item) for item in value]
    return value


def interactive_payload(
    result: PipelineResult,
    cnv: QDNAseqCallReport | None = None,
    evidence_root: Path | None = None,
    *,
    marlin_report: NativeMarlinReport | None = None,
) -> dict[str, Any]:
    if cnv is not None and (
        cnv.sample_id != result.manifest.sample_id
        or cnv.genome_build != result.manifest.assay.genome_build
    ):
        raise ValueError("CNV report sample/build identity differs from the report")
    validate_marlin_identity(result, marlin_report)
    cnv_data = None
    plots: dict[str, dict[str, str | None]] = {}
    if cnv is not None:
        fields = {
            "bin_size_kbp",
            "cellularity",
            "ploidy",
            "fit_error",
            "candidate_count",
            "segment_count",
            "alternatives",
        }
        cnv_data = {
            "primary_fit": cnv.primary_fit.model_dump(include=fields),
            "fits": [fit.model_dump(include=fields) for fit in cnv.fits],
            "chromosomes": [row.model_dump() for row in cnv.chromosome_consensus],
            "warnings": cnv.warnings,
            "limitations": cnv.limitations,
        }
        if evidence_root is not None:
            for fit in cnv.fits:
                plots[str(fit.bin_size_kbp)] = {
                    "copy_number": read_plot(evidence_root, fit.copy_number_plot),
                    "fit": read_plot(evidence_root, fit.fit_plot),
                }
    data = _redact(
        {
            "contract": REPORT_PRESENTATION_VERSION,
            "view": asdict(build_report_view(result)),
            "requested_modules": [module.value for module in result.manifest.analysis.modules],
            "iscn": result.iscn.model_dump(mode="json"),
            "cnv": cnv_data,
            "marlin": marlin_report.model_dump(mode="json") if marlin_report else None,
        }
    )
    # Data URIs are created here, never read from uncontrolled report text.
    data["plots"] = plots
    return dict(data)


def attach_interactive_report(
    document: str,
    result: PipelineResult,
    *,
    cnv: QDNAseqCallReport | None = None,
    evidence_root: Path | None = None,
    marlin_report: NativeMarlinReport | None = None,
) -> str:
    payload = interactive_payload(result, cnv, evidence_root, marlin_report=marlin_report)
    serialized = json.dumps(payload, ensure_ascii=True, allow_nan=False)
    serialized = serialized.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    bundle = Path(__file__).with_name("report_bundle.html").read_text(encoding="utf-8")
    document = document.replace(
        '<meta name="ontseq-report-layout" content="befund-v8-local-1">',
        f'<meta name="ontseq-report-layout" content="{REPORT_PRESENTATION_VERSION}">',
    )
    document = document.replace(
        "<body>", '<body><div id="ontseq-befund-root"></div><div id="ontseq-static-report">', 1
    )
    return document.replace(
        "</body>",
        '</div><script id="ontseq-befund-data" type="application/json">'
        + serialized
        + "</script>"
        + bundle
        + "</body>",
        1,
    )

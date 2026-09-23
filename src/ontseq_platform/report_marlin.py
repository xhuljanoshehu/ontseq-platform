"""Presentation of recorded MARLIN evidence without adding classification decisions."""

from __future__ import annotations

import json
import re

from .marlin_native_contracts import NativeMarlinReport
from .models import PipelineResult
from .report_formatting import cell, redact_paths

ABSENT_MARLIN = "Für diesen Lauf kein MARLIN-Ausführungsnachweis vorhanden."
MARLIN_RESEARCH_NOTE = (
    "UNVALIDATED_RESEARCH · Nur für Forschungszwecke; nicht klinisch validiert. "
    "Modellscores sind keine Erkrankungswahrscheinlichkeiten."
)


def validate_marlin_identity(result: PipelineResult, report: NativeMarlinReport | None) -> None:
    if report is not None and (
        report.sample_id != result.manifest.sample_id
        or report.genome_build != result.manifest.assay.genome_build
        or report.run_id != result.manifest.run_id
    ):
        raise ValueError("MARLIN report sample/build/run identity differs from the report")
    if report is None:
        return
    outcomes = [item for item in result.modules if item.module.value == "marlin"]
    if len(outcomes) != 1 or (outcomes[0].status, outcomes[0].reason, outcomes[0].tools) != (
        report.status,
        report.reason,
        report.tools,
    ):
        raise ValueError("MARLIN report contradicts the normalized module outcome")
    digest = result.provenance.reference_checksums.get("marlin_report")
    if report.status.value in {"COMPLETED", "NO_CALL"}:
        # Byte-level verification belongs to the artifact loader: presentation receives
        # a typed model, not the original serialized bytes (including whitespace).
        if digest is None or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise ValueError("MARLIN concluded report requires artifact checksum provenance")
    elif digest is not None:
        raise ValueError("MARLIN failed/blocked outcome cannot claim a prediction artifact")


def marlin_facts(report: NativeMarlinReport | None) -> list[list[object]]:
    feature = report.feature_summary if report else None
    return [
        ["Status", report.status.value if report else "NOT_RECORDED"],
        ["Reason", report.reason if report else ABSENT_MARLIN],
        ["Decision", report.decision.value if report and report.decision else None],
        ["Observed model CpGs", feature.observed_model_feature_count if feature else None],
        ["Expected model CpGs", feature.expected_feature_count if feature else None],
        ["Observed model CpG fraction", feature.observed_fraction if feature else None],
        ["Explicit NA CpGs", feature.explicit_na_feature_count if feature else None],
        ["Absent model CpGs", feature.absent_feature_count if feature else None],
        ["Leading model class", report.top_class if report else None],
        ["Leading model score", report.top_class_score if report else None],
        ["Raw model score threshold", report.model_score_threshold if report else None],
        ["Model score threshold met", report.model_score_threshold_met if report else None],
        ["Assay assessability", report.assay_assessability if report else None],
        ["Assessment", "Bewertung nicht validiert" if report else None],
        ["Validation", report.validation_status if report else None],
        ["Research use", MARLIN_RESEARCH_NOTE],
    ]


def marlin_html_section(report: NativeMarlinReport | None) -> str:
    section = '<section id="marlin"><h2>Methylierungsklassifikation · MARLIN</h2>'
    if report is None:
        return section + f"<p>{ABSENT_MARLIN}</p></section>"
    section += f"<p>{MARLIN_RESEARCH_NOTE}</p>"
    section += (
        "<p>Bewertung nicht validiert · Die rohe Modellscore-Schwelle ist keine "
        "validierte Grenze für die Beurteilbarkeit dieser Probe. Die Zahl und der Anteil "
        "beobachteter CpGs sind technische Nachweise; eine belastbare "
        "Assay-Beurteilbarkeit ist nicht etabliert.</p>"
    )
    if report.decision and report.decision.value == "UNKNOWN":
        section += (
            "<p>UNKNOWN · Keine hinreichend sichere Klassifikation. "
            "Eine führende Modellklasse ist keine bestätigte Klasse.</p>"
        )
    section += "<table><thead><tr><th>Merkmal</th><th>Wert</th></tr></thead><tbody>"
    section += (
        "".join(
            f"<tr><th>{cell(label)}</th>"
            f"<td>{cell(value) if value is not None else 'nicht verfügbar'}</td></tr>"
            for label, value in marlin_facts(report)
        )
        + "</tbody></table>"
    )
    for title, scores in (
        ("Klassen", report.class_scores),
        ("Familien", report.family_scores),
        ("Linien", report.lineage_scores),
    ):
        if scores:
            section += f"<h3>{title} · Modellscores</h3><table><thead><tr>"
            section += "<th>Rang</th><th>Bezeichnung</th><th>Modellscore</th></tr></thead><tbody>"
            section += (
                "".join(
                    f"<tr><td>{rank}</td><td>{cell(row.label)}</td><td>{row.score:.6f}</td></tr>"
                    for rank, row in enumerate(
                        sorted(scores, key=lambda row: (-row.score, row.label)), 1
                    )
                )
                + "</tbody></table>"
            )
    section += "".join(f"<p>{cell(text)}</p>" for text in [*report.warnings, *report.limitations])
    provenance = {
        "schema_version": report.schema_version,
        "adapter_version": report.adapter_version,
        "run_id": report.run_id,
        "tools": [tool.model_dump(mode="json") for tool in report.tools],
        "parameters": report.parameters,
        "input_fingerprints": {
            key: value.model_dump() for key, value in report.input_fingerprints.items()
        },
        "installation_signature": report.installation_signature,
        "feature_summary": report.feature_summary.model_dump() if report.feature_summary else None,
    }
    section += "<details><summary>MARLIN-Werkzeuge, Parameter und Prüfsummen</summary><pre>"
    section += cell(redact_paths(json.dumps(provenance, ensure_ascii=False, indent=2)))
    return section + "</pre></details></section>"

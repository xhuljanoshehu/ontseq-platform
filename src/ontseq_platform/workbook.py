from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.worksheet import Worksheet

from .models import PipelineResult

HEADER_FILL = PatternFill("solid", fgColor="075985")
HEADER_FONT = Font(color="FFFFFF", bold=True)
WARNING_FILL = PatternFill("solid", fgColor="FEE2E2")


def _style_sheet(sheet: Worksheet) -> None:
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    for cell in sheet[1]:
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(horizontal="center")
    for column in sheet.columns:
        width = min(max(len(str(cell.value or "")) for cell in column) + 2, 60)
        column_index = column[0].column
        if column_index is not None:
            sheet.column_dimensions[get_column_letter(column_index)].width = width


def _write_table(sheet: Worksheet, headers: list[str], rows: Iterable[list[object]]) -> None:
    sheet.append(headers)
    for row in rows:
        sheet.append(row)
    _style_sheet(sheet)


def render_workbook(result: PipelineResult, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    workbook = Workbook()
    active_sheet = workbook.active
    if active_sheet is not None:
        workbook.remove(active_sheet)

    summary = workbook.create_sheet("00_Summary")
    _write_table(
        summary,
        ["Field", "Value"],
        [
            ["Sample ID", result.manifest.sample_id],
            ["Run ID", result.manifest.run_id],
            ["Assay", result.manifest.assay.mode.value],
            ["Genome build", result.manifest.assay.genome_build.value],
            ["QC verdict", result.qc.verdict.value],
            ["Release status", result.release_status.value],
            ["ISCN proposal", result.iscn.notation],
            ["ISCN profile", result.iscn.conformance_profile],
        ],
    )

    qc = workbook.create_sheet("01_QC")
    _write_table(
        qc,
        ["Metric", "Value"],
        [[key, value] for key, value in result.qc.metrics.items()],
    )

    cnv_chr = workbook.create_sheet("02_CNV_Chromosomes")
    chromosome_events = [
        event for event in result.events if "chromosome_" in event.event_type.value
    ]
    _write_table(
        cnv_chr,
        ["event_id", "event_type", "chromosome", "copy_number", "confidence", "reportable"],
        [
            [
                event.event_id,
                event.event_type.value,
                event.primary.chromosome,
                event.copy_number,
                event.confidence,
                event.reportable,
            ]
            for event in chromosome_events
        ],
    )

    cnv_segments = workbook.create_sheet("03_CNV_Segments")
    segment_events = [
        event
        for event in result.events
        if event.event_type.value in {"deletion", "duplication"} and event.copy_number is not None
    ]
    _write_table(
        cnv_segments,
        ["event_id", "type", "chromosome", "start", "end", "band_start", "band_end", "copy_number"],
        [
            [
                event.event_id,
                event.event_type.value,
                event.primary.chromosome,
                event.primary.start,
                event.primary.end,
                event.primary.cytoband_start,
                event.primary.cytoband_end,
                event.copy_number,
            ]
            for event in segment_events
        ],
    )

    sv = workbook.create_sheet("04_SV")
    sv_events = [
        event
        for event in result.events
        if event.event_type.value in {"inversion", "translocation", "insertion"}
        or (event.event_type.value in {"deletion", "duplication"} and event.copy_number is None)
    ]
    _write_table(
        sv,
        [
            "event_id",
            "type",
            "length_bp",
            "locus_1",
            "locus_2",
            "confidence",
            "reportable",
            "evidence",
        ],
        [
            [
                event.event_id,
                event.event_type.value,
                event.length_bp,
                f"{event.primary.chromosome}:{event.primary.start}-{event.primary.end}",
                ""
                if not event.secondary
                else (
                    f"{event.secondary.chromosome}:{event.secondary.start}-{event.secondary.end}"
                ),
                event.confidence,
                event.reportable,
                json.dumps([item.model_dump(mode="json") for item in event.evidence]),
            ]
            for event in sv_events
        ],
    )

    fusions = workbook.create_sheet("05_Fusions")
    fusion_events = [event for event in result.events if event.event_type.value == "fusion"]
    _write_table(
        fusions,
        ["event_id", "genes", "breakpoint_1", "breakpoint_2", "confidence", "reportable"],
        [
            [
                event.event_id,
                "::".join(event.genes),
                f"{event.primary.chromosome}:{event.primary.start}",
                ""
                if not event.secondary
                else f"{event.secondary.chromosome}:{event.secondary.start}",
                event.confidence,
                event.reportable,
            ]
            for event in fusion_events
        ],
    )

    iscn = workbook.create_sheet("06_ISCN")
    _write_table(
        iscn,
        ["Field", "Value"],
        [
            ["Notation", result.iscn.notation],
            ["Edition", result.iscn.standard_edition],
            ["Conformance profile", result.iscn.conformance_profile],
            ["Review status", result.iscn.review_status.value],
            ["Source event IDs", ", ".join(result.iscn.source_event_ids)],
        ],
    )

    warnings = workbook.create_sheet("07_Warnings")
    warning_rows: list[list[object]] = [
        [item] for item in result.warnings + result.qc.warnings + result.iscn.warnings
    ]
    _write_table(warnings, ["Warning"], warning_rows)
    for cell in warnings["A"][1:]:
        cell.fill = WARNING_FILL

    methods = workbook.create_sheet("08_Methods_Versions")
    _write_table(
        methods,
        ["Tool", "Version", "Parameters", "Container digest"],
        [
            [
                tool.name,
                tool.version,
                json.dumps(tool.parameters, sort_keys=True),
                tool.container_digest,
            ]
            for tool in result.provenance.tools
        ],
    )

    run_log = workbook.create_sheet("09_Run_Log")
    _write_table(
        run_log,
        ["Field", "Value"],
        [
            ["Pipeline version", result.provenance.pipeline_version],
            ["Git commit", result.provenance.git_commit],
            ["Created at", result.provenance.created_at.isoformat()],
            [
                "Reference checksums",
                json.dumps(result.provenance.reference_checksums, sort_keys=True),
            ],
        ],
    )

    modules = workbook.create_sheet("10_Module_Status")
    _write_table(
        modules,
        ["Module", "Status", "Reason", "Tools"],
        [
            [
                outcome.module.value,
                outcome.status.value,
                outcome.reason,
                ", ".join(f"{tool.name} {tool.version}" for tool in outcome.tools),
            ]
            for outcome in result.modules
        ],
    )

    workbook.save(output_path)
    return output_path

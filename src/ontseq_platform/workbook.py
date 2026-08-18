from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.worksheet import Worksheet

from .models import EventAnnotation, GenomicEvent, PipelineResult

HEADER_FILL = PatternFill("solid", fgColor="075985")
HEADER_FONT = Font(color="FFFFFF", bold=True)
WARNING_FILL = PatternFill("solid", fgColor="FEE2E2")

#: Printed above the annotation table, not merely beside it. The HTML report can put this
#: in red prose the eye passes on the way to the table; a spreadsheet cannot. A grid with a
#: column headed "assertion" containing "Pathogenic" one cell away from an event identifier
#: is read as a verdict about the sample unless something in the sheet itself says it is not.
ANNOTATION_BANNER = (
    "These rows are classifications of database records, not findings about this sample. "
    "Each assertion is stated in the vocabulary of the source that made it: a germline "
    "classification does not become a statement about a somatic finding by appearing "
    "here. Nothing in this sheet makes a finding reportable, and no column here was used "
    "to compute confidence or reportability anywhere in the pipeline. Read the "
    "scope_alignment and caveats columns before reading the assertion."
)


def _style_header(sheet: Worksheet, row: int) -> None:
    for cell in sheet[row]:
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(horizontal="center")


def _fit_columns(sheet: Worksheet, *, from_row: int) -> None:
    """Size each column to its widest cell, ignoring rows above ``from_row``.

    The annotation sheet's banner spans every column, so measuring from row 1 there would
    size the first column to the length of a paragraph.
    """
    for column in sheet.iter_cols(min_row=from_row):
        width = min(max(len(str(cell.value or "")) for cell in column) + 2, 60)
        index = column[0].column
        if index is not None:
            sheet.column_dimensions[get_column_letter(index)].width = width


def _style_sheet(sheet: Worksheet) -> None:
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    _style_header(sheet, 1)
    _fit_columns(sheet, from_row=1)


def _write_table(sheet: Worksheet, headers: list[str], rows: Iterable[list[object]]) -> None:
    sheet.append(headers)
    for row in rows:
        sheet.append(row)
    _style_sheet(sheet)


def _annotation_count(event: GenomicEvent) -> int:
    """How many database records matched this event.

    A count, deliberately, and named for what it counts. A reviewer who never opens the
    annotation sheet would otherwise not know there was anything there to open, and a
    count cannot be misread as a classification of the finding.
    """
    return len(event.annotations)


def _stars(annotation: EventAnnotation) -> object:
    """The star rating, or the word for not knowing it.

    ``None`` means the review status was not in the vocabulary this build knows, which is
    a different thing from a status that earns no stars. Writing ``0`` for both would erase
    the difference in the one direction that flatters the record.
    """
    return "unknown" if annotation.review_stars is None else annotation.review_stars


def _write_annotations(sheet: Worksheet, result: PipelineResult) -> None:
    """Write the annotation table under a banner that states how to read it.

    The banner occupies row 1 and the header row 2, so the freeze pane and the filter are
    one row lower than on every other sheet. That is the cost of having the reading rule
    travel with the data: a sheet can be exported, copied into another workbook or printed
    on its own, and the caveat has to survive all three.
    """
    headers = [
        "event_id",
        "source_id",
        "source_release",
        "source_sha256",
        "record_id",
        "record_type",
        "assertion",
        "assertion_vocabulary",
        "record_origin",
        "scope_alignment",
        "scope_note",
        "match_type",
        "reciprocal_overlap",
        "review_status",
        "review_stars",
        "genes",
        "conditions",
        "caveats",
    ]
    sheet.append([ANNOTATION_BANNER])
    sheet.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(headers))
    banner = sheet.cell(row=1, column=1)
    banner.fill = WARNING_FILL
    banner.font = Font(bold=True)
    banner.alignment = Alignment(wrap_text=True, vertical="top")
    sheet.row_dimensions[1].height = 60

    sheet.append(headers)
    mismatched: list[int] = []
    for event in result.events:
        for annotation in event.annotations:
            sheet.append(
                [
                    event.event_id,
                    annotation.source_id,
                    annotation.source_release,
                    annotation.source_sha256,
                    annotation.record_id,
                    annotation.record_type,
                    annotation.assertion,
                    annotation.assertion_vocabulary,
                    annotation.record_origin,
                    annotation.scope_alignment,
                    annotation.scope_note,
                    annotation.match_type,
                    annotation.reciprocal_overlap,
                    annotation.review_status,
                    _stars(annotation),
                    ", ".join(annotation.genes),
                    ", ".join(annotation.conditions),
                    " ".join(annotation.caveats),
                ]
            )
            if annotation.scope_alignment == "mismatched":
                mismatched.append(sheet.max_row)

    # A mismatched record is kept and marked, never dropped: a germline assertion matched
    # against a somatic question is information about the match, and hiding it would leave
    # the reviewer unable to see why an assertion looks relevant when it is not.
    for row_index in mismatched:
        for position in range(1, len(headers) + 1):
            sheet.cell(row=row_index, column=position).fill = WARNING_FILL

    sheet.freeze_panes = "A3"
    last_column = get_column_letter(len(headers))
    sheet.auto_filter.ref = f"A2:{last_column}{max(sheet.max_row, 2)}"
    _style_header(sheet, 2)
    _fit_columns(sheet, from_row=2)


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
        [
            "event_id",
            "event_type",
            "chromosome",
            "copy_number",
            "confidence",
            "reportable",
            "db_records_matched",
        ],
        [
            [
                event.event_id,
                event.event_type.value,
                event.primary.chromosome,
                event.copy_number,
                event.confidence,
                event.reportable,
                _annotation_count(event),
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
        [
            "event_id",
            "type",
            "chromosome",
            "start",
            "end",
            "band_start",
            "band_end",
            "copy_number",
            "db_records_matched",
        ],
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
                _annotation_count(event),
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
            "db_records_matched",
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
                _annotation_count(event),
            ]
            for event in sv_events
        ],
    )

    fusions = workbook.create_sheet("05_Fusions")
    fusion_events = [event for event in result.events if event.event_type.value == "fusion"]
    _write_table(
        fusions,
        [
            "event_id",
            "genes",
            "breakpoint_1",
            "breakpoint_2",
            "confidence",
            "reportable",
            "db_records_matched",
        ],
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
                _annotation_count(event),
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

    _write_annotations(workbook.create_sheet("11_Annotations"), result)

    workbook.save(output_path)
    return output_path

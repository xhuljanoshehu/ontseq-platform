"""Typed detail sections retained across the portable-report presentation redesign."""

from __future__ import annotations

from collections.abc import Iterable

from .models import GenomicEvent, PipelineResult, ResolvedResourceContext
from .report_formatting import cell as _cell
from .report_formatting import optional_number
from .report_plots import CoverageBar, coverage_depth_svg
from .reporting import (
    fusion_assessment,
    fusion_review_events,
    gene_pair_label,
    is_structural_variant,
    pathology_label,
    release_state,
    review_priority,
)
from .sv_evidence import sv_review_queue
from .target_coverage import TargetCoverageReport


def _release_label(event: GenomicEvent) -> str:
    state = release_state(event)
    if event.reportable:
        return "REPORTABLE — pipeline flag only; this RUO report is not clinically validated"
    return state + " — release gate not satisfied; not a biological negative result"


def _event_evidence(event: GenomicEvent) -> str:
    return ", ".join(
        f"{item.caller} {item.caller_version} "
        f"(support={item.support_reads}, vaf={item.variant_allele_fraction})"
        for item in event.evidence
    )


def _caller_support(event: GenomicEvent) -> str:
    return ", ".join(
        f"{item.caller}={item.support_reads if item.support_reads is not None else 'n/a'}"
        for item in event.evidence
    )


def _cytobands(event: GenomicEvent) -> str:
    bands = [event.primary.cytoband_start]
    if event.secondary is not None:
        bands.append(event.secondary.cytoband_start)
    return " ↔ ".join(band or "unannotated" for band in bands)


def _local_coverage(event: GenomicEvent) -> str:
    if not event.breakpoint_mean_depths:
        return "not measured"
    return " / ".join(
        "n/a" if depth is None else f"{depth:.1f}×" for depth in event.breakpoint_mean_depths
    )


def _event_loci(event: GenomicEvent) -> tuple[str, str]:
    primary = f"{event.primary.chromosome}:{event.primary.start:,}-{event.primary.end:,}"
    secondary = (
        ""
        if event.secondary is None
        else f"{event.secondary.chromosome}:{event.secondary.start:,}-{event.secondary.end:,}"
    )
    return primary, secondary


def _breakpoint_transcript_context(event: GenomicEvent) -> str:
    parts: list[str] = []
    for annotation in event.breakpoint_annotations:
        preferred = next((item for item in annotation.transcripts if item.preferred), None)
        if preferred is None and annotation.transcripts:
            preferred = annotation.transcripts[0]
        if preferred is None:
            parts.append(f"{annotation.label}: unannotated")
            continue
        location = str(preferred.region)
        if preferred.exon_number is not None:
            location += f" {preferred.exon_number}"
        elif preferred.intron_number is not None:
            location += f" {preferred.intron_number}"
        phase = "" if preferred.cds_phase is None else f", CDS phase {preferred.cds_phase}"
        parts.append(
            f"{annotation.label}: {preferred.gene_name}/{preferred.transcript_id} "
            f"({location}{phase})"
        )
    return "; ".join(parts) or "not annotated"


def _fusion_context(event: GenomicEvent) -> str:
    evidence = event.fusion_evidence
    if evidence is None:
        return "n/a"
    return f"orientation={evidence.orientation or 'unknown'}; frame={evidence.frame_status}"


def _reference_methods(result: PipelineResult) -> tuple[str, str]:
    context = result.reference_context
    if not isinstance(context, ResolvedResourceContext):
        return "<tr><td>Reference context</td><td>legacy_unspecified</td></tr>", ""
    releases = context.resource_releases
    rows = [
        ("Genome assembly", releases.get("reference.genome_fasta", context.genome_build.value)),
        ("ReferenceBundle", f"{context.reference_bundle_id} ({context.reference_bundle_version})"),
        ("BAM dictionary contract", context.reference_dictionary_contract.value),
        ("GENCODE", releases.get("reference.gencode_gtf", "unspecified")),
        ("MANE", releases.get("reference.mane_gff3", "unspecified")),
        ("Cytobands", releases.get("reference.cytobands", "unspecified")),
        (
            "PanelBundle",
            f"{context.panel_bundle_id} ({context.panel_bundle_version})"
            if context.panel_bundle_id is not None
            else "NOT_APPLICABLE",
        ),
        ("KnowledgeBundle", f"{context.knowledge_bundle_id} ({context.knowledge_bundle_version})"),
    ]
    table_rows = "".join(
        f"<tr><td>{_cell(name)}</td><td>{_cell(value)}</td></tr>" for name, value in rows
    )
    checksum_rows = "".join(
        f"<tr><td>{_cell(name)}</td><td><code>{_cell(checksum)}</code></td></tr>"
        for name, checksum in sorted(context.resource_checksums.items())
    )
    return table_rows, checksum_rows


def _panel_resolution_methods(result: PipelineResult) -> str:
    context = result.reference_context
    if not isinstance(context, ResolvedResourceContext) or context.panel_bundle_id is None:
        return ""
    summary = context.panel_resolution
    rows: list[tuple[str, object]]
    if summary is None:
        rows = [("Panel resolution", "Not recorded in this legacy result")]
    elif summary.mapping_status == "native_build_not_required":
        rows = [
            ("Coordinate state", "Native build; coordinate mapping not required"),
            ("Panel build", summary.target_genome_build.value),
            ("Selection intervals", summary.selection_interval_count),
            ("Native analysis ROI intervals", summary.analysis_roi_interval_count),
            (
                "Unresolved/review panel labels",
                ", ".join(summary.unresolved_target_labels) or "No labels recorded",
            ),
        ]
    else:
        rows = [
            ("Coordinate state", summary.mapping_status),
            ("Mapping ID", summary.mapping_id),
            (
                "Source panel",
                f"{summary.source_panel_bundle_id}:{summary.source_panel_resource_id}",
            ),
            (
                "Build mapping",
                f"{summary.source_genome_build.value} to {summary.target_genome_build.value}",
            ),
            ("Mapping method", f"{summary.mapping_tool} ({summary.mapping_method})"),
            (
                "Forward-mapped intervals",
                f"{summary.mapped_interval_count} of {summary.source_interval_count}",
            ),
            ("Final selection intervals", summary.selection_interval_count),
            ("Native analysis ROI intervals", summary.analysis_roi_interval_count),
            (
                "Forward-unmapped panel labels",
                ", ".join(summary.unmapped_target_labels) or "No labels recorded",
            ),
            (
                "Reciprocal-exact intervals",
                f"{summary.reciprocal_exact_interval_count} of "
                f"{summary.roundtrip_mapped_interval_count} roundtrip-mapped",
            ),
            (
                "Roundtrip review-required labels",
                ", ".join(summary.roundtrip_review_required_target_labels) or "No labels recorded",
            ),
            (
                "Unresolved/review panel labels",
                ", ".join(summary.unresolved_target_labels) or "No labels recorded",
            ),
            ("Runtime coordinate mapping", "PROHIBITED"),
        ]
    table_rows = "".join(
        f"<tr><td>{_cell(name)}</td><td>{_cell(value)}</td></tr>" for name, value in rows
    )
    return (
        "<h3>Panel coordinate resolution</h3>"
        "<div class='notice warn'><strong>Research-use panel provenance:</strong> Coordinate "
        "mapping records how the panel design was resolved; mapped does not mean analytically "
        "validated. Unmapped and unresolved/review-required panel targets are not negative "
        "findings for this sample.</div>"
        "<div class='table-wrap'><table><thead><tr><th>Panel field</th><th>Value</th></tr>"
        f"</thead><tbody>{table_rows}</tbody></table></div>"
    )


def _fusion_row(event: GenomicEvent) -> str:
    primary, secondary = _event_loci(event)
    return (
        "<tr>"
        f"<td><strong>{_cell(gene_pair_label(event))}</strong><br>"
        f"<span class='muted'>{_cell(event.event_id)}</span></td>"
        f"<td>{_cell(primary)}<br>{_cell(secondary)}</td>"
        f"<td>{_cell(_cytobands(event))}</td>"
        f"<td>{_cell(_caller_support(event))}</td>"
        f"<td>{_cell(_local_coverage(event))}</td>"
        f"<td>{_cell(event.observability.value)}</td>"
        f"<td>{_cell(_breakpoint_transcript_context(event))}</td>"
        f"<td>{_cell(_fusion_context(event))}</td>"
        f"<td>{_cell(event.known_rearrangement or 'no knowledge match')}</td>"
        f"<td>{_cell(pathology_label(event))}</td>"
        f"<td>{_cell(fusion_assessment(event))}</td>"
        f"<td>{_cell(event.confidence)}</td>"
        f"<td>{_cell(_release_label(event))}</td>"
        "</tr>"
    )


def _review_event_row(event: GenomicEvent) -> str:
    primary, secondary = _event_loci(event)
    return (
        "<tr>"
        f"<td>{_cell(event.event_id)}</td><td>{_cell(event.event_type.value)}</td>"
        f"<td>{_cell(gene_pair_label(event))}</td>"
        f"<td>{_cell(primary)}</td><td>{_cell(secondary)}</td>"
        f"<td>{_cell(_cytobands(event))}</td><td>{_cell(_caller_support(event))}</td>"
        f"<td>{_cell(_local_coverage(event))}</td><td>{_cell(event.observability.value)}</td>"
        f"<td>{_cell(', '.join(event.technical_flags))}</td>"
        f"<td>{_cell(event.known_rearrangement or '')}</td>"
        f"<td>{_cell(pathology_label(event))}</td>"
        f"<td>{_cell(fusion_assessment(event))}</td>"
        f"<td><strong>{_cell(event.confidence)}</strong></td>"
        f"<td>{_cell(_release_label(event))}</td></tr>"
    )


def _full_event_row(event: GenomicEvent) -> str:
    primary, secondary = _event_loci(event)
    return (
        "<tr>"
        f"<td>{_cell(event.event_id)}</td><td>{_cell(event.event_type.value)}</td>"
        f"<td>{_cell(event.length_bp)}</td><td>{_cell(primary)}</td><td>{_cell(secondary)}</td>"
        f"<td>{_cell(gene_pair_label(event))}</td><td>{_cell(event.confidence)}</td>"
        f"<td>{_cell(review_priority(event))}</td><td>{_cell(_release_label(event))}</td>"
        f"<td>{_cell(event.validation_status.value)}</td><td>{_cell(event.observability.value)}</td>"
        f"<td>{_cell(', '.join(event.technical_flags))}</td>"
        f"<td>{_cell(event.known_rearrangement or '')}</td>"
        f"<td>{_cell(pathology_label(event))}</td>"
        f"<td>{_cell(fusion_assessment(event))}</td><td>{_cell(_event_evidence(event))}</td></tr>"
    )


def table(caption: str, headers: list[str], rows: Iterable[list[object]]) -> str:
    heading = "".join(f"<th>{_cell(value)}</th>" for value in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{_cell(value)}</td>" for value in row) + "</tr>" for row in rows
    )
    if not body:
        body = f"<tr><td colspan='{len(headers)}'>No entries recorded.</td></tr>"
    return (
        f"<div class='table-wrap'><table><caption>{_cell(caption)}</caption>"
        f"<thead><tr>{heading}</tr></thead><tbody>{body}</tbody></table></div>"
    )


def resource_details(result: PipelineResult) -> str:
    rows, checksums = _reference_methods(result)
    return (
        "<div class='table-wrap'><table><caption>Resource releases</caption>"
        "<thead><tr><th>Resource</th><th>Release / contract</th></tr></thead>"
        f"<tbody>{rows}</tbody></table></div>"
        + _panel_resolution_methods(result)
        + "<div class='table-wrap'><table><caption>Resource SHA256 provenance</caption>"
        "<thead><tr><th>Resource</th><th>SHA256</th></tr></thead>"
        f"<tbody>{checksums or '<tr><td colspan=2>Not recorded.</td></tr>'}</tbody></table></div>"
    )


def _coverage_figure(report: TargetCoverageReport, label: str) -> str:
    """Bar plot of the report's own per-target means; nothing is recomputed here."""
    bars = [
        CoverageBar(
            region_id=region.region_id,
            chromosome=region.chromosome,
            start=region.start,
            mean_depth=region.mean_depth,
        )
        for region in report.regions
    ]
    reference_depths = [value for value in report.policy.thresholds if value >= 10]
    svg = coverage_depth_svg(
        bars,
        title=f"{label}: per-target mean depth",
        reference_depths=reference_depths,
    )
    if not svg:
        return ""
    labels = "/".join(f"{value}×" for value in reference_depths) or "no"
    caption = (
        f"{_cell(label)}: mean depth per target in genome order; dashed lines mark the "
        f"coverage policy's descriptive {labels} bins. A rendering of the normalized "
        "values only — descriptive technical evidence, not an adequacy or reportability "
        "assessment."
    )
    return (
        "<figure class='plot' style='margin:14px 0'>"
        + svg
        + f"<figcaption style='color:#5e687a;font-size:12px;margin-top:6px'>{caption}"
        "</figcaption></figure>"
    )


def coverage_section(
    target: TargetCoverageReport | None, selection: TargetCoverageReport | None
) -> str:
    pieces = ["<section id='coverage'><h2>Adaptive-sampling target coverage</h2>"]
    if target is None:
        pieces.append(
            "<p>No analysis-ROI coverage sidecar was supplied. Coverage is not assessed.</p>"
        )
    else:
        metrics = target.summary_metrics
        cards = [
            ("Targets assessed", _cell(metrics.get("region_count"))),
            ("Target-weighted mean", optional_number(metrics.get("interval_weighted_mean_depth"))),
            ("Median target mean", optional_number(metrics.get("median_region_mean_depth"))),
            ("Least-covered target", optional_number(metrics.get("minimum_region_mean_depth"))),
            (
                "Target bases ≥20×",
                optional_number(metrics.get("interval_bases_at_20x_fraction"), percent=True),
            ),
            (
                "Buffered selection mean",
                optional_number(
                    selection.summary_metrics.get("interval_weighted_mean_depth")
                    if selection
                    else None
                ),
            ),
        ]
        pieces.append(
            "<div class='identity'>"
            + "".join(
                f"<div><span>{_cell(name)}</span><strong>{value}</strong></div>"
                for name, value in cards
            )
            + "</div>"
        )
    # A selection-only sidecar must not disappear just because analysis-ROI data are absent.
    for report, label in ((target, "Analysis ROI"), (selection, "Buffered selection")):
        if report is None:
            continue
        pieces.append(
            f"<h3>{label}</h3><p>Role: {_cell(report.target_bed_role.value)}; "
            f"design: {_cell(report.target_bed_version)}; status: {_cell(report.status.value)}</p>"
        )
        figure = _coverage_figure(report, label)
        if figure:
            pieces.append(figure)
        pieces.append(
            table(
                label + " coverage",
                [
                    "Target",
                    "Chromosome",
                    "Start (0-based)",
                    "End (exclusive)",
                    "Mean depth",
                    "Bases ≥20×",
                    "Bases ≥30×",
                ],
                [
                    [
                        r.region_id,
                        r.chromosome,
                        r.start,
                        r.end,
                        optional_number(r.mean_depth, places=2),
                        optional_number(r.fraction_at_threshold.get("20x"), percent=True),
                        optional_number(r.fraction_at_threshold.get("30x"), percent=True),
                    ]
                    for r in sorted(
                        report.regions, key=lambda r: (r.mean_depth, r.chromosome, r.start)
                    )
                ],
            )
        )
        pieces.append(
            "<ul>"
            + "".join(f"<li>{_cell(w)}</li>" for w in report.warnings + report.limitations)
            + "</ul>"
        )
    pieces.append(
        "<p class='muted'>Coverage values are descriptive technical evidence, not validated "
        "assay-adequacy thresholds. Missing measurement is not zero coverage or a negative "
        "finding.</p></section>"
    )
    return "".join(pieces)


def sv_details(result: PipelineResult) -> str:
    events = [e for e in result.events if is_structural_variant(e)]
    queue = sv_review_queue(events, limit=max(1, len(events)))
    rows = "".join(_review_event_row(e) for e in queue)
    if not rows:
        rows = (
            "<tr><td colspan='15'>No high/moderate technical-priority SV candidate "
            "was identified.</td></tr>"
        )
    headers = [
        "Event",
        "Type",
        "Genes",
        "Locus 1",
        "Locus 2",
        "Cytobands",
        "Support",
        "Depth",
        "Observability",
        "Technical flags",
        "Known rearrangement",
        "Pathology",
        "Fusion assessment",
        "Confidence",
        "Release state",
    ]
    header = "".join(f"<th>{h}</th>" for h in headers)
    full = "".join(_full_event_row(e) for e in events)
    if not full:
        full = (
            "<tr><td colspan='16'>No structural-variant event was produced. Review module status; "
            "this is not a biological negative result.</td></tr>"
        )
    full_header = "".join(
        f"<th>{h}</th>"
        for h in [
            "Event",
            "Type",
            "Length",
            "Locus 1",
            "Locus 2",
            "Genes",
            "Confidence",
            "Priority",
            "Release state",
            "Validation status",
            "Observability",
            "Flags",
            "Known rearrangement",
            "Pathology",
            "Fusion assessment",
            "Evidence",
        ]
    )
    fusion_rows = "".join(_fusion_row(e) for e in fusion_review_events(events))
    fusion_headers = "".join(
        f"<th>{h}</th>"
        for h in [
            "Genes",
            "Loci",
            "Cytobands",
            "Support",
            "Depth",
            "Observability",
            "Transcript context",
            "Orientation / frame",
            "Known rearrangement",
            "Pathology",
            "Fusion assessment",
            "Confidence",
            "Release state",
        ]
    )
    return (
        "<section id='sv-review'><h2>SV review queue</h2>"
        "<label for='sv-filter'>Filter review queue</label> "
        "<input type='search' id='sv-filter' aria-controls='sv-review-table' "
        "placeholder='Gene, event, locus…'>"
        f"<div class='table-wrap'><table id='sv-review-table'><thead><tr>{header}</tr></thead>"
        f"<tbody>{rows}</tbody></table></div>"
        "<h3>Fusion / rearrangement candidates</h3><p>Neither state asserts an expressed, "
        "in-frame or "
        "functional fusion transcript. Annotation and caller agreement are not "
        "independent validation.</p>"
        f"<div class='table-wrap'><table><thead><tr>{fusion_headers}</tr></thead><tbody>"
        f"{fusion_rows or '<tr><td colspan=13>No annotated candidate available.</td></tr>'}"
        "</tbody></table></div>"
        "</section><section id='technical-appendix'><h2>Technical appendix</h2>"
        "<p>All normalized structural events are retained below, regardless of review priority.</p>"
        f"<div class='table-wrap'><table><thead><tr>{full_header}</tr></thead>"
        f"<tbody>{full}</tbody></table></div>"
        "</section><script>document.getElementById('sv-filter').addEventListener('input', "
        "function() {"
        "const query = this.value.toLowerCase(); "
        "document.querySelectorAll('#sv-review-table tbody tr')"
        ".forEach(row => { row.hidden = !row.textContent.toLowerCase().includes(query); }); "
        "});</script>"
    )


def iscn_details(result: PipelineResult) -> str:
    proposal = result.iscn
    pieces = ["<h3>Technical ISCN proposal — fachzytogenetisch zu prüfen</h3>"]
    pieces.append(
        table(
            "Assessment contract",
            ["Field", "Value"],
            [
                ["Proposal status", proposal.proposal_status.value],
                [
                    "Genome build",
                    proposal.genome_build.value if proposal.genome_build else "legacy_unspecified",
                ],
                ["Clinical release allowed", proposal.clinical_release_allowed],
                ["Requires expert review", proposal.requires_expert_review],
                [
                    "Selection policy",
                    proposal.selection_policy.value
                    if proposal.selection_policy
                    else "legacy_unspecified",
                ],
            ],
        )
    )
    pieces.append(
        table(
            "Event fragments",
            ["Source event", "Fragment", "Event type", "Confidence", "Selection reason"],
            [
                [f.event_id, f.fragment, f.event_type.value, f.confidence, f.selection_reason]
                for f in proposal.event_fragments
            ],
        )
    )
    pieces.append(
        table(
            "Event dispositions",
            ["Source event", "Outcome", "Reason", "Detail"],
            [
                [d.event_id, d.outcome.value, d.reason_code, d.detail]
                for d in proposal.event_dispositions
            ],
        )
    )
    pieces.append(
        table(
            "Assessment blockers",
            ["Reason", "Detail"],
            [[b.reason_code, b.detail] for b in proposal.assessment_blockers],
        )
    )
    pieces.append(
        table(
            "Policy parameters",
            ["Parameter", "Value"],
            [[k, v] for k, v in sorted(proposal.policy_parameters.items())],
        )
    )
    if proposal.resource_provenance is not None:
        pieces.append(
            table(
                "ISCN resource provenance",
                ["Field", "Value"],
                [[k, v] for k, v in proposal.resource_provenance.model_dump(mode="json").items()],
            )
        )
    pieces.append(
        "<ul>"
        + "".join(
            f"<li>{_cell(w)}</li>" for w in proposal.technical_assumptions + proposal.warnings
        )
        + "</ul>"
    )
    return "".join(pieces)

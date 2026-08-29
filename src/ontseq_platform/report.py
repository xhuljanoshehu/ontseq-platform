# ruff: noqa: E501

from __future__ import annotations

import html
import json
from pathlib import Path

from .models import AnalysisModule, GenomicEvent, PipelineResult, ResolvedResourceContext
from .reporting import (
    caller_count,
    fusion_assessment,
    fusion_review_events,
    gene_pair_label,
    is_structural_variant,
    key_findings,
    maximum_support,
    pathology_label,
    release_state,
    review_priority,
)
from .sv_evidence import sv_review_queue
from .target_coverage import TargetCoverageReport


def _cell(value: object) -> str:
    return html.escape("" if value is None else str(value))


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


def _status_pill(value: str, *, css_class: str = "neutral") -> str:
    return f"<span class='pill {css_class}'>{_cell(value)}</span>"


def _key_finding_row(event: GenomicEvent) -> str:
    primary, secondary = _event_loci(event)
    priority = review_priority(event)
    priority_class = "critical" if priority == "HEMATOLOGY_REVIEW" else "review"
    confidence_class = "good" if event.confidence == "high" else "warn"
    return (
        f"<tr class='{priority_class}'>"
        f"<td>{_status_pill(priority, css_class=priority_class)}</td>"
        f"<td><strong>{_cell(gene_pair_label(event))}</strong><br>"
        f"<span class='muted'>{_cell(event.event_id)}</span></td>"
        f"<td>{_cell(event.event_type.value)}</td>"
        f"<td>{_cell(primary)}<br>{_cell(secondary)}</td>"
        f"<td>{_cell(_cytobands(event))}</td>"
        f"<td>{_status_pill(event.confidence, css_class=confidence_class)}</td>"
        f"<td>{caller_count(event)} caller(s); max support {maximum_support(event)}<br>"
        f"<span class='muted'>{_cell(_caller_support(event))}</span></td>"
        f"<td>{_cell(_local_coverage(event))}</td>"
        f"<td>{_cell(pathology_label(event))}</td>"
        f"<td>{_cell(fusion_assessment(event))}</td>"
        f"<td>{_status_pill(release_state(event), css_class='locked')}</td>"
        "</tr>"
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
        f"<td>{_cell(release_state(event))}</td>"
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
        f"<td>{_cell(release_state(event))}</td></tr>"
    )


def _full_event_row(event: GenomicEvent) -> str:
    primary, secondary = _event_loci(event)
    return (
        "<tr>"
        f"<td>{_cell(event.event_id)}</td><td>{_cell(event.event_type.value)}</td>"
        f"<td>{_cell(event.length_bp)}</td><td>{_cell(primary)}</td><td>{_cell(secondary)}</td>"
        f"<td>{_cell(gene_pair_label(event))}</td><td>{_cell(event.confidence)}</td>"
        f"<td>{_cell(review_priority(event))}</td><td>{_cell(release_state(event))}</td>"
        f"<td>{_cell(event.validation_status.value)}</td><td>{_cell(event.observability.value)}</td>"
        f"<td>{_cell(', '.join(event.technical_flags))}</td>"
        f"<td>{_cell(event.known_rearrangement or '')}</td>"
        f"<td>{_cell(pathology_label(event))}</td>"
        f"<td>{_cell(fusion_assessment(event))}</td><td>{_cell(_event_evidence(event))}</td></tr>"
    )


def _coverage_section(
    target: TargetCoverageReport | None,
    selection: TargetCoverageReport | None,
) -> str:
    if target is None:
        return (
            "<section id='coverage'><div class='section-heading'><div><span class='eyebrow'>Assay</span>"
            "<h2>Adaptive-sampling target coverage</h2></div></div>"
            "<div class='notice warn'>No target-coverage sidecar was supplied to the renderer. "
            "Coverage is not assessed in this report.</div></section>"
        )
    metrics = target.summary_metrics
    low_regions = sorted(target.regions, key=lambda item: item.mean_depth)[:10]
    low_rows = "".join(
        "<tr>"
        f"<td>{_cell(region.region_id)}</td><td>{_cell(region.chromosome)}</td>"
        f"<td>{region.mean_depth:.2f}×</td>"
        f"<td>{100 * region.fraction_at_threshold.get('20x', 0):.1f}%</td>"
        f"<td>{100 * region.fraction_at_threshold.get('30x', 0):.1f}%</td></tr>"
        for region in low_regions
    )
    selection_mean = (
        selection.summary_metrics.get("interval_weighted_mean_depth") if selection else None
    )
    selection_mean_text = "n/a" if selection_mean is None else f"{float(selection_mean):.1f}×"
    return f"""
    <section id="coverage">
      <div class="section-heading"><div><span class="eyebrow">Assay</span>
        <h2>Adaptive-sampling target coverage</h2></div>
        <span class="pill neutral">{_cell(target.target_bed_role.value)}</span></div>
      <div class="grid metrics">
        <div class="card"><span>Targets assessed</span><strong>{_cell(metrics.get("region_count"))}</strong></div>
        <div class="card"><span>Target-weighted mean</span><strong>{float(metrics.get("interval_weighted_mean_depth", 0)):.1f}×</strong></div>
        <div class="card"><span>Median target mean</span><strong>{float(metrics.get("median_region_mean_depth", 0)):.1f}×</strong></div>
        <div class="card"><span>Least-covered target</span><strong>{float(metrics.get("minimum_region_mean_depth", 0)):.1f}×</strong></div>
        <div class="card"><span>Target bases ≥20×</span><strong>{100 * float(metrics.get("interval_bases_at_20x_fraction", 0)):.1f}%</strong></div>
        <div class="card"><span>Buffered selection mean</span><strong>{selection_mean_text}</strong></div>
      </div>
      <p class="muted">Coverage values are descriptive technical evidence. The table lists the ten
      least-covered analysis targets; low coverage is not a biological negative result.</p>
      <div class="table-wrap"><table><thead><tr><th>Target</th><th>Chromosome</th>
      <th>Mean depth</th><th>Bases ≥20×</th><th>Bases ≥30×</th></tr></thead>
      <tbody>{low_rows}</tbody></table></div>
    </section>"""


def _module_status(result: PipelineResult, module: AnalysisModule) -> str:
    outcome = next((item for item in result.modules if item.module == module), None)
    return outcome.status.value if outcome is not None else "NOT_RECORDED"


def render_html(
    result: PipelineResult,
    output_path: Path,
    *,
    target_coverage: TargetCoverageReport | None = None,
    selection_coverage: TargetCoverageReport | None = None,
) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    structural_events = [event for event in result.events if is_structural_variant(event)]
    finding_events = key_findings(structural_events)
    fusion_events = fusion_review_events(structural_events)
    review_events = sv_review_queue(structural_events, limit=max(1, len(structural_events)))
    review_ids = {event.event_id for event in review_events}
    finding_ids = {event.event_id for event in finding_events}
    background_events = [
        event
        for event in structural_events
        if event.event_id not in review_ids and event.event_id not in finding_ids
    ][:200]
    html_event_ids = review_ids | finding_ids | {event.event_id for event in background_events}
    html_events = [event for event in structural_events if event.event_id in html_event_ids]

    finding_rows = "".join(_key_finding_row(event) for event in finding_events) or (
        "<tr><td colspan='11'>No hematology knowledge match or high/moderate fusion-evidence "
        "candidate was identified. This is not a biological negative result.</td></tr>"
    )
    fusion_rows = "".join(_fusion_row(event) for event in fusion_events) or (
        "<tr><td colspan='13'>No annotated fusion/rearrangement candidate was available for "
        "assessment. This is not a biological negative result.</td></tr>"
    )
    review_rows = "".join(_review_event_row(event) for event in review_events) or (
        "<tr><td colspan='15'>No high/moderate technical-priority SV candidate was identified. "
        "This is not a biological negative result.</td></tr>"
    )
    event_rows = "".join(_full_event_row(event) for event in html_events) or (
        "<tr><td colspan='16'>No structural-variant event was produced. Review module status; "
        "this is not a biological negative result.</td></tr>"
    )

    preferred_metrics = [
        "number_of_reads",
        "aligned_percent",
        "total_yield_gb",
        "mean_coverage_x",
        "n50_bp",
        "median_length_bp",
    ]
    metric_cards = "".join(
        f"<div class='card'><span>{_cell(key.replace('_', ' ').title())}</span>"
        f"<strong>{_cell(result.qc.metrics[key])}</strong></div>"
        for key in preferred_metrics
        if key in result.qc.metrics
    )
    confidence_counts = {
        confidence: sum(event.confidence == confidence for event in structural_events)
        for confidence in ("high", "moderate", "low", "unclassified")
    }
    raw_sv_count = sum(max(1, len(event.source_event_ids)) for event in structural_events)
    warnings = "".join(
        f"<li>{_cell(item)}</li>"
        for item in result.warnings + result.qc.warnings + result.iscn.warnings
    )
    module_rows = (
        "".join(
            f"<tr><td>{_cell(module.module.value)}</td><td>{_cell(module.status.value)}</td>"
            f"<td>{_cell(module.reason)}</td></tr>"
            for module in result.modules
        )
        or "<tr><td colspan='3'>No module outcomes were recorded.</td></tr>"
    )
    tool_rows = "".join(
        f"<tr><td>{_cell(tool.name)}</td><td>{_cell(tool.version)}</td>"
        f"<td><code>{_cell(json.dumps(tool.parameters, sort_keys=True))}</code></td></tr>"
        for tool in result.provenance.tools
    )
    reference_rows, checksum_rows = _reference_methods(result)
    adaptive_warning = (
        "<div class='notice warn'><strong>Adaptive-sampling CNV caution:</strong> genome-wide "
        "read depth is enrichment-biased. CNV output is exploratory until an assay-matched "
        "normalization and benchmark are available.</div>"
        if result.manifest.assay.mode.value == "adaptive_sampling"
        and any(event.copy_number is not None for event in result.events)
        else ""
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ONTSeq report - {_cell(result.manifest.sample_id)}</title>
  <style>
    :root {{ color-scheme:light; --ink:#10233f; --muted:#62718a; --line:#d9e2ec;
      --brand:#0b557a; --brand-dark:#08374f; --soft:#eef6fa; --page:#f3f6f9;
      --critical:#9f1239; --critical-bg:#fff1f2; --review:#9a4d0a; --review-bg:#fff7ed;
      --good:#166534; --good-bg:#ecfdf3; --locked:#475569; --locked-bg:#f1f5f9; }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; font:15px/1.5 Inter,Segoe UI,system-ui,sans-serif; color:var(--ink);
      background:var(--page); }}
    .banner {{ background:#7f1d1d; color:white; padding:9px 18px; text-align:center;
      font-weight:750; letter-spacing:.045em; font-size:13px; }}
    main {{ max-width:1440px; margin:0 auto; padding:28px clamp(16px,3vw,42px) 64px; }}
    header {{ background:linear-gradient(125deg,var(--brand-dark),var(--brand)); color:white;
      padding:30px clamp(24px,4vw,48px); border-radius:22px; box-shadow:0 14px 40px #08374f24; }}
    header h1 {{ margin:0 0 8px; font-size:clamp(30px,4vw,48px); letter-spacing:-.035em; }}
    header p {{ margin:0; opacity:.88; }}
    nav {{ position:sticky; top:0; z-index:5; display:flex; gap:8px; overflow-x:auto;
      padding:12px 2px; background:#f3f6f9f2; backdrop-filter:blur(8px); }}
    nav a {{ color:var(--brand-dark); text-decoration:none; background:white; border:1px solid var(--line);
      border-radius:999px; padding:7px 12px; white-space:nowrap; font-size:13px; font-weight:650; }}
    section {{ margin-top:18px; background:white; border:1px solid var(--line); border-radius:18px;
      padding:clamp(18px,2.5vw,30px); box-shadow:0 5px 20px #0f172a08; }}
    .section-heading {{ display:flex; align-items:flex-start; justify-content:space-between; gap:16px;
      margin-bottom:16px; }}
    h2 {{ margin:0; font-size:clamp(23px,2.4vw,32px); letter-spacing:-.025em; }}
    h3 {{ margin:24px 0 10px; }}
    .eyebrow {{ color:var(--brand); text-transform:uppercase; font-size:12px; letter-spacing:.12em;
      font-weight:800; }}
    .grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:12px;
      margin:18px 0; }}
    .card {{ background:white; border:1px solid var(--line); border-radius:14px; padding:16px; }}
    .metrics .card {{ background:linear-gradient(180deg,#fff,var(--soft)); }}
    .card span {{ display:block; color:var(--muted); font-size:11px; text-transform:uppercase;
      letter-spacing:.06em; }}
    .card strong {{ display:block; font-size:25px; line-height:1.2; margin-top:7px;
      overflow-wrap:anywhere; }}
    .notice {{ border-radius:12px; padding:13px 15px; margin:14px 0; border:1px solid var(--line); }}
    .notice.info {{ background:var(--soft); color:var(--brand-dark); }}
    .notice.warn {{ background:var(--review-bg); color:#7c2d12; border-color:#fed7aa; }}
    .muted {{ color:var(--muted); font-size:13px; }}
    .pill {{ display:inline-flex; border-radius:999px; padding:4px 9px; font-size:11px;
      font-weight:800; letter-spacing:.035em; white-space:nowrap; }}
    .pill.critical {{ color:var(--critical); background:var(--critical-bg); }}
    .pill.review,.pill.warn {{ color:var(--review); background:var(--review-bg); }}
    .pill.good {{ color:var(--good); background:var(--good-bg); }}
    .pill.locked,.pill.neutral {{ color:var(--locked); background:var(--locked-bg); }}
    .table-wrap {{ overflow:auto; border:1px solid var(--line); border-radius:12px; }}
    table {{ width:100%; border-collapse:separate; border-spacing:0; min-width:880px; }}
    th,td {{ padding:11px 12px; border-bottom:1px solid var(--line); text-align:left;
      vertical-align:top; }}
    th {{ position:sticky; top:0; background:var(--soft); font-size:12px; text-transform:uppercase;
      letter-spacing:.04em; z-index:1; }}
    tbody tr:last-child td {{ border-bottom:0; }} tbody tr:hover {{ background:#f8fbfd; }}
    tr.critical {{ box-shadow:inset 4px 0 var(--critical); }}
    tr.review {{ box-shadow:inset 4px 0 #f59e0b; }}
    .table-filter {{ width:min(480px,100%); padding:10px 12px; border:1px solid var(--line);
      border-radius:9px; margin:0 0 12px; font:inherit; }}
    details {{ margin-top:12px; }} details>summary {{ cursor:pointer; font-weight:750; color:var(--brand); }}
    code {{ font-size:12px; overflow-wrap:anywhere; }}
    .iscn {{ font:700 18px ui-monospace,monospace; color:var(--brand); overflow-wrap:anywhere; }}
    ul {{ padding-left:22px; }} li+li {{ margin-top:7px; }}
    @media (max-width:720px) {{ main {{ padding-inline:12px; }} header {{ border-radius:16px; }}
      section {{ padding:16px; border-radius:14px; }} .section-heading {{ display:block; }} }}
  </style>
</head>
<body>
  <div class="banner">RESEARCH USE ONLY · NOT CLINICALLY VALIDATED · EXPERT REVIEW REQUIRED</div>
  <main>
    <header><span class="eyebrow" style="color:#bae6fd">ONTSeq single-sample report</span>
      <h1>{_cell(result.manifest.sample_id)}</h1>
      <p>Run {_cell(result.manifest.run_id)} · {_cell(result.manifest.assay.mode.value)} ·
      {_cell(result.manifest.assay.genome_build.value)} · pipeline {_cell(result.provenance.pipeline_version)}</p>
    </header>
    <nav aria-label="Report sections"><a href="#overview">Overview</a><a href="#findings">Key findings</a>
      <a href="#coverage">Target coverage</a><a href="#fusions">Fusion assessment</a>
      <a href="#sv-review">SV review</a><a href="#methods">Methods</a></nav>
    <section id="overview"><div class="section-heading"><div><span class="eyebrow">Run overview</span>
      <h2>What needs attention</h2></div>{_status_pill(result.release_status.value, css_class="locked")}</div>
      <div class="grid metrics"><div class="card"><span>QC verdict</span><strong>{_cell(result.qc.verdict.value)}</strong></div>
        <div class="card"><span>Fusion assessment</span><strong>{_cell(_module_status(result, AnalysisModule.FUSION))}</strong></div>
        <div class="card"><span>Key findings</span><strong>{len(finding_events)}</strong></div>
        <div class="card"><span>Technical SV review</span><strong>{len(review_events)}</strong></div>
        <div class="card"><span>All normalized SV</span><strong>{len(structural_events)}</strong></div></div>
      <div class="notice info"><strong>Release state:</strong> BENCHMARK_REQUIRED means that the
      pipeline produced reviewable evidence but no assay-specific analytical release gate has been
      satisfied. It does not mean irrelevant, negative, or absent.</div>{adaptive_warning}</section>
    <section id="findings"><div class="section-heading"><div><span class="eyebrow">Prioritized evidence</span>
        <h2>Key findings for expert review</h2></div><span class="pill neutral">{len(finding_events)} shown</span></div>
      <p class="muted">Knowledge matches and technically high/moderate fusion-evidence candidates.
      Gene-pair matching is order independent. These are candidates, not confirmed fusions.
      Pathology labels are source-database associations, not diagnoses for this sample.</p>
      <div class="table-wrap"><table><thead><tr><th>Review priority</th><th>Finding</th><th>Type</th>
      <th>Loci</th><th>Cytobands</th><th>Technical confidence</th><th>Evidence</th>
      <th>Breakpoint coverage</th><th>Associated pathologies</th><th>Fusion assessment</th><th>Release state</th></tr></thead>
      <tbody>{finding_rows}</tbody></table></div></section>
    {_coverage_section(target_coverage, selection_coverage)}
    <!-- ONTSEQ_CNV_SECTION -->
    <section id="fusions"><div class="section-heading"><div><span class="eyebrow">Breakpoint interpretation</span>
        <h2>Fusion and rearrangement assessment</h2></div><span class="pill neutral">{len(fusion_events)} candidates</span></div>
      <div class="notice info">BREAKPOINT_EVIDENCE means both breakpoints were gene/transcript
      annotated. KNOWLEDGE_MATCH_CANDIDATE additionally matched a locked hematology pattern.
      Neither state asserts a productive transcript or clinical reportability. Associated
      pathology names and DOIDs describe the source record, not the sample.</div>
      <input class="table-filter" data-table="fusion-events" placeholder="Filter gene, locus, band, status…">
      <div class="table-wrap"><table><thead><tr><th>Gene pair</th><th>Loci</th><th>Cytobands</th>
      <th>Caller support</th><th>Coverage</th><th>Observability</th><th>Transcript context</th>
      <th>Orientation/frame</th><th>Knowledge match</th><th>Associated pathologies</th><th>Assessment</th><th>Confidence</th>
      <th>Release state</th></tr></thead><tbody id="fusion-events">{fusion_rows}</tbody></table></div></section>
    <section id="sv-review"><div class="section-heading"><div><span class="eyebrow">Technical evidence</span>
        <h2>SV review queue</h2></div><span class="pill neutral">{len(review_events)} high/moderate</span></div>
      <div class="grid metrics"><div class="card"><span>Normalized caller records</span><strong>{raw_sv_count}</strong></div>
        <div class="card"><span>Consolidated SV</span><strong>{len(structural_events)}</strong></div>
        <div class="card"><span>High</span><strong>{confidence_counts["high"]}</strong></div>
        <div class="card"><span>Moderate</span><strong>{confidence_counts["moderate"]}</strong></div>
        <div class="card"><span>Low/background</span><strong>{confidence_counts["low"]}</strong></div></div>
      <p class="muted">Technical confidence measures caller evidence, not clinical relevance.
      Tumor-only data without a matched normal or population/PON filter can retain germline and
      recurrent technical background.</p>
      <input class="table-filter" data-table="priority-events" placeholder="Filter review queue…">
      <div class="table-wrap"><table><thead><tr><th>ID</th><th>Type</th><th>Gene pair</th>
      <th>Locus 1</th><th>Locus 2</th><th>Cytobands</th><th>Caller support</th><th>Coverage</th>
      <th>Observability</th><th>Context flags</th><th>Knowledge match</th><th>Associated pathologies</th><th>Fusion assessment</th>
      <th>Confidence</th><th>Release state</th></tr></thead>
      <tbody id="priority-events">{review_rows}</tbody></table></div>
      <details><summary>Technical appendix — {len(html_events)} of {len(structural_events)} SV shown in HTML</summary>
        <p class="muted">JSON and XLSX retain every normalized event. HTML includes every key/review
        event plus at most 200 background calls to keep the report usable.</p>
        <input class="table-filter" data-table="all-events" placeholder="Filter technical appendix…">
        <div class="table-wrap"><table><thead><tr><th>ID</th><th>Type</th><th>Length</th>
        <th>Locus 1</th><th>Locus 2</th><th>Gene pair</th><th>Confidence</th><th>Review priority</th>
        <th>Release state</th><th>Validation status</th><th>Observability</th><th>Context flags</th>
        <th>Knowledge match</th><th>Associated pathologies</th><th>Fusion assessment</th><th>Evidence</th></tr></thead>
        <tbody id="all-events">{event_rows}</tbody></table></div></details></section>
    <section><div class="section-heading"><div><span class="eyebrow">Pipeline state</span>
      <h2>Module status</h2></div></div><div class="table-wrap"><table><thead><tr><th>Module</th>
      <th>Status</th><th>Reason</th></tr></thead><tbody>{module_rows}</tbody></table></div></section>
    <section><div class="section-heading"><div><span class="eyebrow">Nomenclature</span>
      <h2>Proposed ISCN notation</h2></div></div><div class="iscn">{_cell(result.iscn.notation)}</div>
      <p>{_cell(result.iscn.standard_edition)} · {_cell(result.iscn.conformance_profile)} ·
      {_cell(result.iscn.review_status.value)}</p></section>
    <section><div class="section-heading"><div><span class="eyebrow">Read metrics</span>
      <h2>Quality control</h2></div></div><div class="grid metrics">{metric_cards}</div></section>
    <section><div class="section-heading"><div><span class="eyebrow">Limitations</span>
      <h2>Warnings and limitations</h2></div></div><ul>{warnings}</ul></section>
    <section id="methods"><div class="section-heading"><div><span class="eyebrow">Provenance</span>
      <h2>Methods and versions</h2></div></div>
      <h3>Reference resources</h3><div class="table-wrap"><table><thead><tr><th>Resource</th>
      <th>Release</th></tr></thead><tbody>{reference_rows}</tbody></table></div>
      <h3>Tools</h3><div class="table-wrap"><table><thead><tr><th>Tool</th><th>Version</th>
      <th>Parameters</th></tr></thead><tbody>{tool_rows}</tbody></table></div>
      <details><summary>Resource SHA256 provenance</summary><div class="table-wrap"><table><thead><tr>
      <th>Resource</th><th>SHA256</th></tr></thead><tbody>{checksum_rows}</tbody></table></div></details></section>
  </main>
  <script>
    for (const input of document.querySelectorAll("input[data-table]")) {{
      input.addEventListener("input", () => {{
        const query = input.value.toLowerCase();
        const body = document.getElementById(input.dataset.table);
        if (!body) return;
        for (const row of body.rows) row.hidden = !row.textContent.toLowerCase().includes(query);
      }});
    }}
  </script>
</body>
</html>
"""
    output_path.write_text(document, encoding="utf-8")
    return output_path

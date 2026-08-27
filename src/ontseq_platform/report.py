from __future__ import annotations

import html
import json
from pathlib import Path

from .models import GenomicEvent, PipelineResult
from .sv_evidence import sv_review_queue


def _cell(value: object) -> str:
    return html.escape("" if value is None else str(value))


def _event_evidence(event: GenomicEvent) -> str:
    return ", ".join(
        f"{item.caller} {item.caller_version} "
        f"(support={item.support_reads}, vaf={item.variant_allele_fraction})"
        for item in event.evidence
    )


def _event_loci(event: GenomicEvent) -> tuple[str, str]:
    primary = f"{event.primary.chromosome}:{event.primary.start}-{event.primary.end}"
    secondary = (
        ""
        if event.secondary is None
        else f"{event.secondary.chromosome}:{event.secondary.start}-{event.secondary.end}"
    )
    return primary, secondary


def _full_event_row(event: GenomicEvent) -> str:
    primary_locus, secondary_locus = _event_loci(event)
    return (
        "<tr>"
        f"<td>{_cell(event.event_id)}</td>"
        f"<td>{_cell(event.event_type.value)}</td>"
        f"<td>{_cell(event.length_bp)}</td>"
        f"<td>{_cell(primary_locus)}</td>"
        f"<td>{_cell(secondary_locus)}</td>"
        f"<td>{_cell(event.primary.cytoband_start)}</td>"
        f"<td>{_cell(', '.join(event.genes))}</td>"
        f"<td>{_cell(event.confidence)}</td>"
        f"<td>{_cell(event.reportable)}</td>"
        f"<td>{_cell(_event_evidence(event))}</td>"
        "</tr>"
    )


def _review_event_row(event: GenomicEvent) -> str:
    primary_locus, secondary_locus = _event_loci(event)
    return (
        "<tr>"
        f"<td>{_cell(event.event_id)}</td>"
        f"<td>{_cell(event.event_type.value)}</td>"
        f"<td>{_cell(primary_locus)}</td>"
        f"<td>{_cell(secondary_locus)}</td>"
        f"<td>{_cell(', '.join(event.genes))}</td>"
        f"<td><strong>{_cell(event.confidence)}</strong></td>"
        f"<td>{_cell(_event_evidence(event))}</td>"
        "</tr>"
    )


def render_html(result: PipelineResult, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    event_rows = [_full_event_row(event) for event in result.events]
    if not event_rows:
        event_rows.append(
            "<tr><td colspan='10'>No events were produced. Review module status; this is not "
            "a biological negative result.</td></tr>"
        )

    review_events = sv_review_queue(result.events, limit=50)
    review_rows = [_review_event_row(event) for event in review_events]
    if not review_rows:
        review_rows.append(
            "<tr><td colspan='7'>No high/moderate technical-priority SV candidate was "
            "identified. This is not a biological negative result.</td></tr>"
        )

    metric_cards = "".join(
        f"<div class='card'><span>{_cell(key.replace('_', ' ').title())}</span>"
        f"<strong>{_cell(value)}</strong></div>"
        for key, value in result.qc.metrics.items()
    )
    warnings = "".join(
        f"<li>{_cell(item)}</li>"
        for item in result.warnings + result.qc.warnings + result.iscn.warnings
    )
    module_rows = "".join(
        f"<tr><td>{_cell(module.module.value)}</td><td>{_cell(module.status.value)}</td>"
        f"<td>{_cell(module.reason)}</td></tr>"
        for module in result.modules
    )
    if not module_rows:
        module_rows = "<tr><td colspan='3'>No module outcomes were recorded.</td></tr>"
    tool_rows = "".join(
        f"<tr><td>{_cell(tool.name)}</td><td>{_cell(tool.version)}</td>"
        f"<td><code>{_cell(json.dumps(tool.parameters, sort_keys=True))}</code></td></tr>"
        for tool in result.provenance.tools
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ONTSeq report - {_cell(result.manifest.sample_id)}</title>
  <style>
    :root {{ color-scheme: light; --ink:#132238; --muted:#5b6778; --line:#dbe3ec;
      --brand:#075985; --soft:#eef7fb; --warn:#9a3412; --priority:#7c2d12; }}
    body {{ margin:0; font:15px/1.5 system-ui,sans-serif; color:var(--ink); background:#f6f8fb; }}
    main {{ max-width:1120px; margin:32px auto; padding:0 24px 48px; }}
    header {{ background:linear-gradient(135deg,#083344,#075985); color:white; padding:28px;
      border-radius:18px; box-shadow:0 12px 35px #0c4a6e33; }}
    .banner {{ background:#7f1d1d; color:white; padding:10px 16px; text-align:center;
      font-weight:800; letter-spacing:.04em; }}
    .grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr));
      gap:12px; margin:18px 0; }}
    .card,section {{ background:white; border:1px solid var(--line); border-radius:14px;
      padding:18px; }}
    .card span {{ display:block; color:var(--muted); font-size:12px; text-transform:uppercase; }}
    .card strong {{ display:block; font-size:24px; margin-top:6px; }}
    section {{ margin-top:18px; overflow-x:auto; }}
    table {{ width:100%; border-collapse:collapse; }}
    th,td {{ padding:10px; border-bottom:1px solid var(--line); text-align:left;
      vertical-align:top; }}
    th {{ background:var(--soft); }} code {{ font-size:12px; }}
    .iscn {{ font:700 20px ui-monospace,monospace; color:var(--brand); overflow-wrap:anywhere; }}
    .warn {{ color:var(--warn); }}
    .priority-note {{ color:var(--priority); font-weight:650; }}
  </style>
</head>
<body>
  <div class="banner">RESEARCH USE ONLY - NOT CLINICALLY VALIDATED</div>
  <main>
    <header>
      <h1>ONTSeq single-sample report</h1>
      <p>Sample {_cell(result.manifest.sample_id)} | Run {_cell(result.manifest.run_id)} |
      {_cell(result.manifest.assay.mode.value)} |
      {_cell(result.manifest.assay.genome_build.value)}</p>
    </header>
    <div class="grid">
      <div class="card"><span>QC verdict</span>
        <strong>{_cell(result.qc.verdict.value)}</strong></div>
      <div class="card"><span>Release status</span>
        <strong>{_cell(result.release_status.value)}</strong></div>
      <div class="card"><span>Review queue</span><strong>{len(review_events)}</strong></div>
      <div class="card"><span>All events</span><strong>{len(result.events)}</strong></div>
      <div class="card"><span>Pipeline</span>
        <strong>{_cell(result.provenance.pipeline_version)}</strong></div>
    </div>
    <section><h2>Module status</h2><table><thead><tr><th>Module</th><th>Status</th>
      <th>Reason</th></tr></thead><tbody>{module_rows}</tbody></table></section>
    <section><h2>Proposed ISCN notation</h2><div class="iscn">{_cell(result.iscn.notation)}</div>
      <p>{_cell(result.iscn.standard_edition)} | {_cell(result.iscn.conformance_profile)} |
      {_cell(result.iscn.review_status.value)}</p></section>
    <section><h2>Quality control</h2><div class="grid">{metric_cards}</div></section>
    <section><h2>SV review queue</h2>
      <p class="priority-note">Automated technical prioritization only — not clinical validation.
      High/moderate candidates are surfaced here for efficient review; every normalized event is
      retained in the complete table below and remains non-reportable pending assay validation.</p>
      <table><thead><tr><th>ID</th><th>Type</th><th>Locus 1</th><th>Locus 2</th>
      <th>Genes</th><th>Technical confidence</th><th>Evidence</th></tr></thead>
      <tbody>{"".join(review_rows)}</tbody></table></section>
    <section><h2>All normalized genomic events</h2><table><thead><tr><th>ID</th><th>Type</th>
      <th>Length (bp)</th><th>Locus 1</th><th>Locus 2</th><th>Band</th><th>Genes</th>
      <th>Confidence</th><th>Reportable</th><th>Evidence</th></tr></thead>
      <tbody>{"".join(event_rows)}</tbody></table></section>
    <section><h2>Warnings and limitations</h2><ul class="warn">{warnings}</ul></section>
    <section><h2>Methods and versions</h2><table><thead><tr><th>Tool</th><th>Version</th>
      <th>Parameters</th></tr></thead><tbody>{tool_rows}</tbody></table></section>
  </main>
</body>
</html>
"""
    output_path.write_text(document, encoding="utf-8")
    return output_path

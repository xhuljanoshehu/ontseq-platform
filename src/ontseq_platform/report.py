from __future__ import annotations

import html
import json
from pathlib import Path

from .models import PipelineResult
from .report_view import AnnotationView, EventView, ReportView, build_report_view


def _cell(value: object) -> str:
    return html.escape("" if value is None else str(value))


def _optional(value: object | None) -> str:
    return "not available" if value is None else _cell(value)


def _fraction(value: float | None) -> str:
    if value is None:
        return "not available"
    return f"{value * 100:.1f}%"


def _metric_name(name: str) -> str:
    return name.replace("_", " ").strip().title()


def _alerts(view: ReportView) -> str:
    if not view.alerts:
        return (
            "<p class='muted'>No FAILED, NO_CALL, QC WARN or QC FAIL alert was derived from "
            "this result contract. This is not a validation or biological-negative claim.</p>"
        )
    return "".join(
        (
            f"<div class='alert alert-{_cell(item.level)}'>"
            f"<strong>{_cell(item.title)}</strong><p>{_cell(item.detail)}</p></div>"
        )
        for item in view.alerts
    )


def _module_strip(view: ReportView) -> str:
    if not view.modules:
        return "<p class='muted'>No module outcomes were recorded.</p>"
    return "".join(
        (
            f"<div class='module-state {_cell(item.css_class)}'>"
            f"<span>{_cell(item.name)}</span><strong>{_cell(item.status.value)}</strong></div>"
        )
        for item in view.modules
    )


def _module_rows(view: ReportView) -> str:
    if not view.modules:
        return "<tr><td colspan='4'>No module outcomes were recorded.</td></tr>"
    rows: list[str] = []
    for item in view.modules:
        reason = _cell(item.reason) or "not recorded"
        rows.append(
            "<tr>"
            f"<td>{_cell(item.name)}</td>"
            f"<td><span class='state-label {_cell(item.css_class)}'>"
            f"{_cell(item.status.value)}</span></td>"
            f"<td>{reason}</td>"
            f"<td>{_cell(item.meaning)}</td>"
            "</tr>"
        )
    return "".join(rows)


def _qc_rows(view: ReportView) -> str:
    if not view.qc_metrics:
        return "<tr><td colspan='2'>No normalized QC metrics were recorded.</td></tr>"
    return "".join(
        f"<tr><td>{_cell(_metric_name(key))}</td><td>{_optional(value)}</td></tr>"
        for key, value in view.qc_metrics
    )


def _failed_gates(view: ReportView) -> str:
    if not view.qc_failed_gates:
        return (
            "<p class='muted'>No failed QC gates were recorded. Metric-specific adequacy must "
            "not be inferred unless the governing QC policy is validated.</p>"
        )
    items = "".join(f"<li>{_cell(item)}</li>" for item in view.qc_failed_gates)
    return f"<div class='gate-failure'><strong>Failed QC gates</strong><ul>{items}</ul></div>"


def _annotation_table(annotations: tuple[AnnotationView, ...]) -> str:
    if not annotations:
        return "<p class='muted'>No knowledge-resource annotations were attached.</p>"
    rows = []
    for item in annotations:
        caveats = "; ".join(item.caveats)
        rows.append(
            "<tr>"
            f"<td>{_cell(item.source_id)} {_cell(item.source_release)}</td>"
            f"<td>{_cell(item.record_id)}</td>"
            f"<td>{_cell(item.assertion)}</td>"
            f"<td>{_cell(item.assertion_vocabulary)}</td>"
            f"<td>{_cell(item.record_origin)}</td>"
            f"<td>{_cell(item.scope_alignment)}</td>"
            f"<td>{_cell(item.scope_note)}</td>"
            f"<td>{_cell(caveats)}</td>"
            "</tr>"
        )
    return (
        "<div class='table-wrap'><table><caption>Knowledge-resource evidence</caption>"
        "<thead><tr><th>Source</th><th>Record</th><th>Assertion</th><th>Vocabulary</th>"
        "<th>Origin</th><th>Scope</th><th>Scope note</th><th>Caveats</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )


def _evidence_table(event: EventView) -> str:
    if not event.evidence:
        return (
            "<p class='muted'>No normalized caller evidence was attached to this event. "
            "Do not infer missing evidence values.</p>"
        )
    rows = []
    for item in event.evidence:
        filters = ", ".join(item.filters) if item.filters else "none recorded"
        rows.append(
            "<tr>"
            f"<td>{_cell(item.caller)}</td>"
            f"<td>{_cell(item.caller_version)}</td>"
            f"<td>{_optional(item.support_reads)}</td>"
            f"<td>{_optional(item.local_coverage)}</td>"
            f"<td>{_fraction(item.variant_allele_fraction)}</td>"
            f"<td>{_optional(item.quality)}</td>"
            f"<td>{_optional(item.supporting_read_strands)}</td>"
            f"<td>{_optional(item.precise)}</td>"
            f"<td>{_cell(filters)}</td>"
            "</tr>"
        )
    return (
        "<div class='table-wrap'><table><caption>Normalized caller evidence</caption>"
        "<thead><tr><th>Caller</th><th>Version</th><th>Support reads</th>"
        "<th>Local coverage</th><th>VAF</th><th>Quality</th><th>Strands</th>"
        "<th>Precise</th><th>Filters</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )


def _event_boundary(event: EventView) -> str:
    text = (
        "Pipeline event label and confidence are analytical metadata, not a clinical "
        "classification. Review caller evidence, observability, validation status and source "
        "provenance separately."
    )
    if event.event_type == "fusion":
        text += (
            " A fusion event label alone does not establish an expressed, in-frame or functional "
            "fusion transcript."
        )
    return text


def _event_card(event: EventView) -> str:
    genes = ", ".join(event.genes) if event.genes else "not available"
    notes = "".join(f"<li>{_cell(item)}</li>" for item in event.notes)
    notes_block = f"<ul>{notes}</ul>" if notes else "<p class='muted'>No event notes recorded.</p>"
    return f"""
    <article class="event-card" id="event-{_cell(event.event_id)}">
      <div class="event-heading">
        <div>
          <span class="eyebrow">Normalized genomic event</span>
          <h3>{_cell(event.event_id)} · {_cell(event.event_type)}</h3>
        </div>
        <span class="reportability">reportable: {_cell(event.reportability_text)}</span>
      </div>
      <dl class="event-grid">
        <div><dt>Locus 1</dt><dd>{_cell(event.primary_locus)}</dd></div>
        <div><dt>Locus 2</dt><dd>{_optional(event.secondary_locus)}</dd></div>
        <div><dt>Cytoband</dt><dd>{_optional(event.cytobands)}</dd></div>
        <div><dt>Length</dt><dd>{_optional(event.length_bp)} bp</dd></div>
        <div><dt>Copy number</dt><dd>{_optional(event.copy_number)}</dd></div>
        <div><dt>Genes</dt><dd>{_cell(genes)}</dd></div>
        <div><dt>Pipeline confidence</dt><dd>{_cell(event.confidence)}</dd></div>
        <div><dt>Evidence records</dt><dd>{len(event.evidence)}</dd></div>
      </dl>
      <div class="boundary"><strong>Interpretation boundary.</strong>
        {_cell(_event_boundary(event))}</div>
      {_evidence_table(event)}
      <h4>Event notes</h4>{notes_block}
      <h4>Knowledge-resource annotations</h4>{_annotation_table(event.annotations)}
    </article>
    """


def _events(view: ReportView) -> str:
    if not view.events:
        return (
            "<div class='empty-state'><strong>No normalized events were produced.</strong>"
            "<p>Review module status and observability. This is not a biological negative "
            "result.</p></div>"
        )
    return "".join(_event_card(item) for item in view.events)


def _warnings(view: ReportView) -> str:
    if not view.warnings:
        return (
            "<p class='muted'>No warning strings were recorded in this result contract. "
            "This does not establish assay adequacy or absence of limitations.</p>"
        )
    return "".join(f"<li>{_cell(item)}</li>" for item in view.warnings)


def _tool_rows(result: PipelineResult) -> str:
    if not result.provenance.tools:
        return "<tr><td colspan='3'>No tool provenance was recorded.</td></tr>"
    return "".join(
        f"<tr><td>{_cell(tool.name)}</td><td>{_cell(tool.version)}</td>"
        f"<td><code>{_cell(json.dumps(tool.parameters, sort_keys=True))}</code></td></tr>"
        for tool in result.provenance.tools
    )


def _checksum_rows(view: ReportView) -> str:
    if not view.reference_checksums:
        return "<tr><td colspan='2'>No reference checksums were recorded.</td></tr>"
    return "".join(
        f"<tr><td>{_cell(name)}</td><td><code>{_cell(value)}</code></td></tr>"
        for name, value in view.reference_checksums
    )


def render_html(result: PipelineResult, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    view = build_report_view(result)
    target_design = view.target_bed_version or "not applicable / not recorded"
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ONTSeq report - {_cell(view.sample_id)}</title>
  <style>
    :root {{
      color-scheme: light;
      --ink:#172033; --muted:#5e687a; --line:#d9dee7; --panel:#ffffff;
      --canvas:#f3f5f8; --accent:#174a6e; --accent-soft:#eaf2f7;
      --critical:#8f1d1d; --critical-soft:#fff0f0; --warning:#8a4b08;
      --warning-soft:#fff7e8; --info:#36566f; --info-soft:#eef5f9;
      --ok:#245c45; --ok-soft:#edf7f1; --neutral:#586174;
      --neutral-soft:#f0f2f5;
    }}
    * {{ box-sizing:border-box; }}
    body {{
      margin:0; font:14px/1.5 Inter,Segoe UI,system-ui,sans-serif; color:var(--ink);
      background:var(--canvas);
    }}
    .ruo {{
      position:sticky; top:0; z-index:20; background:#731c1c; color:white;
      padding:9px 20px; text-align:center; font-weight:800; letter-spacing:.05em;
    }}
    .shell {{ max-width:1460px; margin:0 auto; padding:24px; }}
    .masthead {{
      background:var(--panel); border:1px solid var(--line); border-radius:14px;
      padding:24px;
    }}
    .masthead h1 {{ margin:0 0 6px; font-size:28px; }}
    .eyebrow {{
      color:var(--muted); font-size:11px; font-weight:800; letter-spacing:.08em;
      text-transform:uppercase;
    }}
    .identity {{
      display:grid; grid-template-columns:repeat(6,minmax(120px,1fr)); gap:10px;
      margin-top:20px;
    }}
    .identity div {{ border-top:2px solid var(--line); padding-top:8px; min-width:0; }}
    .identity span {{
      display:block; color:var(--muted); font-size:11px; text-transform:uppercase;
    }}
    .identity strong {{ display:block; margin-top:3px; overflow-wrap:anywhere; }}
    .layout {{
      display:grid; grid-template-columns:220px minmax(0,1fr); gap:18px;
      margin-top:18px;
    }}
    nav {{
      align-self:start; position:sticky; top:58px; background:var(--panel);
      border:1px solid var(--line); border-radius:12px; padding:10px;
    }}
    nav a {{
      display:block; padding:9px 10px; border-radius:8px; color:var(--ink);
      text-decoration:none;
    }}
    nav a:hover, nav a:focus-visible {{ background:var(--accent-soft); outline:none; }}
    main {{ min-width:0; }}
    section {{
      background:var(--panel); border:1px solid var(--line); border-radius:12px;
      padding:20px; margin-bottom:16px;
    }}
    h2 {{ margin:0 0 14px; font-size:20px; }}
    h3 {{ margin:2px 0 0; font-size:17px; }}
    h4 {{ margin:18px 0 8px; font-size:14px; }}
    p {{ margin:6px 0; }}
    .muted {{ color:var(--muted); }}
    .module-strip {{
      display:grid; grid-template-columns:repeat(auto-fit,minmax(125px,1fr));
      gap:8px; margin-top:12px;
    }}
    .module-state {{
      border:1px solid var(--line); border-left-width:5px; border-radius:9px;
      padding:10px;
    }}
    .module-state span {{
      display:block; color:var(--muted); font-size:11px; text-transform:uppercase;
    }}
    .module-state strong {{ display:block; margin-top:3px; }}
    .state-completed {{ border-left-color:var(--ok); background:var(--ok-soft); }}
    .state-no-call {{ border-left-color:var(--warning); background:var(--warning-soft); }}
    .state-failed {{ border-left-color:var(--critical); background:var(--critical-soft); }}
    .state-not-run {{ border-left-color:var(--neutral); background:var(--neutral-soft); }}
    .state-label {{
      display:inline-block; padding:3px 7px; border-radius:999px;
      border:1px solid currentColor; font-size:11px; font-weight:800;
    }}
    .alert {{
      border-left:5px solid; padding:12px 14px; margin:10px 0; border-radius:8px;
    }}
    .alert-critical {{ color:var(--critical); background:var(--critical-soft); }}
    .alert-warning {{ color:var(--warning); background:var(--warning-soft); }}
    .alert-info {{ color:var(--info); background:var(--info-soft); }}
    .table-wrap {{ overflow-x:auto; margin-top:10px; }}
    table {{ width:100%; border-collapse:collapse; min-width:620px; }}
    caption {{ text-align:left; font-weight:800; margin:0 0 8px; }}
    th,td {{
      padding:9px 10px; border-bottom:1px solid var(--line); text-align:left;
      vertical-align:top;
    }}
    th {{ background:#f7f8fa; font-size:12px; }}
    code {{
      font:12px/1.45 ui-monospace,SFMono-Regular,Consolas,monospace;
      overflow-wrap:anywhere;
    }}
    .event-card {{
      border:1px solid var(--line); border-radius:11px; padding:16px;
      margin:14px 0; background:#fcfcfd;
    }}
    .event-heading {{
      display:flex; justify-content:space-between; gap:14px; align-items:flex-start;
    }}
    .reportability {{
      max-width:360px; border:1px solid var(--line); border-radius:8px;
      padding:7px 9px; font-size:12px; font-weight:700; background:white;
    }}
    .event-grid {{
      display:grid; grid-template-columns:repeat(4,minmax(120px,1fr)); gap:10px;
      margin:14px 0;
    }}
    .event-grid div {{ border-top:1px solid var(--line); padding-top:7px; min-width:0; }}
    dt {{ color:var(--muted); font-size:11px; text-transform:uppercase; }}
    dd {{ margin:2px 0 0; overflow-wrap:anywhere; }}
    .boundary {{
      background:var(--info-soft); border-left:4px solid var(--info); padding:10px 12px;
      border-radius:7px;
    }}
    .gate-failure {{
      background:var(--critical-soft); color:var(--critical); border-radius:8px;
      padding:10px 12px; margin-top:10px;
    }}
    .iscn {{
      font:700 19px/1.5 ui-monospace,SFMono-Regular,Consolas,monospace;
      color:var(--accent); overflow-wrap:anywhere;
    }}
    .empty-state {{
      background:var(--info-soft); border:1px solid #bfd1dd; border-radius:8px;
      padding:14px;
    }}
    footer {{ color:var(--muted); font-size:12px; padding:4px 2px 24px; }}
    @media (max-width:1000px) {{
      .identity {{ grid-template-columns:repeat(3,minmax(120px,1fr)); }}
      .layout {{ grid-template-columns:1fr; }}
      nav {{ position:static; display:flex; overflow-x:auto; gap:4px; }}
      nav a {{ white-space:nowrap; }}
      .event-grid {{ grid-template-columns:repeat(2,minmax(120px,1fr)); }}
    }}
    @media (max-width:620px) {{
      .shell {{ padding:12px; }}
      .masthead {{ padding:17px; }}
      .identity {{ grid-template-columns:1fr 1fr; }}
      .event-heading {{ display:block; }}
      .reportability {{ margin-top:9px; max-width:none; }}
      .event-grid {{ grid-template-columns:1fr; }}
      section {{ padding:15px; }}
    }}
    @media print {{
      body {{ background:white; }}
      .ruo {{ position:static; }}
      nav {{ display:none; }}
      .layout {{ display:block; }}
      section,.masthead,.event-card {{ break-inside:avoid; box-shadow:none; }}
    }}
  </style>
</head>
<body>
  <div class="ruo">RESEARCH USE ONLY · NOT CLINICALLY VALIDATED</div>
  <div class="shell">
    <header class="masthead">
      <span class="eyebrow">ONTSeq evidence report</span>
      <h1>Single-sample analytical review</h1>
      <p class="muted">Evidence, execution state and provenance are shown separately from
        interpretation. Missing or non-executed analyses are never displayed as negatives.</p>
      <div class="identity">
        <div><span>Sample</span><strong>{_cell(view.sample_id)}</strong></div>
        <div><span>Run</span><strong>{_cell(view.run_id)}</strong></div>
        <div><span>Assay</span><strong>{_cell(view.assay_mode)}</strong></div>
        <div><span>Genome build</span><strong>{_cell(view.genome_build)}</strong></div>
        <div><span>Reference</span><strong>{_cell(view.reference_id)}</strong></div>
        <div><span>Release state</span><strong>{_cell(view.release_status)}</strong></div>
      </div>
    </header>
    <div class="layout">
      <nav aria-label="Report sections">
        <a href="#overview">Overview</a><a href="#modules">Module status</a>
        <a href="#qc">Quality</a><a href="#events">Events</a>
        <a href="#iscn">ISCN proposal</a><a href="#warnings">Warnings</a>
        <a href="#provenance">Provenance</a>
      </nav>
      <main>
        <section id="overview">
          <h2>1 · Review overview</h2>
          <div class="identity">
            <div><span>QC verdict</span><strong>{_cell(view.qc_verdict)}</strong></div>
            <div><span>Events</span><strong>{len(view.events)}</strong></div>
            <div><span>Analysis profile</span><strong>{_cell(view.analysis_profile)}</strong></div>
            <div><span>Analysis intent</span><strong>{_cell(view.analysis_intent)}</strong></div>
            <div>
              <span>Target design version</span><strong>{_cell(target_design)}</strong>
            </div>
            <div><span>Pipeline</span><strong>{_cell(view.pipeline_version)}</strong></div>
          </div>
          <h3>Execution-state strip</h3>
          <div class="module-strip">{_module_strip(view)}</div>
          <h3 style="margin-top:18px">Interpretation blockers and warnings</h3>
          {_alerts(view)}
        </section>
        <section id="modules">
          <h2>2 · Module execution status</h2>
          <p class="muted">Status is an execution statement, not a biological conclusion.</p>
          <div class="table-wrap"><table><caption>Module outcomes</caption>
            <thead><tr><th>Module</th><th>Status</th><th>Recorded reason</th>
              <th>Meaning</th></tr></thead><tbody>{_module_rows(view)}</tbody></table></div>
        </section>
        <section id="qc">
          <h2>3 · Quality and assay context</h2>
          <p><strong>QC verdict:</strong> {_cell(view.qc_verdict)}</p>
          <p class="muted">Normalized metrics are descriptive unless a validated QC policy
            explicitly defines an adequacy threshold.</p>
          <div class="table-wrap"><table><caption>Normalized QC metrics</caption>
            <thead><tr><th>Metric</th><th>Value</th></tr></thead>
            <tbody>{_qc_rows(view)}</tbody></table></div>{_failed_gates(view)}
        </section>
        <section id="events">
          <h2>4 · Genomic events and evidence</h2>
          <p class="muted">Each normalized event is displayed with caller evidence and an
            explicit interpretation boundary. Missing values remain “not available”.</p>
          {_events(view)}
        </section>
        <section id="iscn">
          <h2>5 · ISCN proposal</h2>
          <p class="muted">This is a proposal generated by an unvalidated conformance subset and
            requires expert review. It is not a released cytogenetic result.</p>
          <div class="iscn">{_cell(result.iscn.notation)}</div>
          <div class="identity">
            <div><span>Edition</span><strong>{_cell(result.iscn.standard_edition)}</strong></div>
            <div>
              <span>Conformance</span>
              <strong>{_cell(result.iscn.conformance_profile)}</strong>
            </div>
            <div>
              <span>Review status</span>
              <strong>{_cell(result.iscn.review_status.value)}</strong>
            </div>
          </div>
        </section>
        <section id="warnings">
          <h2>6 · Warnings and limitations</h2><ul>{_warnings(view)}</ul>
        </section>
        <section id="provenance">
          <h2>7 · Methods and provenance</h2>
          <div class="identity">
            <div>
              <span>Pipeline version</span><strong>{_cell(view.pipeline_version)}</strong>
            </div>
            <div><span>Git commit</span><strong>{_cell(view.git_commit)}</strong></div>
            <div><span>Generated</span><strong>{_cell(view.created_at)}</strong></div>
          </div>
          <div class="table-wrap"><table><caption>Tools and parameters</caption>
            <thead><tr><th>Tool</th><th>Version</th><th>Parameters</th></tr></thead>
            <tbody>{_tool_rows(result)}</tbody></table></div>
          <div class="table-wrap"><table><caption>Reference checksums</caption>
            <thead><tr><th>Resource</th><th>Checksum / lock value</th></tr></thead>
            <tbody>{_checksum_rows(view)}</tbody></table></div>
        </section>
        <footer>ONTSeq portable report · offline/self-contained presentation · RUO.</footer>
      </main>
    </div>
  </div>
</body>
</html>
"""
    output_path.write_text(document, encoding="utf-8")
    return output_path

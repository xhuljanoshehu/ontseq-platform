from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from .cnv.qdnaseq import QDNAseqCallReport
from .marlin_native_contracts import NativeMarlinReport
from .methylation import MODKIT_MODIFICATION_NAMES, MethylationReport
from .models import PipelineResult
from .report_cnv import cnv_html_section
from .report_formatting import cell as _cell
from .report_formatting import tool_parameters
from .report_interactive import attach_interactive_report
from .report_marlin import marlin_html_section, validate_marlin_identity
from .report_methylation import html_rows as methylation_html_rows
from .report_methylation import ordered_regions as ordered_methylation_regions
from .report_methylation import validate_methylation_identity
from .report_plots import (
    MethylationCell,
    ReadLengthBin,
    methylation_heatmap_svg,
    read_length_histogram_svg,
)
from .report_sections import coverage_section, iscn_details, resource_details, sv_details
from .report_style import REPORT_CSS
from .report_summary import status_cards, summary_cards
from .report_view import AnnotationView, EventView, ReportView, build_report_view
from .target_coverage import TargetCoverageReport, validate_report_coverage


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
        row_class = "critical" if item.status.value == "FAILED" else "neutral"
        rows.append(
            f"<tr class='{row_class}'><td>{_cell(item.name)}</td>"
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
        f"<td><code>{_cell(tool_parameters(tool.name, tool.parameters))}</code></td></tr>"
        for tool in result.provenance.tools
    )


def _checksum_rows(view: ReportView) -> str:
    if not view.reference_checksums:
        return "<tr><td colspan='2'>No reference checksums were recorded.</td></tr>"
    return "".join(
        f"<tr><td>{_cell(name)}</td><td><code>{_cell(value)}</code></td></tr>"
        for name, value in view.reference_checksums
    )


def _qc_histogram_figure(bins: Sequence[ReadLengthBin] | None, view: ReportView) -> str:
    """Read-length distribution from the normalized Cramino histogram sidecar."""
    if not bins:
        return ""
    n50 = next((value for key, value in view.qc_metrics if key == "n50_bp"), None)
    svg = read_length_histogram_svg(
        bins,
        title="Read length distribution",
        n50_bp=n50 if isinstance(n50, int) and not isinstance(n50, bool) else None,
    )
    if not svg:
        return ""
    return (
        "<figure class='plot' style='margin:14px 0'>"
        + svg
        + "<figcaption style='color:#5e687a;font-size:12px;margin-top:6px'>Read length "
        "distribution from the normalized Cramino histogram; the dashed line marks N50. "
        "Descriptive technical evidence, not an adequacy assessment.</figcaption></figure>"
    )


def _methylation_table(report: MethylationReport) -> str:
    """Every region row with its counts; below-floor rows say so instead of showing 0%."""
    rows = "".join(
        "<tr>" + "".join(f"<td>{_cell(value)}</td>" for value in row) + "</tr>"
        for row in methylation_html_rows(ordered_methylation_regions(report))
    )
    policy = report.policy
    tagged = report.reads_with_modified_base_tags
    return (
        f"<p class='muted'>Policy {_cell(policy.profile_id)} ({_cell(policy.status)}) · "
        f"{_cell(report.tool.name)} {_cell(report.tool.version)} · pinned call threshold "
        f"{_cell(policy.filter_threshold)} · coverage floor "
        f"{_cell(policy.minimum_valid_coverage)} valid calls per site · region source "
        f"{_cell(report.region_source.value)} · reads with MM tags "
        f"{_cell(tagged if tagged is not None else 'not established')}.</p>"
        f"<details><summary>Region table ({len(report.regions)} row(s))</summary>"
        "<div class='table-wrap'><table><caption>Modified-base fractions by region</caption>"
        "<thead><tr><th>Region</th><th>Locus</th><th>Modification</th>"
        "<th>Sites at floor / sites</th><th>Valid calls</th><th>Modified calls</th>"
        "<th>Mean fraction (call-weighted)</th><th>Median site fraction</th>"
        "<th>Mean valid coverage</th><th>Failed / no-call calls</th></tr></thead>"
        f"<tbody>{rows}</tbody></table></div></details>"
    )


def _methylation_section(report: MethylationReport | None) -> str:
    """Region × modification heatmap and region table from the normalized report."""
    if report is None:
        return ""
    cells = [
        MethylationCell(
            region_id=region.region_id,
            chromosome=region.chromosome,
            start=region.start,
            end=region.end,
            modification_label=MODKIT_MODIFICATION_NAMES[region.modification_code],
            fraction=region.mean_modified_fraction,
            valid_call_count=region.valid_call_count,
            sites_at_minimum_coverage=region.sites_at_minimum_coverage,
            sites_total=region.sites_total,
        )
        for region in report.regions
    ]
    svg = methylation_heatmap_svg(cells, title="Modified-base fractions by region")
    if not svg:
        return ""
    fail_calls = report.summary_metrics.get("fail_call_count", 0)
    nocall_calls = report.summary_metrics.get("nocall_call_count", 0)
    return (
        "<section id='methylation'>"
        "<h2>Modified-base fractions</h2>"
        "<p class='muted'>Call-weighted modified fractions from the normalized "
        "modkit pileup. Grey hatched cells are below the coverage floor "
        "(not measurable), never a measured zero. Descriptive technical "
        "evidence, not a classifier or reportability assessment.</p>"
        f"<p class='muted'>Status {_cell(report.status.value)}; "
        f"failed-threshold calls {fail_calls}; no-call counts {nocall_calls}.</p>"
        "<figure class='plot' style='margin:14px 0'>"
        + svg
        + "<figcaption style='color:#5e687a;font-size:12px;margin-top:6px'>"
        "Rows are modification codes; columns are regions in genome order. "
        "The colour scale is the already-normalized mean modified fraction. "
        "A measured 0% is a pale filled cell; an unmeasurable region is hatched."
        "</figcaption></figure>" + _methylation_table(report) + "</section>"
    )


def render_html(
    result: PipelineResult,
    output_path: Path,
    *,
    target_coverage: TargetCoverageReport | None = None,
    selection_coverage: TargetCoverageReport | None = None,
    qc_histogram: Sequence[ReadLengthBin] | None = None,
    methylation_report: MethylationReport | None = None,
    marlin_report: NativeMarlinReport | None = None,
    cnv_report: QDNAseqCallReport | None = None,
    cnv_evidence_root: Path | None = None,
) -> Path:
    validate_report_coverage(result, target_coverage, selection_coverage)
    validate_marlin_identity(result, marlin_report)
    validate_methylation_identity(result, methylation_report)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    view = build_report_view(result)
    target_design = view.target_bed_version or "not applicable / not recorded"
    metric_cards = "".join(
        f"<div><span>{_cell(_metric_name(key))}</span><strong>{_optional(value)}</strong></div>"
        for key, value in view.qc_metrics
    )
    qc_histogram_figure = _qc_histogram_figure(qc_histogram, view)
    methylation_section = _methylation_section(methylation_report)
    methylation_nav = '<a href="#methylation">Methylation</a>' if methylation_section else ""
    document = f"""<!doctype html>
<html lang="de">
<head>
  <meta charset="utf-8">
  <meta name="ontseq-report-layout" content="befund-v8-local-1">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ONTSeq report - {_cell(view.sample_id)}</title>
  <style>
{REPORT_CSS}
  </style>
</head>
<body>
  <div class="ruo">RESEARCH USE ONLY · NOT CLINICALLY VALIDATED</div>
  <header class="masthead">
    <div class="masthead-inner">
      <div>
      <span class="eyebrow">ONTSeq · Analytischer Befund</span>
      <h1>Auswertung Oxford-Nanopore-Sequenzierung</h1>
      <p class="muted">Maschinell erzeugte Auswertung einer Einzelprobe.
        Ohne fachliche Vidierung nicht freigegeben.</p>
      </div>
      <div class="identity">
        <div><span>Probe</span><strong>{_cell(view.sample_id)}</strong></div>
        <div><span>Lauf</span><strong>{_cell(view.run_id)}</strong></div>
        <div><span>Assay</span><strong>{_cell(view.assay_mode)}</strong></div>
        <div><span>Genom</span><strong>{_cell(view.genome_build)}</strong></div>
        <div><span>Referenz</span><strong>{_cell(view.reference_id)}</strong></div>
        <div><span>Freigabe</span>
          <strong class="release">{_cell(view.release_status)}</strong></div>
      </div>
    </div>
  </header>
  <div class="shell">
    <div class="layout">
      <nav aria-label="Berichtsabschnitte">
        <a href="#overview">Kernbefunde</a><a href="#modules">Ausführungszustand</a>
        <!-- ONTSEQ_CNV_NAV --><a href="#sv-review">Strukturvarianten</a>
        <a href="#events">Ereignisse</a><a href="#qc">Qualitätskontrolle</a>
        <a href="#coverage">Abdeckung</a>{methylation_nav}<a href="#marlin">MARLIN</a>
        <a href="#iscn">ISCN</a><a href="#warnings">Grenzen</a>
        <a href="#provenance">Technische Nachweise</a>
      </nav>
      <main>
        <section id="overview">
          <h2>Kernbefunde · erfasste Ereignisse</h2>
          <p class="muted">{len(view.events)} normalisierte Ereignisse im Ergebnisvertrag.
            Die Karten zeigen analytische Daten; ihre Bewertung erfordert die zugehörigen
            Nachweise und den Ausführungszustand.</p>
          {summary_cards(view)}
          <p class="summary-context">Profil: {_cell(view.analysis_profile)} ·
            Analyseabsicht: {_cell(view.analysis_intent)} ·
            Target design version: {_cell(target_design)} ·
            Pipeline {_cell(view.pipeline_version)}</p>
          {_alerts(view)}
        </section>
        <div class="status-grid" aria-label="Statusübersicht">{status_cards(view)}</div>
        <section id="modules">
          <h2>Ausführungszustand</h2>
          <p class="muted">Status is an execution statement, not a biological conclusion.</p>
          <div class="table-wrap"><table><caption>Module outcomes</caption>
            <thead><tr><th>Module</th><th>Status</th><th>Recorded reason</th>
              <th>Meaning</th></tr></thead><tbody>{_module_rows(view)}</tbody></table></div>
        </section>
        <!-- ONTSEQ_CNV_SECTION -->
        {sv_details(result)}
        <section id="events">
          <h2>Genomische Ereignisse und Nachweise</h2>
          <p class="muted">Each normalized event is displayed with caller evidence and an
            explicit interpretation boundary. Missing values remain “not available”.</p>
          {_events(view)}
        </section>
        <section id="qc">
          <h2>Qualitätskontrolle und Assay-Kontext</h2>
          <p><strong>QC verdict:</strong> {_cell(view.qc_verdict)}</p>
          <p class="muted">Normalized metrics are descriptive unless a validated QC policy
            explicitly defines an adequacy threshold.</p>
          <div class="identity">{metric_cards}</div>
          <div class="table-wrap"><table><caption>Normalized QC metrics</caption>
            <thead><tr><th>Metric</th><th>Value</th></tr></thead>
            <tbody>{_qc_rows(view)}</tbody></table></div>{qc_histogram_figure}{_failed_gates(view)}
        </section>
        {coverage_section(target_coverage, selection_coverage)}
        {methylation_section}
        {marlin_html_section(marlin_report)}
        <section id="iscn">
          <h2>ISCN · Vorschlag zur fachlichen Prüfung</h2>
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
          {iscn_details(result)}
        </section>
        <section id="warnings">
          <h2>Grenzen und Warnungen</h2><ul>{_warnings(view)}</ul>
        </section>
        <section id="provenance">
          <h2>Technische Nachweise und Methoden</h2>
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
          {resource_details(result)}
          <div class="table-wrap"><table><caption>Reference checksums</caption>
            <thead><tr><th>Resource</th><th>Checksum / lock value</th></tr></thead>
            <tbody>{_checksum_rows(view)}</tbody></table></div>
        </section>
        <footer>ONTSeq portable report · offline/self-contained presentation · RUO.
          This self-contained HTML has no CDN or remote runtime dependency.
          Darstellung: befund-v8-local-1. Analyseversion siehe technische Nachweise.</footer>
      </main>
    </div>
  </div>
</body>
</html>
"""
    if cnv_report is not None and cnv_evidence_root is not None:
        document = document.replace(
            "<!-- ONTSEQ_CNV_SECTION -->", cnv_html_section(cnv_evidence_root, cnv_report), 1
        )
        document = document.replace("<!-- ONTSEQ_CNV_NAV -->", '<a href="#cnv">Kopienzahl</a>', 1)
    document = attach_interactive_report(
        document,
        result,
        cnv=cnv_report,
        evidence_root=cnv_evidence_root,
        marlin_report=marlin_report,
    )
    output_path.write_text(document, encoding="utf-8")
    return output_path

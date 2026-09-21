"""Compact report cards derived solely from normalized result data."""

from __future__ import annotations

from .report_formatting import cell
from .report_view import ReportView


def summary_cards(view: ReportView) -> str:
    if not view.events:
        return (
            "<div class='empty-state'>No normalized events were produced. "
            "Dies ist kein biologisch negativer Befund. "
            "Ausführungszustand und Nachweisgrenzen prüfen.</div>"
        )
    cards = []
    for event in view.events[:6]:
        value = (
            f"{event.copy_number:.2f} Kopien"
            if event.copy_number is not None
            else cell(event.event_type)
        )
        cards.append(
            "<article class='summary-card'>"
            f"<strong>{cell(event.primary_locus.split(':', 1)[0])}</strong>"
            f"<div class='summary-value'>{value}</div>"
            f"<p>{cell(event.event_id)} · {cell(event.event_type)}</p>"
            f"<p>{cell(event.primary_locus)}</p>"
            f"<p>Berichtsfähigkeit (Pipeline): {cell(event.reportability_text)}</p>"
            f"<a href='#event-{cell(event.event_id)}'>Nachweise ansehen →</a>"
            "</article>"
        )
    remainder = len(view.events) - len(cards)
    more = (
        f"<p class='muted'>Weitere {remainder} Ereignisse im "
        "<a href='#events'>vollständigen Ereignisabschnitt</a>.</p>"
        if remainder
        else ""
    )
    return "<div class='summary-grid'>" + "".join(cards) + "</div>" + more


def status_cards(view: ReportView) -> str:
    entries = [("Qualitätskontrolle", view.qc_verdict, "Bewertung gemäß Laufvertrag")]
    modules = {item.name: item for item in view.modules}
    for name, label in (("cnv", "Kopienzahl"), ("sv", "Strukturvarianten"), ("iscn", "ISCN")):
        module = modules.get(name)
        entries.append(
            (label, module.status.value, "Ausführungsstatus · Details unten")
            if module
            else (label, "NOT_RECORDED", "Kein Modulstatus aufgezeichnet")
        )
    tones = {
        "PASS": "ok",
        "COMPLETED": "ok",
        "FAIL": "critical",
        "FAILED": "critical",
        "WARN": "warning",
        "NO_CALL": "warning",
    }
    return "".join(
        f"<div class='status-card tone-{tones.get(status, 'neutral')}'>"
        f"<span>{cell(label)}</span><strong>{cell(status)}</strong>"
        f"<p>{cell(note)}</p></div>"
        for label, status, note in entries
    )

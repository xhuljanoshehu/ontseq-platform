from __future__ import annotations

from pathlib import Path

from ontseq_platform.demo import build_demo_result
from ontseq_platform.models import ModuleRunStatus
from ontseq_platform.report import render_html
from ontseq_platform.report_summary import status_cards, summary_cards
from ontseq_platform.report_view import build_report_view


def test_presentation_preserves_normalized_result(tmp_path: Path) -> None:
    result = build_demo_result()
    before = result.model_dump_json()
    document = render_html(result, tmp_path / "report.html").read_text()
    assert result.model_dump_json() == before
    assert document.index('id="overview"') < document.index('id="modules"')
    assert document.index('id="modules"') < document.index("<!-- ONTSEQ_CNV_SECTION -->")
    assert document.index("<!-- ONTSEQ_CNV_SECTION -->") < document.index('id="qc"')
    for event in result.events:
        assert f"href='#event-{event.event_id}'" in document
        assert f'id="event-{event.event_id}"' in document
    assert "RESEARCH USE ONLY" in document
    assert "not a biological negative result" in document


def test_summary_keeps_zero_copy_number_and_escapes_input() -> None:
    result = build_demo_result()
    result.events[0].copy_number = 0
    result.events[0].event_id = '<img src=x onerror="alert(1)">'
    text = summary_cards(build_report_view(result))
    assert "0.00 Kopien" in text
    assert "<img" not in text
    assert "&lt;img" in text


def test_missing_module_is_not_inferred_completed() -> None:
    result = build_demo_result()
    result.modules = []
    text = status_cards(build_report_view(result))
    assert text.count("NOT_RECORDED") == 3
    assert "COMPLETED" not in text


def test_failed_and_no_call_statuses_are_preserved() -> None:
    result = build_demo_result()
    for module in result.modules:
        if module.module.value == "cnv":
            module.status = ModuleRunStatus.FAILED
        if module.module.value == "sv":
            module.status = ModuleRunStatus.NO_CALL
    text = status_cards(build_report_view(result))
    assert "tone-critical" in text and ">FAILED<" in text
    assert "tone-warning" in text and ">NO_CALL<" in text


def test_long_event_list_has_explicit_continuation(tmp_path: Path) -> None:
    result = build_demo_result()
    result.events = [
        result.events[0].model_copy(update={"event_id": f"synthetic-{index}"}) for index in range(8)
    ]
    document = render_html(result, tmp_path / "report.html").read_text()
    assert document.count("class='summary-card'") == 6
    assert "Weitere 2 Ereignisse" in document
    assert 'id="event-synthetic-7"' in document

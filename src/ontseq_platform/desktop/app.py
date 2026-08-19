from __future__ import annotations

import sys
import threading
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import QObject, QThread, QUrl, Signal, Slot
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .. import __version__
from ..models import AssayMode, GenomeBuild
from .backend import (
    DesktopAnalysisRequest,
    DesktopBackend,
    DesktopRunResult,
    DesktopStage,
    DesktopStageStatus,
    ProgressEvent,
    locate_bam_index,
    sanitize_sample_id,
)
from .config import (
    DesktopConfig,
    DesktopReferenceProfile,
    load_desktop_config,
    save_desktop_config,
)


STAGE_LABELS = {
    DesktopStage.INPUT: "Inputprüfung",
    DesktopStage.QC: "Quality Control",
    DesktopStage.CNV: "Copy Number Variants",
    DesktopStage.SV: "Structural Variants",
    DesktopStage.FUSION: "Fusionsevidenz",
    DesktopStage.ISCN: "ISCN-Vorschlag",
    DesktopStage.REPORT: "Report",
}

STATUS_SYMBOLS = {
    DesktopStageStatus.PENDING: "○",
    DesktopStageStatus.RUNNING: "●",
    DesktopStageStatus.PASS: "✓",
    DesktopStageStatus.WARN: "!",
    DesktopStageStatus.FAIL: "×",
    DesktopStageStatus.NOT_RUN: "—",
    DesktopStageStatus.NO_CALL: "∅",
}


class AnalysisWorker(QObject):
    progress = Signal(object)
    log = Signal(str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, backend: DesktopBackend, request: DesktopAnalysisRequest) -> None:
        super().__init__()
        self.backend = backend
        self.request = request
        self.cancel_event = threading.Event()

    @Slot()
    def run(self) -> None:
        try:
            result = self.backend.run(
                self.request,
                on_progress=self.progress.emit,
                on_log=self.log.emit,
                cancel_event=self.cancel_event,
            )
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(result)

    def cancel(self) -> None:
        self.cancel_event.set()


class StageRow(QFrame):
    def __init__(self, stage: DesktopStage) -> None:
        super().__init__()
        self.stage = stage
        layout = QHBoxLayout(self)
        layout.setContentsMargins(12, 8, 12, 8)
        self.symbol = QLabel(STATUS_SYMBOLS[DesktopStageStatus.PENDING])
        self.symbol.setFixedWidth(22)
        self.title = QLabel(STAGE_LABELS[stage])
        self.title.setObjectName("stageTitle")
        self.message = QLabel("Ausstehend")
        self.message.setObjectName("muted")
        self.message.setWordWrap(True)
        self.status = QLabel(DesktopStageStatus.PENDING.value)
        self.status.setFixedWidth(82)
        layout.addWidget(self.symbol)
        layout.addWidget(self.title, 2)
        layout.addWidget(self.message, 5)
        layout.addWidget(self.status)
        self.set_status(DesktopStageStatus.PENDING, "Ausstehend")

    def set_status(self, status: DesktopStageStatus, message: str) -> None:
        self.symbol.setText(STATUS_SYMBOLS[status])
        self.status.setText(status.value)
        self.message.setText(message)
        self.setProperty("runStatus", status.value)
        self.style().unpolish(self)
        self.style().polish(self)


class SettingsDialog(QDialog):
    def __init__(self, config: DesktopConfig, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("ONTSeq Desktop – Einstellungen")
        self.resize(760, 590)
        self._original = config
        root = QVBoxLayout(self)
        tabs = QTabWidget()
        tabs.addTab(self._build_backend_tab(config), "Backend")
        tabs.addTab(self._build_reference_tab(config), "Referenzen")
        root.addWidget(tabs)

        note = QLabel(
            "Diese Einstellungen enthalten nur lokale Tool-/Referenzpfade. Patientendaten "
            "gehören nicht in die Konfiguration."
        )
        note.setWordWrap(True)
        note.setObjectName("muted")
        root.addWidget(note)

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        cancel = QPushButton("Abbrechen")
        cancel.clicked.connect(self.reject)
        save = QPushButton("Speichern")
        save.setObjectName("primaryButton")
        save.clicked.connect(self._save)
        buttons.addWidget(cancel)
        buttons.addWidget(save)
        root.addLayout(buttons)
        self.saved_config: DesktopConfig | None = None

    def _build_backend_tab(self, config: DesktopConfig) -> QWidget:
        widget = QWidget()
        form = QFormLayout(widget)
        self.backend_mode = QComboBox()
        self.backend_mode.addItem("WSL2 (Windows – empfohlen)", "wsl")
        self.backend_mode.addItem("Lokal (Entwicklung/Linux)", "local")
        self.backend_mode.setCurrentIndex(0 if config.backend_mode == "wsl" else 1)
        self.wsl_distribution = QLineEdit(config.wsl_distribution or "")
        self.wsl_distribution.setPlaceholderText("leer = Standard-Distribution")
        self.project_root = QLineEdit(config.wsl_project_root)
        self.project_root.setPlaceholderText("~/ontseq-platform")
        self.output_root = QLineEdit(config.output_root)
        output_row = QWidget()
        output_layout = QHBoxLayout(output_row)
        output_layout.setContentsMargins(0, 0, 0, 0)
        output_layout.addWidget(self.output_root)
        browse = QPushButton("Ordner…")
        browse.clicked.connect(self._browse_output_root)
        output_layout.addWidget(browse)
        self.qc_policy = QLineEdit(config.qc_policy_path)
        self.sniffles_policy = QLineEdit(config.sniffles_policy_path)

        form.addRow("Ausführung:", self.backend_mode)
        form.addRow("WSL-Distribution:", self.wsl_distribution)
        form.addRow("ONTSeq-Projekt in WSL:", self.project_root)
        form.addRow("Ergebnisordner:", output_row)
        form.addRow("QC-Policy:", self.qc_policy)
        form.addRow("Sniffles2-Policy:", self.sniffles_policy)
        return widget

    def _build_reference_tab(self, config: DesktopConfig) -> QWidget:
        widget = QWidget()
        layout = QVBoxLayout(widget)
        existing = {profile.genome_build: profile for profile in config.reference_profiles}
        self.reference_fields: dict[GenomeBuild, dict[str, QLineEdit]] = {}
        for build in (GenomeBuild.GRCH37, GenomeBuild.GRCH38):
            profile = existing.get(build)
            group = QGroupBox(build.value)
            form = QFormLayout(group)
            reference_id = QLineEdit(profile.reference_id if profile else "")
            lock_path = QLineEdit(profile.reference_lock_path if profile else "")
            lock_row = self._path_row(lock_path, "JSON auswählen…")
            bed_path = QLineEdit(
                profile.adaptive_sampling_target_bed_path
                if profile and profile.adaptive_sampling_target_bed_path
                else ""
            )
            bed_row = self._path_row(bed_path, "BED auswählen…")
            bed_version = QLineEdit(
                profile.adaptive_sampling_target_bed_version
                if profile and profile.adaptive_sampling_target_bed_version
                else ""
            )
            form.addRow("Reference ID:", reference_id)
            form.addRow("Reference lock:", lock_row)
            form.addRow("Adaptive-Sampling BED:", bed_row)
            form.addRow("BED-Version:", bed_version)
            layout.addWidget(group)
            self.reference_fields[build] = {
                "reference_id": reference_id,
                "lock_path": lock_path,
                "bed_path": bed_path,
                "bed_version": bed_version,
            }
        layout.addStretch(1)
        return widget

    def _path_row(self, field: QLineEdit, button_text: str) -> QWidget:
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(field)
        button = QPushButton(button_text)
        button.clicked.connect(lambda: self._browse_file(field))
        layout.addWidget(button)
        return row

    def _browse_file(self, target: QLineEdit) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Datei auswählen")
        if path:
            target.setText(path)

    def _browse_output_root(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Ergebnisordner auswählen")
        if path:
            self.output_root.setText(path)

    @Slot()
    def _save(self) -> None:
        profiles: list[DesktopReferenceProfile] = []
        try:
            for build, fields in self.reference_fields.items():
                reference_id = fields["reference_id"].text().strip()
                lock_path = fields["lock_path"].text().strip()
                bed_path = fields["bed_path"].text().strip()
                bed_version = fields["bed_version"].text().strip()
                if not reference_id and not lock_path and not bed_path and not bed_version:
                    continue
                profiles.append(
                    DesktopReferenceProfile(
                        genome_build=build,
                        reference_id=reference_id,
                        reference_lock_path=lock_path,
                        adaptive_sampling_target_bed_path=bed_path or None,
                        adaptive_sampling_target_bed_version=bed_version or None,
                    )
                )
            self.saved_config = DesktopConfig(
                backend_mode=self.backend_mode.currentData(),
                wsl_distribution=self.wsl_distribution.text().strip() or None,
                wsl_project_root=self.project_root.text().strip(),
                output_root=self.output_root.text().strip(),
                qc_policy_path=self.qc_policy.text().strip(),
                sniffles_policy_path=self.sniffles_policy.text().strip(),
                reference_profiles=profiles,
            )
        except Exception as exc:
            QMessageBox.critical(self, "Ungültige Einstellungen", str(exc))
            return
        save_desktop_config(self.saved_config)
        self.accept()


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"ONTSeq Desktop {__version__}")
        self.resize(1120, 790)
        try:
            self.config = load_desktop_config()
        except Exception:
            self.config = DesktopConfig()
        self.thread: QThread | None = None
        self.worker: AnalysisWorker | None = None
        self.last_result: DesktopRunResult | None = None
        self.stage_rows: dict[DesktopStage, StageRow] = {}
        self._build_ui()
        self._apply_styles()

    def _build_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(24, 20, 24, 20)
        root.setSpacing(14)

        header = QHBoxLayout()
        title_block = QVBoxLayout()
        title = QLabel("ONTSeq Desktop")
        title.setObjectName("title")
        subtitle = QLabel("Oxford Nanopore – strukturierte Einzelprobenanalyse")
        subtitle.setObjectName("muted")
        title_block.addWidget(title)
        title_block.addWidget(subtitle)
        header.addLayout(title_block)
        header.addStretch(1)
        system_check = QPushButton("Systemcheck")
        system_check.clicked.connect(self._system_check)
        settings = QPushButton("Einstellungen")
        settings.clicked.connect(self._open_settings)
        header.addWidget(system_check)
        header.addWidget(settings)
        root.addLayout(header)

        banner = QLabel(
            "RESEARCH USE ONLY · nicht klinisch validiert · Ergebnisse erfordern fachliche Prüfung"
        )
        banner.setObjectName("researchBanner")
        banner.setWordWrap(True)
        root.addWidget(banner)

        main_grid = QGridLayout()
        main_grid.setHorizontalSpacing(16)
        root.addLayout(main_grid, 1)
        main_grid.addWidget(self._build_input_panel(), 0, 0)
        main_grid.addWidget(self._build_progress_panel(), 0, 1)
        main_grid.setColumnStretch(0, 4)
        main_grid.setColumnStretch(1, 6)

        self.details_button = QPushButton("Technische Details anzeigen")
        self.details_button.clicked.connect(self._toggle_log)
        root.addWidget(self.details_button)
        self.log_view = QTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMinimumHeight(145)
        self.log_view.setVisible(False)
        root.addWidget(self.log_view)

    def _build_input_panel(self) -> QWidget:
        panel = QFrame()
        panel.setObjectName("card")
        layout = QVBoxLayout(panel)
        heading = QLabel("Analyseauftrag")
        heading.setObjectName("sectionTitle")
        layout.addWidget(heading)

        form = QFormLayout()
        self.sample_id = QLineEdit()
        self.sample_id.setPlaceholderText("pseudonymisierte Sample-ID")
        form.addRow("Sample-ID:", self.sample_id)

        bam_row = QWidget()
        bam_layout = QHBoxLayout(bam_row)
        bam_layout.setContentsMargins(0, 0, 0, 0)
        self.bam_path = QLineEdit()
        self.bam_path.setReadOnly(True)
        self.bam_path.setPlaceholderText("Aligned BAM auswählen")
        bam_layout.addWidget(self.bam_path, 1)
        browse = QPushButton("BAM…")
        browse.clicked.connect(self._browse_bam)
        bam_layout.addWidget(browse)
        form.addRow("BAM:", bam_row)

        self.index_hint = QLabel("BAM-Index wird automatisch erkannt")
        self.index_hint.setObjectName("muted")
        form.addRow("", self.index_hint)

        self.genome_build = QComboBox()
        self.genome_build.addItem("GRCh38", GenomeBuild.GRCH38)
        self.genome_build.addItem("GRCh37", GenomeBuild.GRCH37)
        form.addRow("Referenz:", self.genome_build)

        self.assay_mode = QComboBox()
        self.assay_mode.addItem("Adaptive Sampling", AssayMode.ADAPTIVE_SAMPLING)
        self.assay_mode.addItem("Low-coverage WGS", AssayMode.LOW_COVERAGE_WGS)
        form.addRow("Analyseprofil:", self.assay_mode)
        layout.addLayout(form)

        privacy = QLabel(
            "Nur pseudonymisierte Sample-IDs verwenden. BAM/VCF/Patientendaten werden nicht "
            "an GitHub oder externe Dienste übertragen."
        )
        privacy.setWordWrap(True)
        privacy.setObjectName("muted")
        layout.addWidget(privacy)
        layout.addStretch(1)

        self.start_button = QPushButton("ANALYSE STARTEN")
        self.start_button.setObjectName("primaryButton")
        self.start_button.setMinimumHeight(44)
        self.start_button.clicked.connect(self._start_analysis)
        layout.addWidget(self.start_button)
        self.cancel_button = QPushButton("Analyse abbrechen")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self._cancel_analysis)
        layout.addWidget(self.cancel_button)
        return panel

    def _build_progress_panel(self) -> QWidget:
        panel = QFrame()
        panel.setObjectName("card")
        layout = QVBoxLayout(panel)
        heading = QLabel("Analysefortschritt")
        heading.setObjectName("sectionTitle")
        layout.addWidget(heading)
        for stage in DesktopStage:
            row = StageRow(stage)
            self.stage_rows[stage] = row
            layout.addWidget(row)
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        layout.addWidget(self.progress)
        self.summary = QLabel("Bereit für einen Analyseauftrag")
        self.summary.setWordWrap(True)
        self.summary.setObjectName("muted")
        layout.addWidget(self.summary)

        output_row = QHBoxLayout()
        self.open_report = QPushButton("HTML-Bericht")
        self.open_excel = QPushButton("Excel")
        self.open_folder = QPushButton("Ergebnisordner")
        for button in (self.open_report, self.open_excel, self.open_folder):
            button.setEnabled(False)
            output_row.addWidget(button)
        self.open_report.clicked.connect(lambda: self._open_result_path("report"))
        self.open_excel.clicked.connect(lambda: self._open_result_path("excel"))
        self.open_folder.clicked.connect(lambda: self._open_result_path("folder"))
        layout.addLayout(output_row)
        return panel

    @Slot()
    def _browse_bam(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Aligned BAM auswählen", filter="BAM (*.bam)")
        if not path:
            return
        bam = Path(path)
        self.bam_path.setText(str(bam))
        self.sample_id.setText(sanitize_sample_id(bam.stem))
        index = locate_bam_index(bam)
        if index is None:
            self.index_hint.setText("⚠ Kein .bai neben der BAM gefunden")
        else:
            self.index_hint.setText(f"✓ Index erkannt: {index.name}")

    @Slot()
    def _open_settings(self) -> None:
        dialog = SettingsDialog(self.config, self)
        if dialog.exec() and dialog.saved_config is not None:
            self.config = dialog.saved_config
            self.summary.setText("Einstellungen gespeichert")

    @Slot()
    def _system_check(self) -> None:
        QApplication.setOverrideCursor("WaitCursor")
        try:
            checks = DesktopBackend(self.config).diagnose()
        finally:
            QApplication.restoreOverrideCursor()
        lines = [f"{'✓' if item.ok else '×'} {item.name}: {item.detail}" for item in checks]
        all_ok = bool(checks) and all(item.ok for item in checks)
        box = QMessageBox(self)
        box.setWindowTitle("ONTSeq Desktop – Systemcheck")
        box.setIcon(QMessageBox.Information if all_ok else QMessageBox.Warning)
        box.setText("System bereit" if all_ok else "Konfiguration unvollständig")
        box.setDetailedText("\n".join(lines))
        box.setInformativeText("\n".join(lines[:5]))
        box.exec()

    @Slot()
    def _start_analysis(self) -> None:
        if self.thread is not None:
            return
        bam_text = self.bam_path.text().strip()
        sample_text = self.sample_id.text().strip()
        if not bam_text or not sample_text:
            QMessageBox.warning(self, "Eingaben fehlen", "Bitte BAM-Datei und Sample-ID angeben.")
            return
        sample = sanitize_sample_id(sample_text)
        build = self.genome_build.currentData()
        assay = self.assay_mode.currentData()
        try:
            self.config.reference_for(build)
        except Exception as exc:
            QMessageBox.warning(
                self,
                "Referenz nicht konfiguriert",
                f"{exc}\n\nBitte unter Einstellungen → Referenzen konfigurieren.",
            )
            return

        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        output_dir = Path(self.config.output_root).expanduser() / sample / timestamp
        request = DesktopAnalysisRequest(
            bam_path=Path(bam_text),
            sample_id=sample,
            genome_build=build,
            assay_mode=assay,
            output_dir=output_dir,
        )
        self._reset_run_ui()
        self.start_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.summary.setText("Analyse läuft …")
        self.log_view.clear()

        self.thread = QThread(self)
        self.worker = AnalysisWorker(DesktopBackend(self.config), request)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self._on_progress)
        self.worker.log.connect(self.log_view.append)
        self.worker.finished.connect(self._on_finished)
        self.worker.failed.connect(self._on_failed)
        self.worker.finished.connect(self.thread.quit)
        self.worker.failed.connect(self.thread.quit)
        self.thread.finished.connect(self._thread_finished)
        self.thread.start()

    @Slot()
    def _cancel_analysis(self) -> None:
        if self.worker is not None:
            self.cancel_button.setEnabled(False)
            self.summary.setText("Abbruch angefordert …")
            self.worker.cancel()

    @Slot(object)
    def _on_progress(self, payload: object) -> None:
        if not isinstance(payload, ProgressEvent):
            return
        self.stage_rows[payload.stage].set_status(payload.status, payload.message)
        self.progress.setValue(max(self.progress.value(), payload.percent))
        self.summary.setText(payload.message)

    @Slot(object)
    def _on_finished(self, payload: object) -> None:
        if not isinstance(payload, DesktopRunResult):
            self._on_failed("Backend returned an unexpected result object")
            return
        self.last_result = payload
        self.progress.setValue(100)
        self.summary.setText(
            f"Abgeschlossen · QC {payload.pipeline_result.qc.verdict.value} · "
            f"{len(payload.pipeline_result.events)} SV-Kandidat(en) · "
            "fachliche Prüfung erforderlich"
        )
        self.open_report.setEnabled(True)
        self.open_excel.setEnabled(True)
        self.open_folder.setEnabled(True)
        self.start_button.setEnabled(True)
        self.cancel_button.setEnabled(False)

    @Slot(str)
    def _on_failed(self, message: str) -> None:
        self.summary.setText(f"Analyse fehlgeschlagen: {message}")
        self.start_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        QMessageBox.critical(
            self,
            "ONTSeq Analyse fehlgeschlagen",
            f"{message}\n\nTechnische Details bleiben lokal im Ergebnisordner erhalten.",
        )

    @Slot()
    def _thread_finished(self) -> None:
        if self.thread is not None:
            self.thread.deleteLater()
        self.thread = None
        self.worker = None

    @Slot()
    def _toggle_log(self) -> None:
        visible = not self.log_view.isVisible()
        self.log_view.setVisible(visible)
        self.details_button.setText(
            "Technische Details ausblenden" if visible else "Technische Details anzeigen"
        )

    def _reset_run_ui(self) -> None:
        self.last_result = None
        self.progress.setValue(0)
        for row in self.stage_rows.values():
            row.set_status(DesktopStageStatus.PENDING, "Ausstehend")
        for button in (self.open_report, self.open_excel, self.open_folder):
            button.setEnabled(False)

    def _open_result_path(self, kind: str) -> None:
        result = self.last_result
        if result is None:
            return
        path = {
            "report": result.report_html,
            "excel": result.workbook_xlsx,
            "folder": result.output_dir,
        }[kind]
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))

    def _apply_styles(self) -> None:
        self.setStyleSheet(
            """
            QMainWindow, QWidget { background: #f5f7fa; color: #17202a; font-family: 'Segoe UI'; font-size: 10pt; }
            QLabel#title { font-size: 24pt; font-weight: 700; color: #102a43; }
            QLabel#sectionTitle { font-size: 13pt; font-weight: 700; color: #102a43; padding-bottom: 6px; }
            QLabel#muted { color: #627d98; }
            QLabel#researchBanner { background: #fff3cd; color: #6b4f00; border: 1px solid #e6c65c; border-radius: 7px; padding: 10px; font-weight: 600; }
            QFrame#card { background: white; border: 1px solid #d9e2ec; border-radius: 10px; }
            QFrame[runStatus="RUNNING"] { background: #eef6ff; border: 1px solid #9fc5e8; border-radius: 7px; }
            QFrame[runStatus="PASS"] { background: #effaf3; border: 1px solid #a8d5b5; border-radius: 7px; }
            QFrame[runStatus="WARN"] { background: #fff8e8; border: 1px solid #e6c65c; border-radius: 7px; }
            QFrame[runStatus="FAIL"] { background: #fff0f0; border: 1px solid #e6a3a3; border-radius: 7px; }
            QFrame[runStatus="NO_CALL"], QFrame[runStatus="NOT_RUN"] { background: #f7f8fa; border: 1px solid #d9e2ec; border-radius: 7px; }
            QLabel#stageTitle { font-weight: 600; }
            QLineEdit, QComboBox, QTextEdit { background: white; border: 1px solid #bcccdc; border-radius: 6px; padding: 7px; }
            QPushButton { background: white; border: 1px solid #bcccdc; border-radius: 6px; padding: 7px 12px; }
            QPushButton:hover { background: #edf2f7; }
            QPushButton:disabled { color: #9aa5b1; background: #f1f3f5; }
            QPushButton#primaryButton { background: #0b5fa5; color: white; border: none; font-weight: 700; padding: 9px 14px; }
            QPushButton#primaryButton:hover { background: #084c84; }
            QProgressBar { border: 1px solid #bcccdc; border-radius: 6px; text-align: center; background: #edf2f7; min-height: 20px; }
            QProgressBar::chunk { background: #0b5fa5; border-radius: 5px; }
            QGroupBox { font-weight: 600; border: 1px solid #d9e2ec; border-radius: 8px; margin-top: 10px; padding-top: 12px; }
            QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px; }
            """
        )


def run() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("ONTSeq Desktop")
    app.setOrganizationName("ONTSeq Platform")
    app.setStyle("Fusion")
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(run())

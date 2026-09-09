import { useCallback, useEffect, useRef, useState } from "react";
import { createLiveApi, isTerminalJob, requestedProfile, resolveAvailableProfile, watchRun, workspaceHomeHref } from "./liveApi.js";
import { MethylationOptions, MethylationResults, useMethylationProbe } from "./MethylationOptions.jsx";
import "./liveWorkspace.css";

const NAVIGATION = [["input", "01", "Analyse vorbereiten"], ["execution", "02", "Ausführung"], ["evidence", "03", "Ergebnisse"], ["exports", "04", "Bericht & Dateien"]];
const JOB_LABELS = { running: "Analyse läuft", passed: "Technisch bestanden", failed: "Technisch fehlgeschlagen", error: "Ausführungsfehler" };
const safeValue = (value) => value === null || value === undefined || value === "" ? "Nicht dokumentiert" : String(value);

function Status({ value }) {
  const tone = ["COMPLETED", "PASS", "passed"].includes(value) ? "good"
    : ["FAILED", "FAIL", "failed", "error"].includes(value) ? "bad"
      : ["NO_CALL", "WARN", "WARNING", "running"].includes(value) ? "caution" : "neutral";
  return <span className={`live-status live-status--${tone}`}>{safeValue(value)}</span>;
}

function Notice({ children, error = false }) {
  return <p className={`live-notice${error ? " live-notice--error" : ""}`} role={error ? "alert" : undefined}>{children}</p>;
}

function Facts({ entries }) {
  return <dl className="live-facts">{entries.map(([label, value]) => <div key={label}><dt>{label}</dt><dd>{safeValue(value)}</dd></div>)}</dl>;
}

function StageList({ stages }) {
  if (!stages?.length) return <Notice>Die Analyseschritte werden vorbereitet. Ihr Status erscheint, sobald die Ausführung beginnt.</Notice>;
  return <ol className="live-stages">{stages.map((stage, index) => <li key={stage.stage ?? index}>
    <div className="live-stage-heading"><strong>{stage.title ?? stage.stage}</strong><Status value={stage.status} /></div>
    <p>{stage.reason}</p>
    <small>Adapterprüfung: {safeValue(stage.verification)}{stage.resumed ? " · Vorhandenes Ergebnis wiederverwendet" : ""}{typeof stage.duration_seconds === "number" ? ` · ${stage.duration_seconds.toFixed(1)} s` : ""}</small>
  </li>)}</ol>;
}

function locusText(locus) {
  return locus ? `${locus.chromosome}:${locus.start}–${locus.end}` : "—";
}

function EventList({ events }) {
  const [page, setPage] = useState(0);
  const pageSize = 25;
  const pageCount = Math.max(1, Math.ceil(events.length / pageSize));
  const currentPage = Math.min(page, pageCount - 1);
  if (!events.length) return <Notice>Keine genomischen Ereignisse im Ergebnis. Das schließt eine Veränderung nicht aus: Qualitätskontrolle, Modulstatus und Auswertungsgrenzen sind mit zu bewerten.</Notice>;
  return <>
    <div className="live-event-head"><span>{events.length} technische Ereignisse · keine klinische Freigabe</span><small>Koordinaten unverändert: 0-basiert, Ende exklusiv</small></div>
    <div className="live-event-list">{events.slice(currentPage * pageSize, (currentPage + 1) * pageSize).map((event, index) => <details key={event.event_id ?? index}>
      <summary><span><strong>{event.event_type}</strong><small>{event.event_id}</small></span><span className="live-locus">{locusText(event.primary)}{event.secondary ? <><br />{locusText(event.secondary)}</> : null}</span><span>{event.genes?.join(" · ") || "Keine Genzuordnung"}</span><span className="live-event-expand">Details</span></summary>
      <Facts entries={[["Validierungsstatus", event.validation_status], ["Fusionstatus", event.fusion_status], ["Beurteilbarkeit", event.observability], ["Caller", event.evidence?.map((item) => `${item.caller} ${item.caller_version}`).join(" · ")], ["Technische Flags", event.technical_flags?.join(" · ")]]} />
      <p className="live-subtle">Die vollständigen Originalfelder einschließlich Quellenhinweisen und Einschränkungen:</p>
      <pre>{JSON.stringify(event, null, 2)}</pre>
    </details>)}</div>
    {pageCount > 1 ? <div className="live-pagination"><button disabled={!currentPage} onClick={() => setPage(currentPage - 1)}>Zurück</button><span>Seite {currentPage + 1} / {pageCount}</span><button disabled={currentPage + 1 >= pageCount} onClick={() => setPage(currentPage + 1)}>Weiter</button></div> : null}
  </>;
}

export function LiveWorkspace() {
  const [api] = useState(() => createLiveApi());
  const [desktopProfile] = useState(() => requestedProfile());
  const [brandHref] = useState(() => workspaceHomeHref());
  const [config, setConfig] = useState(null);
  const [configError, setConfigError] = useState("");
  const [profileRequestError, setProfileRequestError] = useState("");
  const [configRevision, setConfigRevision] = useState(0);
  const [connecting, setConnecting] = useState(!api.disabledReason);
  const [listing, setListing] = useState(null);
  const [browseError, setBrowseError] = useState("");
  const [browsing, setBrowsing] = useState(false);
  const [bam, setBam] = useState(null);
  const bamSelectionSequence = useRef(0);
  const [profile, setProfile] = useState(desktopProfile);
  const methylation = useMethylationProbe(api, bam?.posix, bam?.selection_id, profile);
  const [sampleId, setSampleId] = useState("");
  const [job, setJob] = useState(null);
  const [submitting, setSubmitting] = useState(false);
  const [submitError, setSubmitError] = useState("");
  const [submitDetail, setSubmitDetail] = useState("");
  const [pollError, setPollError] = useState("");
  const [pollPaused, setPollPaused] = useState(false);
  const [pollRevision, setPollRevision] = useState(0);
  const [result, setResult] = useState(null);
  const [resultError, setResultError] = useState("");
  const [resultRevision, setResultRevision] = useState(0);
  const [loadingResult, setLoadingResult] = useState(false);
  const [downloadError, setDownloadError] = useState("");
  const [downloading, setDownloading] = useState("");
  const browseController = useRef(null);
  const submitController = useRef(null);
  const downloadController = useRef(null);
  const submitted = useRef(false);
  const running = job?.state === "running";
  const terminal = isTerminalJob(job);

  useEffect(() => {
    if (api.disabledReason) return undefined;
    const controller = new AbortController();
    setConnecting(true);
    setConfigError("");
    api.config(controller.signal).then((value) => {
      if (controller.signal.aborted) return;
      setConfig(value);
      const resolution = resolveAvailableProfile(value.profiles, profile, desktopProfile);
      setProfile(resolution.profile);
      setProfileRequestError(resolution.error);
    }).catch((error) => { if (!controller.signal.aborted) { setConfig(null); setConfigError(error.message); setProfileRequestError(""); } })
      .finally(() => { if (!controller.signal.aborted) setConnecting(false); });
    return () => controller.abort();
  }, [api, configRevision, desktopProfile]);

  const browse = useCallback(async (path = "") => {
    browseController.current?.abort();
    const controller = new AbortController();
    browseController.current = controller;
    setBrowsing(true);
    setBrowseError("");
    try {
      const value = await api.browse(path, controller.signal);
      if (!controller.signal.aborted) setListing(value);
    } catch (error) { if (!controller.signal.aborted) setBrowseError(error.message); }
    finally { if (!controller.signal.aborted) setBrowsing(false); }
  }, [api]);

  useEffect(() => {
    if (config?.roots.length) browse(config.roots[0].posix);
    return () => browseController.current?.abort();
  }, [config?.roots, browse]);

  useEffect(() => () => {
    browseController.current?.abort();
    submitController.current?.abort();
    downloadController.current?.abort();
  }, []);

  useEffect(() => {
    if (!job?.run_id || !running) return undefined;
    const controller = new AbortController();
    setPollPaused(false);
    watchRun(api, job.run_id, {
      signal: controller.signal,
      onSnapshot: (value) => { setJob(value); setPollError(""); },
      onError: (error, stopped) => { setPollError(error.message); setPollPaused(stopped); },
    });
    return () => controller.abort();
  }, [api, job?.run_id, running, pollRevision]);

  useEffect(() => {
    if (!terminal) return undefined;
    const controller = new AbortController();
    setLoadingResult(true);
    setResultError("");
    api.result(job, controller.signal).then((value) => { if (!controller.signal.aborted) setResult(value); })
      .catch((error) => { if (!controller.signal.aborted) { setResult(null); setResultError(error.message); } })
      .finally(() => { if (!controller.signal.aborted) setLoadingResult(false); });
    return () => controller.abort();
  }, [api, job?.run_id, terminal, resultRevision]);

  useEffect(() => { if (terminal) setConfigRevision((value) => value + 1); }, [terminal]);

  async function startRun(event) {
    event.preventDefault();
    if (submitted.current || !bam || !methylation.ready) return;
    submitted.current = true;
    const controller = new AbortController();
    submitController.current = controller;
    setSubmitting(true);
    setSubmitError("");
    setSubmitDetail("");
    try {
      if (!await methylation.confirmBeforeStart(controller.signal)) {
        if (!controller.signal.aborted) setSubmitError("Die BAM-Prüfung hat sich geändert. Bitte die aktualisierte Information prüfen und den Analyseumfang erneut auswählen. Es wurde kein Lauf gestartet.");
        return;
      }
      const value = await api.start({ bam: bam.posix, profile, sample_id: sampleId.trim(), include_methylation: methylation.decision === true }, controller.signal);
      if (controller.signal.aborted) return;
      methylation.clearChoice();
      setResult(null);
      setResultError("");
      setDownloadError("");
      setJob(value);
      document.getElementById("live-execution")?.scrollIntoView({ behavior: "smooth", block: "start" });
    } catch (error) {
      if (!controller.signal.aborted) {
        const referenceMismatch = error.status === 400 && /dictionary|canonical assembly/i.test(error.message);
        setSubmitError(referenceMismatch
          ? "Die BAM passt nicht zu einem vollständigen unterstützten Referenzgenom. Bitte BAM und Referenzprofil prüfen. Es wurde kein Lauf gestartet."
          : "Der Analysestart wurde nicht bestätigt. Bitte zuerst die Eingaben und den Dienststatus prüfen, bevor Sie erneut starten.");
        setSubmitDetail(error.message);
        setConfigRevision((value) => value + 1);
      }
    } finally {
      submitted.current = false;
      if (!controller.signal.aborted) setSubmitting(false);
    }
  }

  async function download(kind) {
    if (!result || loadingResult || resultError) return;
    downloadController.current?.abort();
    const controller = new AbortController();
    downloadController.current = controller;
    setDownloading(kind);
    setDownloadError("");
    try {
      const { blob } = await api.artifact(job, kind, controller.signal);
      if (controller.signal.aborted) return;
      const url = URL.createObjectURL(blob);
      const anchor = document.createElement("a");
      anchor.href = url;
      anchor.download = `ONTSeq-${job.run_id}-${job.sample_id}.${kind}`;
      document.body.append(anchor);
      anchor.click();
      anchor.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    } catch (error) { if (!controller.signal.aborted) setDownloadError(error.message); }
    finally { if (!controller.signal.aborted) setDownloading(""); }
  }

  const disconnected = api.disabledReason || configError || profileRequestError;
  const reference = result?.reference_context;
  const panelResolution = reference?.panel_resolution;
  const boundBuild = result?.manifest?.assay?.genome_build ?? job?.detected_genome_build;
  const disableInputs = !config || connecting || running || submitting || Boolean(config.busy) || Boolean(profileRequestError);
  const validSampleId = /^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$/.test(sampleId.trim());

  return <div className="live-shell">
    <header className="live-header"><a className="live-brand" href={brandHref} aria-label="ONTSeq Live-Arbeitsplatz"><span>ONT</span>Seq</a><span className="live-header-label">Genom &amp; Methylierung</span><span className="live-connection">{connecting ? "Dienst wird geprüft…" : config ? `${config.version} · Lokal` : "Nicht verbunden"}</span><span className="live-ruo">Nur für Forschungszwecke</span></header>
    <nav className="live-outline" aria-label="Arbeitsplatz-Abschnitte"><div className="live-outline-caption">ARBEITSPLATZ</div>{NAVIGATION.map(([id, number, title]) => <a key={id} href={`#live-${id}`}><span>{number}</span>{title}</a>)}<div className="live-outline-note">Ihre Analyse bleibt lokal.<br />Genom und Methylierung<br />in einem Arbeitsplatz.</div></nav>
    <main className="live-canvas">
      <div className="live-page-title"><h1>{job ? job.sample_id : "Neue Analyse"}</h1><p>BAM auswählen, Profil prüfen und Analyseumfang festlegen.</p></div>
      {disconnected ? <Notice error>{disconnected} Die Analysefunktionen bleiben gesperrt.</Notice> : null}
      {configError ? <button className="live-secondary" onClick={() => setConfigRevision((value) => value + 1)} disabled={connecting}>Verbindung erneut prüfen</button> : null}
      <section className="live-section" id="live-input"><div className="live-section-title"><span className="live-step-number">1</span><div><h2>BAM auswählen</h2><p className="live-subtle">Wählen Sie eine BAM-Datei aus dem lokalen Dateisystem.</p></div></div>
        <form onSubmit={startRun}>
          <div className="live-file-browser"><div className="live-file-toolbar"><strong>Lokale Dateien</strong><button type="button" onClick={() => setConfigRevision((value) => value + 1)} disabled={!config || connecting || submitting || running}>Profile aktualisieren</button></div>
            <div className="live-root-buttons">{config?.roots.map((root) => <button type="button" key={root.posix} onClick={() => browse(root.posix)} disabled={disableInputs || browsing} title={root.display}>{root.display}</button>)}</div>
            {listing ? <div className="live-current-directory"><button type="button" onClick={() => browse(listing.parent)} disabled={!listing.parent || disableInputs || browsing}>Übergeordnet</button><code>{listing.display}</code></div> : null}
            {browsing ? <p role="status">Verzeichnis wird gelesen…</p> : null}
            {browseError ? <Notice error>{browseError}</Notice> : null}
            {listing?.unreadable_entries ? <Notice>{listing.unreadable_entries} Einträge konnten nicht gelesen werden.</Notice> : null}
            {listing?.entries?.length ? <ul className="live-files">{listing.entries.map((entry) => <li key={entry.posix}><button type="button" aria-pressed={!entry.directory && bam?.posix === entry.posix} disabled={disableInputs || browsing || (!entry.directory && !entry.indexed)} onClick={() => entry.directory ? browse(entry.posix) : setBam({ ...entry, selection_id: ++bamSelectionSequence.current })}><span className="live-file-kind">{entry.directory ? "DIR" : "BAM"}</span><span>{entry.name}<small>{entry.directory ? "Verzeichnis öffnen" : entry.indexed ? "BAI vorhanden · Integrität wird beim Start geprüft" : "BAI fehlt – nicht startbar"}</small></span>{!entry.directory && typeof entry.size_bytes === "number" ? <small>{(entry.size_bytes / 1024 / 1024).toLocaleString("de-DE", { maximumFractionDigits: 1 })} MiB</small> : null}</button></li>)}</ul> : listing && !browsing ? <p className="live-subtle">Keine BAM-Dateien oder Unterverzeichnisse in diesem Ordner.</p> : null}
          </div>
          {bam ? <p className="live-selected-file"><strong>Ausgewählt:</strong> <code>{bam.display}</code></p> : null}
          <div className="live-setup-step"><div className="live-section-title"><span className="live-step-number">2</span><div><h2>Probe &amp; Analyseprofil</h2><p className="live-subtle">Der Referenzbuild wird vor der Analyse anhand der BAM geprüft.</p></div></div>
          <div className="live-form-grid"><label>Proben-ID<input value={sampleId} onChange={(event) => setSampleId(event.target.value)} placeholder="Pseudonymisierte Proben-ID" disabled={disableInputs} maxLength={64} pattern="[A-Za-z0-9][A-Za-z0-9._-]{2,63}" required autoComplete="off" /><small>3–64 Zeichen: Buchstaben, Ziffern, Punkt, Bindestrich, Unterstrich.</small></label><label>Analyseprofil<select value={profile} onChange={(event) => setProfile(event.target.value)} disabled={disableInputs || !config?.profiles.length} required><option value="">{connecting ? "Profile werden geprüft…" : "Kein Profil ausgewählt"}</option>{config?.profiles.map((id) => <option value={id} key={id}>{id}</option>)}</select><small>Nur lokal auflösbare Profile. Keine automatische Build-Konvertierung.</small></label></div>
          {config && !config.profiles.length ? <Notice>Kein vollständig installiertes Analyseprofil verfügbar. Referenzpaket im lokalen Setup installieren und anschließend die Profile aktualisieren.</Notice> : null}
          </div>
          <MethylationOptions bam={bam} {...methylation} disabled={disableInputs} />
          {config?.busy && !running ? <Notice>Der Dienst führt bereits eine Analyse aus. Ein zweiter Lauf bleibt gesperrt; den Status kannst du über „Profile aktualisieren“ erneut prüfen.</Notice> : null}
          {submitError ? <Notice error>{submitError}</Notice> : null}
          {submitDetail ? <details className="live-details"><summary>Technische Details zur Startanforderung</summary><p>{submitDetail}</p><p>Die Startanforderung wird nicht automatisch wiederholt.</p></details> : null}
          <div className="live-form-footer"><span>{!bam ? "Bitte eine BAM-Datei auswählen." : (methylation.probing || methylation.scanning) ? "BAM-Prüfung läuft…" : !methylation.ready ? "Bitte Analyseumfang auswählen." : methylation.decision === true ? "Genomanalyse + regionale Methylierung" : "Genomanalyse ohne Methylierung"}</span><button className="live-primary" type="submit" disabled={disableInputs || !bam?.indexed || !profile || !validSampleId || !methylation.ready}>{submitting ? "Start wird angefordert…" : running ? "Analyse läuft…" : "Analyse starten"}</button></div>
        </form>
      </section>
      <section className="live-section" id="live-execution"><div className="live-section-title"><span>02</span><h2>Technische Ausführung</h2>{job ? <Status value={job.state} /> : null}</div>
        {job ? <><div className="live-run-heading"><h3>{JOB_LABELS[job.state] ?? job.state}</h3><code>{job.run_id}</code></div><Facts entries={[["Proben-ID", job.sample_id], ["Laufprofil", job.profile], ["Laufgebundener Build", boundBuild], ["Gestartet", job.started_at], ["Beendet", job.finished_at]]} />{job.detail ? <Notice error={job.state === "failed" || job.state === "error"}>{job.detail}</Notice> : null}<StageList stages={job.stages} /><p className="live-subtle">COMPLETED = ausgeführt · NO_CALL = keine belastbare Aussage · FAILED = fehlgeschlagen · NOT_RUN = nicht ausgeführt.</p>{pollError ? <Notice error>{pollError}{pollPaused ? " Die Statusabfrage ist pausiert; der Lauf wird nicht abgebrochen." : " Die Statusabfrage wird erneut versucht."}</Notice> : null}{pollPaused ? <button onClick={() => setPollRevision((value) => value + 1)}>Statusabfrage fortsetzen</button> : null}</> : <Notice>Nach dem Start sehen Sie hier den Fortschritt und den Status der einzelnen Analyseschritte.</Notice>}
      </section>
      <section className="live-section" id="live-evidence"><div className="live-section-title"><span>03</span><h2>Laufgebundene Ergebnisse</h2>{result ? <span className="live-build-tag">{boundBuild}</span> : null}</div>
        {loadingResult ? <p role="status">Originalergebnis wird geladen…</p> : null}
        {resultError ? <><Notice error>Kein darstellbares kanonisches Ergebnis: {resultError}</Notice><button disabled={loadingResult} onClick={() => setResultRevision((value) => value + 1)}>Ergebnis erneut abrufen</button></> : null}
        {result ? <><div className="live-subheading"><h3>Qualitätskontrolle</h3><Status value={result.qc.verdict} /></div><Facts entries={Object.entries(result.qc.metrics ?? {})} />{result.qc.failed_gates?.length ? <Notice error>Fehlgeschlagene QC-Grenzen: {result.qc.failed_gates.join(" · ")}</Notice> : null}{result.qc.warnings?.map((warning, index) => <Notice key={index}>{warning}</Notice>)}<MethylationResults key={job.run_id} result={result} job={job} api={api} Status={Status} /><h3>Genomische Ereignisse</h3><EventList key={job.run_id} events={result.events} />{result.modules?.length ? <details className="live-details"><summary>Modulstatus aus dem PipelineResult</summary><ul>{result.modules.map((module) => <li key={module.module}><strong>{module.module}</strong> <Status value={module.status} /><p>{module.reason}</p></li>)}</ul></details> : null}{result.warnings?.map((warning, index) => <Notice key={index}>{warning}</Notice>)}</> : !loadingResult && !resultError ? <Notice>Nach der Analyse finden Sie hier die Qualitätskontrolle, genomische Ereignisse und die gewählte Methylierungsauswertung.</Notice> : null}
      </section>
      <section className="live-section" id="live-exports"><div className="live-section-title"><span>04</span><h2>Bericht & Originaldateien</h2></div><p className="live-subtle">Berichte und Ergebnisdateien für die fachliche Bewertung herunterladen. Die Berichte sind nicht klinisch freigegeben.</p><div className="live-export-buttons">{[["html", "HTML-Bericht"], ["xlsx", "Excel-Arbeitsmappe"], ["json", "PipelineResult JSON"]].map(([kind, label]) => <button key={kind} disabled={!result || loadingResult || Boolean(resultError) || Boolean(downloading)} onClick={() => download(kind)}>{downloading === kind ? "Wird geladen…" : `${label} herunterladen`}</button>)}</div>{downloadError ? <Notice error>{downloadError}</Notice> : null}{terminal ? <p className="live-subtle">Downloads werden erst nach erfolgreicher Prüfung des PipelineResult freigegeben. Bei einem fehlgeschlagenen Lauf können einzelne Dateien fehlen.</p> : null}<Facts entries={[["Ausgabe-Basisverzeichnis", config?.output_dir], ["Ergebnisordner dieses Laufs", job && config?.output_dir ? `${config.output_dir}/${job.run_id}/${job.sample_id}` : null]]} />{result?.sidecars?.length ? <details className="live-details"><summary>{result.sidecars.length} zusätzliche Artefakte laut PipelineResult</summary><Facts entries={result.sidecars.map((item) => [item.artifact_id, `${item.relative_path} · SHA256 ${item.sha256}`])} /></details> : null}</section>
      <footer className="live-page-footer">Forschungssoftware · Technische Prüfung ersetzt keine fachliche Bewertung, klinische Validierung oder autorisierte Freigabe.</footer>
    </main>
    <aside className="live-inspector" aria-label="Analyseumfang und Herkunft"><div className="live-scope"><h2>Analyseumfang</h2><div className="live-scope-item"><h3>Genomanalyse</h3><p>CNV, strukturelle Varianten und weitere Module gemäß gewähltem Analyseprofil.</p></div><div className="live-scope-item"><h3>Methylierung</h3><p>{(methylation.probing || methylation.scanning) ? "BAM-Prüfung läuft…" : methylation.decision === true ? "Regionale Auswertung ausgewählt" : methylation.decision === false || methylation.probe?.status === "not_detected" ? "Für den nächsten Lauf nicht ausgewählt" : bam ? "Entscheidung offen" : "Nach der BAM-Auswahl verfügbar"}</p></div><details className="live-research"><summary>Mischungsanalyse · Forschung</summary><p>Separater Versuchsablauf mit zwei Quellen: Kalibration, zurückgehaltene Testreads und definierte Mischungen. Nicht Teil einer gewöhnlichen Einzelprobenanalyse.</p><p>Im gemeinsamen ONTSeq-Paket über die Befehle <code>ontseq methylation-mixture</code> verfügbar. Die genaue Befehlshilfe beschreibt die erforderlichen Eingaben.</p></details></div><details className="live-provenance" open={Boolean(result)}><summary>Herkunft &amp; Prüfstatus</summary>{result ? <><p className="live-subtle">Nur dokumentierte Angaben aus dem geladenen Originalergebnis.</p><Facts entries={[["Referenzbuild", result?.manifest?.assay?.genome_build], ["Referenz-ID", result?.manifest?.assay?.reference_id], ["Referenzpaket", reference?.reference_bundle_id], ["Dictionary-Vertrag", reference?.reference_dictionary_contract], ["Panelpaket", reference?.panel_bundle_id], ["Wissenspaket", reference?.knowledge_bundle_id], ["Wissensversion", reference?.knowledge_bundle_version], ["Pipeline-Version", result?.provenance?.pipeline_version], ["Git-Commit", result?.provenance?.git_commit], ["Ergebnis erstellt", result?.provenance?.created_at], ["Eingabe-SHA256 im Manifest", result?.manifest?.input?.sha256], ["Reviewstatus im Ergebnis", result?.release_status]]} />{reference?.panel_bundle_id ? <details className="live-details"><summary>Panel-Koordinatenauflösung</summary><Facts entries={panelResolution ? [["Mappingstatus", panelResolution.mapping_status], ["Mapping-ID", panelResolution.mapping_id], ["Quellpanel", panelResolution.source_panel_bundle_id ? `${panelResolution.source_panel_bundle_id}:${panelResolution.source_panel_resource_id}` : "Native Koordinaten"], ["Build-Pfad", `${panelResolution.source_genome_build} → ${panelResolution.target_genome_build}`], ["Finale Selektionsintervalle", panelResolution.selection_interval_count], ["Native Analyse-ROI-Intervalle", panelResolution.analysis_roi_interval_count], ["Vorwärts gemappte Intervalle", panelResolution.mapped_interval_count == null ? "Nicht erforderlich" : `${panelResolution.mapped_interval_count} von ${panelResolution.source_interval_count}`], ["Vorwärts nicht gemappte Ziele", panelResolution.unmapped_target_labels?.join(" · ") || "Keine Bezeichnungen dokumentiert"], ["Roundtrip-Prüfziele", panelResolution.roundtrip_review_required_target_labels?.join(" · ") || "Keine Bezeichnungen dokumentiert"], ["Ungelöste/zu prüfende Panelziele", panelResolution.unresolved_target_labels?.join(" · ") || "Keine Bezeichnungen dokumentiert"]] : [["Mappingstatus", "In diesem älteren Ergebnis nicht dokumentiert"]]} /><Notice>Koordinaten-Mapping dokumentiert die Auflösung des Paneldesigns; gemappt bedeutet nicht analytisch validiert. Nicht gemappte oder zu prüfende Panelziele sind keine negativen Befunde dieser Probe.</Notice></details> : null}{reference?.resource_checksums ? <details className="live-details"><summary>Ressourcen-Prüfsummen</summary><Facts entries={Object.entries(reference.resource_checksums)} /></details> : null}{result?.provenance?.tools?.length ? <details className="live-details"><summary>Dokumentierte Werkzeuge</summary><Facts entries={result.provenance.tools.map((tool, index) => [`${tool.name} (${index + 1})`, `${tool.version}${tool.container_digest ? ` · ${tool.container_digest}` : " · Kein Container-Digest dokumentiert"}`])} /></details> : null}</> : <p className="live-subtle">Referenzbuild, Versionen und Prüfsummen erscheinen mit dem geprüften Laufergebnis.</p>}</details><Notice>Keine elektronische Signatur oder klinische Freigabe durch diese Ansicht. Ein bestandener technischer Lauf ist kein freigegebener Befund.</Notice></aside>
  </div>;
}

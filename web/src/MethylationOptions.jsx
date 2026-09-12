import { useEffect, useRef, useState } from "react";
import { methylationExplanation, methylationReady, sameMethylationEvidence, selectedMethylationState } from "./liveApi.js";

export function useMethylationProbe(api, bamPath, selectionId, profile) {
  const [probe, setProbe] = useState(null);
  const [choice, setChoice] = useState(null);
  const [revision, setRevision] = useState(0);
  const [scan, setScan] = useState(null);
  const [scanError, setScanError] = useState("");
  const activeScan = useRef(null);
  const generation = useRef(0);
  useEffect(() => {
    const currentGeneration = ++generation.current;
    setProbe(null);
    setChoice(null);
    setScan(null);
    setScanError("");
    if (!bamPath) return undefined;
    const controller = new AbortController();
    api.probeMethylation(bamPath, controller.signal).then((value) => {
      if (!controller.signal.aborted && generation.current === currentGeneration) setProbe({ ...value, selection_id: selectionId });
    }).catch((error) => {
      if (!controller.signal.aborted && generation.current === currentGeneration) setProbe({ bam_path: bamPath, selection_id: selectionId, status: "unknown", reason_code: "transport_error", reason: error.message, complete: false, checked_reads: 0 });
    });
    return () => {
      controller.abort();
      generation.current++;
      const lease = activeScan.current;
      if (lease) {
        lease.cancelled = true;
        lease.controller.abort();
        // Acquisition is deliberately not aborted: its response identifies the job
        // that must be cancelled even when selection changes during the POST.
        if (lease.id) api.cancelMethylationScan(lease.id, lease.bamPath).catch(() => {});
        activeScan.current = null;
      }
    };
  }, [api, bamPath, selectionId, profile, revision]);
  const { probe: current, decision } = selectedMethylationState(probe, choice, bamPath, selectionId);
  const scanning = Boolean(scan && ["starting", "running"].includes(scan.state));
  const applyScan = (snapshot) => {
    setScan(snapshot);
    if (snapshot.result) setProbe({ ...snapshot.result, selection_id: selectionId });
  };
  const startThorough = async () => {
    if (!bamPath || activeScan.current) return;
    const lease = { bamPath, id: null, cancelled: false, stopRequested: false, controller: new AbortController(), generation: generation.current };
    activeScan.current = lease;
    setChoice(null);
    setScanError("");
    setScan({ state: "starting", checked_reads: 0, elapsed_seconds: 0 });
    const isCurrent = () => activeScan.current === lease && generation.current === lease.generation;
    try {
      let snapshot = await api.startMethylationScan(bamPath);
      lease.id = snapshot.scan_id;
      if (lease.cancelled || !isCurrent()) {
        await api.cancelMethylationScan(lease.id, bamPath);
        return;
      }
      if (lease.stopRequested) snapshot = await api.cancelMethylationScan(lease.id, bamPath);
      if (lease.cancelled || !isCurrent()) return;
      applyScan(snapshot);
      while (snapshot.state === "running" && !lease.cancelled && isCurrent()) {
        await new Promise((resolve) => setTimeout(resolve, 500));
        if (lease.cancelled || !isCurrent()) return;
        snapshot = await api.methylationScan(lease.id, bamPath, lease.controller.signal);
        if (isCurrent() && !lease.cancelled) applyScan(snapshot);
      }
    } catch (error) {
      if (isCurrent() && !lease.cancelled) {
        setScanError(error.message);
        setScan({ state: "failed", checked_reads: 0, elapsed_seconds: 0 });
        setProbe({ bam_path: bamPath, selection_id: selectionId, status: "unknown", reason_code: "transport_error", reason: error.message, complete: false, checked_reads: 0 });
      }
      if (lease.id) api.cancelMethylationScan(lease.id, bamPath).catch(() => {});
    } finally {
      if (activeScan.current === lease) activeScan.current = null;
    }
  };
  const cancelThorough = async () => {
    const lease = activeScan.current;
    if (!lease) return;
    setScanError("");
    if (!lease.id) {
      // The acquisition handler will cancel once the server returns its identifier.
      lease.stopRequested = true;
      return;
    }
    try {
      const snapshot = await api.cancelMethylationScan(lease.id, lease.bamPath);
      if (activeScan.current !== lease || generation.current !== lease.generation) return;
      lease.stopRequested = true;
      applyScan(snapshot);
      if (snapshot.state !== "running") {
        lease.cancelled = true;
        lease.controller.abort();
        if (!snapshot.result) setProbe({ bam_path: bamPath, selection_id: selectionId, status: "unknown", reason_code: "cancelled", reason: "Prüfung abgebrochen.", checked_reads: snapshot.checked_reads, complete: false });
      }
    } catch (error) {
      if (activeScan.current === lease) setScanError(`Abbruch nicht bestätigt: ${error.message}`);
    }
  };
  return {
    probe: current, decision, probing: Boolean(bamPath && !current), scan, scanning, scanError,
    startThorough, cancelThorough,
    ready: !scanning && methylationReady(current, bamPath, decision),
    choose: (value) => setChoice({ bam_path: bamPath, selection_id: selectionId, value }),
    clearChoice: () => setChoice(null),
    confirmBeforeStart: async (signal) => {
      if (activeScan.current || scanning) return false;
      let fresh;
      try {
        fresh = await api.probeMethylation(bamPath, signal, true);
      } catch (error) {
        if (signal?.aborted) throw error;
        fresh = { bam_path: bamPath, status: "unknown", reason: error.message, complete: false, checked_reads: 0 };
      }
      if (signal?.aborted) return false;
      setProbe({ ...fresh, selection_id: selectionId });
      if (!sameMethylationEvidence(current, fresh)) {
        setChoice(null);
        return false;
      }
      return methylationReady(fresh, bamPath, decision);
    },
    retry: () => { setProbe(null); setChoice(null); setRevision((value) => value + 1); },
  };
}

export function MethylationOptions({ bam, probe, probing, decision, choose, retry, disabled, scan, scanning, scanError, startThorough, cancelThorough }) {
  return <fieldset className="live-methylation" disabled={disabled}>
    <legend><span className="live-step-number">3</span>Methylierung mitbeurteilen?</legend>
    {!bam ? <p className="live-subtle">Nach der Dateiauswahl wird geprüft, ob die BAM auswertbare Methylierungsinformationen enthält.</p>
      : scanning ? <div className="live-detection" role="status"><strong>Gründliche BAM-Prüfung läuft</strong><span>{scan.checked_reads.toLocaleString("de-DE")} Reads geprüft · {Math.floor(scan.elapsed_seconds)} Sekunden. Die Prüfung kann mehrere Minuten dauern und endet bei einem passenden Fund vorzeitig.</span><button type="button" onClick={cancelThorough}>Prüfung abbrechen</button></div>
      : probing ? <p className="live-detection" role="status">Methylierungsinformationen werden geprüft…</p>
        : <>
          <div className={`live-detection live-detection--${probe?.status}`} role="status">
            <strong>{probe?.status === "detected" ? "Methylierungsinformationen erkannt" : probe?.status === "not_detected" ? "Keine Methylierungsinformationen erkannt" : "Methylierungsstatus noch unklar"}</strong>
            <span>{probe?.status === "detected" ? probe.methylation_available === true ? "Regionale Methylierung kann zusätzlich ausgewertet werden. Sie entscheiden vor dem Start." : "Methylierungsinformationen sind vorhanden. Die zusätzliche Auswertung ist in dieser Installation noch nicht verfügbar." : probe?.status === "not_detected" ? "Die Genomanalyse kann ohne Methylierung gestartet werden." : "Die Prüfung konnte die Verfügbarkeit nicht sicher feststellen. Sie können erneut prüfen oder bewusst ohne Methylierung starten."}</span>
          </div>
          <p className="live-subtle">{methylationExplanation(probe)}</p>
          {probe?.status === "detected" ? <>
            <p className="live-option-question">Regionale Methylierung zusätzlich auswerten?</p>
            <div className="live-choice-row">{[[true, "Ja, Methylierung ergänzen"], [false, "Nein, nur Genomanalyse"]].map(([value, label]) => <label key={String(value)}><input type="radio" name="methylation" disabled={value && probe.methylation_available !== true} checked={decision === value} onChange={() => choose(value)} />{label}</label>)}</div>
            {probe.methylation_available !== true ? <p className="live-notice">Für die zusätzliche Auswertung muss modkit{probe.expected_modkit_version ? ` ${probe.expected_modkit_version}` : ""} mit der passenden Auswertungspolicy eingerichtet sein. Die Genomanalyse bleibt nach Auswahl von „Nein“ verfügbar.</p> : null}
            <p className="live-subtle">Experimentelle Auswertung. Erkannte Tags belegen weder ausreichende Qualität noch klinische Validierung.</p>
          </> : probe?.status === "unknown" ? <div className="live-choice-row"><button type="button" onClick={retry}>Schnellprüfung wiederholen</button><button type="button" onClick={startThorough}>Gründlicher prüfen</button><label><input type="checkbox" checked={decision === false} onChange={(event) => choose(event.target.checked ? false : null)} />Ohne Methylierung fortfahren</label></div> : null}
          <details className="live-details"><summary>Details der BAM-Prüfung</summary><p>{probe?.reason}</p>{probe?.methylation_unavailable_reason ? <p>{probe.methylation_unavailable_reason}</p> : null}<p>{probe?.checked_reads?.toLocaleString("de-DE")} Reads geprüft · {probe?.complete ? "Datei vollständig geprüft" : "Begrenzte Vorprüfung"}. Die vollständige Eingabeprüfung folgt beim Analysestart.</p></details>
        </>}
    {scanError ? <p className="live-notice live-notice--error" role="alert">{scanError}</p> : null}
  </fieldset>;
}

export function MethylationResults({ result, job, api, Status }) {
  const outcome = result.modules?.find((module) => module.module === "methylation");
  const artifacts = result.sidecars?.filter((item) => /methylation|modkit/i.test(`${item.artifact_id} ${item.relative_path}`)) ?? [];
  const [report, setReport] = useState(null);
  const [error, setError] = useState("");
  const [revision, setRevision] = useState(0);
  const [page, setPage] = useState(0);
  const completed = ["COMPLETED", "NO_CALL"].includes(outcome?.status);
  useEffect(() => {
    if (!completed) return undefined;
    const controller = new AbortController();
    setReport(null);
    setError("");
    setPage(0);
    api.methylation({ ...job, detected_genome_build: result.manifest.assay.genome_build }, controller.signal)
      .then((value) => { if (!controller.signal.aborted) setReport(value); })
      .catch((failure) => { if (!controller.signal.aborted) setError(failure.message); });
    return () => controller.abort();
  }, [api, job.run_id, completed, revision, result.manifest.assay.genome_build]);
  const pageCount = Math.max(1, Math.ceil((report?.regions.length ?? 0) / 25));
  const format = (value, fraction = false) => typeof value === "number" && Number.isFinite(value)
    ? `${(fraction ? value * 100 : value).toLocaleString("de-DE", { maximumFractionDigits: 2 })}${fraction ? " %" : ""}` : "Keine Aussage";
  return <div className="live-methylation-result"><div className="live-subheading"><h3>Regionale Methylierung</h3><Status value={outcome?.status ?? "NOT_RUN"} /></div>
    <p className="live-subtle">{outcome?.reason ?? "Für diesen Lauf ist keine Methylierungsauswertung dokumentiert."}</p>
    {completed && !report && !error ? <p role="status">Regionale Ergebnisse werden geladen…</p> : null}
    {error ? <><p className="live-notice live-notice--error" role="alert">{error}</p><button type="button" onClick={() => setRevision((value) => value + 1)}>Methylierungsergebnis erneut laden</button></> : null}
    {report ? <>
      <div className="live-subheading"><span>Regionale Auswertung</span><Status value={report.status} /></div>
      <p className="live-subtle">Methylierungsanteil = modifizierte Aufrufe / gültige Aufrufe über Sites mit ausreichender Abdeckung. Kein Zellanteil.</p>
      {report.regions.length ? <div className="live-region-table-wrap"><table className="live-region-table"><thead><tr><th scope="col">Region</th><th scope="col">Modifikation</th><th scope="col">Methylierungsanteil</th><th scope="col">Ø Abdeckung</th><th scope="col">Sites ≥ Mindestabdeckung</th></tr></thead><tbody>{report.regions.slice(page * 25, (page + 1) * 25).map((region, index) => <tr key={`${region.region_id}-${region.modification_code}-${index}`}><td>{region.region_id}<small>{region.start == null ? region.chromosome : `${region.chromosome}:${region.start}–${region.end}`}</small></td><td>{region.modification_name}</td><td>{format(region.mean_modified_fraction, true)}</td><td>{format(region.mean_valid_coverage)}</td><td>{format(region.sites_at_minimum_coverage)} / {format(region.sites_total)}</td></tr>)}</tbody></table></div> : <p className="live-subtle">Keine regionalen Werte im Bericht vorhanden.</p>}
      {pageCount > 1 ? <div className="live-pagination"><button disabled={!page} onClick={() => setPage(page - 1)}>Vorherige Regionen</button><span>{page + 1} / {pageCount}</span><button disabled={page + 1 >= pageCount} onClick={() => setPage(page + 1)}>Weitere Regionen</button></div> : null}
      <p className="live-subtle">Koordinaten: 0-basiert, Ende exklusiv. Abdeckungsregeln gemäß versionierter Methylierungspolicy.</p>
      <details className="live-details"><summary>Auswertungsregeln &amp; Werkzeug</summary><p>Policy: {report.policy?.profile_id ?? "Nicht dokumentiert"} · Mindestabdeckung: {format(report.policy?.minimum_valid_coverage)} · {report.tool?.name ?? "modkit"} {report.tool?.version ?? "Version nicht dokumentiert"}</p></details>
      {[...(report.warnings ?? []), ...(report.limitations ?? [])].map((message, index) => <p className="live-notice" key={index}>{message}</p>)}
    </> : null}
    {outcome?.tools?.length ? <p className="live-subtle">{outcome.tools.map((tool) => `${tool.name} ${tool.version}`).join(" · ")}</p> : null}
    {artifacts.length ? <details className="live-details"><summary>Regionale Ergebnisdateien ({artifacts.length})</summary><p>Die Originaldateien liegen im lokalen Laufordner. Pfad und Prüfsumme sind im Ergebnis dokumentiert.</p><ul>{artifacts.map((item) => <li key={item.artifact_id}><code>{item.relative_path}</code><small>SHA256 {item.sha256}</small></li>)}</ul></details> : null}
    <p className="live-subtle">Forschungsauswertung; keine validierte Tumor- oder Blastenzellfraktion.</p>
  </div>;
}

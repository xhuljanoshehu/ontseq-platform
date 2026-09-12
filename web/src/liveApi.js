// The live workspace has no dependency on the prototype's seed data or report state.
// Session credentials remain in memory and are sent only to this page's own origin.
export const TOKEN_HEADER = "X-ONTSeq-Token";
const LOOPBACK_HOSTS = new Set(["127.0.0.1", "localhost", "[::1]", "::1"]);
const IDENTIFIER = /^[A-Za-z0-9][A-Za-z0-9._-]{2,63}$/;
const TERMINAL_STATES = new Set(["passed", "failed", "error"]);

export class LiveApiError extends Error {
  constructor(message, status = 0) {
    super(message);
    this.name = "LiveApiError";
    this.status = status;
  }
}

export function getLiveSession(location, document) {
  if (!location || !["http:", "https:"].includes(location.protocol) || !LOOPBACK_HOSTS.has(location.hostname)) {
    return { disabledReason: "Der Live-Arbeitsplatz muss vom lokalen ONTSeq-Dienst geöffnet werden. Eine HTML-Datei oder externe Website kann keine Analyse starten." };
  }
  const token = document?.querySelector('meta[name="ontseq-session-token"]')?.content?.trim();
  // The server replaces its exact placeholder throughout the bundled HTML. Keep that
  // literal out of the client or it would also rewrite this guard to reject real tokens.
  if (!token || /^__ONTSEQ_[A-Z]+__$/.test(token)) {
    return { disabledReason: "Keine lokale Sitzung vorhanden. Diese Ansicht ist nicht mit dem Analyse-Dienst verbunden." };
  }
  return { token, disabledReason: null };
}

export function requestedProfile(location = globalThis.location) {
  try {
    const value = new URLSearchParams(location?.search ?? "").get("profile")?.trim() ?? "";
    return IDENTIFIER.test(value) ? value : "";
  } catch {
    return "";
  }
}

export function resolveAvailableProfile(profiles, current = "", requested = "") {
  const available = Array.isArray(profiles)
    ? profiles.filter((profile) => typeof profile === "string" && profile.length > 0)
    : [];
  // A Desktop handoff is an identity contract, not a preference. Never replace an
  // unavailable requested profile with a different build/profile without telling the user.
  if (requested && !available.includes(requested)) {
    return {
      profile: "",
      error: `Das vom Desktop angeforderte Analyseprofil ${requested} ist in dieser Dienstinstanz nicht verfügbar. Bitte das zugehörige Referenzpaket in „System einrichten“ installieren oder reparieren und den Live-Arbeitsplatz danach neu öffnen.`,
    };
  }
  if (available.includes(current)) return { profile: current, error: "" };
  if (available.includes(requested)) return { profile: requested, error: "" };
  return { profile: available[0] ?? "", error: "" };
}

export function workspaceHomeHref(location = globalThis.location) {
  const parameters = new URLSearchParams({ runtime: "local" });
  const profile = requestedProfile(location);
  if (profile) parameters.set("profile", profile);
  return `?${parameters.toString()}`;
}

function assertIdentifier(value) {
  if (!IDENTIFIER.test(value ?? "")) throw new LiveApiError("Ungültige Lauf- oder Proben-ID.");
  return value;
}

export function isTerminalJob(job) {
  return TERMINAL_STATES.has(job?.state);
}

export function validateMethylationProbe(value, bamPath) {
  if (!value || !["detected", "not_detected", "unknown"].includes(value.status)
    || typeof value.reason !== "string" || value.bam_path !== bamPath
    || !Number.isInteger(value.checked_reads) || value.checked_reads < 0
    || typeof value.complete !== "boolean"
    || (value.status !== "unknown" && (typeof value.bam_identity !== "string" || value.bam_identity_kind !== "stat-fingerprint-v1"))
    || (value.status === "not_detected" && !value.complete)) {
    throw new LiveApiError("Die Methylierungsprüfung hat keine eindeutig zur BAM passende Antwort geliefert.");
  }
  return value;
}

export function validateMethylationScan(value, bamPath, scanId = null) {
  if (!value || !/^[a-f0-9]{32}$/.test(value.scan_id ?? "")
    || (scanId && value.scan_id !== scanId) || value.bam_path !== bamPath
    || !["running", "completed", "cancelled", "failed"].includes(value.state)
    || !Number.isInteger(value.checked_reads) || value.checked_reads < 0
    || !Number.isFinite(value.elapsed_seconds) || value.elapsed_seconds < 0
    || (value.state === "completed" && !value.result)) {
    throw new LiveApiError("Die gründliche Prüfung hat keine eindeutig zur BAM passende Statusantwort geliefert.");
  }
  if (value.result) validateMethylationProbe(value.result, bamPath);
  return value;
}

export function methylationExplanation(probe) {
  const messages = {
    detected_5mc: "Informationen zur Cytosin-Methylierung (5mC) erkannt. Die Qualität wird erst in der Auswertung beurteilt.",
    complete_no_tags: "Die gesamte BAM wurde fehlerfrei gelesen. Es wurden keine Methylierungsinformationen gefunden.",
    sample_incomplete: "In der begrenzten Stichprobe wurden bisher keine passenden Informationen gefunden. Die Datei ist nicht vollständig geprüft.",
    timeout: "Das Zeitlimit der Vorprüfung wurde erreicht. Das Vorhandensein von Methylierungsinformationen bleibt offen.",
    cancelled: "Die Prüfung wurde abgebrochen. Es liegt keine vollständige Aussage zur Datei vor.",
    reader_unavailable: "Der direkte BAM-Leser ist in dieser Laufzeit nicht verfügbar. Bitte die ONTSeq-Laufzeit prüfen.",
    input_unavailable: "Die ausgewählte BAM konnte nicht geöffnet werden. Bitte Dateizugriff und Pfad prüfen.",
    file_changed: "Die BAM wurde während der Prüfung geändert. Bitte die abgeschlossene Datei erneut auswählen.",
    record_limit: "Ein sehr großer Datensatz hat die Speichergrenze der Vorprüfung erreicht. Das Ergebnis bleibt offen.",
    invalid_tags: "Modifikationsinformationen sind unvollständig oder widersprüchlich und konnten nicht sicher gelesen werden.",
    unsupported_modification: "Es wurden Basenmodifikationen gefunden, die diese 5mC-Auswertung nicht unterstützt.",
    read_error: "Die BAM konnte nicht fehlerfrei gelesen werden. Das ist kein Nachweis für fehlende Methylierungsinformationen.",
    worker_failed: "Die Vorprüfung konnte technisch nicht abgeschlossen werden. Bitte die Laufzeit prüfen.",
    busy: "Eine andere BAM-Prüfung läuft bereits. Bitte deren Abschluss abwarten.",
    transport_error: "Der lokale Dienst hat keine bestätigte Antwort geliefert. Bitte Verbindung und Dienststatus prüfen.",
  };
  return messages[probe?.reason_code] ?? probe?.reason ?? "Die Vorprüfung konnte keine sichere Aussage liefern.";
}

export function methylationReady(probe, bamPath, decision) {
  if (!bamPath || probe?.bam_path !== bamPath) return false;
  return probe.status === "not_detected" || (probe.status === "detected" && (decision === false || (decision === true && probe.methylation_available === true)))
    || (probe.status === "unknown" && decision === false);
}

export function sameMethylationEvidence(previous, current) {
  if (!previous || !current) return false;
  if (previous.status === "detected" && (!previous.bam_identity || !current.bam_identity)) return false;
  return previous.bam_path === current.bam_path && previous.status === current.status
    && previous.methylation_available === current.methylation_available
    && previous.expected_modkit_version === current.expected_modkit_version
    && previous.bam_identity_kind === current.bam_identity_kind
    && JSON.stringify(previous.bam_identity) === JSON.stringify(current.bam_identity);
}

export function selectedMethylationState(probe, choice, bamPath, selectionId) {
  return {
    probe: probe && probe.bam_path === bamPath && probe.selection_id === selectionId ? probe : null,
    decision: choice && choice.bam_path === bamPath && choice.selection_id === selectionId ? choice.value : null,
  };
}

export function validateMethylationReport(value, job) {
  if (!value || value.run_id !== job.run_id || value.sample_id !== job.sample_id
    || !job.detected_genome_build || value.genome_build !== job.detected_genome_build
    || !Array.isArray(value.regions) || typeof value.status !== "string") {
    throw new LiveApiError("Methylierungsergebnis zurückgewiesen: Lauf-, Proben- oder Build-Bindung fehlt oder ist widersprüchlich.");
  }
  return value;
}

export function validatePipelineResult(result, job) {
  const manifest = result?.manifest;
  if (!manifest || !result.qc || !result.provenance || !Array.isArray(result.events)) {
    throw new LiveApiError("Der Dienst hat kein lesbares PipelineResult geliefert.");
  }
  if (manifest.run_id !== job.run_id || manifest.sample_id !== job.sample_id) {
    throw new LiveApiError("Ergebnis zurückgewiesen: Lauf- oder Probenbindung stimmt nicht überein.");
  }
  const build = manifest.assay?.genome_build;
  if (!build || (job.detected_genome_build && build !== job.detected_genome_build)
    || (result.reference_context?.genome_build && result.reference_context.genome_build !== build)) {
    throw new LiveApiError("Ergebnis zurückgewiesen: widersprüchliche Referenzbuilds.");
  }
  if (!job.profile || !result.reference_context?.profile_id
    || result.reference_context.profile_id !== job.profile) {
    throw new LiveApiError("Ergebnis zurückgewiesen: das Analyseprofil stimmt nicht mit dem Lauf überein.");
  }
  return result;
}

export function createLiveApi({ location = globalThis.location, document = globalThis.document, fetchImpl = globalThis.fetch } = {}) {
  const session = getLiveSession(location, document);

  async function request(path, { method = "GET", payload, signal, timeout = 30000, blob = false } = {}) {
    if (session.disabledReason) throw new LiveApiError(session.disabledReason);
    // Do not turn this client into an arbitrary URL or credential-forwarding adapter.
    if (!path.startsWith("/api/") || path.includes("\\") || path.includes("#")) throw new LiveApiError("Ungültiger API-Pfad.");
    const controller = new AbortController();
    let timedOut = false;
    const abort = () => controller.abort();
    if (signal?.aborted) abort();
    else signal?.addEventListener("abort", abort, { once: true });
    const timer = setTimeout(() => { timedOut = true; abort(); }, timeout);
    try {
      const response = await fetchImpl(path, {
        method,
        headers: { [TOKEN_HEADER]: session.token, ...(payload ? { "Content-Type": "application/json" } : {}) },
        body: payload ? JSON.stringify(payload) : undefined,
        signal: controller.signal,
        credentials: "same-origin",
        mode: "same-origin",
        cache: "no-store",
        redirect: "error",
      });
      if (!response.ok) {
        const error = await response.json().catch(() => ({}));
        const message = String(error.error || `Anfrage fehlgeschlagen (HTTP ${response.status}).`).replaceAll(session.token, "[Sitzung]").slice(0, 1500);
        throw new LiveApiError(message, response.status);
      }
      if (blob) return { blob: await response.blob(), disposition: response.headers.get("Content-Disposition") };
      if (!response.headers.get("Content-Type")?.includes("application/json")) throw new LiveApiError("Kein Analyse-Dienst an dieser Adresse: JSON-Antwort fehlt.");
      return await response.json();
    } catch (error) {
      if (timedOut) throw new LiveApiError("Zeitüberschreitung. Der Dienst kann weiterarbeiten; ein angeforderter Lauf wird dadurch nicht abgebrochen.");
      throw error;
    } finally {
      clearTimeout(timer);
      signal?.removeEventListener("abort", abort);
    }
  }

  return {
    disabledReason: session.disabledReason,
    async config(signal) {
      const value = await request("/api/config", { signal });
      if (!Array.isArray(value.profiles) || !Array.isArray(value.roots) || typeof value.version !== "string") {
        throw new LiveApiError("Die Konfiguration des lokalen Dienstes ist unvollständig.");
      }
      return { ...value, profiles: value.profiles.filter((profile) => typeof profile === "string" && profile.length > 0) };
    },
    browse: (path = "", signal) => request(`/api/browse?${new URLSearchParams({ path })}`, { signal }),
    async probeMethylation(bam_path, signal, forceRefresh = false) {
      if (!bam_path) throw new LiveApiError("Bitte zuerst eine BAM-Datei auswählen.");
      return validateMethylationProbe(await request("/api/methylation/probe", {
        method: "POST", payload: { bam_path, ...(forceRefresh ? { force_refresh: true } : {}) }, signal, timeout: 30000,
      }), bam_path);
    },
    async startMethylationScan(bam_path, signal) {
      if (!bam_path) throw new LiveApiError("Bitte zuerst eine BAM-Datei auswählen.");
      return validateMethylationScan(await request("/api/methylation/scans", {
        method: "POST", payload: { bam_path }, signal, timeout: 30000,
      }), bam_path);
    },
    async methylationScan(scanId, bamPath, signal) {
      if (!/^[a-f0-9]{32}$/.test(scanId ?? "")) throw new LiveApiError("Ungültige Kennung der BAM-Prüfung.");
      return validateMethylationScan(await request(`/api/methylation/scans/${scanId}`, { signal }), bamPath, scanId);
    },
    async cancelMethylationScan(scanId, bamPath, signal) {
      if (!/^[a-f0-9]{32}$/.test(scanId ?? "")) throw new LiveApiError("Ungültige Kennung der BAM-Prüfung.");
      return validateMethylationScan(await request(`/api/methylation/scans/${scanId}/cancel`, {
        method: "POST", payload: {}, signal,
      }), bamPath, scanId);
    },
    async start({ bam, profile, sample_id, include_methylation = false }, signal) {
      if (!bam || !profile) throw new LiveApiError("Bitte BAM und installiertes Analyseprofil auswählen.");
      if (typeof include_methylation !== "boolean") throw new LiveApiError("Bitte den Methylierungsumfang ausdrücklich auswählen.");
      assertIdentifier(sample_id);
      // POST is deliberately never retried: a lost response can still have started a run.
      return request("/api/runs", { method: "POST", payload: { bam, profile, sample_id, include_methylation }, signal, timeout: 120000 });
    },
    async run(runId, signal) {
      const value = await request(`/api/runs/${encodeURIComponent(assertIdentifier(runId))}`, { signal });
      if (value.run_id !== runId) throw new LiveApiError("Die Statusantwort gehört zu einem anderen Lauf.");
      return value;
    },
    async result(job, signal) {
      const params = new URLSearchParams({ run_id: assertIdentifier(job.run_id), sample_id: assertIdentifier(job.sample_id) });
      return validatePipelineResult(await request(`/api/results?${params}`, { signal }), job);
    },
    async methylation(job, signal) {
      const params = new URLSearchParams({ run_id: assertIdentifier(job.run_id), sample_id: assertIdentifier(job.sample_id) });
      return validateMethylationReport(await request(`/api/methylation?${params}`, { signal }), job);
    },
    artifact(job, kind, signal) {
      if (!["html", "xlsx", "json"].includes(kind)) throw new LiveApiError("Unbekannter Exporttyp.");
      const params = new URLSearchParams({ run_id: assertIdentifier(job.run_id), sample_id: assertIdentifier(job.sample_id), kind });
      return request(`/api/artifacts?${params}`, { signal, blob: true, timeout: 120000 });
    },
  };
}

function pause(milliseconds, signal) {
  return new Promise((resolve, reject) => {
    const abort = () => { clearTimeout(timer); reject(new DOMException("Aborted", "AbortError")); };
    const timer = setTimeout(() => { signal?.removeEventListener("abort", abort); resolve(); }, milliseconds);
    if (signal?.aborted) abort();
    else signal?.addEventListener("abort", abort, { once: true });
  });
}

export async function watchRun(api, runId, { signal, onSnapshot, onError, interval = 2000, maxFailures = 5, maxPolls = 21600 } = {}) {
  let failures = 0;
  for (let count = 0; count < maxPolls && !signal?.aborted; count += 1) {
    try {
      const job = await api.run(runId, signal);
      if (signal?.aborted) return;
      failures = 0;
      onSnapshot(job);
      if (isTerminalJob(job)) return;
    } catch (error) {
      if (signal?.aborted || error.name === "AbortError") return;
      failures += 1;
      const stopped = failures >= maxFailures || [401, 403, 404].includes(error.status);
      onError(error, stopped);
      if (stopped) return;
    }
    try { await pause(Math.min(interval * 2 ** failures, 30000), signal); } catch { return; }
  }
  if (!signal?.aborted) onError(new LiveApiError("Die automatische Statusabfrage wurde nach ihrem Zeitlimit pausiert. Die Analyse wird nicht abgebrochen."), true);
}

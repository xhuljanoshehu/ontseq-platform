import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";
import { createLiveApi, getLiveSession, isTerminalJob, LiveApiError, methylationExplanation, methylationReady, requestedProfile, resolveAvailableProfile, sameMethylationEvidence, selectedMethylationState, TOKEN_HEADER, validateMethylationProbe, validateMethylationReport, validateMethylationScan, validatePipelineResult, watchRun, workspaceHomeHref } from "../src/liveApi.js";

const location = { protocol: "http:", hostname: "127.0.0.1" };
const document = { querySelector: () => ({ content: "unit-test-session-token" }) };
const job = { run_id: "run-001", sample_id: "sample-001", detected_genome_build: "GRCh37", profile: "AML_LCWGS_GRCh37" };
const json = (value, status = 200) => new Response(JSON.stringify(value), { status, headers: { "Content-Type": "application/json" } });
const result = () => ({ manifest: { run_id: job.run_id, sample_id: job.sample_id, assay: { genome_build: "GRCh37" } }, qc: { verdict: "WARN", metrics: {} }, provenance: { pipeline_version: "0.6.2" }, reference_context: { genome_build: "GRCh37", profile_id: job.profile }, events: [] });

test("file, foreign host, missing and placeholder sessions cannot access the live API", async () => {
  for (const context of [
    { location: { protocol: "file:", hostname: "" }, document },
    { location: { protocol: "https:", hostname: "example.com" }, document },
    { location, document: { querySelector: () => null } },
    { location, document: { querySelector: () => ({ content: "__ONTSEQ_TOKEN__" }) } },
  ]) {
    let calls = 0;
    const api = createLiveApi({ ...context, fetchImpl: () => { calls += 1; } });
    assert.ok(api.disabledReason);
    await assert.rejects(api.config(), LiveApiError);
    assert.equal(calls, 0);
  }
  assert.equal(getLiveSession(location, document).disabledReason, null);
});

test("requests stay same-origin and supply only the session header, never a URL credential", async () => {
  const calls = [];
  const api = createLiveApi({ location, document, fetchImpl: async (...args) => { calls.push(args); return json({ version: "0.6.2", roots: [], profiles: ["AML_LCWGS_GRCh37", "AML_LCWGS_GRCh38"] }); } });
  const config = await api.config();
  assert.equal(config.profiles.length, 2);
  const [path, options] = calls[0];
  assert.equal(path, "/api/config");
  assert.equal(options.headers[TOKEN_HEADER], "unit-test-session-token");
  assert.equal(options.mode, "same-origin");
  assert.equal(options.credentials, "same-origin");
  assert.equal(options.redirect, "error");
  assert.equal(options.cache, "no-store");
  assert.equal(path.includes("token"), false);
});

test("Desktop profile handoff is accepted only when the service advertises that exact profile", () => {
  const canonical25 = "AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25";
  const profiles = ["AML_LCWGS_GRCh38", canonical25];
  assert.equal(requestedProfile({ search: `?profile=${canonical25}` }), canonical25);
  assert.deepEqual(resolveAvailableProfile(profiles, "", canonical25), { profile: canonical25, error: "" });
  assert.deepEqual(resolveAvailableProfile(profiles, "AML_LCWGS_GRCh38", canonical25), { profile: "AML_LCWGS_GRCh38", error: "" });
  assert.equal(requestedProfile({ search: "?profile=../../not-allowed" }), "");
  assert.equal(requestedProfile({ search: "?profile=" }), "");
});

test("an unavailable Desktop-requested profile fails closed instead of selecting another build", () => {
  const unavailable = resolveAvailableProfile(
    ["AML_LCWGS_GRCh38"], "AML_LCWGS_GRCh38", "AML_LCWGS_GRCh37");
  assert.equal(unavailable.profile, "");
  assert.match(unavailable.error, /AML_LCWGS_GRCh37/);
  assert.match(unavailable.error, /nicht verfügbar/);
  assert.deepEqual(
    resolveAvailableProfile(["AML_LCWGS_GRCh38"], "", ""),
    { profile: "AML_LCWGS_GRCh38", error: "" },
  );
});

test("workspace brand navigation preserves the Desktop profile handoff", () => {
  const canonical25 = "AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25";
  assert.equal(
    workspaceHomeHref({ search: `?runtime=local&profile=${canonical25}` }),
    `?runtime=local&profile=${canonical25}`,
  );
  assert.equal(workspaceHomeHref({ search: "" }), "?runtime=local");
});

test("start uses installed profile + BAM contract and never retries an uncertain POST", async () => {
  const calls = [];
  const api = createLiveApi({ location, document, fetchImpl: async (path, options) => { calls.push([path, options]); throw new TypeError("connection interrupted"); } });
  await assert.rejects(api.start({ bam: "/approved/sample.bam", profile: job.profile, sample_id: job.sample_id }), /connection interrupted/);
  assert.equal(calls.length, 1);
  assert.equal(calls[0][0], "/api/runs");
  assert.deepEqual(JSON.parse(calls[0][1].body), { bam: "/approved/sample.bam", profile: job.profile, sample_id: job.sample_id, include_methylation: false });
  assert.equal(calls[0][1].method, "POST");
});

const probe = (status = "detected", bam_path = "/approved/sample.bam") => ({ status, bam_path, reason: "Synthetic test probe", checked_reads: 10, complete: status === "not_detected", methylation_available: true, bam_identity: "a".repeat(64), bam_identity_kind: "stat-fingerprint-v1" });

test("a detected BAM requires explicit opt-in or opt-out and a new BAM cannot reuse the decision", () => {
  assert.equal(methylationReady(probe(), probe().bam_path, null), false);
  assert.equal(methylationReady(probe(), probe().bam_path, true), true);
  assert.equal(methylationReady(probe(), probe().bam_path, false), true);
  assert.equal(methylationReady({ ...probe(), methylation_available: false }, probe().bam_path, true), false);
  assert.equal(methylationReady({ ...probe(), methylation_available: false }, probe().bam_path, false), true);
  assert.equal(methylationReady(probe(), "/approved/other.bam", true), false);
  assert.equal(methylationReady(null, probe().bam_path, true), false);
  assert.equal(methylationReady(probe("unknown"), probe().bam_path, null), false);
  assert.equal(methylationReady(probe("unknown"), probe().bam_path, true), false);
  assert.equal(methylationReady(probe("unknown"), probe().bam_path, false), true);
  assert.equal(methylationReady(probe("not_detected"), probe().bam_path, null), true);
});

test("probe responses must belong to the selected BAM and partial negatives remain unresolved", async () => {
  const calls = [];
  const api = createLiveApi({ location, document, fetchImpl: async (path, options) => { calls.push([path, options]); return json(probe()); } });
  assert.deepEqual(await api.probeMethylation(probe().bam_path), probe());
  assert.equal(calls[0][0], "/api/methylation/probe");
  assert.deepEqual(JSON.parse(calls[0][1].body), { bam_path: probe().bam_path });
  assert.equal(calls[0][1].headers[TOKEN_HEADER], "unit-test-session-token");
  for (const value of [{ ...probe(), bam_path: "/other.bam" }, { ...probe(), status: "not_detected", complete: false }, { ...probe(), checked_reads: -1 }, { ...probe(), status: "guess" }]) {
    assert.throws(() => validateMethylationProbe(value, probe().bam_path), /BAM passende/);
  }
});

test("methylation opt-in is serialized explicitly and truthy strings never enable it", async () => {
  const calls = [];
  const api = createLiveApi({ location, document, fetchImpl: async (path, options) => { calls.push([path, options]); return json(job); } });
  await api.start({ bam: probe().bam_path, profile: job.profile, sample_id: job.sample_id, include_methylation: true });
  assert.equal(JSON.parse(calls[0][1].body).include_methylation, true);
  await assert.rejects(api.start({ bam: probe().bam_path, profile: job.profile, sample_id: job.sample_id, include_methylation: "false" }), /ausdrücklich/);
  assert.equal(calls.length, 1);
});

test("thorough scan states are bound to their BAM and server-issued identifier", () => {
  const snapshot = { scan_id: "a".repeat(32), bam_path: probe().bam_path, state: "running", checked_reads: 42, elapsed_seconds: 1.5, result: null };
  assert.equal(validateMethylationScan(snapshot, probe().bam_path, snapshot.scan_id), snapshot);
  for (const invalid of [
    { ...snapshot, scan_id: "../other" }, { ...snapshot, bam_path: "/other.bam" },
    { ...snapshot, state: "completed" }, { ...snapshot, checked_reads: -1 },
    { ...snapshot, elapsed_seconds: -1 }, { ...snapshot, result: { ...probe(), bam_path: "/other.bam" } },
  ]) assert.throws(() => validateMethylationScan(invalid, probe().bam_path, snapshot.scan_id), /BAM passende/);
  assert.throws(() => validateMethylationScan(snapshot, probe().bam_path, "b".repeat(32)), /BAM passende/);
});

test("scan start, polling and cancellation use the authenticated local job contract", async () => {
  const calls = [];
  const snapshot = { scan_id: "a".repeat(32), bam_path: probe().bam_path, state: "running", checked_reads: 0, elapsed_seconds: 0, result: null };
  const api = createLiveApi({ location, document, fetchImpl: async (path, options) => { calls.push([path, options]); return json(snapshot); } });
  await api.startMethylationScan(snapshot.bam_path);
  await api.methylationScan(snapshot.scan_id, snapshot.bam_path);
  await api.cancelMethylationScan(snapshot.scan_id, snapshot.bam_path);
  assert.deepEqual(calls.map(([path]) => path), ["/api/methylation/scans", `/api/methylation/scans/${snapshot.scan_id}`, `/api/methylation/scans/${snapshot.scan_id}/cancel`]);
  assert.deepEqual(JSON.parse(calls[0][1].body), { bam_path: snapshot.bam_path });
  assert.deepEqual(JSON.parse(calls[2][1].body), {});
  for (const [, options] of calls) assert.equal(options.headers[TOKEN_HEADER], "unit-test-session-token");
  await assert.rejects(api.cancelMethylationScan("../other", snapshot.bam_path), /Ungültige Kennung/);
  assert.equal(calls.length, 3);
});

test("pre-start methylation requests bypass display cache explicitly", async () => {
  const calls = [];
  const api = createLiveApi({ location, document, fetchImpl: async (path, options) => { calls.push(JSON.parse(options.body)); return json(probe()); } });
  await api.probeMethylation(probe().bam_path, undefined, true);
  assert.equal(calls[0].force_refresh, true);
});

test("incomplete discovery and unsupported tags are explained without asserting absence", () => {
  assert.match(methylationExplanation({ reason_code: "sample_incomplete" }), /nicht vollständig/);
  assert.match(methylationExplanation({ reason_code: "unsupported_modification" }), /nicht unterstützt/);
  assert.match(methylationExplanation({ reason_code: "reader_unavailable" }), /BAM-Leser/);
  assert.match(methylationExplanation({ reason_code: "file_changed" }), /geändert/);
  assert.match(methylationExplanation({ reason_code: "complete_no_tags" }), /gesamte BAM.*fehlerfrei/);
});

test("an empty initial BAM selection renders without dereferencing absent probe or choice", () => {
  assert.deepEqual(selectedMethylationState(null, null, undefined, undefined), { probe: null, decision: null });
  assert.deepEqual(selectedMethylationState(null, null, null, null), { probe: null, decision: null });
  assert.equal(methylationReady(null, undefined, null), false);
});

test("reselecting the same BAM invalidates both prior detection and the prior decision", () => {
  const oldProbe = { ...probe(), selection_id: 1 };
  const choice = { bam_path: probe().bam_path, selection_id: 1, value: true };
  assert.deepEqual(selectedMethylationState(oldProbe, choice, probe().bam_path, 2), { probe: null, decision: null });
  const currentProbe = { ...probe(), selection_id: 2 };
  const current = selectedMethylationState(currentProbe, choice, probe().bam_path, 2);
  assert.equal(current.probe, currentProbe);
  assert.equal(current.decision, null);
  assert.equal(methylationReady(current.probe, probe().bam_path, current.decision), false);
  assert.equal(selectedMethylationState(currentProbe, null, probe().bam_path, 2).decision, null);
});

test("fresh preflight invalidates consent when same-path BAM metadata or tool availability changes", () => {
  const initial = { ...probe(), bam_identity: "synthetic-metadata-fingerprint-v1", expected_modkit_version: "0.4.1" };
  assert.equal(sameMethylationEvidence(initial, { ...initial }), true);
  assert.equal(sameMethylationEvidence(initial, { ...initial, bam_identity: "synthetic-metadata-fingerprint-v2" }), false);
  assert.equal(sameMethylationEvidence(initial, { ...initial, status: "unknown" }), false);
  assert.equal(sameMethylationEvidence(initial, { ...initial, methylation_available: false }), false);
  assert.equal(sameMethylationEvidence(initial, { ...initial, bam_identity: null }), false);
  assert.equal(sameMethylationEvidence(null, initial), false);
});

test("regional methylation responses bind run, sample and reference independently", async () => {
  const report = { run_id: job.run_id, sample_id: job.sample_id, genome_build: "GRCh37", status: "COMPLETED", regions: [] };
  const paths = [];
  const api = createLiveApi({ location, document, fetchImpl: async (path) => { paths.push(path); return json(report); } });
  assert.deepEqual(await api.methylation(job), report);
  assert.equal(paths[0], "/api/methylation?run_id=run-001&sample_id=sample-001");
  for (const value of [{ ...report, sample_id: "other" }, { ...report, run_id: "other" }, { ...report, genome_build: "GRCh38" }]) {
    assert.throws(() => validateMethylationReport(value, job), /Bindung/);
  }
  const noCall = { ...report, status: "NO_CALL", regions: [{ region_id: "chr1", mean_modified_fraction: null, sites_total: 5, sites_at_minimum_coverage: 0 }] };
  assert.equal(validateMethylationReport(noCall, job), noCall);
});

test("invalid identifiers are refused before a request can name another run", async () => {
  let calls = 0;
  const api = createLiveApi({ location, document, fetchImpl: () => { calls += 1; } });
  await assert.rejects(api.run("../run"), /Ungültige/);
  await assert.rejects(api.start({ bam: "/sample.bam", profile: job.profile, sample_id: "../outside" }), /Ungültige/);
  assert.throws(() => api.artifact({ ...job, run_id: "other/run" }, "json"), /Ungültige/);
  assert.throws(() => api.artifact(job, "../../secret"), /Unbekannter/);
  assert.equal(calls, 0);
});

test("result identity and both reference bindings fail closed", () => {
  const original = result();
  assert.equal(validatePipelineResult(original, job), original);
  assert.throws(() => validatePipelineResult({ ...result(), manifest: { ...result().manifest, sample_id: "different" } }, job), /Probenbindung/);
  assert.throws(() => validatePipelineResult({ ...result(), manifest: { ...result().manifest, assay: { genome_build: "GRCh38" } } }, job), /Referenzbuilds/);
  assert.throws(() => validatePipelineResult({ ...result(), reference_context: { genome_build: "GRCh38" } }, job), /Referenzbuilds/);
  assert.throws(() => validatePipelineResult({ ...result(), reference_context: { genome_build: "GRCh37", profile_id: "other-profile" } }, job), /Analyseprofil/);
  assert.throws(() => validatePipelineResult({}, job), /PipelineResult/);
  assert.equal(original.events.length, 0);
});

test("result validation rejects a missing profile binding", () => {
  assert.throws(() => validatePipelineResult({ ...result(), reference_context: { genome_build: "GRCh37" } }, job), /Analyseprofil/);
  assert.throws(() => validatePipelineResult({ ...result(), reference_context: undefined }, job), /Analyseprofil/);
});

test("results and exports are fetched with explicit run + sample, not arbitrary paths", async () => {
  const paths = [];
  const api = createLiveApi({ location, document, fetchImpl: async (path) => { paths.push(path); return path.includes("artifacts") ? new Response("original backend report") : json(result()); } });
  const loaded = await api.result(job);
  assert.deepEqual(loaded, result());
  const exported = await api.artifact(job, "html");
  assert.equal(await exported.blob.text(), "original backend report");
  assert.equal(paths[0], "/api/results?run_id=run-001&sample_id=sample-001");
  assert.equal(paths[1], "/api/artifacts?run_id=run-001&sample_id=sample-001&kind=html");
});

test("wrong response identity, wrong content type and server error remain explicit failures", async () => {
  const wrongRun = createLiveApi({ location, document, fetchImpl: async () => json({ run_id: "different" }) });
  await assert.rejects(wrongRun.run(job.run_id), /anderen Lauf/);
  const htmlServer = createLiveApi({ location, document, fetchImpl: async () => new Response("<!doctype html>", { headers: { "Content-Type": "text/html" } }) });
  await assert.rejects(htmlServer.config(), /JSON-Antwort fehlt/);
  const denied = createLiveApi({ location, document, fetchImpl: async () => json({ error: "missing unit-test-session-token" }, 401) });
  await assert.rejects(denied.config(), (error) => error.status === 401 && !error.message.includes("unit-test-session-token"));
});

test("watch preserves all four stage states and stops on terminal result", async () => {
  let calls = 0;
  const snapshots = [];
  const stages = ["COMPLETED", "NO_CALL", "FAILED", "NOT_RUN"].map((status) => ({ status }));
  await watchRun({ run: async () => { calls += 1; return { ...job, state: calls === 1 ? "running" : "failed", stages }; } }, job.run_id, { onSnapshot: (value) => snapshots.push(value), onError: assert.fail, interval: 0 });
  assert.equal(calls, 2);
  assert.deepEqual(snapshots[1].stages.map((stage) => stage.status), ["COMPLETED", "NO_CALL", "FAILED", "NOT_RUN"]);
  assert.equal(isTerminalJob(snapshots[1]), true);
  assert.equal(isTerminalJob(snapshots[0]), false);
});

test("polling has bounded error retries and immediate authentication stop", async () => {
  let calls = 0;
  const errors = [];
  await watchRun({ run: async () => { calls += 1; throw new LiveApiError("offline"); } }, job.run_id, { onSnapshot: assert.fail, onError: (_, stopped) => errors.push(stopped), interval: 0, maxFailures: 3 });
  assert.equal(calls, 3);
  assert.deepEqual(errors, [false, false, true]);
  calls = 0;
  await watchRun({ run: async () => { calls += 1; throw new LiveApiError("session expired", 401); } }, job.run_id, { onSnapshot: assert.fail, onError: (_, stopped) => assert.equal(stopped, true), interval: 0 });
  assert.equal(calls, 1);
});

test("aborting or unmounting stops polling without claiming to cancel the pipeline", async () => {
  const controller = new AbortController();
  let calls = 0;
  await watchRun({ run: async () => { calls += 1; return { ...job, state: "running" }; } }, job.run_id, { signal: controller.signal, onSnapshot: () => controller.abort(), onError: assert.fail, interval: 0 });
  assert.equal(calls, 1);
  assert.equal(controller.signal.aborted, true);
});

test("a perpetually running service has a finite polling budget without changing the run", async () => {
  let calls = 0;
  const errors = [];
  await watchRun({ run: async () => { calls += 1; return { ...job, state: "running" }; } }, job.run_id, {
    onSnapshot: (snapshot) => assert.equal(snapshot.state, "running"),
    onError: (error, stopped) => errors.push([error.message, stopped]), interval: 0, maxPolls: 2,
  });
  assert.equal(calls, 2);
  assert.equal(errors.length, 1);
  assert.equal(errors[0][1], true);
  assert.match(errors[0][0], /nicht abgebrochen/);
});

test("live UI is isolated from demo evidence, mutable builds, persisted tokens and review mutations", async () => {
  const ui = await readFile(new URL("../src/LiveWorkspace.jsx", import.meta.url), "utf8");
  const client = await readFile(new URL("../src/liveApi.js", import.meta.url), "utf8");
  const main = await readFile(new URL("../src/main.jsx", import.meta.url), "utf8");
  assert.doesNotMatch(ui, /from ["']\.\/data|from ["']\.\/reportState|from ["']\.\/hematologyKnowledge|cnvDemoData|evidenceRows|runRecord/);
  assert.doesNotMatch(`${ui}\n${client}`, /localStorage|sessionStorage|\/api\/review|dangerouslySetInnerHTML/);
  // Server-side token injection is a global text replacement in a single-file build.
  assert.equal(client.includes("__ONTSEQ_TOKEN__"), false);
  assert.doesNotMatch(ui, /<iframe|window\.open|target=["']_blank/);
  assert.match(ui, /if \(!result \|\| loadingResult \|\| resultError\) return/);
  assert.match(ui, /disabled=\{!result \|\| loadingResult \|\| Boolean\(resultError\)/);
  assert.match(ui, /const disableInputs = !config \|\| connecting \|\| running/);
  assert.match(main, /from "\.\/LiveWorkspace\.jsx"/);
  assert.doesNotMatch(main, /App\.jsx|data\.js/);
  assert.match(ui, /Keine elektronische Signatur/);
  assert.match(ui, /Keine automatische Build-Konvertierung/);
  assert.match(ui, /Nur für Forschungszwecke/);
  assert.match(ui, /Panel-Koordinatenauflösung/);
  assert.match(ui, /gemappt bedeutet nicht analytisch validiert/);
  assert.match(ui, /keine negativen Befunde dieser Probe/);
});

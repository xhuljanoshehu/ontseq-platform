# ONTSeq Desktop ↔ local backend contract (v0.1)

Status: **engineering contract / research only**. This document fixes what the Windows shell may rely on without making any clinical claim.

## Boundary

The desktop application is an operator surface. It must not implement analytical rules or reinterpret pipeline outcomes. The authoritative records remain the backend result contracts and the run envelope.

Transport is HTTP on loopback only (`127.0.0.1`, default port `8765`). Every API request requires the per-process `X-ONTSeq-Token`. The service must not be exposed on a LAN.

## Bootstrap

The desktop starts the backend in WSL and polls `/` until the service is ready. The served page contains the same short-lived token used by the existing browser surface. After obtaining it, the desktop calls `GET /api/config` and refuses to continue if startup failed.

This bootstrap mechanism is adequate for the local engineering prototype because it has the same security boundary as the existing same-origin browser UI. A future installer/service architecture may replace it with an authenticated local IPC bootstrap without changing analytical contracts.

## Start

`POST /api/runs`

Request fields used by the GRCh38 profile-backed Desktop:

```json
{
  "bam": "P:\\Lab\\sample.bam",
  "sample_id": "SAMPLE_001",
  "profile": "AML_LCWGS_GRCh38",
  "genome_build": "GRCh38",
  "assay": "lcwgs",
  "target_bed": null,
  "target_bed_version": null
}
```

`profile` is authoritative for resource and dictionary-contract resolution. The published values
are `AML_LCWGS_GRCh38`, `AML_AS_111_GRCh38`, `AML_LCWGS_GRCh38_CANONICAL25` and
`AML_AS_111_GRCh38_CANONICAL25`; all are build-isolated to GRCh38. The two unsuffixed profiles
retain the `exact_full` Primary-Assembly contract. The two Canonical-25 profiles require exactly
`chr1`-`chr22`, `chrX`, `chrY`, `chrM`. They resolve the same installed GRCh38 bundles, with no
additional reference download, liftover or fallback. The
`genome_build`, `assay`, `target_bed` and `target_bed_version` fields remain in the transport
contract for one compatibility version. New profile runs derive their build and assay from the
profile and leave the two explicit target-BED fields null.
Desktop also omits `run_id`; Core derives the canonical `<sample>-<UTC timestamp>` value and
returns it in the HTTP 202 job. Desktop uses that returned ID for polling and output links.

The service is started with `--resource-root`, taken from `resourceRootWsl` (default
`~/.local/share/ontseq/resources`, below the WSL user's `$HOME`). Core/CLI keeps its separate
`/opt/ontseq` default. The service must resolve the profile's pinned bundle IDs and must not fall back to an
explicit path or another build.

Expected response: HTTP 202 with a run job whose initial `state` is `running`.

`GET /api/config` advertises only profiles whose complete pinned Reference, Knowledge and optional
Panel context passes the fast local manifest/presence/declared-size check. Starting a profile run
uses the same bounded preflight and does not hash every multi-gigabyte resource. Full SHA256
verification remains an explicit `ontseq references validate` operation.

## Status

`GET /api/runs/<run_id>`

Job state is a transport/orchestration state. Stage status is an analytical-execution state and must not be collapsed into job state.

The current backend stage vocabulary is preserved verbatim:

- `COMPLETED` — the stage ran to a defined conclusion and produced its expected contract.
- `NO_CALL` — the stage ran but could not make a call under its rules. **Not a negative biological result.**
- `FAILED` — the stage was attempted and failed.
- `NOT_RUN` — the stage was not attempted, including non-applicability or blocking dependencies. **Not a negative biological result.**

The desktop also polls `<output>/<run>/<sample>/provenance/run.json`. That file is atomically rewritten by the runner after each stage, so it provides live stage progress before the in-memory HTTP job receives its final stage list. The desktop displays those values; it does not derive new scientific statuses.

When the job response contains `detected_genome_build`, Desktop displays it. Independently,
as soon as `provenance/run.json` exists, Desktop reads its top-level `genome_build` and labels it
as provenance-backed. Before that point the UI shows the unambiguous GRCh38 profile expectation,
not a claim that the BAM has already passed dictionary validation.

## Results

The desktop opens files from the run envelope rather than copying them into a second result store:

```text
<output>/<run>/<sample>/
  provenance/run.json
  normalized/<sample>.result.json
  reports/<sample>.report.html
  reports/<sample>.results.xlsx
  release/release.json
  release/checksums.sha256
```

A missing report file is shown as unavailable, never manufactured from partial state.

## Cancellation gate

Desktop v0.1 does **not** send a cancellation request because the backend currently executes many tools through blocking subprocess calls. A UI-only cancel would be dishonest: it could stop polling while the analysis continued, or kill a process without giving the runner a chance to record what happened.

The button may be enabled only after all of the following exist and are tested:

1. an explicit `POST /api/runs/<run_id>/cancel` contract;
2. backend job states that distinguish `cancel_requested` and `cancelled` from failure;
3. cancellation propagation into the command-runner/orchestrator boundary;
4. atomic cleanup or quarantine of any tool output interrupted mid-write;
5. every downstream stage recorded as `NOT_RUN` with a cancellation reason;
6. resume behaviour proven after cancellation and process restart;
7. real-tool failure-injection CI.

## Filesystem rule

The desktop passes the selected BAM's parent directory as the service `--allow-root`. Windows drive-letter paths are mapped to `/mnt/<drive>/...`. This is only valid when that drive is visible inside WSL. UNC paths are not silently invented; they require a deliberate mount.

## Versioning

Breaking changes to request/response fields or stage semantics require a new documented desktop API version or a compatibility layer. UI text may be translated, but machine values must remain stable and traceable to the run report.

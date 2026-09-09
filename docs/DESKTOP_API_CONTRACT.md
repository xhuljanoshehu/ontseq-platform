# ONTSeq Desktop ↔ local backend contract (v0.1)

Status: **engineering contract / research only**. This document fixes what the Windows shell may rely on without making any clinical claim.

## Boundary

The desktop application is an operator surface. It must not implement analytical rules or reinterpret pipeline outcomes. The authoritative records remain the backend result contracts and the run envelope.

Transport is HTTP on loopback only (`127.0.0.1`, default port `8765`). Every API request requires the per-process `X-ONTSeq-Token`. The service must not be exposed on a LAN. The configured port is a preference: if another local service already owns it, Desktop selects a free ephemeral loopback port for the new process rather than adopting or stopping the existing service.

## Bootstrap

The desktop starts the backend in WSL with a fresh 32-character lowercase hexadecimal
`--instance-id` and polls `/` until the service is ready. The served page contains the same
short-lived token used by the existing browser surface. After obtaining it, the desktop calls
`GET /api/config` and commits the client only if `instance_id`, `resource_root`, `output_dir`, the
single allowed input root and the exact requested profile all match the launch expectation. A
foreign listener or simultaneous bind winner is rejected rather than adopted; Desktop retries a
detected collision on a fresh port and does not retain the failed candidate. The selected profile
may be carried to `/workspace` as a non-secret `profile` query parameter; the workspace accepts it
only when it occurs verbatim in `config.profiles`. Session credentials are never put in the URL.

This bootstrap mechanism is adequate for the local engineering prototype because it has the same security boundary as the existing same-origin browser UI. A future installer/service architecture may replace it with an authenticated local IPC bootstrap without changing analytical contracts.

## Start

### Optional methylation in 0.7.1

`POST /api/methylation/probe` accepts `{"bam_path": "<selected local BAM>", "force_refresh": false}`.
It uses the same allowed-root boundary as intake. The response contains the canonical
`bam_path`, `status` (`detected`, `not_detected`, `unknown`), `reason`, `checked_reads`,
`complete`, and a versioned v2 `probe_version`. Additional fields are `reason_code`,
`elapsed_seconds`, `scan_mode`, `reader` and `detected_modifications`. Only successful EOF
without candidate tags yields `not_detected`; a bounded sample, unsupported modification
or incomplete tag pair stays `unknown`. Display results may be cached briefly; clients set
`force_refresh=true` immediately before starting a run. A positive cached BAM-record offset
is an internal witness: that record is freshly read and checked, never exposed in the API.
No alignments, read names or sequences are returned or retained.

`POST /api/methylation/scans` with `{"bam_path": "<selected local BAM>"}` starts an optional
thorough scan and returns HTTP 202. `GET /api/methylation/scans/<scan_id>` polls progress;
`POST /api/methylation/scans/<scan_id>/cancel` with `{}` requests cancellation. All three
use the same authentication and path boundaries. A snapshot binds `scan_id` (32 lowercase
hex characters), `bam_path`, `state` (`running`, `completed`, `cancelled`, `failed`),
`checked_reads`, `elapsed_seconds` and `result` (null while running, otherwise a probe
response when available). Only one discovery operation runs at a time; a conflicting scan
start returns HTTP 409. Cancelled/failed scans never become negative results. Clients cancel
owned scans and invalidate pending responses when the BAM/profile selection changes.

`bam_identity` is a SHA256 fingerprint of device, inode, size and nanosecond mtime;
`bam_identity_kind=stat-fingerprint-v1` explicitly distinguishes it from a content hash.
The file is stat-checked before and after probing, including positive detection.
Changes invalidate the probe. This convenience guard does not replace intake's full
input checksum or eliminate every possible concurrent filesystem change.

Both this response and `/api/config` advertise `methylation_available`,
`methylation_unavailable_reason` and `expected_modkit_version`. The service validates
the exact tool version and policy before accepting a requested additional analysis.
Desktop can name a separately installed executable through optional
`modkitExecutableWsl`; it is passed as one `--modkit` argument, never shell code.

Detected annotations require an explicit yes/no choice, initially unset. Unknown status
requires explicit confirmation to continue without methylation. Each file selection,
including the same path, resets the choice. The workspace re-probes before submission
and rejects changed evidence; subsequent runs require a fresh choice. Desktop probes
again immediately before each start and presents Yes/No/Cancel. Tool absence disables
inclusion while leaving an explicitly selected genome-only run available.

`POST /api/runs`

Request fields used by the profile-backed Desktop (GRCh38 example):

```json
{
  "bam": "P:\\Lab\\sample.bam",
  "sample_id": "SAMPLE_001",
  "profile": "AML_LCWGS_GRCh38",
  "genome_build": "GRCh38",
  "assay": "lcwgs",
  "include_methylation": false,
  "target_bed": null,
  "target_bed_version": null
}
```

`profile` is authoritative for resource and dictionary-contract resolution. Published values are
the four GRCh38 profiles `AML_LCWGS_GRCh38`, `AML_AS_111_GRCh38`,
`AML_LCWGS_GRCh38_CANONICAL25`, `AML_AS_111_GRCh38_CANONICAL25` and the four GRCh37 profiles
`AML_LCWGS_GRCh37`, `AML_AS_111_GRCh37`, `AML_LCWGS_GRCh37_UCSC_HG19_CANONICAL25` and
`AML_AS_111_GRCh37_UCSC_HG19_CANONICAL25`. Every profile is build-isolated.
Unsuffixed profiles retain their complete exact-full contract. Canonical-25 profiles require
exactly `chr1`-`chr22`, `chrX`, `chrY`, `chrM` with the selected build's pinned lengths; the
GRCh37 UCSC-hg19 contract fixes `chrM=16571`. There is no additional reference download,
run-time liftover, cross-build or dictionary fallback. Both GRCh37 Adaptive-Sampling profiles
pin only `AML_AS_111_GRCh37_v1`; the resource records a controlled build-time projection of
110/111 source intervals, 109 reciprocal-exact intervals and 108 native GENCODE-19 analysis
ROIs. `ACACA` remains roundtrip-review-required. It never invokes coordinate
conversion while processing a BAM. The
`genome_build`, `assay`, `target_bed` and `target_bed_version` fields remain in the transport
contract for one compatibility version. New profile runs derive their build and assay from the
profile and leave the two explicit target-BED fields null.
Desktop also omits `run_id`; Core derives the canonical `<sample>-<UTC timestamp>` value and
returns it in the HTTP 202 job. Desktop uses that returned ID for polling and output links.

The service is started with `--resource-root`, taken from `resourceRootWsl` (default
`~/.local/share/ontseq/resources-v0.7.1`, below the WSL user's `$HOME`). The default follows
the Desktop/Core release; explicitly saved resource paths remain selected. Core/CLI keeps its separate
`/opt/ontseq` default. The service must resolve the profile's pinned bundle IDs and must not fall back to an
explicit path or another build.

Expected response: HTTP 202 with a run job whose initial `state` is `running`.

`include_methylation` is a strict JSON boolean, defaults to false and changes the
manifest's requested modules. String truthiness never enables it. lcWGS uses the
versioned chromosome-summary policy; Adaptive Sampling selects the target-BED policy
and the profile's locked analysis ROIs. No analytical threshold is chosen by the UI.

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
as provenance-backed. Before that point the UI shows the unambiguous selected-profile build
expectation, not a claim that the BAM has already passed dictionary validation.

## Results

The workspace retrieves regional methylation through
`GET /api/methylation?run_id=<run>&sample_id=<sample>`. The response is the typed
`MethylationReport` plus `run_id`, for `COMPLETED` and `NO_CALL` outcomes. It is bound
to the current run envelope, completed assembly and methylation stage checksums,
sample, reference build, tool identity and bedMethyl provenance. Failed/resumed runs
cannot mix an older summary with newly generated regional values. `null` fractions
remain unavailable values; they are never converted to zero.

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

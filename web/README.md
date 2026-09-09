# Local live workspace

This is the source of the real service workspace at `/workspace`. It does not import the
separate UX demo, synthetic findings, chart mockups, or external network assets. The existing
service page at `/` is unchanged.

## Reproducible build

Prerequisites: Node.js 20 or newer and npm. From this directory:

```sh
npm ci
npm test
npm run build
npm run check:build
```

`npm ci` installs the integrity-pinned dependency graph in `package-lock.json`. React,
React DOM and esbuild have exact direct versions. Subsequent builds and normal workspace
use do not require network access. `scripts/build_live_workspace.mjs` checks installed
versions, the input module graph, product-version consistency, and absence of external assets.
It creates the deterministic, self-contained
`src/ontseq_platform/service/workspace.html`, already included by Python's `*.html` package
data rule. `check:build` rebuilds in memory and fails if the checked-in HTML differs.

The HTML contains the runtime libraries' license notices and exactly one session-token
placeholder. It contains no actual session credential. When serving `/workspace`, the local
service replaces that marker with its fresh per-process token. Do not publish the authenticated
page, expose this service on the network, or use a public hosting provider for this workflow.

### Offline bootstrap with existing dependencies

For a controlled initial build without another download, set `ONTSEQ_WEB_NODE_MODULES` to an
existing `node_modules` directory containing the exact lockfile versions, then run:

```sh
node scripts/build_live_workspace.mjs
node scripts/build_live_workspace.mjs --check
```

Those commands are run from the repository root. The override changes only package resolution;
it cannot bypass the version checks. Normal repeatable development should use `npm ci` under
`web/` and omit the override.

## Runtime boundaries

- Same-origin loopback API only; the token stays in memory, never local/session storage or URLs.
- A static `file://` page or page without a service-issued token remains disabled.
- BAMs must come from server-authorized directories; only locally available profiles are offered.
- Starting a run sends `bam`, `sample_id` and `profile`; POST requests are never automatically retried.
- Polling is bounded, stops on terminal results or repeated failures, and aborts on unmount.
  Stopping polling does not cancel the pipeline.
- Results must match the original run, sample, profile and reference build. A new profile selection
  never changes a prior result. Missing results and an actual empty event list remain distinct.
- Artifact downloads are enabled only after result validation. HTML is downloaded, not executed
  in an iframe or injected into the live page. The original backend export is not reformatted or signed.
- No review mutations, electronic signature, clinical release or invented execution/provenance fields.

The Node tests verify these client contracts. Rendered-browser checks, real reference provisioning,
backend boundary tests and analytical/clinical validation are separate gates; a passing web build
does not claim any of them.

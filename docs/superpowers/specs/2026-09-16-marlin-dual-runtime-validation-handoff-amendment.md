# MARLIN Dual-Runtime Validation Handoff Amendment

## Status

Corrective amendment to `2026-09-15-marlin-dual-runtime-compatibility-design.md`.

## Reason

The dual-runtime comparison qualifies a live candidate runtime against a frozen reference runtime, but ONTSeq's existing `marlin-validate` safety contract intentionally requires the supplied `MarlinRuntimeCompatibilityProfile.reference_runtime_lock_id` to equal the artifact lock used for biological classification.

A reference profile is therefore not directly substituted for the candidate's own validation profile. Doing so would weaken or bypass the existing same-lock fail-closed boundary.

## Required handoff

The controlled sequence is:

```text
reference live probe
-> reference execution lock
-> reference fixture/profile/runtime identity freeze
-> candidate live probe
-> candidate execution lock over identical artifact bytes
-> marlin-compare-runtimes
-> require dual-runtime verdict PASS
-> candidate marlin-freeze-runtime under the candidate execution lock
-> candidate fixture/profile self-baseline
-> marlin-validate with candidate lock + candidate profile + candidate fixture
-> AL_001, then GSE280090
```

The candidate freeze after a dual-runtime PASS is not a second biological qualification. It creates the same-lock reproducibility profile required by the existing validation harness. The deterministic fixture generator and predeclared engineering tolerances remain unchanged.

## Fail-closed rules

- A dual-runtime FAIL blocks biological validation.
- `marlin-validate` must continue rejecting a reference profile when the supplied artifact lock is the candidate execution lock.
- Candidate validation uses a profile whose `reference_runtime_lock_id` equals the candidate lock ID and whose backend matches the candidate lock.
- The candidate fixture must have the same deterministic feature-vector identity as the reference fixture; no biological data are involved.
- Thresholds and runtime tolerances remain fixed prospectively and are not retuned from AL_001/GSE280090 outcomes.
- Technical PASS remains distinct from analytical and clinical validity.

## Repository boundary

Reference/candidate fixture files, profiles, runtime identities, dual-runtime reports, model bytes and biological validation payloads remain controlled local artifacts and are not committed to Git.

# Changelog

All notable changes to this research software are recorded here. The project has no clinically
validated release.

## Unreleased

- Fail closed in the standard SV lane unless the cuteSV executable resolves to an
  exact reviewed script body and the portable `python-shebang-v1` launch contract
  before any version probe or analysis command is allowed to run. Missing,
  unresolved or unknown executable identities are rejected; record the shebang SHA,
  launch contract and observed interpreter token while retaining full-source drift
  checks and temporary qualified staging.
- Validation impact: this is technical executable/launch qualification only. It does
  not change caller thresholds, assay rules, reference resources or reportability.
  Synthetic unknown/missing-executable regressions and the real identical-body,
  non-Python-shebang, worker-failure and positive-control paths are covered in CI.

- cuteSV worker-failure hardening (#91/#93): the exact known 2.1.3 entry-script
  copy now retrieves every signature-extraction task result and re-raises
  clustering exceptions. Both stock and the prior integer-only script upgrade
  to `ontseq-cutesv-2.1.3-worker-errors-v2`; the installed source is unchanged.
- Reject reviewed 2.1.3/2.1.4 `LocalError:` diagnostics before VCF normalization or
  atomic promotion. Record this marker guard as `not_qualified` for other versions.
- Validation impact: partial SV results after failed workers are no longer accepted.
  No successful-call thresholds or reference resources changed. Synthetic real-worker
  injection, positive controls, integer boundaries and cleanup are tested; this is
  technical failure-propagation qualification, not analytical/clinical validation.

- Fail closed when an Adaptive Sampling manifest declares a target BED but the
  supplied methylation policy sets `region_source=chromosome`. Reject the
  combination in preflight (`methylation.regions`) and again in
  `run_methylation()` before modkit or samtools is invoked, so no chromosome-wide
  pileup runs under a target-BED provenance fingerprint it did not constrain.
  Standard target-based Adaptive Sampling profiles are unchanged. Closes #87.
- Validation impact: closes a configuration-safety/provenance gap reachable only
  through a custom/direct policy that bypasses the standard profiles; no default
  profile, threshold, caller, or clinical-reportability rule changes.

- cuteSV now writes large signature/pickle intermediates to a managed user-cache
  directory (`~/.cache/ontseq/cutesv`, with `ONTSEQ_CUTESV_SCRATCH_ROOT` override)
  while retaining destination-local VCF staging and atomic promotion.
  A synthetic WSL comparison reproduced ENOMEM on the Windows mount and completed
  from native Linux storage. Worker count reduction alone was insufficient; a
  subsequent RAM-backed `/tmp`
  attempt exhausted its small filesystem, so the default now uses the user cache.

- Correct floating work-interval boundaries in the exact known cuteSV 2.1.3 entry
  script using a temporary, hash-verified copy. Integer chunks prevent read loss
  where pysam truncates a float fetch bound but cuteSV compares read starts to the
  untruncated bound. Preserve the installed tool and interpreter; record original
  and executed hashes plus `ontseq-cutesv-2.1.3-integer-chunks-v1`. Unknown script
  bodies are not rewritten or identified as qualified. Bind identity to resume and
  reject executable changes between planning and completion.

- Limit cuteSV to one worker by default and expose `--cutesv-threads` independently
  of other tools in run/analyze/serve/watch. Record the effective count in the SV
  plan and cuteSV provenance. This reduces concurrent signature-list memory; it
  does not guarantee completion for every input or change analytical thresholds.
- Carry the actual SV stage outcome into both result assemblers and their resume
  signatures. A later caller failure or a deselected stage cannot promote leftover
  caller files to current events, a completed SV module or fusion evidence.
  Existing archived reports are not silently rewritten.

- Identify the methylation executable by SHA-256 and include its identity in result
  provenance and resume signatures. Allow the independent-MM-group correction only
  for explicitly qualified builds of modkit PR #709 at source commit
  `9fb9aea763aa1b78ac1172fa6d6734e7955ccc70`; the unchanged `0.6.4` version banner alone
  does not qualify a binary. Reject changes between planning, preflight and completion.
- Retain the stock-tool independent-group guard and discard zero-exit partial results
  when modkit reports failed records. Add literal synthetic CpG position/count oracles,
  strand/clipping controls and full adapter acceptance/refusal checks.
- Validation impact: a qualified scanner can change previously lost or misassigned
  5mC/5hmC calls. Confidence/coverage thresholds and clinical release rules are unchanged;
  synthetic implementation checks establish neither patient accuracy nor clinical validity.

- Restrict the legacy GATK research run directory to mode `0700` on POSIX and launch native child processes with umask `077`, preventing local disclosure of BAM-header/sample and variant artifacts. No caller parameters, thresholds, result semantics, or clinical-reportability rules change.

- Add opt-in GATK 4.6.2.0 research adapters, separate candidate/config/result schemas,
  tumor-only and paired modes, native artifact provenance and synthetic contract tests.
- Repair strict typing and exact-version checks; reject suffixed or extended GATK builds.
  Include both adapter entry points in focused CI without weakening full-project gates.
- Validation impact: adapters remain separately invoked, non-reportable and not
  real-GATK-qualified or analytically validated. No canonical runner, CNV/SV policy,
  reference bundle, package version or automatic clinical release changes.

- Integrate the owner-selected complete Befund design as `befund-interactive-v1`:
  offline React report, real run-bound evidence, original CNV plots, model selection,
  synchronized chromosome/dilution views, explicit hypothetical VAF anchoring,
  module states, methylation, ISCN and complete technical evidence. Package the
  reproducible report bundle in wheels and stamp the report resume signature.

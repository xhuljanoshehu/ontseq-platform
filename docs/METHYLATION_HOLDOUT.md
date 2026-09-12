# Prospective two-source read-holdout recovery

**Research use only. This study evaluates recovery within two fixed biological
sources; it does not establish donor-independent or clinical validation.**

The two-source lane is the immediate experiment for new, compatible public or
institutionally released ONT methylation inputs. It supports Nanopolish and the
Dorado/modBAM/modkit adapters. It measures whether a prospectively registered
estimator recovers the fraction of retained whole read groups from source A on
read groups held out from calibration.

This study needs **two biological sources**, not four independently sampled donors.
`source_a` and `source_b` each retain their actual accession, version, biological
sample identity, provenance and input SHA-256. A source can share a donor with the
other source, and donor identity can remain unknown when source identity and
technical provenance are sufficient. Neither case establishes donor independence.
The optional stronger four-sample study remains described in
[METHYLATION_VALIDATION.md](METHYLATION_VALIDATION.md).

## Frozen experiment

Before inspecting the new test outcomes, register:

1. The two-source cohort and source distinction evidence. A/B must be genuinely
   declared distinct biological sources with distinct biological sample identifiers.
   Reusing identical input bytes or differently exported partitions of the same
   parent modBAM does not establish two source identities.
2. Complete comparable technical metadata: platform, flow cell, library chemistry,
   reference build/SHA-256, basecaller/model, modification caller/model, modification
   semantics and the adapter/filter policy. Dorado 5mC-only and 5mC/5hmC models must
   preserve their actual probability semantics; modkit provenance includes its
   version and parent modBAM SHA-256.
3. `PairedHoldoutSplitPolicy`: the calibration fraction, separate explicit A/B split
   seeds and the `sorted_read_identifiers_seeded_shuffle_v1` algorithm. The calibration
   fraction must match the same value in the registered estimator policy.
4. The versioned matrix, acceptance criteria and exact software/runtime fingerprint.
   The existing matrix `configs/methylation/independent_validation.v1.json` contains
   nine fractions × three budgets × ten replicates = **270 planned levels** and can
   be reused unchanged for this separately identified study. A reused matrix does
   not reuse an old registration or make the result donor-independent.

Registration records an operator declaration that outcomes were unseen and locks
the complete cohort, split, matrix and implementation by canonical JSON SHA-256.
It does not provide an authenticated external timestamp. Source holdout outcomes
already evaluated or used to adjust acceptance criteria cannot subsequently be
presented as a new prospective test of those choices. Metadata and file-integrity
checks before registration do not constitute recovery-outcome evaluation.

## Read partitions and marker isolation

For each source, sort read identifiers, shuffle once using its registered seed,
and take `round(total_read_groups * calibration_fraction)` groups for calibration,
clamped to leave at least one group in each pool. All calls belonging to a read
group stay together. Source read counts must equal calibration plus test counts;
pool sizes must reproduce the registered split rule.

The four resulting keys `calibration_a`, `calibration_b`, `test_a` and `test_b`
refer only to **computational read pools from two sources**. They must not be
reported as four biological specimens or four donors.

Marker discovery uses the calibration pools exclusively. Marker loci and calibration
counts are then fixed. Every in-silico mixture draws whole read groups only from
the held-out test pools, using the registered grid seeds and a constant requested
read-group budget at each depth. Insufficient pool size yields explicit NO_CALL.
No readgroup budget is genomic coverage in x or an independent CpG count.

The evidence records two source fingerprints/counts, two calibration counts, two
test counts and four pool-selection SHA-256 values. `split_manifest_sha256` binds
these values to the registration digest. Source read identifiers are not exported.
The runner establishes disjointness from the actual input read groups; reloading
summary JSON verifies content consistency, not independent biological authenticity.
A SHA-256 content lock is not a digital signature.

Replaying evidence JSON recomputes numerical results from the stored counts and
checks the registration and pool-manifest consistency. It does not reconstruct
read selections from a selection digest or prove that the supplied counts came
from the declared files. That stronger check requires re-executing the registered
runner against the checksummed local inputs. In particular, an edited per-level
selection digest alone is not independently authenticated by numerical replay.

## Reports and decisions

`PairedHoldoutRecoveryReport` always contains:

```json
{
  "validation_scope": "held_out_read_group_recovery_within_two_sources",
  "donor_independent_validation": "NOT_EVALUATED",
  "estimand": "fraction_of_retained_whole_read_groups_from_source_a",
  "research_only": true,
  "sensitive_output": true
}
```

Its top-level `decision` has the following meaning:

- **PASS:** source eligibility, prior registration, complete requested grid and all
  versioned depth/fraction recovery criteria pass within these held-out read pools.
- **FAIL:** an eligible complete experiment produces numerical estimates but at
  least one prespecified error, interval-coverage or NO_CALL criterion fails.
- **NO_CALL:** evidence is missing, registration/provenance is insufficient, the
  requested grid is incomplete, read budget cannot be supplied, or no levels are
  quantitatively evaluable. NO_CALL is not zero signal or a negative biological result.

Bias, MAE, RMSE, conditional interval coverage and mean interval width use completed
levels only. Errors are estimate minus the **realised** source-A read-group fraction.
The NO_CALL denominator includes every registered level, including missing results.
Reports retain depth-specific and fraction-specific metrics so that cancellation of
positive/negative bias or selective abstention cannot be hidden by one mean value.

The numerical estimator and replay checks are shared with the stronger independent
validation lane: marker counts, seeds, fractions, read budgets, constrained/raw
points, fit values, status, uncertainty draws, conditional intervals and derived
metrics are recomputed or checked. Contradictory evidence is an invalid input and
raises an error. The separate holdout contract prevents its PASS from being promoted
by changing the label to donor-independent validation.

## Local APIs and storage

The versioned split template is `configs/methylation/paired_holdout_split.v1.json`.
It reserves 25% of each source's retained whole reads using two fixed, distinct split
seeds. The existing 270-cell candidate matrix can be shared across study designs;
the registration and report determine the scientific scope.

```powershell
ontseq methylation-holdout-plan `
  --matrix configs/methylation/independent_validation.v1.json `
  --split-policy configs/methylation/paired_holdout_split.v1.json `
  --output results/methylation-holdout/readiness.json

ontseq methylation-holdout-fingerprint `
  --inputs results/methylation-holdout/local-inputs.json `
  --output results/methylation-holdout/input-fingerprints.json

ontseq methylation-holdout-register `
  --matrix configs/methylation/independent_validation.v1.json `
  --split-policy configs/methylation/paired_holdout_split.v1.json `
  --cohort results/methylation-holdout/cohort.json `
  --registration-id paired-source-study-001 `
  --confirm-test-outcomes-unseen `
  --output results/methylation-holdout/registration.json

ontseq methylation-holdout-run `
  --registration results/methylation-holdout/registration.json `
  --inputs results/methylation-holdout/local-inputs.json `
  --output-dir results/methylation-holdout/run-001
```

The local input file binds `samples.source_a` and `samples.source_b`, each using
`calls_path`, `input_format` and `metadata_path`. SAM/BAM or modkit TSV also require
`reference_path` and `adapter_policy_path`. These are two genuine source records,
not four duplicated specimen declarations. Fingerprints are computed from local
bytes; biological source identity and provenance must come from the approved run
record. Full input examples are generated only in synthetic tests; schemas describe
the real contracts without inventing missing biological metadata.

The run writes JSON, CSV and HTML with suffix `.methylation-holdout`. The CSV carries
the study scope and `donor_independent_validation` on every level, and the HTML states
the same limitation alongside its accuracy and NO_CALL metrics. Output files are
never overwritten. `methylation-holdout-evaluate --registration ... --evidence ...
--output-dir ...` replays a stored evidence object; omit `--evidence` to document
missing evidence as NO_CALL. No CLI command downloads inputs or installs software.

The models and evaluator are in `src/ontseq_platform/methylation_holdout.py`:

- `PairedHoldoutCohort`, containing exactly `samples.source_a` and `samples.source_b`,
  `source_distinction_evidence` and `access_basis`.
- `PairedHoldoutSplitPolicy` and `preregister_paired_holdout(...)`, creating the full
  source/split/matrix/code lock before the test.
- `PairedHoldoutEvidence` and `evaluate_paired_holdout(...)`, checking submitted
  evidence and producing a report. Without evidence the result is explicit NO_CALL.
- `holdout_split_manifest_sha256(...)`, binding source hashes, pool counts and pool
  selection digests to that registration.

`execute_registered_paired_holdout(...)` in
`src/ontseq_platform/methylation_holdout_runner.py` binds local adapter inputs,
checks the registration, constructs the actual read partitions, discovers markers
and executes the mixtures. It uses the same adapter contracts as the four-source
runner. Registered source metadata and implementation must match local execution.

Store registrations, local path bindings, evidence, reports and all genomic inputs
in the ignored `results/` tree or an institutionally approved local data directory.
No real source data, read identifiers, donor metadata, credentials or analysis results
belong in Git. The repository contains only contracts, study policy and synthetic tests.

## Validation impact and next evidence

### Prospective targeted FMR1 pilot

`configs/methylation/fmr1_targeted_holdout.v1.json` is a separate candidate matrix
for the two native iPSC sources in [GSE218679](METHYLATION_ADDITIONAL_DATA_SEARCH.md).
It fixes the same nine fractions and ten replicate seeds at budgets of 20, 50 and
100 retained read groups, producing 270 planned levels. These smaller budgets are
fixed before genomic acquisition or methylation outcome inspection. Archived read
counts only rule out infeasible larger budgets; they are not usable read-group counts.
Admission still requires intact usable pools and complete compatible provenance.

The original 1000-read candidate's accuracy limits and all estimator gates, including
the minimum 20 markers, are retained. The selected data may contain only 19 promoter
CpGs; no marker threshold is reduced to accommodate this possibility. This study is
restricted to the deposited FMR1 targeting and filtering domain. Its rationale and
matrix ID accompany reports so that any recovery cannot be presented as genome-wide
or donor-independent performance. The matrix is a prospective template, not a
registered study or a claim of biological success. Generate its readiness record
with the same `methylation-holdout-plan` command and this matrix path.

### Scope of the implementation change

This additional lane removes the unnecessary four-donor prerequisite for the requested
two-source feasibility/recovery experiment. It preserves the stricter independent
lane and keeps the biological scope of each explicit. Synthetic tests cover
reproducible whole-read splitting, marker isolation from test biology, numerical
recovery and deliberate reversal, pool-count/digest changes, incomplete grids,
low marker counts, provenance and timestamp inconsistency, and attempts to change
the scientific scope or aggregate decision. A synthetic PASS verifies software
behaviour and is not evidence of recovery on public biological samples.

Read partitions share the same biological sources; repeated mixtures reuse finite
test pools. Neither operation creates independent donors or independent physical
dilution preparations. The current nominal 95% intervals remain conditional on fixed
calibration and marker choices and do not propagate calibration, donor, marker-selection
or within-read/cross-marker correlation uncertainty. These remain later studies.

Successful recovery here permits studying a conversion to DNA mass using independent
reference measurements; it does not establish that conversion. Cell fraction,
tumour purity and subclone abundance require their own truth definitions and evidence.
SNV, CNV and ploidy remain separate layers. A 10–20% detection limit remains a
hypothesis, irrespective of the result of this holdout experiment.

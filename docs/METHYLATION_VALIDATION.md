# Prospective independent methylation recovery

**Research use only. No analytical or clinical validation has been established.**

For the immediate two-source experiment using disjoint read groups, see
[`METHYLATION_HOLDOUT.md`](METHYLATION_HOLDOUT.md). This document describes the
additional donor-independent study; its four-source requirement is not imposed on
the paired-source holdout design.

This lane evaluates whether the existing source-A methylation coefficient recovers the
known **fraction of retained whole read groups** in an independent in-silico experiment.
It does not convert the coefficient into DNA mass, cell fraction, tumour purity, blast
fraction or a subclone fraction. It does not establish a detection limit. The previous
GSE210260 feasibility result remains a failed technical recovery experiment; it is not
independent evidence for this matrix.

The prospective candidate is
`configs/methylation/independent_validation.v1.json`. The file is a study template;
registration also requires a fixed four-sample cohort, exact input SHA-256 values and
the executable software/runtime fingerprint **before examining test outcomes**.

## Locked matrix and candidate acceptance criteria

The matrix contains fractions 0, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95 and 1.00;
budgets of 1,000, 5,000 and 10,000 retained whole read groups; and ten explicit
replicate seeds. This is 270 planned grid cells. A read-group budget is not genomic
coverage in x, a count of CpGs or a count of independent biological observations.
Each mixture samples without replacement within its source pool. Mixtures and depths
can reuse the same finite pools, so their errors are dependent.

All error limits below use absolute fraction units, e.g. 0.05 means five percentage
points. The limits are engineering acceptance **candidates**, selected before the
next independent experiment; they are not validated assay thresholds.

| Retained read groups | Maximum absolute bias | Maximum MAE | Maximum RMSE | Minimum empirical coverage of nominal 95% intervals | Maximum unconditional NO_CALL |
| --- | --- | --- | --- | --- | --- |
| 1,000 | 0.05 | 0.10 | 0.12 | 0.90 | 0.10 |
| 5,000 | 0.03 | 0.07 | 0.09 | 0.90 | 0.10 |
| 10,000 | 0.03 | 0.05 | 0.07 | 0.90 | 0.10 |

Every fraction at every depth must also meet a maximum absolute bias of 0.10, 0.08
or 0.07 respectively, empirical interval coverage of at least 0.80 and a NO_CALL rate
of at most 0.20. These strata prevent errors at the endpoints, selective abstention
and positive/negative bias cancellation from being hidden by aggregate statistics.
With only ten resamples per fraction the coverage step is coarse; these criteria do
not constitute a confidence statement about population coverage.

The estimator policy fixes marker selection, counts, uncertainty draws, fit and
abstention before execution. Its historical `calibration_read_group_fraction` field
is unused in this independent lane: all retained reads in the two calibration
samples are used for marker discovery, and none of those samples enters a test mix.
The historical policy's aggregate technical recovery thresholds are also unused;
the independent matrix acceptance limits are authoritative.

## Eligibility and independence

Four explicit roles are required: `calibration_a`, `calibration_b`, `test_a` and
`test_b`. Each role carries an accession, dataset version, biological sample and
donor identifiers, biological source, cell type, material type, provenance reference,
input SHA-256 and full technical metadata. Store these records locally as sensitive
analysis inputs. Public metadata references can be catalogued separately; no
patient identifiers, read identifiers, genomic files or results belong in Git.

Calibration/test samples and input fingerprints must be different. Calibration and
test donor sets must not overlap. Calibration and test within A or B must describe
the same biological source and cell type; A and B must be distinct sources. All four
roles require matching declared platform, flow cell, library kit, basecaller/model,
caller/version, reference build/hash, modification semantics and adapter policy.
Unknown fields block independent recovery. Dorado requires its modification model
and complete cytosine modification set; modkit requires the version and parent
modBAM fingerprint. Reusing a parent BAM through differently filtered exports does
not create independent samples and is rejected as overlapping provenance.

This conservative donor criterion defines the scope of this study. Separate cell
line aliquots, technical library replicates and random read splits can support
technical feasibility but do not satisfy donor-independent biological recovery.
Declared metadata and sample labels cannot themselves prove biological identity;
their external provenance remains reviewable evidence. Four donors also do not
establish generalisation across a donor population.

## Registration and local execution

The Python contracts and evaluator live in
`src/ontseq_platform/methylation_validation.py`; local file execution is in
`src/ontseq_platform/methylation_validation_runner.py`. Both the Nanopolish and
Dorado/modkit adapters feed the same estimator through typed read-call sources.
The old two-file Nanopolish experiment remains available as a technical lane.

The independent commands are also callable using
`python -m ontseq_platform.methylation_validation_cli` followed by the command below.
All output paths should remain in the ignored `results/` tree or an institutionally
approved local data directory. Existing output files are not overwritten.

```powershell
ontseq methylation-validation-plan `
  --matrix configs/methylation/independent_validation.v1.json `
  --output results/methylation-validation/readiness.json

ontseq methylation-validation-fingerprint `
  --inputs results/methylation-validation/local-inputs.json `
  --output results/methylation-validation/input-fingerprints.json

ontseq methylation-validation-register `
  --matrix configs/methylation/independent_validation.v1.json `
  --cohort results/methylation-validation/cohort.json `
  --registration-id independent-study-001 `
  --confirm-test-outcomes-unseen `
  --output results/methylation-validation/registration.json

ontseq methylation-validation-run `
  --registration results/methylation-validation/registration.json `
  --inputs results/methylation-validation/local-inputs.json `
  --output-dir results/methylation-validation/run-001
```

`local-inputs.json` binds each of the four roles under `samples` to `calls_path`,
`input_format` (`nanopolish`, `modkit` or `modbam`) and `metadata_path`. The modification
adapters additionally require `reference_path` and `adapter_policy_path`. These local
bindings contain no inferred biological assignments; the cohort must provide those.

Registration locks canonical JSON digests for matrix and cohort plus the software
source/runtime digest. JSON whitespace and key order do not change these digests.
Source files and the implementation are rechecked during execution. Changing the
policy or choosing a different cohort requires a new prospectively registered study;
re-labelling an already observed result is not preregistration. The timestamp records
an operator declaration, **not a cryptographically authenticated external timestamp**.
Content integrity is not proof of biological provenance or proof that outcomes were
unseen. A controlled external registry and reviewer record are future governance work.

Numerical evidence replay uses stored marker/observation counts. It cannot recover
read identifiers from calibration or per-level selection digests and does not
authenticate those selections against the raw inputs. Re-execute the registered
runner on the checksummed local files to verify the actual seeded read choices;
changing a selection digest alone need not change the numerical replay decision.

To evaluate a stored evidence contract, use `methylation-validation-evaluate` with
`--registration`, `--evidence` and `--output-dir`. Omitting `--evidence` creates an
explicit NO_CALL report with all planned cells missing. The `plan` command records
NO_CALL/NOT_REGISTERED when no cohort or experiment exists.

## Decision and integrity rules

`decision` is the top-level scientific result:

- **PASS:** eligible, prospectively declared, complete grid and every locked depth and
  fraction stratum meets the candidate acceptance limits. This is only the defined
  donor-independent in-silico retained-read recovery scope.
- **FAIL:** an eligible complete experiment has quantitative results but at least one
  accuracy, interval or unconditional abstention criterion fails. Completing an
  estimator calculation does not override this decision.
- **NO_CALL:** missing or ineligible biological evidence, no prior registration,
  mismatched execution provenance, incomplete grid or no quantitative results.
  NO_CALL means insufficient evidence, never zero signal or a negative biological result.

Bias, MAE, RMSE, interval coverage and mean interval width use **completed levels only**.
The unconditional NO_CALL denominator contains every preregistered cell: explicit
NO_CALL plus missing cells count as unavailable. Reports show expected, observed,
completed, explicit NO_CALL and missing counts together. No omitted cell can improve
the denominator. There is no post-hoc exclusion or best-depth selection.

Saved reference markers, marker counts, fractions, budgets, seeds, constrained/raw
estimates, model fit, conditional intervals, statuses and aggregate metrics are
checked on reload. Numerical values and the conditional Monte Carlo intervals are
recomputed from counts, the registered policy and deterministic seeds. Duplicate or
unregistered grid cells and contradictory numerical content are invalid inputs;
they raise a validation error rather than produce a scientific estimate. Checksums
detect inconsistent content; they are not digital signatures or independent
authentication of the submitted source data.

## Validation impact and remaining studies

This change adds a biological eligibility boundary and a new prospective decision
layer. It does not change the existing estimator into an analytically validated
method. Synthetic tests cover successful deterministic software recovery, reversed
biology with cancelling aggregate bias, missing grid cells, low marker count,
post-hoc declarations, donor reuse, parent-BAM reuse, altered hashes/seeds/intervals,
and derived-result tampering. Synthetic PASS is a software regression result only.

The current conditional interval does not propagate uncertainty from calibration,
between-donor variability, marker discovery or correlated methylation calls. Ten
seeded resamples are not ten donors or ten independent physical dilution experiments.
Next studies must lock a donor/replicate design and assess hierarchical or clustered
uncertainty using calibration resampling and held-out biological donors; altered
uncertainty methods require a new code and study lock.

Only after successful independent recovery should orthogonal reference measurements
test a conversion from retained-read coefficient to DNA mass or another quantity.
Read length, capture efficiency, methylation state, DNA yield, CNV and ploidy can all
affect such a conversion. SNV, CNV, ploidy and subclones remain distinct evidence
layers. A proposed 10–20% detection limit remains a hypothesis requiring its own
predefined truth material, error criteria and independent dilution experiment.

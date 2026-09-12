# Paired-source methylation mixtures and a methylation-based source coefficient

For the prospective four-sample validation matrix, use
[`METHYLATION_VALIDATION.md`](METHYLATION_VALIDATION.md). The original two-file Nanopolish
experiment described below remains available. Independent Dorado/modBAM and modkit per-read
input handling is documented in [`MODBAM_ADAPTER.md`](MODBAM_ADAPTER.md); these formats are not
silently re-labelled as Nanopolish markers.

## Status and intended quantity

Research use only. This path is a deterministic technical experiment for two declared,
biologically distinct methylation sources. It has not been analytically or clinically
validated.

The reported value \(\hat f\) is a methylation-based source-A mixture coefficient and is
calibrated against the **known fraction of held-out read groups drawn from source A**. A read group
here is one Nanopolish `read_name` together with all methylation calls assigned to that read.
Without an additional, independently validated biological calibration, \(\hat f\) is not a
tumour fraction, blast fraction, fraction of cells or fraction of DNA mass. Calling one input
"tumour" and the other "normal" does not change that meaning.

This distinction is the central result contract. The lane answers whether the software can
recover known read-group mixtures of these two sources. It does not yet answer how much tumour
is present in an unknown specimen.

## Inputs

The current Nanopolish adapter is `nanopolish-call-table-v2` with security profile
`bounded-tsv-gzip-v1`. `maximum_input_bytes_per_source` independently limits both
compressed file bytes and decoded bytes, including blank lines. Before CSV parsing,
headers are limited to 16 KiB, physical data lines to 64 KiB and columns to 64.
C0/C1 control characters are rejected in read identifiers and selection hashes.
Source summaries record `adapter_version`, `input_security_profile` and the actual
`decoded_size_bytes`. These guards do not establish biological eligibility or provide
an antivirus certification. Existing v1 artifacts are retained as historical records;
they are not relabelled v2 or made valid under a new software lock automatically.

The first adapter accepts two Nanopolish `call-methylation` tables as uncompressed `.tsv` or
gzip-compressed `.tsv.gz`. Both inputs must:

- contain the complete Nanopolish call-table columns, including coordinates, strand,
  `read_name`, likelihoods, motif count and sequence;
- have different source IDs and different SHA-256 file contents;
- be explicitly declared by the operator as biologically distinct and use the same declared
  genome build;
- contain at least two read groups, so calibration and mixture data can be kept separate.

The declared genome build is a contract supplied by the operator. A call table alone cannot
prove which reference was used upstream. `GRCh38` is a coordinate reference and does not supply
a healthy methylation pattern. Each source in this experiment is measured data from a separate
input table.

The mandatory `--confirm-biologically-distinct-sources` flag records that declaration. Different
IDs and file hashes prove only that two files differ; the software cannot infer whether they came
from different donors, specimens or biological materials. That relationship remains auditable
run metadata.

Only canonical chromosomes `chr1` through `chr22`, `chrX` and `chrY` enter the calculation;
names without the `chr` prefix are normalized. Other contigs are excluded and counted. For a
row containing `num_motifs = q`, the event threshold is
`log_likelihood_ratio_threshold * q`. A log-likelihood ratio at least that value is methylated,
a value at or below its negative is unmethylated, and an intermediate value is ambiguous. The
multiplication matters because one Nanopolish row can represent more than one motif. Ambiguous
calls are reported but do not enter beta fractions.

No human sequencing table belongs in Git. Public or institutional source data remain in the
controlled external reference bundle; only synthetic fixtures may be committed. The report
stores source IDs, file size and SHA-256 fingerprints, never source paths or read names.

## Upstream provenance is required

Comparable interpretation requires more provenance than the call-table columns provide. The CLI
accepts a typed metadata sidecar for each source through `--source-a-metadata` and
`--source-b-metadata`. Each sidecar records:

- the exact reference FASTA identifier, declared genome build and SHA-256 checksum;
- caller name (restricted to `nanopolish` for this adapter) and version;
- sequencing platform, flow-cell product code and library kit;
- basecaller name, version and model;
- the modification context/model. Source identity and relationship remain separate controlled
  run metadata.

The structured report embeds those declared fields under each source's `upstream` object, along
with the call-table SHA-256, operator-declared build, ONTSeq software version and Git commit.
Known values must agree between the two sources or the run fails closed. Omitted metadata fields
become the literal `unknown` (and an omitted reference checksum becomes `null`) and produce a
reviewer-visible incomplete-provenance warning; an unspecified Git commit is `UNKNOWN`.

These fields are declared provenance, not facts inferred from the TSV. The call table still
cannot prove the reference or upstream command. Retain complete caller/basecaller commands and
any additional protocol context in the approved run record because the current sidecar has no
fields for them. A technically executable result with required provenance marked unknown must
not be compared across experiments or given biological meaning.

## Experimental separation

For each source, the program sorts the read names and performs a seeded split:

1. the **calibration pool** defines source-specific methylation at each marker;
2. the disjoint **held-out pool** supplies reads for the in-silico mixtures.

No read group used to define a marker is used in a mixture level. Every requested fraction uses
the same total read-group budget. Source A and source B read groups are sampled whole, with
deterministic seeds. The report records the target fraction, the realized fraction after integer
rounding and a SHA-256 selection digest. The digest makes a draw auditable without disclosing
the read names.

Replicates are repeated seeded subsamples from the same two held-out pools. They measure the
stability of this computational sampling procedure; they are not independent biological
replicates.

## Exact markers and the model

A reference marker must match between the two calibration pools by the exact tuple

```
(canonical chromosome, strand, Nanopolish start, Nanopolish end, number of motifs,
 SHA-256 of the normalized Nanopolish sequence context)
```

The sequence context is upper-cased before hashing; the raw sequence is not exported. This keeps
rows with different sequence contexts from being treated as the same exact marker. Broad regions,
nearest-neighbour matching and coordinate lift-over are not used. Each source
must meet `minimum_reference_valid_calls`, and the absolute contrast must satisfy
`minimum_absolute_delta_beta`. Markers are ranked by absolute contrast before the optional
maximum is applied. This is an absolute beta-difference filter, not a fold-change rule.

Beta values are used to select markers with source contrast:

\[
\beta_j = \frac{m_j}{m_j+u_j}, \qquad
\Delta\beta_j = \beta_{A,j}-\beta_{B,j}.
\]

The fitted quantity itself uses **methylated (M) and unmethylated (U) call rates per total read
group**, as separate channels. For source or mixture \(s\), marker \(j\), and channel
\(c\in\{M,U\}\):

\[
r_{s,j,c}=\frac{n_{s,j,c}}{N_s}, \qquad
r_{mix,j,c}=r_{B,j,c}+f\left(r_{A,j,c}-r_{B,j,c}\right),
\]

where \(N_s\) is the total number of read groups in that calibration or mixture pool. Therefore
an ambiguous call or the absence of a call affects the denominator rather than being silently
renormalized away. This call-rate formulation carries marker callability into the fit.

The unconstrained weighted least-squares coefficient is

\[
f^*=\left(
\frac{\sum_{j,c} w_{j,c}d_{j,c}(r_{mix,j,c}-r_{B,j,c})}
     {\sum_{j,c} w_{j,c}d_{j,c}^2}
\right), \qquad
d_{j,c}=r_{A,j,c}-r_{B,j,c}, \qquad
\hat f=\operatorname{clip}_{[0,1]}(f^*).
\]

Each channel is weighted by the inverse sum of three Jeffreys-binomial posterior variances:

\[
w_{j,c}=\left[V_J(n_{A,j,c},N_A)+V_J(n_{B,j,c},N_B)
                 +V_J(n_{mix,j,c},N_{mix})\right]^{-1},
\]

\[
V_J(k,N)=\frac{(k+\tfrac12)(N-k+\tfrac12)}{(N+1)^2(N+2)}.
\]

The reported estimate is restricted to zero through one, while the unconstrained coefficient
remains visible for model diagnostics. Accuracy summaries compare the estimate with the
**realized** source-A read-group fraction rather than only the requested target.

The fit diagnostic is a variance-standardized, dimensionless RMSE across the retained M and U
channel residuals:

\[
\operatorname{RMSE}_{std}=\sqrt{
  \frac{\sum_{j,c} w_{j,c}\left[r_{mix,j,c}-\hat r_{mix,j,c}\right]^2}{K-1}
},
\]

where \(K\) is the number of channel residuals. It is not an RMSE on the beta scale and must not
be compared with an absolute methylation-percentage error. Because M and U channels from one
marker are correlated, this technical disagreement heuristic is not a calibrated chi-square
statistic.

## Conditional uncertainty

For each usable marker, the interval model divides calibration read groups into four mutually
exclusive categories:

1. methylated (M);
2. unmethylated (U);
3. ambiguous;
4. missing, derived as `total_read_groups - M - U - ambiguous`.

The source-A and source-B four-state rates stay fixed at their locked calibration values. For the
fitted coefficient \(\hat f\), the expected mixture distribution is calculated for all four states
as \(p_{mix}=p_B+\hat f(p_A-p_B)\). Each Monte Carlo iteration draws the mixture rates from a
Dirichlet distribution centred exactly on this fitted distribution, with concentration equal to
the mixture read-group count, and refits the coefficient with the point fit's fixed
inverse-variance weights. A state with fitted probability zero remains zero. The confidence
level, number of draws and minimum successful draw fraction are policy fields.

This is a conditional, parametric mixture-pool interval under the fitted linear model. Keeping the
calibration lock fixed avoids the errors-in-variables displacement produced by independently
perturbing both noisy endpoints. The interval does not propagate uncertainty of the calibration
sources, model misspecification, donor variation, cell composition, library preparation,
platform transfer or correlation of CpGs on the same molecule. It is therefore neither a full
frequentist coverage guarantee nor a clinical confidence interval. Calibration-source variation
must be assessed with independent sources and experimental replicates.

## Fail-closed and `NO_CALL`

Each mixture level is `COMPLETED` only when all locked gates pass. It becomes `NO_CALL`, without
a point estimate or interval, when any of these conditions applies:

- fewer than `minimum_markers` exact calibration markers pass coverage and contrast filters;
- too few mixture markers reach `minimum_mixture_valid_calls`;
- the unconstrained coefficient lies farther than `maximum_extrapolation` outside the interval
  from zero to one;
- variance-standardized model-fit RMSE exceeds
  `maximum_standardized_model_fit_rmse`;
- too few uncertainty draws can be estimated;
- the resulting interval is wider than `maximum_confidence_interval_width`.

`NO_CALL` means that the model declined to quantify this level. It is not zero source-A signal
and not a negative biological result. The aggregate report is `COMPLETED` only when every level
is quantifiable, `PARTIAL` when completed and no-call levels coexist, and `NO_CALL` when no level
can be quantified. Accuracy metrics cover completed levels only; `quantifiable_fraction` and
`grid_fully_quantified` make that scope visible.

Quantifiability and known-fraction recovery are separate decisions. Once enough of the requested
grid is quantifiable, the report compares the estimates with the known realized fractions and
emits `technical_recovery_assessment`:

- `PASS` when the locked technical limits for absolute mean bias, MAE, RMSE and empirical
  confidence-interval coverage all pass;
- `FAIL` when at least one of those limits is exceeded;
- `NOT_EVALUABLE` when too little of the grid was quantifiable to assess recovery.

This assessment is an engineering regression control. Its thresholds are unvalidated and it is
not an analytical or clinical validation verdict. In particular, a report may be `COMPLETED`
because every level produced a number while its recovery assessment is `FAIL` because those
numbers do not recover the known mixtures accurately.

## Running the experiment

Use two call tables that have been generated under a compatible upstream protocol and keep them
outside the repository:

```bash
ontseq methylation-mixture \
  --source-a-calls source_a.tsv.gz \
  --source-b-calls source_b.tsv.gz \
  --source-a-id SOURCE_A \
  --source-b-id SOURCE_B \
  --confirm-biologically-distinct-sources \
  --source-a-metadata source_a.metadata.json \
  --source-b-metadata source_b.metadata.json \
  --source-a-genome-build GRCh38 \
  --source-b-genome-build GRCh38 \
  --analysis-id METH_MIX_001 \
  --policy configs/methylation/paired_source_nanopolish.technical.yaml \
  --output-dir results/methylation-mixture/METH_MIX_001
```

The output directory contains a schema-validated JSON report, a level table in CSV and a
self-contained HTML view. The JSON is the structured source of truth. The report includes input
fingerprints, policy, algorithm version, exact reference markers, known and estimated fractions,
conditional intervals, standardized fit RMSE, seeded selection digests, warnings and limitations.
It also records the empirical interval coverage and the separate technical recovery assessment.
The HTML and CSV display requested and realized read-group fractions separately. On reload, the
typed JSON contract recomputes completed point estimates, fit statistics and seeded intervals
from its locked marker counts, verifies that calibration sizes and level seeds follow the locked
policy, and rejects inconsistent values even when the report status is `NO_CALL`.

Every output is a **sensitive derived genomic artifact** and the JSON declares
`sensitive_output: true`. In particular, the JSON contains exact marker coordinates and source-
and mixture-specific call counts. CSV and HTML contain fewer details but still expose derived
sample signals. Store all three only in the approved local/controlled result area; do not commit,
upload or publish them. The omission of source paths and read names reduces disclosure but does
not make an output anonymous or non-sensitive.

## Technical policy

`configs/methylation/paired_source_nanopolish.technical.yaml` contains executable engineering
defaults:

| Field | Technical default | Meaning |
| --- | ---: | --- |
| `log_likelihood_ratio_threshold` | 2.5 | Symmetric per-motif threshold; the row threshold is this value times `num_motifs` |
| `source_a_fractions` | 0, 0.05, 0.1, 0.2, 0.5, 1 | Requested source-A read-group fractions |
| `replicates` | 3 | Seeded subsamples per level |
| `calibration_read_group_fraction` | 0.25 | Fraction reserved separately in each source |
| `total_mixture_read_groups` | `null` | Derive the largest common fixed budget from the held-out pools |
| `minimum_reference_valid_calls` | 5 | Per source and exact calibration marker |
| `minimum_mixture_valid_calls` | 3 | Per marker in one mixture level |
| `minimum_absolute_delta_beta` | 0.2 | Minimum absolute source contrast |
| `minimum_markers` / `maximum_markers` | 20 / 2000 | Model marker-count gates |
| `uncertainty_draws` | 500 | Seeded four-state conditional Dirichlet draws around the fitted mixture; calibration fixed |
| `confidence_level` | 0.95 | Central interval probability |
| `minimum_successful_uncertainty_fraction` | 0.8 | Draw-completion gate |
| `maximum_confidence_interval_width` | 0.5 | Width gate before `NO_CALL` |
| `maximum_standardized_model_fit_rmse` | 3.0 | Dimensionless variance-standardized residual gate before `NO_CALL` |
| `maximum_extrapolation` | 0.1 | Maximum tolerated distance of the unconstrained coefficient outside 0–1 |
| `maximum_input_bytes_per_source` | 100,000,000 | Fail-closed on the on-disk input size before parsing |
| `maximum_rows_per_source` | 2,000,000 | Fail-closed when the non-empty input-row count exceeds this value |
| `minimum_quantifiable_fraction_for_recovery` | 0.8 | Minimum completed share of the requested grid before recovery is assessed |
| `maximum_absolute_mean_bias` | 0.1 | Technical recovery limit for absolute mean bias |
| `maximum_mean_absolute_error` | 0.15 | Technical recovery limit for MAE |
| `maximum_root_mean_squared_error` | 0.2 | Technical recovery limit for RMSE |
| `minimum_confidence_interval_coverage_fraction` | 0.8 | Minimum empirical coverage of known fractions among completed levels |

These numbers make the workflow executable and reviewable. None is a validated detection limit,
reportability threshold or clinical acceptance criterion. Deriving the read budget separately
for every source pair also makes experiments differ in depth; pin `total_mixture_read_groups`
before comparing experiments.

## SMALL public-data smoke scope

A run labelled `SMALL` public-data smoke may establish only that the Nanopolish adapter, input
contract, deterministic split, M/U call-rate fit, conditional Dirichlet interval and
JSON/CSV/HTML writers can
execute on a small public input subset. It is a technical integration check, not validation. It
cannot establish tumour purity, blast fraction, biological source deconvolution, a detection
limit, marker generalisation, clinical accuracy or full-scale resource behaviour. Public source
files and derived outputs also remain outside Git and under their applicable access and licence
terms.

## Runtime and memory boundary

Plain and gzip-compressed tables are opened and read line by line, but the current parser retains
every accepted canonical call grouped by `read_name` in RAM so it can perform disjoint splitting
and whole-read resampling. Memory therefore grows with the number of canonical calls and read
groups. The policy fails closed before parsing when the on-disk file exceeds
`maximum_input_bytes_per_source` and during parsing when `maximum_rows_per_source` is exceeded.
For `.gz` input the byte limit applies to the compressed file, so it is not a direct peak-RAM
bound. `maximum_markers` is applied only after parsing and does not cap this peak.

The byte and row limits are technical safety caps, not a validated capacity claim or measured RAM
preflight. There is currently no bounded-memory end-to-end stream or disk-backed spill. Large
whole-genome call tables require explicit resource sizing; a successful `SMALL` smoke run says
nothing about production-scale memory. A future large-input adapter should use an externally
sorted or indexed, disk-backed two-pass design while preserving whole-read selection.

## Relationship to the other lanes

This command is standalone. It consumes Nanopolish call tables and does not use the canonical
modkit stage's chromosome- or target-region summaries. Those summaries are useful descriptive
outputs, but they cannot provide the exact shared markers or per-read grouping needed here.

The paired-source equation assumes that each retained marker follows a linear mixture of two
stable endpoints. Copy number, ploidy, allelic imbalance and structural variants can change the
number or identity of contributing molecules and can violate that assumption. The current model
does not correct them and must not be described as CNV-independent normalization. CNV-aware or
haplotype-aware inference belongs in a later, separately validated extension.

## What is needed before biological interpretation

The current implementation can test reproducibility and recovery of known read-group fractions
from computational mixtures. Promotion beyond that requires a pre-registered study with independent donors and
orthogonal truth, wet-lab dilution series, controlled genome builds and upstream callers,
coverage and input-mass strata, repeatability and reproducibility, interference studies, CNV and
ploidy strata, predefined bias/precision/no-call criteria and held-out external validation.

An observed 10 or 20 percent performance boundary is a result to measure under those conditions,
not a default detection limit. Translation from \(\hat f\) to tumour or blast fraction needs its
own biological model and validation.

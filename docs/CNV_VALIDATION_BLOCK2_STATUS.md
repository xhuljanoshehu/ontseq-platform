# CNV validation Block 2 status

Research Use Only.

Implementation status: complete.

Acceptance rule: Block 2 is a technical PASS only when the exact current implementation head has successful same-head CI, MARLIN runtime smoke, Desktop bundle contract and Desktop CI evidence recorded in PR #79. Workflow run IDs are recorded in the PR conversation so recording them does not mutate the verified repository head.

Block 2 implements Full Evidence + Traceability only. It does not calculate validation performance metrics.

Implemented on the isolated `feat/cnv-validation-program` branch:

- caller-native full-evidence records for run summaries, fits, bins, segments, events and chromosome summaries;
- explicit preservation of source coverage, tumour/blast fraction, bin size, caller version, dependency versions, parameter snapshots, execution identity and reference/input fingerprints;
- selected and alternative caller fits retained as separate records;
- 100/500/1000-kbp evidence remains independently addressable rather than being collapsed into a preferred resolution;
- recursive rejection of non-finite values and non-JSON caller metadata;
- safe relative native-artifact paths plus SHA-256 and byte-size fingerprints;
- explicit `OBSERVED`, `FAILED`, `NO_CALL`, `NOT_ASSESSABLE` and orthogonally supported biological-negative states;
- contribution states that preserve excluded, secondary and exploratory evidence rather than deleting it;
- derived normalized CNV events with mandatory source full-evidence IDs;
- fail-closed identity binding between normalized events and source evidence across specimen, caller/runtime, build, data basis, reference, input, coverage, tumour fraction, bin size and replicate context;
- a sealed evidence manifest with `retain_all_evidence=true`, duplicate-address rejection, dangling-reference rejection, orphan-artifact rejection and canonical SHA-256 content locking;
- aggregate traceability records that store numerator, denominator, excluded and normalized-event memberships without storing or computing the future metric value;
- numerator membership constrained to the denominator, with exclusion memberships kept disjoint;
- traceability verification from future aggregate/question/stratum identity through normalized events and full evidence to native caller artifacts;
- duplicate analytical trace-address rejection and tamper-evident traceability-index locking;
- versioned JSON schemas for the complete evidence manifest and traceability index.

Evidence-retention rule:

- every evaluated caller-native value is retained when emitted into the validation evidence layer;
- normalization is an additive derived layer and never replaces caller-native evidence;
- exclusion from a primary metric changes contribution state only and never removes the underlying record;
- a future aggregate is only a view over retained evidence; the aggregate cannot become the sole surviving representation;
- every normalized event and aggregate membership remains auditable to explicit source IDs and fingerprinted native artifacts.

Out of scope for Block 2:

- real biological or patient genomic data;
- QDNAseq/ACE, ichorCNA or Spectre execution changes;
- caller result-adapter implementation beyond the generic evidence contracts;
- TP/FP/FN assignment logic;
- sensitivity, specificity, precision, F1, LoD or reproducibility calculations;
- production caller selection or clinical reportability decisions;
- merge or version bump.

Package version remains `0.8.2`. PR #79 remains Draft and unmerged. A technical PASS establishes only the engineering evidence/traceability contract, not analytical or clinical validity.

# Supplemental declarative routing contract

Merge scope (2026-09-20): this is a separate, planning-only contract. It does not replace multicaller_contracts.MultiCallerPlan, the canonical lane planner, caller runtimes, or their qualification registry. Its existing_adapter flags describe this supplemental route, not platform-wide installation status. All its lanes remain NOT_RUN and execution_enabled=false. No implicit conversion into the canonical execution plan is provided.

# ONTSeq multi-caller integration design

Date: 2026-09-19. Scope: Research Use Only; human review required.
Baseline: `b6e0338b6f042239a1399e3438e6c2bab4e48b08` on
`feat/cnv-validation-program` (Draft PR #79). Package remains 0.8.2.

## Intent and scope

Integrate the important, complementary CNV and SV methods into one traceable system,
not an unqualified "run every caller" button. Preserve the existing production graph,
all caller-native evidence, and Blocks 1–4. No merge, release, patient data or real
biological execution is authorized by this engineering change.

The first independently reviewable deliverable is **5A: a typed catalog and a
fail-closed, content-addressed research planner**. It does not execute external
programs. Catalog membership, a compatible input declaration, a runtime adapter,
a successful runtime smoke test and intended-use analytical validation are different
milestones. In particular, a planned lane never means an executed or negative result.

## Chosen architecture

Use a common intake / identity / assessability boundary, then parallel complementary
analytical lanes, then the existing full-evidence and normalization layers.
A single winner or majority-vote consensus is explicitly excluded. Adding a second
independent monolithic pipeline would duplicate reference, QC and reporting rules;
running every program unconditionally would violate their different input assumptions.

| Lane | Role | First integration scope |
| --- | --- | --- |
| QDNAseq + ACE | depth, segments, fitted total CN and alternative purity/ploidy fits | reuse existing runtime and Block-4 evidence adapter; retain 100/500/1000 kbp |
| ONT-Spectre | long-read depth / SNV-informed CNV comparison | distinct provider `ont_spectre`; require compatible coverage bins and SNV evidence |
| ichorCNA | cfDNA/ULP-WGS CN and tumour-fraction research | separate cfDNA route; no automatic equivalence with adaptive-sampling or ONT tumour WGS |
| Sniffles2 | SV breakpoint candidates | reuse existing adapter; explicit germline/mosaic policy, no automatic somatic assertion |
| cuteSV | alternative SV breakpoint candidates | reuse existing module; retain disagreement and source-event identities |
| SAVANA | tumour SV plus CN; selected/alternative purity/ploidy and allelic CN | paired and tumour-only routes separate; CN requires SNP information |
| Severus | tumour/normal and complex-SV breakpoint graph | paired somatic route; unpaired calls stay candidates, not verified somatic calls |
| Wakhan | haplotype-specific CN, purity/ploidy, LOH-related evidence | phased SNP / haplotagged route; chosen SV-guided vs change-point segmentation explicit |

The original Spectre and ONT-Spectre fork are not interchangeable. A future second
Spectre provider must carry its own version, adapter and policy identities and a
shared method-family tag; it is not a second independent biological experiment.

## Routing rules

Inputs declare sample identity, reference fingerprint/build, platform, assay basis,
coverage and artifact fingerprints. Unknown coverage stays unknown. No missing tumour
fraction is replaced by ACE/SAVANA/Wakhan/ichor fitted purity, and no ploidy of two is
silently assigned. CNV coverage and breakpoint-support coverage are distinct validation
questions. A global depth check is a necessary engineering gate, not a detection guarantee.

Every selected caller requires a prospectively locked, matching assay/platform/build
policy with an explicit minimum coverage. Those are **study parameters**, not universal
clinical recommendations. Policies also bind the exact reference and assessability-mask
hash and the coverage measurement definition. Measured tumour fraction requires method
and timepoint; these declarations are not independently verified by the planner.
Paired somatic lanes additionally require a distinct normal sample, normal BAM identity
and a separately registered normal-depth criterion and measurement definition.

Whole-genome read-depth CNV lanes on adaptive-sampling data are blocked in this initial
planner. Unlocking a future AS route requires a separately validated off-/on-target
selection, callable territory and bias-correction design, not a generic force flag.
SV-only AS research remains possible with its own locked policy and mask; absence of
an event outside assessable territory is not a negative call.

ONT-Spectre requires SNV and coverage-bin artifacts from the same sample/reference and
an exact bin-size match. SAVANA CN requests require SNP or heterozygous-count evidence.
Wakhan's allele-specific route requires phased-SNV and haplotagged-BAM declarations;
SV-guided segmentation additionally requires explicit breakpoint evidence. These are
chosen integration preconditions, not assertions that every upstream optional mode
has the same requirements. Source artifacts must never be substituted across samples.
Normal-derived SNPs are allowed only when explicitly bound to the declared normal BAM.

## State and identity contract

A plan retains a row for every catalog caller, including disabled and blocked callers.
It records input blockers independently from adapter availability. Existing adapters
are not re-executed by this new planner; new callers remain `ADAPTER_PENDING` when
declared inputs are eligible, or `BLOCKED` with `existing_adapter=false` otherwise.
Every plan row has `execution_status=NOT_RUN`. No plan may claim `NO_CALL`, `FAILED`
or `COMPLETED`, since these are execution outcomes, not planning outcomes.

Each request locks caller version, runtime digest, adapter digest, parameter-policy
digest and study-policy digest. Unpinned `latest`, moving branch names and wildcard
versions are rejected. Catalog and plan hashes are canonical and contain no runtime
clock. Reordering requests or input records cannot change a plan; changing a sample,
reference, coverage, request, input digest or registered threshold must change it.
Revalidated inputs prevent unchecked model copies from bypassing contracts.
A study-registration digest is a declared reference in 5A, not verification of the
external registration content. Native execution must resolve and verify that content.

## Evidence and reporting boundary (subsequent adapter blocks)

All native outputs, fit alternatives, exclusions and raw quantitative values remain
available locally with fingerprints. Normalized events are additive. Event identifiers
must be namespaced by study, specimen, caller/version, run, resolution and native ID.
Caller-specific CN semantics, minor/major CN, BAF, purity and ploidy remain explicit;
copy-neutral LOH is not encoded as a deletion. A BND is not proof of an expressed fusion.

Reuse Block-2 evidence manifests and Block-3 matching/denominator logic. Do not create a
parallel metric implementation. Secondary resolutions require their own normalized
analytical lanes before per-resolution performance is claimed. Benchmark all callers
against the same prospectively locked orthogonal truth and assessability masks; caller
concordance is never used as truth.

Maintain a dependency graph for evidence: SAVANA CN can depend on SAVANA breakpoints;
Wakhan can depend on Severus breakpoints; SNV/phase and depth tracks can be shared.
Such dependence prevents interpreting multiple outputs as independent confirmations.
No vote count is converted into a probability, diagnosis, winner or automatic release.
Conflicts remain visible and traceable for human review.

## Native-adapter delivery sequence and acceptance

5A delivers only the catalog, declared-input checks and deterministic plan, with tests.
5B implements a pinned ONT-Spectre runtime + native import + Block-2/3 bridge.
5C implements a pinned SAVANA paired/tumour-only adapter; preserve all ranked CN fits.
5D implements Severus and then Wakhan, preserving breakpoint / phase dependencies.
5E implements the isolated ichorCNA cfDNA path and resource identities.
5F connects the reviewed adapters to existing engineering CLI/desktop controls and
reports. Standard/default execution stays unchanged until separately approved.

Every new native adapter requires: exact upstream version/API and license/resource
review; failing tests before code; safe local paths and atomic non-overwriting outputs;
exit/status/timeout handling; native-to-normalized identity checks; build/sample mismatch
and non-finite-value tests; real-tool smoke on permitted synthetic material; full
same-head repository CI. Biological performance requires the separate real-validation
block and does not follow from synthetic tests.

## Primary sources checked for this design

- https://github.com/tgac-vumc/ACE
- https://github.com/nanoporetech/ont-spectre
- https://github.com/broadinstitute/ichorCNA
- https://github.com/fritzsedlazeck/Sniffles
- https://github.com/tjiangHIT/cuteSV
- https://github.com/cortes-ciriano-lab/savana
- https://github.com/KolmogorovLab/Severus
- https://github.com/KolmogorovLab/Wakhan

These references establish candidate scope, not installed versions, clinical validity,
or permission to redistribute software/resources. Versions and license acceptance must
be locked during each native-adapter block. No third-party code/data is vendored here.

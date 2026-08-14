# Architecture decision log

## ADR-001: Keep computation local by default

**Decision:** Raw and derived human genomic data remain on approved institutional compute.

**Reason:** Genomic files are sensitive and potentially identifying. Cloud orchestration is
metadata-only unless explicitly approved.

## ADR-002: Use Snakemake for the initial workflow

**Decision:** Use a file-oriented Snakemake DAG for the first validated implementation, with
every scientific tool isolated behind an adapter.

**Reason:** Snakemake supports explicit file dependencies, isolated software environments,
container execution, provenance and testable modular rules. The decision is independent of any
specific caller or prior implementation.

**Revisit when:** The operational target requires native cloud orchestration, a centrally
maintained Nextflow platform, or evidence shows that another engine materially improves
reproducibility or supportability. Scientific adapters and contracts must remain portable.

## ADR-003: Use one normalized event model

**Decision:** HTML, Excel, JSON and ISCN proposals are rendered from the same validated
event model.

**Reason:** Parallel report logic causes silent discrepancies. One source of truth makes each
displayed result traceable to caller evidence.

## ADR-004: Treat ISCN as a proposal until validated

**Decision:** No automated notation can be released without expert review.

**Reason:** Syntax alone does not establish semantic correctness, clonality, uncertainty or
clinical validity.

## ADR-005: Treat the thesis as context, not specification

**Decision:** Lea Evers' master's thesis remains a historical comparison source only. Public
primary literature, official documentation, controlled benchmarks and the intended-use
statement govern architecture and tool selection.

**Reason:** A single exploratory dataset cannot establish general performance, clinical
thresholds or a durable software architecture. Preserving traceability is useful, but inheriting
its tools or parameters without independent validation would create avoidable anchoring bias.

## ADR-006: Select callers through assay-specific benchmarks

**Decision:** No CNV, SV or fusion caller is a production default until it passes a predeclared
benchmark for the relevant assay mode, coverage, tumor/blast fraction, reference build and event
class. Candidate lists are not rankings.

**Reason:** Published results are not interchangeable across germline, tumor-normal, tumor-only,
low-coverage WGS and adaptive-sampling settings. Caller agreement is an evidence feature, not a
truth label.

**Initial candidates:** Compare ichorCNA, QDNAseq + ACE and Spectre for copy number; compare
Sniffles2 with at least one cancer-aware or independent long-read SV method. Evaluate SAVANA and
Severus only in data regimes matching their assumptions. See `docs/EVIDENCE_BASE.md`.

## ADR-007: Integrate Sniffles2 first as conservative candidate evidence

**Decision:** Pin Sniffles2 v2.8.0 in the executable environment and normalize its PASS records,
but force every event to `unclassified`, `reportable: false` under the technical policy.

**Reason:** Sniffles2 has strong primary evidence and a stable long-read VCF interface, making it
useful for exercising the complete adapter and reporting boundary. Published cross-platform
performance does not constitute AML, tumor-only, low-coverage or adaptive-sampling validation.

**Privacy and interpretation constraints:** Request symbolic alleles, never request read-name
output, count rejected records, map an empty accepted set to `NO_CALL`, and never infer a fusion or
ISCN assertion directly from a BND record.

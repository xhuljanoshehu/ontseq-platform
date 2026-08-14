# Evidence base and tool-selection record

**Status:** living scoping review  
**Last searched:** 2026-08-14  
**Scope:** single-sample Oxford Nanopore analysis for hematologic malignancies, with emphasis
on low-coverage whole-genome copy number, adaptive-sampling structural variants/fusions,
reproducible reporting and expert-reviewed ISCN proposals.

This document is the scientific decision record for candidate selection. It is not a systematic
review, a clinical claim or a substitute for local analytical validation. Lea Evers' thesis is a
context source only; it neither outranks independent evidence nor fixes the implementation.

## Review method

The search prioritized peer-reviewed primary studies, benchmark studies, professional guidance,
public reference-material programs and official workflow/tool documentation. Search themes
combined `nanopore`, `long-read`, `acute leukemia` or `AML`, `low coverage`, `copy number`,
`structural variant`, `fusion`, `adaptive sampling`, `benchmark`, `validation` and `ISCN`.

Sources were included when they informed at least one of these decisions:

1. intended assay and sample model;
2. candidate algorithm or parameter range;
3. QC, no-call or reportability requirement;
4. benchmark dataset or scoring method;
5. clinical validation, reproducibility or nomenclature control.

For each production dependency, the repository must retain the exact software version, official
manual, parameter rationale, reference bundle, benchmark result and known limitations. A paper
can nominate a candidate; it cannot promote it to clinical use.

## Evidence levels

| Level | Meaning | Permitted use |
| --- | --- | --- |
| A | Peer-reviewed ONT study in AML/acute leukemia with clinically characterized specimens | Direct design input; local validation still required |
| B | Peer-reviewed long-read cancer or caller benchmark | Candidate selection and failure-mode design |
| C | Transferable low-coverage or clinical-pipeline evidence from another platform/disease | Hypothesis and benchmark design only |
| D | Official documentation, standard, preprint or draft public benchmark | Implementation reference; never a clinical performance claim |

Evidence level describes applicability to this project, not general publication quality.

## Evidence matrix

| Source | Level | Design and principal result | Applicability and limitation | Repository decision |
| --- | --- | --- | --- | --- |
| Heuser et al., *Blood Advances* (2026), [doi:10.1182/bloodadvances.2026019960](https://doi.org/10.1182/bloodadvances.2026019960) | A | Low-coverage long-read karyotyping in 100 diagnostic AML samples; reports high agreement, inter-laboratory reproducibility and a roughly 34-hour median turnaround time. | Closest disease and intended workflow evidence. Publication-level performance must be decomposed by event class, coverage, purity and truth method before reuse. | Treat as the primary external design reference for an AML lcWGS lane; reproduce only against a locked local truth set. |
| Geyer et al., *Leukemia* (2025), [doi:10.1038/s41375-025-02565-y](https://doi.org/10.1038/s41375-025-02565-y) | A | Adaptive-sampling ONT in 57 pediatric acute leukemias; reported 100% specificity and 96% sensitivity for genomic subtype, with both gross-karyotype misses occurring below 30% blasts. | Strong support for rapid integrated CNV/fusion analysis. Pediatric mix, custom analysis and high blast fractions limit direct transfer to adult AML. | Make tumor/blast fraction a mandatory manifest/QC field; validate explicit low-purity no-call thresholds. |
| Salmon et al., *Leukemia* (2026), [doi:10.1038/s41375-025-02801-5](https://doi.org/10.1038/s41375-025-02801-5) | A | Adaptive sampling of 240 genes in 20 hematologic cases detected all 12 known tyrosine-kinase fusions and resolved novel/complex partners. | Supports targeted long-read fusion resolution. Small selected cohort; a breakpoint outside the target design was missed and off-target CNV was exploratory. | Version the target BED as an assay component; report insufficient partner/breakpoint coverage as `NO_CALL`, never as negative. |
| Smolka et al., *Nature Biotechnology* (2024), [doi:10.1038/s41587-023-02024-y](https://doi.org/10.1038/s41587-023-02024-y) | B | Sniffles2 benchmarked long-read SV calling across ONT/HiFi, SV classes and 5-50x coverage, including mosaic calling. | Strong general-purpose SV evidence, but not an AML tumor-only clinical validation. | Retain Sniffles2 as a conservative candidate and provenance-rich evidence source, not a validated truth generator. |
| Aydin et al., *Scientific Reports* (2025), [doi:10.1038/s41598-025-92750-x](https://doi.org/10.1038/s41598-025-92750-x) | B | Compared eight long-read SV callers and caller combinations on paired cancer cell lines; higher-support combinations improved precision. | Motivates consensus evidence. Only two cell lines and no systematic coverage/purity study prevent universal thresholds. | Preserve every caller's support separately; benchmark evidence tiers rather than encode “two callers equals true.” |
| Elrick et al., *Nature Methods* (2025), [doi:10.1038/s41592-025-02708-0](https://doi.org/10.1038/s41592-025-02708-0) | B | SAVANA jointly analyzes somatic SVs, copy number, purity and ploidy from long-read cancer data and supports tumor-only analysis. | Cancer-aware and relevant for complex genomes, but evaluated at substantially higher coverage than the proposed lcWGS lane. | Add a research adapter for matched high-coverage or tumor-only evaluation; do not use as the default lcWGS caller. |
| Keskus et al., *Nature Biotechnology* (2025), [doi:10.1038/s41587-025-02618-8](https://doi.org/10.1038/s41587-025-02618-8) | B | Severus uses phased breakpoint graphs for somatic SV and complex rearrangement detection in long-read cancer genomes. | Useful independent cancer-aware method; tumor-normal assumptions and coverage must match the assay. | Evaluate as a secondary research candidate, especially when a matched normal is available. |
| Wang et al., *Briefings in Bioinformatics* (2025), [doi:10.1093/bib/bbaf514](https://doi.org/10.1093/bib/bbaf514) | C | Benchmarked six CNV methods from 0.1-10x; ichorCNA showed the strongest overall precision/runtime balance, with performance dependent on purity and coverage. | Directly relevant to low-pass CNV design, but largely based on short-read/downsampled data rather than ONT AML. | Add ichorCNA to the CNV benchmark beside QDNAseq + ACE and Spectre; predeclare depth/purity strata and do not select a default yet. |
| ONT EPI2ME, [`wf-human-variation`](https://epi2me.nanoporetech.com/epi2me-docs/workflows/wf-human-variation/) | D | Official workflow uses Sniffles2 for SV, Spectre or QDNAseq for CNV and accepts a single BAM. Documentation states a 20x minimum and recommends over 30x for the complete workflow. | Excellent compatibility and implementation reference; its full-workflow coverage assumptions do not fit a roughly 3x lcWGS assay. | Use as an interoperability baseline, not as a drop-in low-coverage solution. Keep a separate lcWGS profile. |
| NIST, [Cancer Genome in a Bottle](https://www.nist.gov/programs-projects/cancer-genome-bottle) | D | Public HG008/HG009 tumor-normal long-read resources; HG008 draft somatic SV/CNV benchmark v0.5 includes subclonal SV options. | Reproducible public technical benchmark, but draft and pancreatic rather than hematologic disease. | Use version-pinned HG008 for somatic engineering tests and HG002 for germline SV checks; retain AML validation as a separate gate. |
| Jennings et al., *Journal of Molecular Diagnostics* (2017), [doi:10.1016/j.jmoldx.2017.01.011](https://doi.org/10.1016/j.jmoldx.2017.01.011) | C | AMP/CAP consensus guidance for validating oncology NGS panels. | General clinical validation framework; not ONT-, CNV- or SV-specific. | Structure intended use, accuracy, precision, LoD, interference and reference-material studies by variant class. |
| Roy et al., *Journal of Molecular Diagnostics* (2018), [doi:10.1016/j.jmoldx.2017.11.003](https://doi.org/10.1016/j.jmoldx.2017.11.003) | C | AMP/CAP recommendations for validation and change control of NGS bioinformatics pipelines. | Written primarily around small variants, but its software-validation principles transfer. | Lock versions and fixtures; require validation-impact assessment and proportionate revalidation for every biological-output change. |
| Hastings, Moore and Chia, [ISCN 2024](https://karger.com/books/book/6011/ISCN-2024An-International-System-for-Human) and [2026 erratum](https://doi.org/10.1159/000549238) | D | Current controlled nomenclature reference plus published corrections. | Copyrighted standard requires authorized access; coordinate conversion alone cannot establish semantic validity. | Emit only a traceable `ISCN_PROPOSAL`; validate with authorized cases, current errata and cytogenetic expert review. |
| Snakemake, [deployment documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) | D | Official guidance for standardized workflow layout, environments, containers, testing and archival. | Supports reproducible engineering but does not validate scientific algorithms. | Keep Snakemake initially while making all adapters and result contracts workflow-engine neutral. |

## Locked implementation interfaces

| Dependency | Locked interface | Safety-relevant decision |
| --- | --- | --- |
| Sniffles2 | v2.8.0; official [CLI source](https://github.com/fritzsedlazeck/Sniffles/blob/master/src/sniffles/config.py), [VCF source](https://github.com/fritzsedlazeck/Sniffles/blob/master/src/sniffles/vcf.py) and [Bioconda recipe](https://bioconda.github.io/recipes/sniffles/README.html) | Explicit support/size/MAPQ parameters; symbolic PASS records; no read-name output; unvalidated candidates only. |
| samtools | v1.24 in `workflow/envs/aligned_bam.yaml` | Fail-closed quickcheck/header/dictionary/index gate before a scientific caller. |
| Cramino | v1.3.0 in `workflow/envs/aligned_bam.yaml` | Normalize aggregate QC only; do not copy source path or read-level data. |

An interface lock establishes reproducible execution, not analytical validity. Patch and minor
upgrades that can alter normalized events require regression results and validation-impact review.

## Working conclusions

### Assay lanes

- **Primary implementation lane:** one aligned BAM per run, with separate low-coverage WGS and
  adaptive-sampling profiles.
- **Later lane:** POD5/basecalling-to-report after the aligned-BAM path is reproducible.
- **Research lane:** high-coverage and/or matched tumor-normal somatic SV methods.

### Copy number

Benchmark ichorCNA, QDNAseq + ACE and Spectre on the same version-pinned datasets. Stratify by
coverage, tumor/blast fraction, event size and event class. Required metrics include precision,
recall/sensitivity, breakpoint tolerance, copy-number error, false-positive burden, no-call rate,
runtime and memory. No method or ACE penalty is currently selected.

### Structural variants and fusions

Use Sniffles2 as one conservative evidence source and compare it with at least one independent,
cancer-aware method. Store per-caller support, mapping quality, strand/orientation, breakpoint
confidence, local/partner coverage and target-design observability. Consensus raises or lowers an
evidence tier; it never establishes truth by itself. SAVANA and Severus remain research candidates
until coverage and input assumptions match the intended assay.

### QC and no-call policy

Tumor/blast fraction, global depth, per-target depth, both fusion-partner coverage, read length,
mapping quality, target-BED version and reference build are mandatory reportability inputs. A
module that was not run, a locus outside observable target space and a sample below a validated
gate must remain distinguishable from a biological negative.

### ISCN

The pipeline may create an evidence-linked ISCN proposal only. It must expose every source event,
assumption, uncertainty and unsupported construct. Clinical release requires authorized ISCN
materials, current errata, an expert-reviewed conformance corpus and human sign-off.

## Benchmark promotion gate

A candidate can be promoted only when all items are committed or recorded in the controlled
validation system:

- intended use and reportable event classes;
- immutable dataset identifiers, truth versions and inclusion/exclusion rules;
- software/container/reference/target checksums;
- predefined coverage and tumor/blast-fraction strata;
- predefined metrics and minimum acceptance thresholds;
- raw and normalized result retention with failure/no-call accounting;
- independent review of errors and clinically important discordances;
- validation-impact decision and rollback plan.

## Maintenance

Review this evidence base at least quarterly, before promoting a caller, and before changing a
caller, model, reference genome, target BED or reportability threshold. Add new sources with an
explicit applicability statement; do not silently replace prior evidence. Superseded decisions
remain in `docs/DECISIONS.md` and `CHANGELOG.md`.

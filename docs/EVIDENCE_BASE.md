# Evidence base and tool-selection record

**Status:** living scoping review  
**Last searched:** 2026-08-31  
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

## Candidate survey 2026-08-31: small variants, methylation and clinical annotation

A second search pass, prompted by the observation that the report can name *where* something
changed but not *what* it means, and that three whole classes of evidence the instrument
already produces are currently discarded. Candidates were additionally screened for whether
they are **packaged and pinnable**, because an unpackaged tool is a materially larger
provenance commitment than a conda pin.

### A coverage figure in this document needs correcting

The EPI2ME row above reasons against "a roughly 3x lcWGS assay". The measured local run
`260611_RAD114_AS_S700` produced **8.8x genome-wide off-target coverage**, not 3x. Several
conclusions below depend on which figure is right, so the assay's real off-target depth
should be confirmed across runs and this document corrected. Published adaptive-sampling work
reports comparable magnitudes — one clinical CNV study reports 28.4x on-target against 5.3x
off-target and explicitly uses the off-target background for non-targeted regions — so 8.8x
is not anomalous.

### Evidence matrix additions

| Source | Level | Design and principal result | Applicability and limitation | Repository decision |
| --- | --- | --- | --- | --- |
| Steinicke et al., *Nature Genetics* (2025), [doi:10.1038/s41588-025-02321-z](https://doi.org/10.1038/s41588-025-02321-z) | A | MARLIN: neural-network acute-leukemia classification from **sparse** nanopore methylation. Reference cohort n=2,540, 38 methylation classes. Retrospective nanopore concordance 25/26; real-time 5/5, typically within 2 h of sample receipt. | The most directly applicable new evidence found: same disease, same platform, and it consumes low-coverage genome-wide data rather than target coverage. Classes and cohort are not this laboratory's. | Evaluate as an **independent diagnostic axis** on data already generated. It is a classification complement, never a substitute for genetic findings, and its output must be a proposal under review like ISCN. |
| Shah et al., *medRxiv* preprint (2026), [doi:10.64898/2026.06.16.26355747](https://doi.org/10.64898/2026.06.16.26355747) | D | MOSAIC. Reports that on held-out **low-purity** specimens MARLIN was concordant in 7/10 and ALMA in 5/10; MOSAIC concordant in every case, including one at 1.4% blasts. | Preprint, not peer-reviewed, single group. But it characterises the failure mode of the method in the row above, which matters more than the headline accuracy. | Treat low blast fraction as the predeclared failure mode of any methylation classifier. Blast fraction is already a mandatory manifest field; make it a hard stratum for this lane and a `NO_CALL` trigger, not a footnote. |
| Zheng et al., *Nature Communications* (2025), [doi:10.1038/s41467-025-64547-z](https://doi.org/10.1038/s41467-025-64547-z) | B | ClairS-TO: deep-learning **tumor-only** somatic small-variant calling for long reads; ensemble of two networks trained for opposite tasks. Reported to outperform DeepSomatic, Mutect2, Octopus and Pisces across ONT, PacBio and Illumina. 50x ONT WGS in ~214 min on 24 cores. | The tumor-only design matches routine AML diagnostics, which usually has no matched normal. Not AML-specific, and the reported evaluation is at coverage far above this assay's off-target depth. | Highest-value candidate for the missing small-variant lane, scoped **to the panel**, where measured depth is roughly 80x (1.37 Gb over a 17.03 Mb design). Do not run it genome-wide at 8.8x. |
| Blätte et al., *Leukemia* (2019), [doi:10.1038/s41375-019-0483-z](https://doi.org/10.1038/s41375-019-0483-z) | C | getITD: dedicated FLT3-ITD detection and MRD monitoring, Needleman-Wunsch based, reporting ITD length, insertion site and VAF per clone. | Developed for **short-read** amplicon/hybrid-capture data. Transfer to 500 bp+ ONT reads is unestablished. | Do not adopt on the assumption that it transfers. Benchmark it against the SV lane on ONT reads before treating it as an FLT3-ITD method. |
| Multiclonal FLT3-ITD profiling on MinION, *OncoTargets and Therapy* (2025), [tandfonline.com/doi/full/10.2147/OTT.S526628](https://www.tandfonline.com/doi/full/10.2147/OTT.S526628) | B | ONT MinION FLT3-ITD profiling with a tailored clustering approach for subclonal detection; notes that general-purpose SV callers including Sniffles underrepresent minor clones. | Directly on-platform and on-target for a key AML marker, but the method is described rather than distributed as a maintained package. | Records that Sniffles2 alone is expected to **underrepresent FLT3-ITD subclones**. Treat low-VAF ITD as a known limitation of the current SV lane rather than an absence. |
| Somatic SV caller comparison in lung cancer, *BMC Genomics* (2024), [doi:10.1186/s12864-024-10792-3](https://doi.org/10.1186/s12864-024-10792-3) | B | ONT somatic SV benchmark: with minimap2, SAVANA 79.5% and Severus 79.25% recall, nanomonsv 72.5%. Runtimes Sniffles2 ~16 min, nanomonsv ~2 h, SAVANA ~4.8 h, Severus ~7.3 h. | Lung cancer, not AML, and tumor/normal designs. The runtime spread is a real operational constraint for a diagnostic turnaround. | Strengthens the existing SAVANA and Severus rows with concrete recall and cost figures. Note that both assume a matched normal, which this workflow does not have. |
| Geoffroy et al., AnnotSV, *Nucleic Acids Research* (2021), [academic.oup.com/nar/article/49/W1/W21/6281473](https://academic.oup.com/nar/article/49/W1/W21/6281473) | D | SV/CNV annotation from 20+ sources — genes, haploinsufficiency, triplosensitivity, known pathogenic and benign regions, regulatory elements — plus an ACMG/ClinGen-compliant 5-class ranking module. | **The ranking is a germline vocabulary.** ACMG/ClinGen SV classification answers a constitutional question; AML asks a somatic one. This is precisely the confusion ADR-022 already records for ClinVar. | Adopt the **annotation content**; do not adopt the ranking as a somatic classification. If the ranking is surfaced at all it must be labelled as germline-vocabulary, exactly as ClinVar assertions are. |
| Nakken et al., PCGR, *Bioinformatics* (2018), [academic.oup.com/bioinformatics/article/34/10/1778/4764004](https://academic.oup.com/bioinformatics/article/34/10/1778/4764004) | D | Somatic variant interpretation report engine: SNVs/indels, CNAs and fusions, classified by **oncogenicity and actionability**. Extends VEP annotations via vcfanno. Python/R, distributed via Docker. | The somatic counterpart to the germline vocabulary above, and it consumes exactly the event classes this pipeline produces. Not packaged on Bioconda. | Evaluate as the interpretation layer once a small-variant lane exists. Its vocabulary — oncogenicity, actionability — is the one this project should be attaching, and it is not ACMG germline. |

### Packaging audit against the live Bioconda index

Checked against `noarch` and `linux-64` repodata on 2026-08-31. Version shown is the newest
present. This determines cost, not merit: an unpackaged tool is not disqualified, it is a
larger provenance commitment.

| Purpose | Packaged and pinnable | Not on Bioconda |
| --- | --- | --- |
| Small variants, germline | `clair3` 2.0.2, `deepvariant` 1.10.0, `medaka` 2.2.2, `longshot`, `nanocaller` | — |
| Small variants, somatic tumor-only | — | **ClairS-TO**, ClairS |
| Somatic SV | `severus` 1.7, `nanomonsv` 0.9.0, `savana` 1.3.8 | — |
| SV comparators | `cutesv` 2.1.4, `svim` 2.0.0, `dysgu` 1.9.0, `delly` 2.6.0 | — |
| CNV comparators | `cnvkit`, `r-ichorcna` 0.5.1, `cnvpytor`, `wisecondorx` | **Spectre** |
| Methylation | `ont-modkit` 0.6.4, `methylartist` 1.5.4 | MARLIN/MOSAIC models |
| Phasing | `whatshap` 2.8, `longphase` 2.0.2 | — |
| Annotation | `ensembl-vep` 116.1, `annotsv` 3.5.10, `vcfanno` 0.3.9 | ClassifyCNV |
| Nomenclature | `hgvs` 1.5.7 | — |
| Knowledge access | `civicpy` 5.4.0 | OncoKB annotator |
| Clinical report | — | **PCGR** |
| FLT3-ITD / repeats | `getitd` 1.5.17, `straglr`, `trgt`, `tandem-genotypes` | — |

Note that `Spectre`, referenced in the EPI2ME row above, is **not** available on Bioconda.
Any Spectre comparison therefore carries the same packaging cost as ClairS-TO and should not
be treated as the cheaper option.

### Working conclusions from this pass

**The largest gap is small variants, not another CNV caller.** The pipeline cannot see NPM1,
FLT3-ITD, CEBPA, TP53 or the myelodysplasia-related genes, which is most of what decides AML
classification and risk. Ten of the twenty-four criteria drafted in
`GUIDELINE_CRITERIA_DRAFT_v0` are unevaluable for exactly this reason. Adding a fifth CNV
caller does not move that.

**Refined 2026-09-01, and the refinement matters.** The sentence above conflates two gaps
that have different owners. Checking the criteria bundle against the shipped panel BED --
now automated in `ontseq_platform.panel_reachability` and reported in
[`PANEL_REACHABILITY.md`](PANEL_REACHABILITY.md) -- separates them:

* **A software gap, for three of the seven small-variant criteria.** NPM1 (43.1 kb), FLT3
  (117.3 kb) and RUNX1 (281.5 kb) are already panel targets and already sequenced at roughly
  80x on every run. Nothing is missing from the assay for these; no caller is wired in. Note
  that FLT3-ITD is a tandem duplication rather than a small variant, so a SNV/indel caller
  alone does not deliver it.
* **A design gap, for the other four.** CEBPA and TP53 are not targeted at all, and eight of
  the nine myelodysplasia-related genes are absent -- only RUNX1 is present. At 8.8x
  off-target, no caller makes these criteria evaluable. Ten genes would have to be added to
  the design: ASXL1, BCOR, CEBPA, EZH2, SF3B1, SRSF2, STAG2, TP53, U2AF1, ZRSR2.

The partially covered myelodysplasia-related criterion is the one to watch. Evaluated
against the single gene that happens to be targeted it can report no myelodysplasia-related
mutation having never looked at the other eight, and that error runs towards favourable.
`panel_reachability` therefore reports partial coverage as unreachable rather than folding it
in with the reachable criteria.

Extending an adaptive-sampling design costs no reagent per added region; it costs depth,
because the same yield is divided across a larger enriched territory. That is a laboratory
decision about assay design, and nothing in this repository makes it.

**Methylation is evidence already being paid for and discarded.** The instrument emits 5mC
from the same reads, and the published classifier consumes sparse genome-wide data rather
than deep target coverage — which is the shape of the off-target fraction this assay already
produces. No additional sequencing, no additional run.

**Depth should be spent where it exists.** Roughly 80x on-target supports variant calling;
8.8x off-target supports copy number and methylation and supports neither confident somatic
SNV calling nor a complete karyotype. Any lane added should declare which depth regime it
belongs to.

**Two germline vocabularies are now on the shortlist.** AnnotSV's ACMG/ClinGen ranking joins
ClinVar as a source whose *content* is useful and whose *classification* answers a different
question than AML asks. The existing ADR-022 boundary applies unchanged.

## Survey 2026-09-01: external validation, the pore-split architecture and the tool stack

Two sources were supplied directly (Heuser et al.; Schoenung et al.) and the rest were found by
searching from them. Full texts could not be retrieved: `biorxiv.org`, `ashpublications.org` and
`europepmc.org` are blocked by this environment's network policy, so every row below is graded
from abstract-level content only. Read the methods before acting on any of it.

### Evidence matrix additions

| Source | Level | What it shows | Limitation | What it changes here |
| --- | --- | --- | --- | --- |
| Heuser et al., *Blood Advances* (2026), [doi:10.1182/bloodadvances.2026019960](https://doi.org/10.1182/bloodadvances.2026019960) | A | ONT lcWGS karyotyping in 100 AML samples (50 retrospective adverse-risk, 50 prospective de novo): 93% sensitivity, specificity and accuracy; complex-karyotype AUC 0.971; cross-laboratory reproducibility R=0.99; lcWGS complex karyotype predicted shorter OS and RFS; ~34 h to bioinformatics. Clone-size estimates correlated only moderately with conventional cytogenetics (R=0.54). | Pure lcWGS, not adaptive sampling, so the read distribution differs from this assay. Clone-size correlation is the weakest reported metric. | The external benchmark this repository lacked. Its 93% sensitivity and 0.971 AUC are the obvious acceptance criteria for the CNV benchmark gate. R=0.54 is independent empirical support for refusing to report clonal fraction as a quantitative result. |
| Capilla-Guerra et al., *Blood* (2025), ASH abstract, "Rapid diagnosis of acute leukemia with integrated epigenetic and genetic profiling" | D | Splits the flow cell: **90% of pores to adaptive sampling** over 274 genes, **10% to conventional non-adaptive sequencing** for genome-wide methylation classification and copy number. Median 133x over ROIs (75-192x); NPM1 105x, TP53 113x, FLT3 163x. 27/31 clinically reported SNVs/indels recovered, including NPM1 p.W288Cfs*12, FLT3-ITD, TP53 and IDH1; four fusions; copy number concordant including del(5q)/-7/del(12p)/del(17p). Methylation classification confident in 10/10 within two hours. Stack: Dorado, minimap2, **Clair3**, QDNAseq, Sniffles. | Conference abstract, n=10, retrospective. | The most architecturally decisive row. It answers the depth-allocation question with a pore split rather than a reliance on off-target reads, and its tool stack is this repository's stack plus Clair3. |
| Steinicke et al., *Nature Genetics* (2025), "Rapid epigenomic classification of acute leukemia" (MARLIN) | A | Reference cohort n=2,540, 38 methylation classes; neural network classifying from **sparse** nanopore methylation; 25/26 retrospective cases concordant; 5/5 real-time cases, typically within 2 h of sample receipt. | Sparse-profile classification; rare classes underrepresented. | Confirms methylation classification as a first-class output of the same reads, not a separate assay. |
| Steinicke et al., *Blood* (2024), ASH abstract | D | The coverage figure the peer-reviewed paper's abstract omits: confident predictions in 18/19 nanopore samples at **~3x genome coverage**, and good performance from only 3% of the input data (r=0.93). | Conference abstract. | Answers the open question directly. This assay's measured 8.8x off-target is roughly threefold above the depth reported as sufficient, so the methylation lane is viable on reads already produced. |
| Achterberg et al., *Blood* (2025), ASH abstract (Lamprey) | D | Reference atlas >5,400 arrays, 38 subtypes spanning AML (18), BCP-ALL (10), T-ALL, low-risk MDS, MPN, JMML, CMML, BPDCN, B-PLL, CLL and two control classes. Trained on **simulated nanopore-style sparse reads**; micro-F1 0.96 on hold-out; 47/52 correct on a retrospective adaptive-sampling cohort. Confidently wrong on near-haploid/low-hypodiploid ALL, PAX5-altered ALL and KMT2A-PTD AML. | Conference abstract; named failure modes are confidently wrong, not merely uncertain. | Validated on adaptive-sampling data specifically. The enumerated high-confidence errors are a ready-made no-call list. |
| Schoenung et al., bioRxiv (2026), DOI 10.64898/2026.07.02.735835 | D | Hierarchical classification from whole-genome nanopore sequencing: 5,420 training samples, 21 entities including healthy controls, then 44 epitypes, then integration of genetic data. Reports diagnosis-defining alterations missed by standard-of-care work-up. | Preprint, not peer reviewed. Licence is `cc_no`, so reuse permission is unstated. Coverage requirement not given in the abstract. | The broadest entity coverage found, and the closest match to a single-assay ambition. Coverage requirement and model availability must be read from the methods before this can be planned against. |
| Marchi et al., *Nature Communications* (2025) | A | Acute Leukemia Methylome Atlas over 3,314 samples and 11 harmonized cohorts; a genome-wide prognostic model and a **targeted 38-CpG panel** both predict five-year survival; specimen-to-result nanopore protocol. | Rare karyotypes limited by training data. | Methylation carries prognostic and not only diagnostic signal, and a 38-CpG panel is a far smaller target than a genome-wide profile. |
| crossNN, *Nature Cancer* (2025) | B | Explainable neural framework classifying tumours from **sparse methylomes across platforms and coverage depths**; pan-cancer model over 170+ tumour types validated in >5,000 tumours including nanopore; 97.8% pan-cancer precision. | Not haematology-specific. | Cross-platform and coverage-agnostic by design, which is the property a classifier needs to survive a change of chemistry or depth. |
| Abel et al., *Journal of Molecular Diagnostics* (2025) | A | Tumour-only targeted long-read analysis of 26 AML/MDS samples at mean 52x: SNVs >96% recall and 91% precision; **indels 66% recall and 42% precision**, worst where few phased reads are available; copy number 95% accurate; all recurrent structural variants detected with no false positives. | n=26; targeted analysis of WGS rather than adaptive sampling. | The number that quantifies the NPM1 problem. A 4 bp insertion is exactly the class with 42% precision, so an NPM1 call from this assay needs orthogonal confirmation before it is reported. |
| Aganezov et al., AACR (2026) abstract, Oxford Nanopore | D | Tumour-only adaptive sampling profiled as 100-200x on-target with 8-15 kb reads and 5-15x off-target at ~500 bp; genome-wide copy number derived **from the off-target reads**; SNV recall high to 0.05 allele fraction and precision high to 0.10-0.20; structural variants recovered with only one breakend on target. Benchmarked on COLO829 and an in-silico synthetic genome across coverage and purity ladders. | Vendor abstract; cell lines, not patients. | The read profile matches this laboratory's measured run almost exactly (80x on-target, 8.8x and 536 bp off-target). Establishes COLO829/COLO829BL and in-silico purity ladders as the benchmark design, and confirms copy number should be derived from off-target reads. |
| Kato et al., *Journal of Clinical Oncology* (2024), ASCO abstract | D | Adaptive sampling of 466 genes on GridION in 28 paediatric leukaemias, mean on-target depth 21x. Against short-read WGS: **60.9% of SNVs, 17.6% of small indels, 89.2% of structural variants**, with poor efficiency at low allele fraction. Chromosome-level and focal copy number both detectable from off-target reads. An NPM1 frameshift outside known hotspots was missed. | Conference abstract; 21x is well below this assay's on-target depth. | The failure mode at shallow on-target depth, and a second independent warning that small indels are the weak class. |
| Martinez-Serra et al., *OncoTargets and Therapy* (2025) | B | A dedicated clustering pipeline for FLT3-ITD on MinION resolves multiclonal ITD architecture and duplications as short as 15 bp. **Sniffles failed to call several biologically validated ITDs** that the clustering approach found. | Small cohort; bespoke pipeline. | Direct evidence that this repository's Sniffles-based structural-variant lane will not deliver FLT3-ITD. A dedicated ITD method is required, not a caller swap. |
| ClairS-TO, *Nature Communications* (2025), [github.com/HKU-BAL/ClairS-TO](https://github.com/HKU-BAL/ClairS-TO) | B | Long-read **tumour-only** somatic small-variant calling; ensemble of two networks trained for opposite tasks; outperforms DeepSomatic and smrest on ONT and PacBio, and Mutect2/Octopus/Pisces on short reads; evaluated across coverage, allele fraction and purity. | Not AML-specific; not packaged on Bioconda. | Remains the leading tumour-only candidate. Note that Capilla-Guerra used plain Clair3, which is packaged, so the cheaper first step differs from the best endpoint. |
| ClairS, *Nature Methods* (2026) | B | Tumour-normal long-read somatic calling. On ONT Q20+ HCC1395/HCC1395BL at 50/25x: F1 89.83% for SNVs and 73.38% for indels, improving to 96.19% and 79.67% with real cell lines added to training. Read phasing identified as the key mechanism at low allele fraction. | Requires a matched normal, which routine AML diagnostics usually lacks. | Quantifies the ceiling a matched normal would buy, and confirms indels trail SNVs by roughly fifteen F1 points even with one. |
| SAVANA, *Nature Methods* (2025) | B | Somatic structural variants and copy-number aberrations at single-haplotype resolution **and estimation of tumour purity and ploidy**, with or without a germline control; 99 tumour-normal pairs; 13- and 82-fold higher specificity than the next two methods. | Whole-genome design; behaviour on adaptive-sampling coverage unknown. | A published long-read implementation of exactly the purity-and-ploidy problem `ontseq_platform.quantitation` models arithmetically, and it does not require a matched normal. |
| Severus, *Nature Biotechnology* (2025) | B | Breakpoint-graph somatic structural-variant caller for long reads; supports unbalanced karyotypes and complex multi-break patterns; highest F1 against Sniffles2, nanomonsv and SAVANA on a multi-technology cell-line panel; found cryptic rearrangements missed by standard panels in paediatric leukaemia. | Uses a matched normal. | A candidate for the structural-variant benchmark, and evidence that Sniffles2 is not the ceiling. |
| DeepSomatic, bioRxiv (2024) | D | Somatic SNV and indel calling for short and long reads with tumour-normal, **tumour-only** and FFPE modes. Publishes five matched tumour-normal cell-line pairs sequenced on Illumina, PacBio HiFi and ONT **with benchmark variant sets**. | Preprint. | The openly available truth data matters more here than the caller: it is material for the in-silico dilution series without new patient samples. |
| Nakamura et al., *npj Genomic Medicine* (2024) | B | Adaptive-sampling workflow over 33 genomes; SNV accuracy comparable to short reads; complex structural variants resolved. **Off-target reads, normally discarded, genotyped common SNPs genome-wide** well enough to compute a polygenic risk score; allele-specific promoter hypermethylation detected. | Germline cancer predisposition, not somatic leukaemia. | Establishes that the off-target fraction supports genome-wide SNP genotyping, which is the input a B-allele-frequency route to tumour fraction would need. |
| Furtado et al., Research Square (2026) | D | Nanopore WGS at 6.7x and adaptive sampling at 8.5x detected 94% of clinically relevant copy-number alterations in Wilms tumour with no false positives, but the depth was **insufficient for methylation assessment at the 11p15 imprinting control regions**. | Preprint; xenografts; n=15. | The counterweight to the methylation optimism above. Genome-wide sparse classification is not the same task as locus-specific methylation, and this assay's depth may support the first and not the second. |
| Hansen et al., *Journal of Molecular Diagnostics* (2023) | B | ONT whole-genome cytogenomics in mantle cell lymphoma: ~99% copy-number reproducibility between replicates at 100 kb resolution and 98% concordance with Illumina, from 1.5-7.5 million long reads. | Cell line; lymphoma. | A concrete resolution and read-count target for copy-number reproducibility. |
| Yang et al., *Genome Biology* (2025) | B | Benchmark of six adaptive-sampling tools across enrichment and depletion tasks; 1.50-4.86-fold coverage enrichment; basecalling plus minimap2 alignment the most accurate read-classification strategy. | Not cancer-specific. | Relevant if the selection strategy itself is ever revisited; the enrichment factors bound what a panel change can deliver. |
| Geoffrion et al. (2025), nf-core-oncoseq | D | Adaptive-sampling whole-genome workflow for paediatric oncology over 31 samples, unifying genomic, structural and epigenomic detection, with an **open-source pipeline**. Reports clonal alterations confidently supported within the first sequencing day. | Preprint; paediatric solid tumours. | An existing open pipeline covering the same assay shape; worth reading before building further lanes. |
| ROBIN, *Neuro-Oncology* (2025) | C | Single nanopore assay giving intraoperative methylome classification plus next-day SNV, CNV and structural-variant profiling; 50 prospective cases; 90% concordance with the final integrated diagnosis. | CNS tumours. | The staged-output pattern -- fast classification first, comprehensive profile later -- transfers directly to an AML workflow. |
| Kuschel et al., *Neuropathology and Applied Neurobiology* (2022) | C | Shallow nanopore methylation classification across 46 brain-tumour types: **1,000 random CpG features sufficed** for high-confidence classification; cross-laboratory concordance 10/11; 100% specificity in validation at a calibrated confidence threshold; median 21.1 h. | CNS tumours. | The methodological precedent for calibrated confidence scoring and an explicit unclassifiable outcome, which is the shape a no-call policy for a classifier needs. |

### Working conclusions from this pass

**The depth-allocation problem has a published answer, and it is not a software change.**
Capilla-Guerra et al. devote 90% of pores to adaptive sampling and 10% to conventional
sequencing, so the genome-wide lane is uniform rather than being reconstructed from off-target
reads. That removes the enrichment bias from copy number and methylation at the cost of a
tenth of the on-target yield. It is a sequencing-setup decision for the laboratory, not an
analysis decision, and it should be evaluated against the current off-target-only approach.

**The missing caller may be smaller than assumed.** The reference stack is Dorado, minimap2,
Clair3, QDNAseq and Sniffles. This repository already pins four of those five. Clair3 is
packaged on Bioconda and is a germline caller pressed into service; ClairS-TO is the
tumour-only design and is not packaged. That is a cheap first step and a better endpoint, and
they are not the same step.

**Small indels are the weak class, twice measured.** 66% recall and 42% precision at 52x, and
17.6% of small indels at 21x. NPM1's canonical alteration is a 4 bp insertion. An NPM1 call
from this assay is a hypothesis requiring orthogonal confirmation, not a reportable result.

**FLT3-ITD needs a dedicated method.** Sniffles demonstrably misses validated ITDs that a
clustering approach recovers, including duplications as short as 15 bp. No swap among general
structural-variant callers fixes this.

**Methylation classification is plausible at this assay's off-target depth, with a caveat.**
Confident classification is reported at ~3x genome coverage and from 1,000 CpG features, well
below the 8.8x measured here. But locus-specific methylation at 8.5x was reported as
insufficient in Wilms tumour. Genome-wide sparse classification and locus-specific methylation
are different questions and must not be promised together.

**Purity and ploidy are solved elsewhere and can be borrowed.** SAVANA estimates both from
long reads without a matched normal. `ontseq_platform.quantitation` carries the arithmetic and
its detection limits; it does not need to carry the estimator.

**Benchmark material exists and needs no patients.** COLO829/COLO829BL, HCC1395/HCC1395BL and
the five DeepSomatic cell-line pairs are published with truth sets, and in-silico purity and
coverage ladders are the established design. This is the in-silico dilution series, already
standard practice rather than a novel proposal.

**Acceptance criteria are now available.** 93% sensitivity and 0.971 complex-karyotype AUC
against conventional karyotyping; ~99% copy-number reproducibility at 100 kb between
replicates and R=0.99 across laboratories. These are external numbers to be measured against,
not targets this project set for itself.

## Locked implementation interfaces

| Dependency | Locked interface | Safety-relevant decision |
| --- | --- | --- |
| Sniffles2 | v2.8.0; official [CLI source](https://github.com/fritzsedlazeck/Sniffles/blob/master/src/sniffles/config.py), [VCF source](https://github.com/fritzsedlazeck/Sniffles/blob/master/src/sniffles/vcf.py) and [Bioconda recipe](https://bioconda.github.io/recipes/sniffles/README.html) | Explicit support/size/MAPQ parameters; symbolic candidates with PASS-only normalization; no read-name output; unvalidated candidates only. |
| cuteSV | v2.1.3; official [source](https://github.com/tjiangHIT/cuteSV) and [Bioconda recipe](https://bioconda.github.io/recipes/cutesv/README.html) | Independent candidate evidence with locked clustering/support parameters; caller agreement raises technical priority but never creates truth or reportability. |
| SV annotation resources | GENCODE v19 for GRCh37, GENCODE v50 for GRCh38 and build-matched UCSC cytobands | Original and normalized resources are locally checksum locked. Build disagreement fails closed; annotation cannot be inferred from an expected fusion name. |
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

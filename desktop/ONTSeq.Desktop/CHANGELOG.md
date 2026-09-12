# ONTSeq Desktop changelog

## 0.8.1-engineering — 2026-09-09

- Restore the Core report API, coverage, structural-variant, ISCN and resource details.
- Persist failed report/adapter stages and invalidate stale release markers on a retry.
- Synchronize current application/versioned runtime/resource defaults; preserve saved
  custom installations and independently versioned reference, panel and tool assets.
- Engineering repair candidate only; no new clinical accuracy claim.

## 0.7.1-engineering — 2026-09-09

- Distinguish preparation errors from a running BAM reader failure. Report the exact
  unavailable drive, directory, reference or runtime component instead of a combined
  prerequisite message. Do not advise reinstalling references for drive access failures.
- Explain the methylation probe outcome with a concrete German reason and progress.
- Add optional thorough BAM inspection with asynchronous status and cancellation.
- Invalidate stale results on selection changes and refresh the actual input before
  requesting explicit methylation opt-in. Keep full negatives distinct from uncertainty.
- Coordinate Desktop, Core and packaged local workspace version 0.7.1.

## 0.7.0-engineering — 2026-09-08

- Unified genome and optional methylation analysis on the current reference/profile core.
- Explicit modified-base inspection and user choice before a run, with unknown and
  unavailable states kept distinct from absent methylation information.
- Guided local workspace and coordinated 0.7.0 Core/Desktop identity.
- Research use only; no automatic clinical release or new analytical validation claim.

## Unreleased

- Show each stage's measured duration (persisted `duration_seconds` provenance) next
  to its status in the run list and in the timeline tooltips; absent durations stay
  empty instead of showing a placeholder. Presentation only; no backend, contract
  or provenance change.
- Show run progress additionally as a horizontal stage timeline in the Desktop panel:
  one node per stage with German status caption and tooltip. NO_CALL stays neutral
  amber, NOT_RUN neutral gray, only FAILED is red, and stages without a report yet
  stay hollow. Presentation only; stage semantics, provenance and backend status are
  unchanged.
- Default new resource installations to `~/.local/share/ontseq/resources-v0.7.1`, derived
  from the same release constant as the Desktop label and runtime installer. Keep saved
  custom roots selected; relocating an existing installation requires an explicit local
  migration. Reference/panel/knowledge versions and checksums retain their own identities.

## 0.6.2-engineering

- add an exact UCSC-hg19 Canonical-25 Adaptive-Sampling path alongside native full GRCh37.p13;
  bind each BAM dictionary to its own pinned FASTA/FAI with no runtime conversion or fallback
- expose native full-GENCODE-19 and UCSC-hg19-Canonical-25 Adaptive-Sampling choices for
  GRCh37, each bound to the dedicated `AML_AS_111_GRCh37_v1` panel
- label the GRCh37 choices as `110/111 mapped`; record 109 reciprocal-exact intervals, use 108
  native GENCODE-19 analysis ROIs and keep ACACA roundtrip review plus the CT45A2/IGH/GPR128
  mapping/ROI gaps visible
- validate the forward and reciprocal mapping artifacts and display source/final/ROI counts in
  the Desktop-backed HTML, XLSX and Live Workspace reports
- package Desktop and Core together as 0.6.2, require a version-bound isolated runtime prefix
  and keep existing 0.6.0/0.6.1 runtime and resource roots unchanged

## 0.6.1-engineering

- expose the regular Core report containing a traceable, partial CNV-derived ISCN proposal;
  keep expert review mandatory and automatic clinical release disabled
- prevent a second Desktop instance from connecting to an older service on the same configured
  port; select a free loopback port without stopping or adopting the existing process
- authenticate each spawned service by a fresh 128-bit launch ID plus exact resource, output and
  input-root scope; retry a detected bind race without publishing or reusing the failed candidate
- bind Live-Workspace navigation to the exact selected profile and fail closed unless the
  backend advertises that profile in `/api/config`; reject results whose profile binding is absent

## 0.6.0-engineering

- Add the native GRCh37.p13 lcWGS profile alongside the four isolated GRCh38 profiles.
- Add a separate UCSC hg19 Canonical-25 lcWGS choice with visible `chrM=16571` contract;
  keep it isolated from the complete native GENCODE 19 dictionary.
- Install/repair references by build and check the selected family without requiring the other.
- Open the authenticated local live workspace backed by the same canonical pipeline.
- Preserve unsigned RUO status, explicit failures and the prohibition on silent liftover.

## 0.5.3-engineering

- package `HEMATOLOGY_v3` with source-attributed pathology/DOID associations and the redesigned HTML/XLSX review reports
- show fusion assessment and reviewable rearrangement evidence without treating a missing
  analytical benchmark as a biological negative result
- include Adaptive-Sampling target coverage in persisted reports and retain cuteSV provenance
- keep all four GRCh38 dictionary-contract profiles and their input gates unchanged

## 0.5.2-engineering

- add explicit `AML_LCWGS_GRCh38_CANONICAL25` and
  `AML_AS_111_GRCh38_CANONICAL25` choices for BAMs whose ordered dictionary is exactly
  `chr1`-`chr22`, `chrX`, `chrY`, `chrM`
- keep `AML_LCWGS_GRCh38` and `AML_AS_111_GRCh38` on the existing `exact_full`
  Primary-Assembly dictionary contract
- derive both dictionary contracts from the same installed, checksum-pinned GRCh38 bundles;
  no second multi-GiB download, liftover or silent profile fallback is introduced
- keep the existing full profiles and default selection unchanged; scientific policies and
  interpretation thresholds are unchanged

## 0.5.1-engineering

- resolve every packaged Core configuration independently of the application working directory
  and pass absolute cuteSV, consensus and evidence-policy paths to profile services
- include cuteSV 2.1.3 in the WSL runtime and require the complete SV/CNV policy and executable
  contract during installation and before starting a Desktop profile
- reject duplicated packaged configuration trees and smoke-test the relocated service from an
  unrelated directory; scientific thresholds and interpretation policies are unchanged

## 0.5.0-engineering

- add `resourceRootWsl` with user-writable `~/.local/share/ontseq/resources` Desktop default,
  while Core retains `/opt/ontseq`; bridge bundle status, install and repair to the backend CLI
- replace independent build/assay selection with the build-isolated `AML_LCWGS_GRCh38` and
  `AML_AS_111_GRCh38` profiles; lcWGS is the default
- show the GRCh38 profile expectation before execution and the backend/provenance-detected build
  when it becomes available
- retain explicit Reference-Lock and Adaptive-Sampling BED settings as legacy-readable fields
  for one compatibility release without mixing them into profile-backed runs
- use fast presence/size/pin checks at interactive start, advertise only locally resolvable
  profiles, and let Core return the canonical sample-plus-UTC run ID; full SHA256 validation
  remains available through `ontseq references validate`
- make **Reparieren** restore the full GRCh38 profile-resource family (Reference, HEMATOLOGY,
  AML AS panel and both profiles) through the backend transaction; no manual bundle deletion is
  required

## 0.4.1-engineering

- advance Core, Desktop, WSL runtime directory and Windows bundle identity together to 0.4.1
- integrate the versioned Sniffles2/cuteSV consensus, annotation, observability and AML-priority
  evidence layer into reviewer reports
- retain the unsigned Research Use Only boundary; technical priority remains distinct from
  analytical or clinical validation and all automated SVs remain non-reportable

## 0.4.0-engineering

- advance Core, Desktop, WSL runtime directory and Windows bundle identity together to 0.4.0
- allow an operator to remove a configured Adaptive Sampling BED without deleting the
  original source file
- decouple the Research Use Only disclaimer from the measurable QC verdict while preserving
  the disclaimer and reportability boundaries in reports
- harden the local service against review-path traversal, concurrent job claims, unsafe BAM
  symlinks, unreadable directory entries and untrusted loopback Host headers
- preflight Mosdepth, target BED, target-coverage policy and selected component policies
  before an Adaptive Sampling run takes a run lock
- enforce the current release number across Desktop code, runtime installation, CI artifact
  naming and operator documentation with the repository version-consistency gate
- retain the unsigned Research Use Only boundary; no analytical or clinical validation is
  claimed by this engineering release

## 0.3.5-engineering

- add the pinned open-source `QDNAseq.hg38` 1.2.0 annotation package to the bundled runtime, locked to upstream commit `cf7c07e39de0ac64a9c38cb030cba4626e2aae83`
- real-tool test the canonical QDNAseq+ACE CNV path on a deterministic whole-genome-shaped GRCh38 fixture at 100/500/1000 kbp
- require recovery of the synthetic chromosome 7 loss at CN~1 and chromosome 8 gain at CN~3 and verify content-addressed resume
- request CNV for GRCh38 Desktop runs as well as GRCh37 after the GRCh38 real-tool gate passed
- preserve exact hg38 annotation source identity in the packed runtime and verify it after relocation
- record the pinned hg38 source commit in the Windows bundle build information and verify it before artifact upload
- keep the existing QDNAseq 1.42.0 / ACE 1.24.0 multi-bin profile and the Research Use Only boundary
- do not claim cohort-level analytical or clinical validation; both GRCh37 and GRCh38 still require validation on controlled real/reference samples

## 0.3.4-engineering

- combine Core 0.3.4 with the working Desktop/reference-integrity hardening from v0.2.2
- report Desktop and assembly version 0.3.4 and install the bundled runtime under `runtime-v0.3.4`
- require canonical GRCh37/GRCh38 reference identity during setup and again before each Desktop run
- create content-addressed Reference-Locks from the selected FAI SHA-256 and publish new locks atomically
- preserve the previously active reference lock if a newly generated lock fails validation
- bind aligned-BAM intake resume identity to BAM, declared BAI, reference lock, pipeline version and git commit
- verify the explicitly declared BAI with `samtools idxstats -X` rather than allowing implicit index discovery
- pass Core 0.3.4 target-coverage policy and component selection through the Desktop service path
- embed the exact runtime git commit in the packed Linux backend
- package a relocatable Linux/R runtime and test it on stock Ubuntu without host Python or R
- build and test the self-contained Windows x64 WPF shell before creating the engineering bundle
- retain the RUO boundary; GRCh38 QDNAseq/ACE CNV and cooperative cancellation remain deliberately unavailable

## 0.1.3-engineering

- support both common BAM index names: `sample.bam.bai` and `sample.bai`
- prefer `sample.bam.bai` deterministically when both files exist
- validate the selected index again inside the local service boundary
- add Windows and backend regression tests for index selection

## 0.1.0-engineering

- initial .NET 10 WPF operator shell
- BAM selection and sample-ID suggestion
- GRCh37/GRCh38 and lcWGS/Adaptive Sampling profiles
- automatic WSL backend preflight and launch
- local API start/status integration
- live stage status from persisted run provenance
- HTML/XLSX/result-folder launch buttons
- RUO boundary and disabled cancellation pending backend support

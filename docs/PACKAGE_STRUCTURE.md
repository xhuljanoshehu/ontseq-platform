# Package structure

This page maps `src/ontseq_platform/` by responsibility and states the dependency rules
between the parts. The rules marked **enforced** are checked by `tests/test_architecture.py`
on every test run; changing one means editing that test in review, with the reason next to
the allowlist entry.

For the scientific boundaries (evidence before interpretation, fail closed, no automatic
clinical release) read [`ARCHITECTURE.md`](ARCHITECTURE.md) first. This page is only about
where code lives and what may import what.

## Layers

| Layer | Responsibility | Modules |
| --- | --- | --- |
| Contracts | Versioned, validated data shapes shared by everything else | `models.py`, `iscn_syntax.py`, `schemas/` (exported by `scripts/export_schemas.py`), `io.py`, `coordinates.py` |
| Inputs, references and resources | Decide what a run is allowed to read: BAM intake, reference locks, bundles, panels, knowledge, caches | `bam_intake.py`, `bam_resolution.py`, `reference.py`, `reference_catalog.py`, `resource_registry.py`, `resource_bootstrap.py`, `resource_commands.py`, `annotation_cache.py`, `panel_bundle.py`, `panel_compiler.py`, `panel_reachability.py`, `model_lock.py`, `path_safety.py`, `sidecars.py` |
| Tool adapters | One external tool each: build a shell-free command, lock the version, normalize output into a typed report | `align.py`, `basecall.py`, `qc.py`, `target_coverage.py`, `sniffles.py`, `cutesv.py` (+`cutesv_build.py`), `methylation.py` (+`modkit_build.py`, `methylation_probe.py`, `_methylation_worker.py`), `cnv/qdnaseq.py`, `cnv/ichorcna.py`, `cnv/spectre.py`, `sv/severus.py`, `tumor/savana.py`, `tumor/wakhan.py`, `gatk_mutect2.py`, `gatk_vcf.py`, `gatk_adapter/`, `clair3_vcf.py`, `marlin_native.py` (+`_marlin_native_worker.py`) |
| Evidence layers | Combine and annotate normalized evidence without asserting clinical meaning | `sv_consensus.py`, `sv_annotation.py`, `breakends.py`, `breakpoint_annotation.py`, `sv_observability.py`, `sv_evidence.py`, `aml_rearrangements.py`, `fusion_evidence.py`, `guideline_criteria.py`, `cnv/cytoband.py`, `iscn.py`, `quantitation.py`, `small_variants.py`, `mvp.py` (result assembly) |
| Execution | Plan, run, record and resume one sample | `pipeline/stages.py` (declared graph), `pipeline/context.py` (stage contract), `pipeline/runner.py` (orchestrator), `pipeline/envelope.py`, `pipeline/state.py`, `pipeline/lock.py`, `pipeline/components.py`, `pipeline/input_digest.py`, `pipeline/checks.py`, `pipeline/review.py`, `pipeline/watch.py`, `pipeline/panel_lock.py`, `coverage_artifacts.py`, `execution.py` |
| Lanes | Stage implementations that own a whole analysis branch and contribute to assembly/report | `pipeline/marlin.py` (MARLIN), `cnv/lane.py` (QDNAseq/ACE copy number) |
| Reporting | Render validated contracts; never drive or query a run | `report.py`, `report_*.py`, `reporting.py`, `workbook.py`, `report_bundle.html`, `web/` (live workspace and Befund view sources) |
| Surfaces and operations | Turn configuration into runs and show their state | `entrypoint.py`, `runtime_cli.py`, `cli.py`, `marlin_cli.py`, `methylation_*_cli.py`, `profile_analysis.py`, `service/`, `watchfolder.py`, `preflight.py`, `status.py`, `review.py`, `system_smoke.py`, `smoke.py`, `demo.py`, `desktop/` (Windows WPF, no bioinformatics) |
| Research and validation programs | Separately executable studies with their own locks and contracts | `benchmark.py`, `dilution.py`, `cnv_validation_*.py`, `methylation_mixture.py`, `methylation_lanes/` (haplotype-resolved lane), `modbam.py`, `methylation_holdout*.py`, `methylation_validation*.py`, `marlin_*.py` (legacy R bridge, runtime freeze/compare, validation), `multicaller_*.py`, `tumor_inputs.py` |

## Dependency rules

1. **One declared stage graph (enforced).** `pipeline/stages.py` declares every stage and
   `SPEC_BY_STAGE` is read-only. No module assigns to, deletes from or updates
   `SPEC_BY_STAGE` or `IMPLEMENTATIONS`, and no production module uses `global`. What
   differs between runs belongs in `RunConfiguration`.
2. **Lanes depend on the contract, not on the orchestrator (enforced).** Only composition
   roots import `pipeline.runner` when they load: `runtime_cli`, `profile_analysis`,
   `service.app`, `service.befund`, `watchfolder`, `system_smoke`. Lanes import
   `pipeline.context`; type-only imports go under `TYPE_CHECKING`.
3. **Declarative modules stay dependency-free (enforced).** `pipeline/stages.py` and
   `pipeline/envelope.py` import nothing from the package and no pydantic, so the graph and
   the envelope can be tested without the rest of the system.
4. **Contracts do not reach into adapters (enforced).** `models.py` may import only
   `iscn_syntax.py`.
5. **Renderers never run anything (enforced).** Reporting modules import no execution,
   service or command-line module. A report is a function of validated contracts.
6. **One way in (enforced).** Only `runtime_cli` starts the local service; only entry points
   (`entrypoint`, `__main__` modules) import command-line modules; production code never
   imports tests.
7. **Current evidence only (reviewed, exercised by tests).** A later stage reads a lane's
   output only through `current_artifact()` / the lane's `load_current_*()` reader: an
   artifact the producing stage recorded in this run that still verifies byte for byte.
   Readers of archived envelopes outside a run (`load_run_coverage()`, `render`, the Befund
   view) are the only file-name readers.

## Where new code goes

- **A new external tool:** an adapter module (or a module in a domain subpackage such as
  `cnv/`, `sv/`, `tumor/`) with a versioned policy model, a version lock, shell-free argv, a
  typed report and `NOT_RUN`/`FAILED`/`NO_CALL` semantics. Real-tool tests are opt-in and run
  in a dedicated workflow.
- **A new stage of the canonical run:** a `StageSpec` in `pipeline/stages.py`, an entry in
  `IMPLEMENTATIONS`, configuration in `RunConfiguration`. If the stage owns a whole branch
  (its own assembly and report contributions), put it in a lane module that depends on
  `pipeline/context.py`.
- **A research study:** its own module, CLI and lock contracts. It must not change the
  canonical run's thresholds or reporting boundary.
- **Presentation:** a renderer function fed by typed reports, with escaping and
  determinism tests.

## Known structural debt

- `models.py` holds most contracts in one module (about 2,000 lines). Splitting it is
  worthwhile but touches nearly every import; do it in its own change.
- Methylation code is flat at the package root (`methylation*.py`, `modbam.py`,
  `marlin_*.py`). New methylation research lanes go into the `methylation_lanes/`
  subpackage (the haplotype lane is the first) rather than extend the root; moving the
  existing modules needs compatibility re-exports and should not collide with open work on
  the MARLIN files.
- `methylation_lanes/haplotype.py` reuses underscore helpers of `methylation.py` (bedMethyl
  sites, region assignment and summaries) so both lanes count identically. When a second
  lane needs them, promote them to a public module instead of widening the private use.
- `pipeline/runner.py` still implements the core stages inline; further lanes can move out
  the same way `cnv/lane.py` did.

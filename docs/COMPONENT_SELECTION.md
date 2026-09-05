# Choosing components per run

A run can now say which program runs each stage and at exactly which version, without any
code change. The selection is a file; running the same sample under two selections and
comparing the envelopes is the intended way to answer "does this caller version change the
result".

## The shape of a selection

```yaml
schema_version: 0.1.0
selection_id: legacy-sniffles-2.4-comparison-v1
status: technical_defaults_only
components:
  sv:
    provider: sniffles2
    version: 2.4.0
    policy: configs/sv/sniffles2.legacy_2_4.technical.yaml
  cnv:
    provider: qdnaseq_ace
    enabled: false
```

Selectable stages are those that run an external tool: `basecall`, `align`, `intake`, `qc`,
`target_coverage`, `cnv`, `sv`. Assemble, report and release are deliberately not
selectable — a run must not be able to switch off its own report or its checksum bundle.

Two files ship with the repository:

* `configs/components/default.yaml` — the toolchain CI exercises.
* `configs/components/legacy_sniffles_2.4.yaml` — the structural-variant component of the
  historical pipeline, for comparison runs.

## Using one

```bash
ontseq run sample.manifest.json \
  --reference-lock /approved/references/reference.lock.json \
  --run-id RUN_001 \
  --components configs/components/default.yaml

# switch a stage off for one run without editing any file
ontseq run sample.manifest.json ... --without cnv --without sv
```

The resolved selection is printed before the first stage and written to
`provenance/components.json` inside the envelope before anything executes, so an
interrupted run still records what it was asked to use.

## What the selection actually enforces

**A pinned version must be the installed version.** Every stage already probes the tool it
is about to run. That probe is now compared against the pin, in one place, before the stage
executes. A mismatch fails the stage naming both numbers:

```
sv: the selection pins sniffles2 to 2.4.0, but the installed tool reports 2.8.0.
Refusing to run a selected component at an unselected version
```

This is the difference between a reproducible run and a run that used whatever happened to
be on the machine.

**A deselected stage is recorded, not omitted.** It becomes `NOT_RUN` with the selection
named as the reason, and the record says in words that a stage which was switched off is
not a negative finding.

**An unpinned component warns.** A component with no version still runs, but the run report
carries a warning that it is reproducible only against the toolchain that was installed.

**A stage that ran no tool is never failed by a pin.** Target coverage on an lcWGS run
probes nothing, so a `mosdepth` pin does not fail it.

**Changing a component re-runs the stage.** Tool versions are part of the resume signature,
so switching Sniffles 2.8.0 to 2.4.0 in the same envelope re-executes structural-variant
calling instead of resuming the previous result. Two versions cannot end up mixed in one
envelope.

## An honest limitation

The CNV lane is still installed by registering an implementation into the runner rather
than by being a first-class member of the execution graph. That is deliberate — the core
must not depend on a specific copy-number caller — and the selection gates it, so a run
that deselects CNV does not get it.

What registration may do is now bounded. It supplies the CNV stage and *contributes* to
assemble and report; it can no longer replace a stage it does not own, and it can no
longer rewrite `SPEC_BY_STAGE`. Those two powers cost three defects before they were
removed: the SV consensus silently missing from every result, a resume signature that
ignored artifacts the stage read, and a preflight that described a different run than the
one it was clearing. See `docs/PIPELINE_EXECUTION.md`, "Extensions contribute; they do not
replace".

What remains is narrower: a reader of `stages.py` alone still cannot tell that a CNV
adapter exists, because only registration reveals it.

## Comparing two callers

```bash
ontseq run sample.manifest.json ... --run-id CMP_28 \
  --components configs/components/default.yaml
ontseq run sample.manifest.json ... --run-id CMP_24 \
  --components configs/components/legacy_sniffles_2.4.yaml
```

Two envelopes, each recording its own selection, each with its own checksummed release
bundle. Neither is the reference for the other; agreement between them is agreement, not
truth.

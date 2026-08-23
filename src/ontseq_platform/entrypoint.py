"""Dispatch to whichever CLI owns the requested command, and say that both exist.

The commands are split across two parsers: ``runtime_cli`` owns execution and operations,
``cli`` owns the scientific single-step adapters. Dispatching on the first argument keeps
the execution core from importing the whole scientific surface, which is worth having.

What it cost, before this module printed anything of its own, was discoverability: ``ontseq
--help`` reached only one parser, so ``run``, ``preflight``, ``status``, ``watch``,
``serve`` and ``review`` — the commands an operator needs first — were invisible unless you
already knew to type them. So a bare invocation lists both groups and then hands over.
"""

from __future__ import annotations

import sys

from .runtime_cli import RUNTIME_COMMANDS

_SCIENTIFIC_COMMANDS = (
    ("demo", "Generate a synthetic JSON/HTML/Excel report"),
    ("validate-manifest", "Validate a sample manifest against the schema"),
    ("validate-result", "Validate a pipeline result against the schema"),
    ("render", "Render HTML and Excel from an existing result JSON"),
    ("reference-lock", "Create a versioned reference lock from a FASTA index"),
    ("inspect-bam", "Run the aligned-BAM integrity and reference gate"),
    ("qc-cramino", "Run Cramino and normalize descriptive BAM QC"),
    ("qc-target-coverage", "Run Mosdepth over a target design"),
    ("call-sniffles", "Run Sniffles2 and normalize candidate SV evidence"),
    ("local-smoke", "Exercise the real toolchain on generated synthetic alignments"),
    ("system-smoke", "Full installed-system self-test including QDNAseq/ACE"),
    ("benchmark", "Score a benchmark case"),
    ("assemble-aligned-mvp", "Assemble a result from existing adapter outputs"),
    ("annotate", "Attach knowledge-base records to a result's events"),
    ("cnv-evaluate", "Score a CNV call set against a truth set"),
    ("cnv-aggregate", "Pool CNV evaluations into a stratified summary"),
    ("cnv-compare-methods", "Compare two CNV methods on shared truth events"),
    ("cnv-karyotype-truth", "Convert an ISCN karyotype into a CNV truth set"),
    ("cnv-demo-benchmark", "Run the synthetic CNV benchmark end to end"),
)

_RUNTIME_COMMANDS = (
    ("run", "Execute one sample into a resumable run envelope"),
    ("preflight", "Check every run precondition without creating output"),
    ("status", "Summarize run envelopes"),
    ("watch", "Process ready sample directories in a drop folder"),
    ("serve", "Run the loopback-only local operator service"),
    ("review", "Record or inspect review state"),
    ("model-lock", "Fingerprint a downloaded Dorado model directory"),
    ("align-fixture", "Generate a synthetic real-alignment fixture"),
)


def _overview() -> str:
    lines = [
        "usage: ontseq <command> [options]",
        "",
        "ONTSeq platform. Research use only; no output is a clinical result.",
        "",
        "Execution and operations:",
    ]
    lines += [f"  {name:<22}{summary}" for name, summary in _RUNTIME_COMMANDS]
    lines += ["", "Single-step adapters and reporting:"]
    lines += [f"  {name:<22}{summary}" for name, summary in _SCIENTIFIC_COMMANDS]
    lines += [
        "",
        "Run 'ontseq <command> --help' for the options of one command.",
        "Start with 'ontseq preflight' before any run against real data.",
    ]
    return "\n".join(lines)


def main() -> None:
    """Dispatch execution commands without coupling the legacy/scientific CLI to runtime code."""
    command = sys.argv[1] if len(sys.argv) > 1 else None
    if command is None or command in {"-h", "--help", "help"}:
        print(_overview())
        return
    if command in RUNTIME_COMMANDS:
        from .runtime_cli import main as runtime_main

        runtime_main()
    else:
        from .cli import main as legacy_main

        legacy_main()

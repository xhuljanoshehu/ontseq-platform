"""One rule, because there is one execution path.

This workflow used to re-implement the pipeline as five Snakemake rules, each invoking a
per-stage CLI command. That made Snakemake a *second* way to run a sample: it produced flat
files instead of a run envelope, recorded no per-artifact tool versions, had no release
bundle, and resumed on file timestamps rather than on content. Two paths mean two
behaviours, and only one of them was ever proven in CI.

ADR-013 settles this: a scheduler is a caller of ``ontseq run``, not an alternative to it.
What Snakemake is genuinely good for here — cluster submission, resource declarations,
fanning many samples out through an executor plugin — is preserved. What it was doing
badly, namely re-deriving the stage graph, is not.

Resume is deliberately left to the runner rather than expressed as per-stage Snakemake
outputs. The runner compares content hashes, parameters and tool versions; Snakemake
compares mtimes, and an mtime comparison would happily accept an artifact produced under
different parameters.
"""


rule ontseq_run:
    """Execute the whole pipeline for one sample into a resumable run envelope."""
    input:
        manifest=config["manifest"],
        reference_lock=config["reference_lock"],
        qc_policy=config["qc_policy"],
        sniffles_policy=config["sniffles_policy"],
    output:
        run_report=config["outputs"]["run_report"],
        release=config["outputs"]["release"],
        checksums=config["outputs"]["checksums"],
    params:
        output_dir=config["output_dir"],
        run_id=config["run_id"],
        threads=config["threads"],
        git_commit=config["git_commit"],
        samtools=config["tools"]["samtools"],
        cramino=config["tools"]["cramino"],
        sniffles=config["tools"]["sniffles"],
        minimap2=config["tools"]["minimap2"],
        # Only needed when the run starts from POD5 or an unaligned BAM; an aligned-BAM
        # run must not be forced to name a FASTA it will never read.
        reference_fasta=(
            f"--reference-fasta {config['reference_fasta']}"
            if config.get("reference_fasta")
            else ""
        ),
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform run {input.manifest:q} "
        "--reference-lock {input.reference_lock:q} "
        "--qc-policy {input.qc_policy:q} "
        "--sniffles-policy {input.sniffles_policy:q} "
        "--output-dir {params.output_dir:q} "
        "--run-id {params.run_id:q} "
        "--threads {params.threads} "
        "--git-commit {params.git_commit:q} "
        "--samtools {params.samtools:q} "
        "--cramino {params.cramino:q} "
        "--sniffles {params.sniffles:q} "
        "--minimap2 {params.minimap2:q} "
        "{params.reference_fasta}"

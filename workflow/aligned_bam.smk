configfile: "workflow/config/aligned_bam.example.yaml"

include: "rules/aligned_bam.smk"

rule all:
    input:
        config["outputs"]["run_report"],
        config["outputs"]["release"],
        config["outputs"]["checksums"],

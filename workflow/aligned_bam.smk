configfile: "workflow/config/aligned_bam.example.yaml"

include: "rules/aligned_bam.smk"

rule all:
    input:
        config["outputs"]["intake"],
        config["outputs"]["qc"],
        config["outputs"]["sniffles_vcf"],
        config["outputs"]["sniffles"],
        config["outputs"]["result"],
        config["outputs"]["html"],
        config["outputs"]["xlsx"],

configfile: "workflow/config/benchmark.yaml"

include: "rules/benchmarking.smk"

rule all:
    input:
        config["outputs"]["cnv"],
        config["outputs"]["sv"],
        config["outputs"]["cnv_state"],
        config["outputs"]["cnv_demo_aggregate"],

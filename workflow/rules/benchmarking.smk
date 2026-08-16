rule benchmark_synthetic_cnv:
    input:
        config["cases"]["cnv"],
    output:
        config["outputs"]["cnv"],
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform benchmark {input:q} --output {output:q}"


rule benchmark_synthetic_sv:
    input:
        config["cases"]["sv"],
    output:
        config["outputs"]["sv"],
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform benchmark {input:q} --output {output:q}"


rule evaluate_synthetic_cnv_case:
    """Score one locked CNV case with the state-based comparator."""
    input:
        config["cases"]["cnv_state"],
    output:
        config["outputs"]["cnv_state"],
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform cnv-evaluate {input:q} --output {output:q}"


rule cnv_demo_benchmark:
    """Simulate, call, evaluate and compare two caller configurations across strata."""
    output:
        config["outputs"]["cnv_demo_comparison"],
    params:
        output_dir=config["outputs"]["cnv_demo_dir"],
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform cnv-demo-benchmark"
        " --output-dir {params.output_dir:q}"

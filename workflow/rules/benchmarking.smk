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

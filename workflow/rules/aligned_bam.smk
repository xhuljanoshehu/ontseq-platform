rule inspect_aligned_bam:
    """Fail closed on BAM/index/reference incompatibility before scientific analysis."""
    input:
        manifest=config["manifest"],
        reference_lock=config["reference_lock"],
    output:
        config["outputs"]["intake"],
    params:
        samtools=config["tools"]["samtools"],
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform inspect-bam {input.manifest:q} "
        "--reference-lock {input.reference_lock:q} --samtools {params.samtools:q} "
        "--output {output:q}"


rule cramino_qc:
    """Normalize descriptive long-read QC without inventing clinical thresholds."""
    input:
        manifest=config["manifest"],
        intake=config["outputs"]["intake"],
        policy=config["qc_policy"],
    output:
        config["outputs"]["qc"],
    params:
        cramino=config["tools"]["cramino"],
        threads=config["threads"],
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform qc-cramino {input.manifest:q} "
        "--policy {input.policy:q} --cramino {params.cramino:q} "
        "--threads {params.threads} --output {output:q}"


rule sniffles2_candidates:
    """Call and normalize non-reportable Sniffles2 SV candidates."""
    input:
        manifest=config["manifest"],
        intake=config["outputs"]["intake"],
        policy=config["sniffles_policy"],
    output:
        vcf=config["outputs"]["sniffles_vcf"],
        report=config["outputs"]["sniffles"],
    params:
        sniffles=config["tools"]["sniffles"],
        threads=config["threads"],
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform call-sniffles {input.manifest:q} "
        "--intake {input.intake:q} --policy {input.policy:q} "
        "--sniffles {params.sniffles:q} --threads {params.threads} "
        "--vcf {output.vcf:q} --output {output.report:q}"


rule assemble_aligned_bam_mvp:
    """Create a result with candidate-only SV evidence and explicit unrun modules."""
    input:
        manifest=config["manifest"],
        intake=config["outputs"]["intake"],
        qc=config["outputs"]["qc"],
        sniffles=config["outputs"]["sniffles"],
    output:
        config["outputs"]["result"],
    params:
        git_commit=config["git_commit"],
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform assemble-aligned-mvp {input.manifest:q} "
        "--intake {input.intake:q} --qc {input.qc:q} --sniffles {input.sniffles:q} "
        "--git-commit {params.git_commit:q} "
        "--output {output:q}"


rule render_aligned_bam_mvp:
    """Render self-contained reviewer artifacts from the structured result."""
    input:
        result=config["outputs"]["result"],
    output:
        html=config["outputs"]["html"],
        xlsx=config["outputs"]["xlsx"],
    params:
        output_dir=config["output_dir"],
    conda:
        "../envs/aligned_bam.yaml"
    shell:
        "PYTHONPATH=src python -m ontseq_platform render {input.result:q} "
        "--output-dir {params.output_dir:q}"

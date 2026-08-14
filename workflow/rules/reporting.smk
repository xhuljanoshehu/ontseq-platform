rule generate_synthetic_demo:
    """Generate the contract-complete synthetic demo bundle."""
    output:
        json=config["outputs"]["json"],
        html=config["outputs"]["html"],
        xlsx=config["outputs"]["xlsx"],
    shell:
        "python -m ontseq_platform demo --output-dir {config[output_dir]}"

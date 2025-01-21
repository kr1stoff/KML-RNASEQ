rule gene_quantitaion_boxplot:
    input:
        filtered_tpm=rules.filter_tpm.output[0],
        metadata=config["metadata"],
    output:
        directory("gene_quantitaion/boxplot"),
    log:
        "logs/gene_quantitaion/boxplot.log",
    benchmark:
        "logs/gene_quantitaion/boxplot.benchmark"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/gene_quantitaion_boxplot.py"

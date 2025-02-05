rule gene_quant_boxplot:
    input:
        filtered_tpm=rules.filter_tpm.output[0],
        metadata=config["metadata"],
    output:
        directory("gene_quant/boxplot"),
    log:
        "logs/gene_quant/boxplot.log",
    benchmark:
        "logs/gene_quant/boxplot.benchmark"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/gene_quant_boxplot.py"


rule gene_quant_violinplot:
    input:
        filtered_tpm=rules.filter_tpm.output[0],
        metadata=config["metadata"],
    output:
        directory("gene_quant/violinplot"),
    log:
        "logs/gene_quant/violinplot.log",
    benchmark:
        "logs/gene_quant/violinplot.benchmark"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/gene_quant_violinplot.py"


rule intersample_correlation:
    input:
        rules.filter_tpm.output[0],
    output:
        table="gene_quant/intersample_correlation/output.tsv",
        fig="gene_quant/intersample_correlation/output.png",
    log:
        "logs/gene_quant/intersample_correlation.log",
    benchmark:
        "logs/gene_quant/intersample_correlation.benchmark"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/intersample_correlation.py"


rule gene_quant_pca:
    input:
        filtered_tpm=rules.filter_tpm.output[0],
        metadata=config["metadata"],
    output:
        directory("gene_quant/pca"),
    log:
        "logs/gene_quant/pca.log",
    benchmark:
        "logs/gene_quant/pca.benchmark"
    conda:
        config["conda"]["R"]
    script:
        "../scripts/gene_quant_pca.R"

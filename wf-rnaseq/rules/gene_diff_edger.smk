rule gene_diff_edger:
    input:
        filt_gene_count=rules.filter_gene_count.output[0],
        metadata=config["metadata"],
    output:
        directory("gene_diff/edgeR/ALL"),
    log:
        "logs/gene_diff/gene_diff_edger.log",
    benchmark:
        "logs/gene_diff/gene_diff_edger.benchmark"
    conda:
        config["conda"]["R"]
    script:
        "../scripts/gene_diff_edger.R"

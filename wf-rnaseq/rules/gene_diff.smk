rule gene_diff_deseq2:
    input:
        filt_gene_count=rules.filter_gene_count.output[0],
        metadata=config["metadata"],
    output:
        directory("gene_diff/DESeq2"),
    log:
        "logs/gene_diff/gene_diff_deseq2.log",
    benchmark:
        "logs/gene_diff/gene_diff_deseq2.benchmark"
    conda:
        config["conda"]["R"]
    script:
        "../scripts/gene_diff_deseq2.R"

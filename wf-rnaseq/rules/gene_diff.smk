rule gene_diff_deseq2:
    input:
        filt_gene_count=rules.filter_gene_count.output[0],
        metadata=config["metadata"],
    output:
        # 1. 输出不同分组方案下的所有基因 logFC 等统计值
        # 2. 输出 1 对应的火山图
        # 3. 输出不同分组方案的均一化基因表达数据表, 用于热图等下游分析
        directory("gene_diff/DESeq2/DEG"),
    log:
        "logs/gene_diff/gene_diff_deseq2.log",
    benchmark:
        "logs/gene_diff/gene_diff_deseq2.benchmark"
    conda:
        config["conda"]["R"]
    script:
        "../scripts/gene_diff_deseq2.R"


rule gene_diff_edger:
    input:
        filt_gene_count=rules.filter_gene_count.output[0],
        metadata=config["metadata"],
    output:
        directory("gene_diff/edgeR/DEG"),
    log:
        "logs/gene_diff/gene_diff_edger.log",
    benchmark:
        "logs/gene_diff/gene_diff_edger.benchmark"
    conda:
        config["conda"]["R"]
    script:
        "../scripts/gene_diff_edger.R"

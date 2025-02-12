rule gene_diff_deseq2:
    input:
        filt_gene_count=rules.filter_gene_count.output[0],
        metadata=config["metadata"],
    output:
        # 1. 输出不同分组方案下的所有基因 logFC 等统计值
        # 2. 输出 1 对应的火山图
        # 3. 输出不同分组方案的均一化基因表达数据表, 用于热图等下游分析
        directory("gene_diff/DESeq2/ALL"),
    log:
        "logs/gene_diff/gene_diff_deseq2.log",
    benchmark:
        "logs/gene_diff/gene_diff_deseq2.benchmark"
    conda:
        config["conda"]["R"]
    script:
        "../scripts/gene_diff_deseq2.R"


rule deseq2_get_deg:
    input:
        rules.gene_diff_deseq2.output[0],
    output:
        directory("gene_diff/DESeq2/DEG"),
    log:
        "logs/gene_diff/deseq2_get_deg.log",
    benchmark:
        "logs/gene_diff/deseq2_get_deg.benchmark"
    conda:
        config["conda"]["R"]
    script:
        "../scripts/deseq2_get_deg.R"


rule deseq2_pheatmap:
    input:
        alldir=rules.gene_diff_deseq2.output[0],
        degdir=rules.deseq2_get_deg.output[0],
        metadata=config["metadata"],
    output:
        directory("gene_diff/DESeq2/pheatmap"),
    log:
        "logs/gene_diff/deseq2_pheatmap.log",
    benchmark:
        "logs/gene_diff/deseq2_pheatmap.benchmark"
    conda:
        config["conda"]["Rpheatmap"]
    script:
        "../scripts/deseq2_pheatmap.R"

# * log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# library
library(DESeq2)
library(ggplot2)

# * io
genefile <- snakemake@input[["filt_gene_count"]]
metafile <- snakemake@input[["metadata"]]
outdir <- snakemake@output[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# * main
gene <- read.delim(genefile, row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.delim(metafile, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

# 定义函数来统计每个分组的 deseq2
deseq2_stat <- function(group_column) {
    group <- metadata[colnames(gene), group_column]
    # 指定分组因子顺序
    coldata <- data.frame(group = factor(group, levels = unique(group)))
    # DESeq2 默认流程
    # 第一步，构建 DESeqDataSet 对象
    dds <- DESeqDataSetFromMatrix(countData = gene, colData = coldata, design = ~group)
    # 获取归一化的基因表达矩阵
    vsd <- assay(vst(dds, blind = FALSE))
    # 输出
    # write.table(vsd, 'norm_matrix.txt', sep = '\t', col.names = NA, quote = FALSE)
    # 第二步，差异分析，详见 ?DESeq 和 ?results
    # 标准方法
    dds <- DESeq(dds, parallel = TRUE) # parallel = TRUE 将启用多线程模式
    suppressMessages(dds)
    res <- results(dds, contrast = c("group", unique(group)), pAdjustMethod = "fdr", alpha = 0.05)
    # an alternate analysis: likelihood ratio test
    ddsLRT <- DESeq(dds, test = "LRT", reduced = ~1)
    suppressMessages(ddsLRT)
    resLRT <- results(ddsLRT, contrast = c("group", unique(group)), pAdjustMethod = "fdr", alpha = 0.05)
    # 可以先按校正和 p 值由小到大排个序，方便查看
    deseq_res <- as.data.frame(res[order(res$padj), ])
    # * 输出
    deseq_res$gene_id <- rownames(deseq_res)
    write.table(deseq_res[c(7, 1:6)], paste0(outdir, "/", group_column, ".txt"), row.names = FALSE, sep = "\t", quote = FALSE)

    # ggplot2 差异火山图
    deseq_res <- read.delim(paste0(outdir, "/", group_column, ".txt"), sep = "\t")
    # 例如这里根据 |log2FC| >= 1 & FDR p-value < 0.05 定义“差异”
    deseq_res[which(deseq_res$padj %in% NA), "sig"] <- "NS"
    deseq_res[which(deseq_res$log2FoldChange >= 1 & deseq_res$padj < 0.05), "sig"] <- "UP (p.adj < 0.05, log2FC >= 1)"
    deseq_res[which(deseq_res$log2FoldChange <= -1 & deseq_res$padj < 0.05), "sig"] <- "DOWN (p.adj < 0.05, log2FC <= -1)"
    deseq_res[which(abs(deseq_res$log2FoldChange) < 1 | deseq_res$padj >= 0.05), "sig"] <- "NS"
    # 纵轴为显著性 p 值
    volcano_p <- ggplot(deseq_res, aes(log2FoldChange, -log(padj, 10))) +
        geom_point(aes(color = sig), alpha = 0.6, size = 1) +
        scale_color_manual(values = c("blue2", "gray30", "red2")) +
        theme(panel.grid = element_blank(), panel.background = element_rect(color = "black", fill = "transparent"), legend.position = c(0.4, 0.9)) +
        theme(legend.title = element_blank(), legend.key = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent")) +
        geom_vline(xintercept = c(-1, 1), color = "gray", size = 0.25) +
        geom_hline(yintercept = -log(0.05, 10), color = "gray", size = 0.25) +
        labs(x = "log2 Fold Change", y = "-log10 p-value", color = NA) +
        xlim(-5, 5)
    # * 输出
    ggsave(paste0(outdir, "/", group_column, ".png"), volcano_p, width = 5, height = 6)
}

# 迭代元数据中的每个分组列，并生成 PCA 图
group_columns <- colnames(metadata)
stats_res <- lapply(group_columns, deseq2_stat)

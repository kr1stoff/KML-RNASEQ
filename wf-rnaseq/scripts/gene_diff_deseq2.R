# DESeq2 差异分析
# ! deepseek 说 DESeq2 更适合用与中或大样本量, 分析速度不如 edgeR
# DESeq2 默认流程, 然后根据 meta 拆分不同分组方案
# 在各个分组方案下进行分组分析
# 分组两两比对分析差异基因
# 细化可视化
# 输出:
# 1. B27_pos_vs_neg.tsv   :   分组方案组间基因差异表格
# 2. B27_neg_vs_pos.png   :   分组方案组间基因差异火山图


# * log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# library
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)

# * io
genefile <- snakemake@input[["filt_gene_count"]]
metafile <- snakemake@input[["metadata"]]
outdir <- snakemake@output[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


######################################## Function ########################################
#' 获取分组方案两两组合
#' group: 分组列表
#' 返回: 分组方案两两组合列表
get_pairwise_combinations <- function(group) {
  # 应对 3 组及以上情况, 如果不清楚 case control 就都做一遍
  categories <- unique(group)
  pairwise_combinations <- expand.grid(category1 = categories, category2 = categories)
  pairwise_combinations <- pairwise_combinations[pairwise_combinations$category1 != pairwise_combinations$category2, ]
  return(pairwise_combinations)
}


#' 输出差异表格添加 symbol 列. 删除转化不出 symbol 的行
#' group_column: 分组列名
#' cat1: 分组1
#' cat2: 分组2
#' deseq_res: 差异分析结果
#' outfile: 输出文件名
write_symbol_table <- function(group_column,
                               cat1,
                               cat2,
                               deseq_res,
                               outfile) {
  # * ENSG 转 symbol 也输出
  deseq_res$symbol <- mapIds(
    org.Hs.eg.db,
    keys = deseq_res$gene_id,
    column = "SYMBOL",
    keytype = "ENSEMBL"
  )
  notisna_symb_data <- deseq_res[!is.na(deseq_res$symbol), ]
  write.table(
    notisna_symb_data[c(7:8, 1:6)],
    outfile,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}


#' 输出差异分析结果
#' pairwise_combinations: 分组方案两两组合列表
#' group_column: 分组列名
#' dds2: DESeq2 对象
#' 返回: 差异分析结果列表
write_deseq_res <- function(pairwise_combinations,
                            group_column,
                            dds2) {
  # 迭代两两组合输出差异表格
  lapp_deseq_res <- lapply(1:nrow(pairwise_combinations), function(i) {
    cat1 <- as.character(pairwise_combinations$category1[i])
    cat2 <- as.character(pairwise_combinations$category2[i])
    outfile <- sprintf("%s/%s.%s_vs_%s.tsv", outdir, group_column, cat1, cat2)
    res <- results(
      dds2,
      contrast = c("group", cat1, cat2),
      pAdjustMethod = "fdr",
      alpha = 0.05
    )
    # 可以先按校正和 p 值由小到大排个序，方便查看
    deseq_res <- as.data.frame(res[order(res$padj), ])

    deseq_res$gene_id <- rownames(deseq_res)
    # * # 输出. 添加 symbol 列，从 ENSG 转换
    write_symbol_table(group_column, cat1, cat2, deseq_res, outfile)
    return(outfile)
  })
  return(lapp_deseq_res)
}


#' 输出差异分析结果火山图
#' intab: 差异分析结果表格
#' 输出: 差异分析结果火山图
deseq2_volcano <- function(intab) {
  # ggplot2 差异火山图
  out_fig <- gsub("tsv", "volcano.png", intab)
  deseq_res <- read.delim(intab, sep = "\t")
  # 例如这里根据 |log2FC| >= 1 & FDR p-value < 0.05 定义“差异”
  deseq_res[which(deseq_res$padj %in% NA), "sig"] <- "NS"
  deseq_res[which(deseq_res$log2FoldChange >= 1 &
    deseq_res$padj < 0.05), "sig"] <- "UP (p.adj < 0.05, log2FC >= 1)"
  deseq_res[which(deseq_res$log2FoldChange <= -1 &
    deseq_res$padj < 0.05), "sig"] <- "DOWN (p.adj < 0.05, log2FC <= -1)"
  deseq_res[which(abs(deseq_res$log2FoldChange) < 1 |
    deseq_res$padj >= 0.05), "sig"] <- "NS"
  # 纵轴为显著性 p 值
  volcano_p <- ggplot(deseq_res, aes(log2FoldChange, -log(padj, 10))) +
    geom_point(aes(color = sig), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("blue2", "gray30", "red2")) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      legend.position = c(0.26, 0.92)
    ) +
    theme(
      legend.title = element_blank(),
      legend.key = element_rect(fill = "transparent"),
      legend.background = element_rect(fill = "transparent")
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      color = "gray",
      size = 0.25
    ) +
    geom_hline(
      yintercept = -log(0.05, 10),
      color = "gray",
      size = 0.25
    ) +
    labs(x = "log2 Fold Change", y = "-log10 p-value", color = NA) +
    xlim(-5, 5)
  # 输出
  ggsave(out_fig, volcano_p, width = 12, height = 9)
}


#' 输出归一化的基因表达矩阵
#' group_column: 分组方案名
#' 输出: 不同分组方案下归一化的基因表达矩阵
write_norm_matrix_file <- function(group_column) {
  # 归一化前 ENSG 转 symbol, 重复的 symbol 取均值
  # ensg 转 symbols
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = rownames(gene),
    column = "SYMBOL",
    keytype = "ENSEMBL"
  )
  notisna_symbs <- symbols[!is.na(symbols)]
  gene_count_symbols <- gene[names(notisna_symbs), ]
  gene_count_symbols$symbol <- notisna_symbs
  # 相同基因取均值
  gene_count_mean <- aggregate(. ~ symbol, mean, data = gene_count_symbols)
  rownames(gene_count_mean) <- gene_count_mean$symbol
  # 删除 symbol 列
  gene_count_mean <- gene_count_mean[, -1]
  # 转成整数格式
  gene_count_mean_int <- data.frame(lapply(gene_count_mean, as.integer))
  rownames(gene_count_mean_int) <- rownames(gene_count_mean)

  # 获取归一化的基因表达矩阵, ENSG 编号转成基因名称
  group <- metadata[colnames(gene_count_mean_int), group_column]
  coldata <- data.frame(group = factor(group, levels = unique(group)))
  # 构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(
    countData = gene_count_mean_int,
    colData = coldata,
    design = ~group
  )
  # 输出归一化的基因表达矩阵
  vsd <- assay(vst(dds, blind = FALSE))
  norm_matrix_file <- sprintf("%s/%s.norm_matrix.tsv", outdir, group_column)
  write.table(
    vsd,
    norm_matrix_file,
    sep = "\t",
    col.names = NA,
    quote = FALSE
  )
}

# * 主函数, deseq2 统计每个分组方案
#' group_column: 分组方案名
deseq2_stat <- function(group_column) {
  group <- metadata[colnames(gene), group_column]
  coldata <- data.frame(group = factor(group, levels = unique(group)))
  # DESeq2 默认流程
  # 第一步，构建 DESeqDataSet 对象
  dds <- DESeqDataSetFromMatrix(
    countData = gene,
    colData = coldata,
    design = ~group
  )
  # 输出归一化的基因表达矩阵
  vsd <- assay(vst(dds, blind = FALSE))
  # 第二步，差异分析，详见 ?DESeq 和 ?results
  dds2 <- DESeq(dds, parallel = FALSE)
  pairwise_combinations <- get_pairwise_combinations(group)
  lapp_deseq_res <- write_deseq_res(pairwise_combinations, group_column, dds2)
  # 火山图
  volcano_res <- sapply(unlist(lapp_deseq_res), deseq2_volcano)
}
######################################## Function ########################################


# * main
gene <- read.delim(genefile, row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.delim(metafile, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
group_columns <- colnames(metadata)
stat_res <- sapply(group_columns, deseq2_stat)
# ! 热图. pheatmap 和当前其他包版本不兼容, 其他脚本画
write_norm_res <- sapply(group_columns, write_norm_matrix_file)

# * log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


# library
library(tidyverse)
library(pheatmap)
library(org.Hs.eg.db)

# * io
in_dir <- snakemake@input[["deseq2_deg_dir"]]
meta_file <- snakemake@input[["metadata"]]
out_dir <- snakemake@output[[1]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

######################################## Function ########################################
#' 获取 top20 差异基因
#' scheme: 分组方案
#' return: top20 基因数组
get_top20_diff_genes <- function(scheme) {
    myprn <- paste0(scheme, ".*_vs_.*\\.tsv$")
    diff_file <- list.files(
        path = in_dir,
        pattern = myprn,
        full.names = TRUE
    )[1]
    diff_data <- read.delim(diff_file, nrows = 100)
    symbols <- mapIds(
        org.Hs.eg.db,
        keys = diff_data$gene_id,
        column = "SYMBOL",
        keytype = "ENSEMBL"
    )
    # 映射为NA的丢弃
    isna_symbols <- symbols[!is.na(symbols)]
    return(isna_symbols[1:20])
}


#' 画热图
#' top20_genes_df: top20 基因均一化表达矩阵
#' scheme: 分组方案
#' 无返回值
pheatmap_by_scheme <- function(top20_genes_df, scheme) {
    # pheatmap 分组
    mygrp <- metadata[rownames(top20_genes_df), scheme] %>% as.factor()
    df_anno_row <- data.frame(SampleGroup = mygrp)
    rownames(df_anno_row) <- rownames(top20_genes_df)
    # pheatmap 热图
    out_fig <- sprintf("%s/%s.pheatmap.png", out_dir, scheme)
    png(
        file = out_fig,
        width = 15,
        height = 9,
        units = "in",
        res = 300
    )
    pheatmap(
        top20_genes_df,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        fontsize_col = 8,
        angle_col = 90,
        # main = 'cgMLST Sample-Allele Heatmap',
        annotation_row = df_anno_row,
        annotation_legend = FALSE,
        annotation_names_row = FALSE,
        legend = FALSE
    )
    dev.off()
}


#' 主函数
#' norm_file: 均一化表达文件
#' 无返回值
analysis_by_scheme_files <- function(norm_file) {
    # 当前方案名称 'B27'
    scheme <- strsplit(basename(norm_file), "\\.")[[1]][1]
    # 读均一化表达数据
    df <- read.delim(norm_file, row.names = 1)
    # row 样本, col 基因
    t_df <- t(df)
    # top20基因
    top20_genes <- get_top20_diff_genes(scheme)
    top20_genes_df <- t_df[, top20_genes]
    # 画热图
    pheatmap_by_scheme(top20_genes_df, scheme)
}
######################################## Function ########################################


# * main
# 元数据
metadata <- read.delim(
    meta_file,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
# 不同分组方案生成的均一化表格
files <- list.files(
    path = in_dir,
    pattern = "\\.norm_matrix\\.tsv$",
    full.names = TRUE
)
analysis_res <- sapply(files, analysis_by_scheme_files)

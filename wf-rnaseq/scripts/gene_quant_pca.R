# * log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# library
library(tidyverse)
library(ggrepel)
library(dplyr)

# * io
genefile <- snakemake@input[["filtered_tpm"]]
metafile <- snakemake@input[["metadata"]]
outdir <- snakemake@output[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# * main
# 读取基因表达数据和元数据
genedata <- read.delim(genefile, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE) %>%
    t() %>%
    as.data.frame()

metadata <- read.delim(metafile, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

# 定义函数来生成每个分组的 PCA 图
generate_pca_plot <- function(group_column) {
    genedata$group <- metadata[rownames(genedata), group_column]
    pca <- prcomp(genedata[, -ncol(genedata)], scale. = TRUE)
    pca_var <- pca$sdev^2
    pca_var_per <- round(pca_var / sum(pca_var) * 100, 2)

    # [250211] 孟博, 图上加样本名
    pca_table <- as.data.frame(pca$x) %>%
        mutate(group = genedata$group, sample_name = rownames(genedata))

    p <- ggplot(aes(PC1, PC2, color = group), data = pca_table) +
        geom_point(size = 3, aes(shape = group)) +
        geom_text_repel(aes(label = sample_name), size = 3) +
        theme_bw() +
        labs(
            x = paste("PC1(", pca_var_per[1], "%)", sep = ""),
            y = paste("PC2(", pca_var_per[2], "%)", sep = ""),
            color = group_column,
            shape = group_column
        ) +
        stat_ellipse(level = 0.68) +
        theme(
            axis.text = element_text(size = 11, color = "black"),
            axis.ticks = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5)
        )

    ggsave(paste0(outdir, "/", group_column, ".png"), dpi = 300)
}

# 迭代元数据中的每个分组列，并生成 PCA 图
group_columns <- colnames(metadata)
pca_plots <- lapply(group_columns, generate_pca_plot)

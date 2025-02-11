# * log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# library
library(edgeR)

# * io
genefile <- snakemake@input[["filt_gene_count"]]
metafile <- snakemake@input[["metadata"]]
outdir <- snakemake@output[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# * 后面有需要再优化这一堆代码
main <- function(group_column) {
    # 分组方案 B27, Strand ...
    # group_column <- "B27"
    group_list <- metadata[colnames(targets), group_column]
    # 应对 3 组及以上情况, 如果不清楚 case control 就都做一遍
    categories <- unique(group_list)
    pairwise_combinations <- expand.grid(category1 = categories, category2 = categories)
    pairwise_combinations <- pairwise_combinations[pairwise_combinations$category1 != pairwise_combinations$category2, ]
    # 迭代两两组合输出差异表格
    lapp_deseq_res <- lapply(1:nrow(pairwise_combinations), function(i) {
        # i = 1
        cat1 <- as.character(pairwise_combinations$category1[i])
        cat2 <- as.character(pairwise_combinations$category2[i])
        glmLRT_file <- sprintf(
            "%s/%s_%s_vs_%s.glmLRT.tsv",
            outdir,
            group_column,
            cat1,
            cat2
        )

        # design 分组关系
        samples <- rownames(metadata[metadata[[group_column]] %in% c(cat1, cat2), ])
        now_targets <- targets[, samples]
        now_group <- metadata[metadata[[group_column]] %in% c(cat1, cat2), group_column]
        now_group <- factor(now_group)
        now_group <- relevel(now_group, cat1) # 将对照组的因子设置为1
        design <- model.matrix(~ 0 + now_group)
        rownames(design) <- colnames(now_targets)
        colnames(design) <- levels(now_group)

        # （1）构建 DGEList 对象
        dgelist <- DGEList(counts = now_targets, group = now_group)
        # （2）过滤 low count 数据, 使用 CPM 标准化（推荐）
        keep <- rowSums(cpm(dgelist) > 1) >= 2
        dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
        # （3）标准化，TMM 标准化
        cnf <- calcNormFactors(dgelist, method = "TMM")
        norm_counts <- cpm(cnf)
        norm_counts_file <- sprintf("%s/%s.edgeR_TMM.tsv", outdir, group_column)
        # ! 没有转symbol. 没有先合并 symbol 再均一化, 如果后面优化
        # 输出 TMM 标准化数据表格
        write.table(
            norm_counts,
            norm_counts_file,
            sep = "\t",
            col.names = NA,
            quote = FALSE
        )
        # PCA 样本无监督聚类
        # plotMDS(cnf, col = rep(c('red', 'blue'), each = 5), dim = c(1, 2))

        # （4）估算离散值
        dge <- estimateDisp(cnf, design, robust = TRUE) # 估算离散值
        # （5）差异分析
        # negative binomial generalized log-linear model 拟合
        fit <- glmFit(dge, design, robust = TRUE) # 拟合模型
        lrt <- glmLRT(fit) # 统计检验
        # 输出主要结果
        write.table(
            topTags(lrt, n = nrow(dgelist$counts)),
            glmLRT_file,
            sep = "\t",
            col.names = NA,
            quote = FALSE
        )
    })
}

# * main
# 读取数据，设置分组
targets <- read.delim(
    genefile,
    row.names = 1,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
metadata <- read.delim(
    metafile,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
group_columns <- colnames(metadata)
stat_res <- sapply(group_columns, main)

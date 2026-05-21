# * log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


# * library
library(dplyr)

# * io
indir <- snakemake@input[[1]]
outdir <- snakemake@output[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# 输出参数到 log 文件
cat("=== Input/Output Parameters ===\n")
cat("indir:", indir, "\n")
cat("outdir:", outdir, "\n")
cat("===============================\n")

# 不同分组方案生成的均一化表格
files <- list.files(
    path = indir,
    pattern = ".*_vs_.*\\.tsv$",
    full.names = TRUE
)

for (dsfile in files) {
    outfile <- sprintf("%s/%s", outdir, basename(dsfile))
    # ! 过滤
    # 1. symbol NA
    # 2. padj NA
    # 3. |log2FoldChange| < 1
    # 4. padj > 0.05
    data <- read.delim(dsfile)
    filt_data <- data %>%
        filter(
            !is.na(symbol),
            !is.na(padj),
            abs(log2FoldChange) > 1,
            padj < 0.05
        )
    sort_filt_data <- filt_data %>%
        arrange(desc(abs(log2FoldChange)), padj)
    write.table(
        sort_filt_data,
        outfile,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    # 输出中间结果到 log 文件
    cat("=== Intermediate Results ===\n")
    cat("dsfile:", dsfile, "\n")
    cat("data dimensions:", dim(data), "\n")
    cat("filt_data dimensions:", dim(filt_data), "\n")
    cat("sort_filt_data dimensions:", dim(sort_filt_data), "\n")
    cat("outfile:", outfile, "\n")
    cat("===============================\n")
}

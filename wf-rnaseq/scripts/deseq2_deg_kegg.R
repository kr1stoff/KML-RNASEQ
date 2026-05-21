# * log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


# * library
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pathview)


# IO
indir <- snakemake@input[[1]]
outdir <- snakemake@output[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# 输出参数到 log 文件
cat("=== Input/Output Parameters ===\n")
cat("indir:", indir, "\n")
cat("outdir:", outdir, "\n")
cat("===============================\n")


######################################## Function ########################################
enrich_kegg_pipe <- function() {
    # 不同分组方案生成的均一化表格
    files <- list.files(
        path = indir,
        pattern = ".*_vs_.*\\.tsv$",
        full.names = TRUE
    )
    for (dsfile in files) {
        outprfx <- sprintf("%s/%s", outdir, gsub(".tsv", "", basename(dsfile)))
        geneFC <- generate_geneFC(dsfile)
        # ! 如果没有差异基因，跳过该方案
        if (length(geneFC) == 0) {
            cat("Warning: 分组: ", dsfile, " 没有差异基因，跳过该方案。\n")
            next
        }
        # 20260521 报错, 添加日志
        tryCatch({
            enrich_kegg_res <- analyza_enrich_kegg(names(geneFC), outprfx)
        }, error = function(e) {
            cat("Error in analyza_enrich_kegg:\n")
            cat("  dsfile:", dsfile, "\n")
            cat("  geneFC dimensions:", dim(geneFC), "\n")
            cat("  entrezids:", paste(names(geneFC), collapse = ", "), "\n")
            cat(e$message, "\n")
            stop(e$message)
        })
        # ! 没有显著的富集通路就跳过
        if (dim(enrich_kegg_res)[1] == 0) next
        kegg_plot(enrich_kegg_res, outprfx)
        kegg_pathway_plot(geneFC, enrich_kegg_res$ID[1], outprfx)
    }
}


generate_geneFC <- function(dsfile) {
    data <- read.delim(dsfile)
    rownames(data) <- data$gene_id
    cnvrt_ids <- bitr(data$gene_id,
        fromType = "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
    )
    uniq_cnvrt_ids <- cnvrt_ids[!duplicated(cnvrt_ids$ENSEMBL), ]
    genes <- uniq_cnvrt_ids$ENTREZID
    geneFC <- data[uniq_cnvrt_ids$ENSEMBL, ]$log2FoldChange
    names(geneFC) <- genes
    return(geneFC)
}


analyza_enrich_kegg <- function(entrezids, outprfx) {
    enrich_kegg_res <- enrichKEGG(
        gene = entrezids,
        organism = "hsa",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    write.table(
        enrich_kegg_res@result,
        file = paste0(outprfx, ".KEGG.tsv"),
        sep = "\t",
        quote = F,
        row.names = F
    )
    return(enrich_kegg_res)
}


kegg_plot <- function(enrich_kegg_res, outprfx) {
    # 绘制条形图
    barplot_file <- paste0(outprfx, ".barplot.png")
    png(
        file = barplot_file,
        width = 15,
        height = 15,
        units = "in",
        res = 300
    )
    print(barplot(
        enrich_kegg_res,
        drop = TRUE,
        showCategory = 20,
        main = "30min"
    ))
    dev.off()
    # 绘制点图
    bubbleplot_file <- paste0(outprfx, ".bubble.png")
    png(
        file = bubbleplot_file,
        width = 15,
        height = 15,
        units = "in",
        res = 300
    )
    print(dotplot(enrich_kegg_res))
    dev.off()
}


kegg_pathway_plot <- function(geneFC, pwid, outprfx) {
    # 绘制通路图, 画富集最高的通路
    orgn_wd <- getwd()
    setwd(sprintf("%s/%s", orgn_wd, outdir))
    pathview(
        gene.data = geneFC,
        pathway.id = pwid,
        species = "hsa",
        out.suffix = basename(outprfx)
    )
    setwd(orgn_wd)
}
######################################## Function ########################################

# * main
enrich_kegg_pipe()

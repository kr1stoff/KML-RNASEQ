# * log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


# * library
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


# IO
indir <- snakemake@input[[1]]
outdir <- snakemake@output[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


######################################## Function ########################################
enrich_go_pipe <- function() {
    # 不同分组方案生成的均一化表格
    files <- list.files(
        path = indir,
        pattern = ".*_vs_.*\\.tsv$",
        full.names = TRUE
    )
    for (dsfile in files) {
        outprfx <- sprintf("%s/%s", outdir, gsub(".tsv", "", basename(dsfile)))
        enrich_go_res <- analyza_enrich_go(dsfile, outprfx)
        tempres <- go_plot(enrich_go_res, outprfx)
    }
}

go_plot <- function(enrich_go_res, outprfx) {
    # 绘制条形图
    barplot_file <- paste0(outprfx, ".barplot.png")
    png(
        file = barplot_file,
        width = 15,
        height = 15,
        units = "in",
        res = 300
    )
    print(
        barplot(
            enrich_go_res,
            drop = TRUE,
            showCategory = 10,
            split = "ONTOLOGY"
        ) + facet_grid(ONTOLOGY ~ ., scale = "free")
    )
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
    print(
        dotplot(enrich_go_res, showCategory = 10, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free")
    )
    dev.off()
}

# enrichGO 富集分析
analyza_enrich_go <- function(dsfile, outprfx) {
    data <- read.delim(dsfile)
    converted_ids <- bitr(
        data$gene_id,
        fromType = "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
    )
    enrich_go_res <- enrichGO(
        gene = converted_ids$ENTREZID,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "ALL",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = FALSE
    )
    write.table(
        enrich_go_res,
        file = paste0(outprfx, ".GO.tsv"),
        sep = "\t",
        quote = F,
        row.names = F
    )
    return(enrich_go_res)
}
######################################## Function ########################################


# * main
enrich_go_pipe()

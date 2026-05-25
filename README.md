# KML-RNASEQ

RNA-seq analysis pipeline

## 命令行

poetry 更新到版本 2 之后使用 -C 指定工作目录, 程序的输入文件和输出文件不能用相对路径了, 因为这个相对路径是相对于工作目录, 不是当前目录

### python

```bash
~/miniforge3/envs/python3.12/bin/python main.py \
  --sample_table /data/mengxf/GitHub/KML-RNASEQ/tests/250115.input.tsv \
  --metadata /data/mengxf/GitHub/KML-RNASEQ/tests/250115.metadata.tsv \
  --work_dir /data/mengxf/Project/KML260520-RNASEQ-UPGRADE/result/260520
```

### snakemake

```bash
snakemake --cores 32 --use-conda --rerun-incomplete --scheduler \
  --snakefile /data/mengxf/GitHub/KML-RNASEQ/wf-rnaseq/Snakefile \
  --configfile .temp/snakemake.yaml
```

## 更新

### [v1.0.1] 20260525 修复 BUG

- qualimap rnaseq 部分, 指定临时目录, 避免因为存储不空不够导致流程中断
- gene_diff_deseq2.R, data.frame  check.names = FALSE 避免列名里面的 - 改成 . 导致后续分析错误
- 部分 R 脚本添加 tryCatch 错误处理

### [v1.0.0] 20260521 初始版本

- DESeq2 分组无差异时 (即没有基因 padj < 0.05), deseq2_pheatmap.R top20 gene 部分为空, 会引发报错. 如果没有差异基因, 则跳过该方案.
- (情况同上) KEGG&GO 富集分析时, 如果没有差异基因, 则跳过该方案.

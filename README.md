# KML-RNASEQ

RNA-seq analysis pipeline

## 命令行

poetry 更新到版本 2 之后使用 -C 指定工作目录, 程序的输入文件和输出文件不能用相对路径了, 因为这个相对路径是相对于工作目录, 不是当前目录

```bash
~/miniforge3/envs/python3.12/bin/python main.py \
  --sample_table /data/mengxf/GitHub/KML-RNASEQ/tests/250115.input.tsv \
  --metadata /data/mengxf/GitHub/KML-RNASEQ/tests/250115.metadata.tsv \
  --work_dir /data/mengxf/Project/KML260520-RNASEQ-UPGRADE/result/260520
```

## 更新

### [v1.0.0] 20260521 初始版本

- DESeq2 分组无差异时 (即没有基因 padj < 0.05), deseq2_pheatmap.R top20 gene 部分为空, 会引发报错. 如果没有差异基因, 则跳过该方案.
- (情况同上) KEGG&GO 富集分析时, 如果没有差异基因, 则跳过该方案.

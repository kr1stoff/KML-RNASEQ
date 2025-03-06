# KML-RNASEQ

RNA-seq analysis pipeline

## 命令行

poetry 更新到版本 2 之后使用 -C 指定工作目录, 程序的输入文件和输出文件不能用相对路径了, 因为这个相对路径是相对于工作目录, 不是当前目录

```bash
mamba run -n python3.12 \
  poetry -C /data/mengxf/GitHub/KML-RNASEQ run python /data/mengxf/GitHub/KML-RNASEQ/main.py \
  -s /data/mengxf/Project/KML250122_rnaseq_ZiYan/work/250306_strand_spec/input.tsv \
  -m /data/mengxf/Project/KML250122_rnaseq_ZiYan/work/250306_strand_spec/metadata.tsv \
  -w /data/mengxf/Project/KML250122_rnaseq_ZiYan/result/250306
```

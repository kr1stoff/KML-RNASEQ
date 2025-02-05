# 根据 featureCounts 结果整理基因表达表格，基因名，基因长度，各个样本基因表达量（reas 数）

import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

# * IO
feature_counts = snakemake.input[0]
output_file = snakemake.output.gene_count
output_withNTC = snakemake.output.withNTC

# * main
df = pd.read_csv(feature_counts, sep='\t', comment='#')
df.drop(columns=['Chr', 'Start', 'End', 'Strand'], inplace=True)
newheader = [col.split('/')[1].replace('.bam', '') if '/' in col else col for col in df.columns]
df.columns = newheader
df.to_csv(output_withNTC, sep='\t', index=False)
# ! 删除 NTC 样本数据
df.drop(columns=df.filter(like='NTC').columns, inplace=True)
df.to_csv(output_file, sep='\t', index=False)

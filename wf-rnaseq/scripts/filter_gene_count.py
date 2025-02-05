# 删除全是零的 gene count 条目，删除 Length 列

import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


def filter_tpm_table(df, filtered_file, threshold=0):
    if threshold == 0:
        # 至少1个
        filtered_df = df[(df != 0).sum(axis=1) >= 1]
    else:
        threshold = int(threshold * df.shape[1])
        filtered_df = df[(df != 0).sum(axis=1) >= threshold]
    filtered_df.to_csv(filtered_file, sep='\t')


# * io
gene_count_file = snakemake.input[0]
filtered_gene_count_file = snakemake.output[0]

# * main
df = pd.read_csv(gene_count_file, sep='\t', index_col=0)
df.drop('Length', axis=1, inplace=True)
filter_tpm_table(df, filtered_gene_count_file, threshold=0)

# 过滤 TPM 表，根据 “检出率(count不为0的比例)小于0.25的基因” 规则进行

import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

tpm_file = snakemake.input[0]
filtered_tpm_file_0 = snakemake.output[0]

df = pd.read_csv(tpm_file, sep='\t', index_col=0)


def filter_tpm_table(df, filtered_tpm_file, threshold=0.25):
    """
    过滤 TPM 表，根据 “检出率(count不为0的比例)小于0.25的基因” 规则进行
    """
    if threshold == 0:
        # 至少1个
        filtered_df = df[(df != 0).sum(axis=1) >= 1]
    else:
        threshold = int(threshold * df.shape[1])
        filtered_df = df[(df != 0).sum(axis=1) >= threshold]
    filtered_df.to_csv(filtered_tpm_file, sep='\t')


filter_tpm_table(df, filtered_tpm_file_0, threshold=0)

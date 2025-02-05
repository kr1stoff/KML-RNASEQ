# 计算 RPKM 和 TPM

import pandas as pd
import numpy as np
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")


def calc_rpkm(df):
    """
    计算 RPKM 值
    :param df: 包含基因名称、长度和样本计数的 DataFrame
    """
    sample_name = df.columns[2:]
    sample_reads = df.loc[:, sample_name].copy()
    gene_len = df.loc[:, ['Length']]
    total_reads = np.sum(sample_reads.values, axis=0).reshape(1, -1)
    rate = sample_reads.values / gene_len.values
    rpkm = rate / total_reads * 1e9
    return pd.DataFrame(data=rpkm, columns=sample_name, index=df['Geneid'])


def calc_tpm(rpkm_df):
    """
    计算 TPM 值
    :param rpkm_df: 包含基因名称和样本 RPKM 值的 DataFrame
    """
    sample_name = rpkm_df.columns
    total_rpkm = np.sum(rpkm_df.values, axis=0).reshape(1, -1)
    tpm = rpkm_df.values / total_rpkm * 1e6
    return pd.DataFrame(data=tpm, columns=sample_name, index=rpkm_df.index)


# main
gene_count = snakemake.input[0]  # 'result/250120/feature_counts/gene_count.tsv'
rpkm_out = snakemake.output.rpkm
tpm_out = snakemake.output.tpm

df = pd.read_csv(gene_count, sep='\t')
rpkm_df = calc_rpkm(df)
tpm_df = calc_tpm(rpkm_df)

rpkm_df.to_csv(rpkm_out, sep='\t')
tpm_df.to_csv(tpm_out, sep='\t')

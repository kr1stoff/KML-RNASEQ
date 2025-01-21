# 过滤 TPM 表，根据 “检出率(count不为0的比例)小于0.25的基因” 规则进行

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")

# tpm_file = Path('.').resolve().parents[1].joinpath('result/250120/feature_counts/tpm.tsv')
tpm_file = snakemake.input[0]
filtered_tpm_file = snakemake.output[0]

df = pd.read_csv(tpm_file, sep='\t', index_col=0)
# ! 筛选至少 25% 样本中表达量非零的基因
threshold = 0.25
non_zero_threshold = int(threshold * df.shape[1])
filtered_df = df[(df != 0).sum(axis=1) >= non_zero_threshold]

filtered_df.to_csv(filtered_tpm_file, sep='\t')

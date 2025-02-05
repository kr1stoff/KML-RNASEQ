# 使用筛选后的 TPM 表格和 metadata 画基因表达量的箱线图，metadata 可以有多种分组方式

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")


# IO
filt_tpm_file = snakemake.input.filtered_tpm
metadata_file = snakemake.input.metadata
outdir = snakemake.output[0]
Path(outdir).mkdir(parents=True, exist_ok=True)

# main
filt_df = pd.read_csv(filt_tpm_file, sep='\t', index_col=0)
meta_df = pd.read_csv(metadata_file, sep='\t', index_col=0)

for i in range(meta_df.shape[1]):
    sample_group_dict = meta_df.iloc[:, i].to_dict()
    filt_long_df = filt_df.melt(var_name='Sample', value_name='TPM')
    filt_long_df['TPM'] = np.log10(filt_long_df['TPM']+1)
    filt_long_df['Group'] = filt_long_df['Sample'].apply(lambda x: sample_group_dict[x])
    filt_long_df = filt_long_df.sort_values(by='Group')

    # 画图
    sns.set_theme(style="ticks")
    f, ax = plt.subplots(figsize=(7, 6))
    sns.violinplot(data=filt_long_df, x="Sample", y="TPM", hue="Group", fill=False, palette="Pastel1")
    ax.set(ylabel="", xlabel="")
    plt.xticks(rotation=90)
    ax.set_ylabel("log10(TPM+1)")
    sns.despine()
    plt.savefig(outdir + f'/{meta_df.columns[i]}.png', dpi=300, bbox_inches='tight')

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# import scipy
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

filt_tpm_file = snakemake.input[0]
out_table = snakemake.output.table
out_fig = snakemake.output.fig

filt_df = pd.read_csv(filt_tpm_file, sep='\t', index_col=0)
correlation_matrix = filt_df.corr(method='spearman')
correlation_matrix.to_csv(out_table, sep='\t', index=True)

# 绘制 Spearman 相关性系数矩阵的热图
plt.figure(figsize=(8, 6))
sns.clustermap(correlation_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1, linewidths=.5)
plt.title('Spearman Correlation Heatmap')
plt.savefig(out_fig, dpi=300, bbox_inches='tight')

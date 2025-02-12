import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")
count_file = snakemake.input.genecount
metafile = snakemake.input.metadata
output_file = snakemake.output[0]

with open(count_file) as f:
    header = f.readline().strip().split('\t')
    # * 样本名对应列数 {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    name_pos = {h[1]: h[0]+1 for h in enumerate(header)}


with open(metafile) as f, open(output_file, 'w') as g:
    g.write("#Sample\tCondition\tPath\tColumn\n")
    # * 跳过表头
    next(f)
    for line in f:
        lns = line.strip().split('\t')
        name, group = lns[:2]
        path = Path(count_file).resolve()
        column = name_pos[name]
        g.write(f"{name}\t{group}\t{path}\t{column}\n")

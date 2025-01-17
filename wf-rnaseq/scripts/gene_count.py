import sys

sys.stderr = open(snakemake.log[0], "w")
feature_counts = snakemake.input[0]
output_file = snakemake.output[0]


with open(feature_counts) as f, open(output_file, 'w') as g:
    # 跳过注释
    next(f)
    # Geneid Chr Start End Strand Length mapped/SRR23955793.bam mapped/SRR23955884.bam mapped/SRR23955887.bam ...
    header = f.readline().strip().split('\t')
    rawnames = header[6:]
    names = [rn.split('/')[-1].replace('.bam', '') for rn in rawnames]
    g.write('Gene\t' + '\t'.join(names) + '\n')

    for line in f:
        lns = line.strip().split('\t')
        gene = lns[0]
        counts = lns[6:]
        g.write(gene + '\t' + '\t'.join(counts) + '\n')
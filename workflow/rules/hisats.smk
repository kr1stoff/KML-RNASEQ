"""
# --rg-id A1 --rg SM:A1 --un-conc ./
hisat2 -p 16 -t \
    -x /data/mengxf/Database/genome/hg38/hg38.fa \
    -1 trimmed/SRR23955793.1.fastq \
    -2 trimmed/SRR23955793.2.fastq \
    --dta | samtools sort -@ 4 -o SRR23955793.bam
"""

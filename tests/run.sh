cd /data/mengxf/Project/KML250113_RNAseq_pipeline/work/250113

# index 从 HISAT2 网站下载
# https://daehwankimlab.github.io/hisat2/download/#h-sapiens
mamba run -n rnaseq hisat2 -t -p 16 --dta -x /data/mengxf/Database/HISAT2/hg38/genome -1 SRR23955793.1.fastq -2 SRR23955793.2.fastq |
    samtools view -hb |
    samtools sort -@ 16 -o SRR23955793.bam

# mamba run -n rnaseq qualimap bamqc -nt 16 -c -gff /data/mengxf/Database/reference/hg38/hg38.ensGene.gtf -bam SRR23955793.bam -outdir qualimap_bamqc

# ! 链特异参数 -p,--sequencing-protocol <arg>
# * 超内存需要 --java-mem-size=30G 参数控制
mamba run -n rnaseq qualimap rnaseq --java-mem-size=30G -pe -bam SRR23955793.bam -gtf /data/mengxf/Database/reference/hg38/hg38.ensGene.gtf -outdir qualimap_rnaseq

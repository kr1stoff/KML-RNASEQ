rule qualimap:
    input:
        bam=rules.hisat2_align.output,
        gtf=config["database"]["gtf"],
    output:
        directory("qc/qualimap_rnaseq/{sample}"),
    conda:
        config["conda"]["rnaseq"]
    log:
        "logs/qualimap/rna-seq/{sample}.log",
    benchmark:
        "logs/qualimap/rna-seq/{sample}.benchmark"
    params:
        # TODO 链特异参数 -p,--sequencing-protocol <arg>
        # * 超内存需要 --java-mem-size=30G 参数控制
        extra="--java-mem-size=30G -pe",
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/qualimap/rnaseq"


# TODO 合并质控结果到一个表格

rule qualimap_rnaseq:
    input:
        bam=rules.hisat2_align.output,
        gtf=config["database"]["gtf"],
    output:
        directory("qc/qualimap_rnaseq/{sample}"),
    conda:
        config["conda"]["qualimap"]
    log:
        "logs/qualimap/rna-seq/{sample}.log",
    benchmark:
        "logs/qualimap/rna-seq/{sample}.benchmark"
    params:
        # todo 链特异参数 -p,--sequencing-protocol <arg>
        # 临时目录需要设置, SH DRAGEN 服务器根目录下 /tmp 太小; 超内存需要 --java-mem-size=30G 参数控制
        java_opts=f"-Djava.io.tmpdir={config['resource']['tmp_dir']}",
        extra=f"-pe --java-mem-size=30G",
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/qualimap/rnaseq"


# todo 合并质控结果到一个表格

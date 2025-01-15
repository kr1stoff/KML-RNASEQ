rule hisat2_align:
    input:
        reads=rules.fastp.output.trimmed,
        idx=multiext(
            config["database"]["hisat2_index"],
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
        ),
    output:
        "mapped/{sample}.bam",
    conda:
        config["conda"]["rnaseq"]
    log:
        "logs/hisat2_align/{sample}.log",
    benchmark:
        "logs/hisat2_align/{sample}.benchmark"
    params:
        # * --dta 为下游组装作准备
        extra="-t --dta",
    threads: config["threads"]["high"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/hisat2/align"

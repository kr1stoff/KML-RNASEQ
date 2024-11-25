rule copy_fq:
    input:
        lambda wildcards: df_sample.loc[wildcards.sample, "fq1"],
    output:
        "rawdata/{sample}.1.fq.gz",
    log:
        "logs/copy_fq_{sample}.log",
    run:
        if input[0].endswith(".gz"):
            shell("cp {input} {output} 2> {log}")
        else:
            shell("gzip -c {input} > {output} 2> {log}")


use rule copy_fq as copy_fq2 with:
    input:
        lambda wildcards: df_sample.loc[wildcards.sample, "fq2"],
    output:
        "rawdata/{sample}.2.fq.gz",
    log:
        "logs/copy_fq_{sample}.log",

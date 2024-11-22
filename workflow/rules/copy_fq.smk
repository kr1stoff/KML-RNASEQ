rule copy_fq:
    input:
        get_copy_fq_input,
    output:
        "rawdata/{sample}.1.fq.gz",
        "rawdata/{sample}.2.fq.gz",
    log:
        "logs/copy_fq_{sample}.log",
    run:
        if input[0].endwith(".gz"):
            shell("cp {input[0]} {output[0]} && cp {input[1]} {output[1]} 2> {log}")
        else:
            shell(
                "gzip -c {input[0]} > {output[0]} && gzip -c {input[1]} > {output[1]} 2> {log}"
            )

rule fastp:
    input:
        sample=[rules.copy_fq.output, rules.copy_fq2.output],
    output:
        trimmed=["trimmed/{sample}.1.fastq", "trimmed/{sample}.2.fastq"],
        html="trimmed/{sample}.html",
        json="trimmed/{sample}.json",
    conda:
        config["conda"]["basic"]
    log:
        "logs/fastp/{sample}.log",
    benchmark:
        "logs/fastp/{sample}.bm"
    params:
        extra="-q 15 -u 40 -l 15 --cut_right --cut_window_size 4 --cut_mean_quality 20 --correction",
    threads: config["threads"]["low"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/fastp"


rule qc_stat:
    input:
        expand("trimmed/{sample}.json", sample=df_sample.index),
    output:
        "trimmed/fastp.stats.tsv",
    benchmark:
        "logs/fastp/qc_stat.bm"
    log:
        "logs/fastp/qc_stat.log",
    script:
        "../scripts/fastp_all_samples_qc.py"

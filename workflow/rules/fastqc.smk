rule fastqc:
    input:
        rules.copy_fq.output,
    output:
        html="qc/fastqc/{sample}.1.html",
        zip="qc/fastqc/{sample}.1_fastqc.zip",
    params:
        extra="--quiet",
    benchmark:
        "logs/fastqc/{sample}.1.bm"
    log:
        "logs/fastqc/{sample}.1.log",
    resources:
        mem_mb=1024,
    conda:
        config["conda"]["basic"]
    threads: config["threads"]["low"]
    # ! 需要在 wrapper 的环境安装 snakemake-wrapper-utils.
    wrapper:
        # * workflow.basedir 是 workflow 的根目录 (来自源码)
        f"file:{workflow.basedir}/wrappers/bio/fastqc"


use rule fastqc as fastqc2 with:
    input:
        rules.copy_fq2.output,
    output:
        html="qc/fastqc/{sample}.2.html",
        zip="qc/fastqc/{sample}.2_fastqc.zip",
    benchmark:
        "logs/fastqc/{sample}.2.bm"
    log:
        "logs/fastqc/{sample}.2.log",


rule multiqc:
    input:
        expand(
            "qc/fastqc/{sample}.{pe}_fastqc.zip", sample=df_sample.index, pe=["1", "2"]
        ),
    output:
        "qc/multiqc.html",
        directory("qc/multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    benchmark:
        "logs/multiqc.bm"
    conda:
        config["conda"]["basic"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/multiqc"

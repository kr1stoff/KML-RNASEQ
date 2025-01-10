rule fastqc:
    input:
        ".rawdata/{sample}_{pe}.fastq.gz",
    output:
        html="qc/fastqc/{sample}_{pe}.html",
        zip="qc/fastqc/{sample}_{pe}_fastqc.zip",
    params:
        extra="--quiet",
    benchmark:
        "logs/fastqc/{sample}_{pe}.bm"
    log:
        "logs/fastqc/{sample}_{pe}.log",
    resources:
        mem_mb=1024,
    conda:
        config["conda"]["basic"]
    threads: config["threads"]["low"]
    # ! 需要在 wrapper 的环境安装 snakemake-wrapper-utils.
    wrapper:
        # * workflow.basedir 是 workflow 的根目录 (来自源码)
        f"file:{workflow.basedir}/wrappers/bio/fastqc"


rule multiqc:
    input:
        expand(
            "qc/fastqc/{sample}_{pe}_fastqc.zip",
            pe=["1", "2"],
            sample=config["samples"],
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

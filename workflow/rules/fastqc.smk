rule fastqc1:
    input:
        get_fastqc_input1,
    output:
        html="qc/fastqc/{sample}_1.html",
        zip="qc/fastqc/{sample}_1_fastqc.zip",
    params:
        extra="--quiet",
    benchmark:
        "logs/fastqc/{sample}_1.bm"
    log:
        "logs/fastqc/{sample}_1.log",
    resources:
        mem_mb=1024,
    conda:
        config["conda"]["basic"]
    threads: config["threads"]["low"]
    # ! 需要在 wrapper 的环境安装 snakemake-wrapper-utils.
    wrapper:
        # * workflow.basedir 是 workflow 的根目录 (来自源码)
        f"file:{workflow.basedir}/wrappers/bio/fastqc"


userule fastqc1 as fastqc2:
    input:
        get_fastqc_input2,
    output:
        html="qc/fastqc/{sample}_2.html",
        zip="qc/fastqc/{sample}_2_fastqc.zip",
    benchmark:
        "logs/fastqc/{sample}_2.bm"
    log:
        "logs/fastqc/{sample}_2.log",


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

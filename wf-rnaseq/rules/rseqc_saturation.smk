from pathlib import Path


rule rseqc_saturation:
    input:
        bam=rules.hisat2_align.output,
        bed=config["database"]["rseqc_bed"],
    output:
        "qc/RPKM_saturation/{sample}.saturation.pdf",
    log:
        "logs/rseqc_saturation/{sample}.log",
    benchmark:
        "logs/rseqc_saturation/{sample}.benchmark"
    run:
        outdir = str(Path(output[0]).parent)
        out_prefix = output[0].replace(".saturation.pdf", "")
        env = config["conda"]["rnaseq"]
        activate = config["conda"]["activate"]
        shell(
            """
            mkdir -p {outdir}
            source {activate} {env}
            RPKM_saturation.py -r {input.bed} -i {input.bam} -o {out_prefix}
            """
        )

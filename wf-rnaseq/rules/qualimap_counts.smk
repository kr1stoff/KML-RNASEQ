rule qualimap_counts_input:
    input:
        genecount=rules.gene_count.output[0],
        metadata=config["metadata"],
    output:
        "qc/qualimap_counts/input.txt",
    log:
        "logs/qualimap_counts/input.log",
    benchmark:
        "logs/qualimap_counts/input.benchmark"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/qualimap_counts_input.py"


rule qualimap_counts:
    input:
        rules.qualimap_counts_input.output,
    output:
        directory("qc/qualimap_counts/output"),
    log:
        "logs/qualimap_counts/output.log",
    benchmark:
        "logs/qualimap_counts/output.benchmark"
    conda:
        config["conda"]["rnaseq"]
    params:
        "-c",
    shell:
        "qualimap counts -d {input} {params} -outdir {output} 2> {log}"

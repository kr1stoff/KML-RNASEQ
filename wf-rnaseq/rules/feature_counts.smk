rule feature_counts:
    input:
        # list of sam or bam files
        # samples="{sample}.bam",
        samples=expand("mapped/{sample}.bam", sample=config["samples"]),
        annotation=config["database"]["gtf"],
    output:
        multiext(
            "feature_counts/all",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    threads: config["threads"]["high"]
    params:
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path="",  # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -J -p -t exon -g gene_id",
    log:
        "logs/feature_counts/all.log",
    benchmark:
        "logs/feature_counts/all.benchmark"
    conda:
        config["conda"]["rnaseq"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/subread/featurecounts"

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


rule gene_count:
    input:
        "feature_counts/all.featureCounts",
    output:
        gene_count="feature_counts/gene_count.tsv",
        withNTC="feature_counts/gene_count_withNTC.tsv",
    log:
        "logs/feature_counts/gene_count.log",
    benchmark:
        "logs/feature_counts/gene_count.benchmark"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/gene_count.py"


rule filter_gene_count:
    input:
        "feature_counts/gene_count.tsv",
    output:
        "feature_counts/gene_count_filtered.tsv",
    log:
        "logs/feature_counts/filter_gene_count.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/filter_gene_count.py"

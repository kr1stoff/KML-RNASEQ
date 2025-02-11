rule calculate_RPKM_TPM:
    input:
        rules.gene_count.output.gene_count,
    output:
        rpkm="feature_counts/rpkm.tsv",
        tpm="feature_counts/tpm.tsv",
    log:
        "logs/feature_counts/rpkm_tpm.log",
    benchmark:
        "logs/feature_counts/rpkm_tpm.benchmark"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/calc_rpkm_tpm.py"


rule filter_tpm:
    input:
        rules.calculate_RPKM_TPM.output.tpm,
    output:
        "feature_counts/tpm_filter.tsv",
    log:
        "logs/feature_counts/filter_tpm.log",
    benchmark:
        "logs/feature_counts/filter_tpm.benchmark"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/filter_tpm_table.py"

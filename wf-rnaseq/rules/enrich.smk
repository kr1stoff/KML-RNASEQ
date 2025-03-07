rule deseq2_deg_go:
    input:
        rules.deseq2_get_deg.output[0],
    output:
        directory("enrich/GO"),
    log:
        "logs/enrich/deseq2_deg_go.log",
    benchmark:
        "logs/enrich/deseq2_deg_go.benchmark"
    conda:
        config["conda"]["R"]
    script:
        "../scripts/deseq2_deg_go.R"


rule deseq_deg_kegg:
    input:
        rules.deseq2_get_deg.output[0],
    output:
        directory("enrich/KEGG"),
    log:
        "logs/enrich/deseq2_deg_kegg.log",
    benchmark:
        "logs/enrich/deseq2_deg_kegg.benchmark"
    conda:
        config["conda"]["R"]
    script:
        "../scripts/deseq2_deg_kegg.R"

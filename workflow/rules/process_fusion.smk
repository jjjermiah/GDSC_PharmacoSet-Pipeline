
rule download_FUSION:
    input:
        gene_fusions = HTTP.remote(molecularProfiles['fusion']['gene_fusions']['url']),
    output:
        gene_fusions = "rawdata/fusion/Fusions_20230725.zip"
    shell:
        "mv {input.gene_fusions} {output.gene_fusions}"

rule preprocess_FUSION:
    input:
        gene_fusions = rules.download_FUSION.output.gene_fusions,
        metadata = rules.preprocess_METADATA.output.metadata
    output:
        preprocessed_fusions = "rawdata/fusion/preprocessed_fusions.qs"
    log:
        "logs/fusion/preprocess_FUSION.log"
    script:
        "../scripts/fusion/preprocess_FUSION.R"

rule make_FUSION_SE:
    input:
        preprocessed_fusions = rules.preprocess_FUSION.output.preprocessed_fusions,
    output:
        fusion_se = "procdata/fusion/fusion_se.qs"
    log:
        "logs/fusion/make_FUSION_SE.log"
    script:
        "../scripts/fusion/make_FUSION_SE.R"
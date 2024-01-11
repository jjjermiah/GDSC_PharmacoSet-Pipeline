
rule download_FUSION:
    input:
        gene_fusions = HTTP.remote(molecularProfiles['fusion']['gene_fusions']['url']),
    output:
        gene_fusions = "rawdata/fusion/Fusions_20230725.zip"
    shell:
        """
        mv {input.gene_fusions} {output.gene_fusions}
        """

rule preprocess_FUSION:
    input:
        gene_fusions = rules.download_FUSION.output.gene_fusions
    output:
        preprocessed_fusions = "rawdata/fusion/preprocessed_fusions.qs"
    script:
        "../scripts/fusion/preprocess_FUSION.R"


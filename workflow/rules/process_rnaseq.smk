################################################################################
## RNA-SEQ

rnaseq_conda_env = "../envs/rnaseq.yaml"
rnaseq = molecularProfiles['rnaseq']

rule download_RNASEQ:
    input:
        processed = HTTP.remote(rnaseq['processed']['url'])
    output:
        processed = "rawdata/rnaseq/rnaseq_all_20220624.zip",
    shell:
        """
        mv {input.processed} {output.processed} 
        """

rule preprocess_RNASEQ:
    input:
        processed = rules.download_RNASEQ.output.processed,
        metadata = rules.preprocess_METADATA.output.metadata
    output:
        preprocessed = "procdata/rnaseq/preprocessed_rnaseq.qs",
    log:
        "logs/rnaseq/preprocess_RNASEQ.log"
    threads:
        3
    conda:
        rnaseq_conda_env
    script:
        "../scripts/rnaseq/preprocess_RNASEQ.R"

rule make_RNASEQ_SE:
    input:
        preprocessed = rules.preprocess_RNASEQ.output.preprocessed,
    output:
        rnaseq_se = "results/data/rnaseq/rnaseq_SE.qs"
    log:     
        "logs/rnaseq/make_RNASEQ_SE.log"
    threads:
        3
    conda:
        rnaseq_conda_env
    script:
        "../scripts/rnaseq/make_RNASEQ_SE.R"
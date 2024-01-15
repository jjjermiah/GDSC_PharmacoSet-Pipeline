################################################################################
## RNA-SEQ

rule download_RNASEQ:
    input:
        all = HTTP.remote(molecularProfiles['rnaseq']['processed']['all']['url'])
    output:
        all = "rawdata/rnaseq/rnaseq_all_20220624.zip",
    shell:
        """
        mv {input.all} {output.all} 
        """

rule preprocess_RNASEQ:
    input:
        all = rules.download_RNASEQ.output.all,
        metadata = rules.preprocess_METADATA.output.metadata
    output:
        preprocessed = "procdata/rnaseq/preprocessed_rnaseq.qs",
    log:
        "logs/rnaseq/preprocess_RNASEQ.log"
    threads:
        3
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
    script:
        "../scripts/rnaseq/make_RNASEQ_SE.R"
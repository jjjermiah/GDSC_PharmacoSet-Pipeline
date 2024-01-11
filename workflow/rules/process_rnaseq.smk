from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

################################################################################
## RNA-SEQ
# 
# note: all, sanger, and broad are included here but only 'all' will be used
# This is because the other two are subsets of the 'all' dataset and does not 
# include as much data

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
        metadata = "procdata/metadata.qs"
    output:
        preprocessed = "procdata/rnaseq/preprocessed_rnaseq.qs",
    log:
        "logs/rnaseq/preprocess_RNASEQ.log"
    threads:
        3
    script:
        "../scripts/rnaseq/process_RNASEQ.R"

rule make_RNASEQ_SE:
    input:
        preprocessed = rules.preprocess_RNASEQ.output.preprocessed,
        metadata = "procdata/metadata.qs"
    output:
        rnaseq_se = "procdata/rnaseq/rnaseq_SE.qs"
    log:     
        "logs/rnaseq/make_RNASEQ_SE.log"
    threads:
        3
    script:
        "../scripts/rnaseq/make_RNASEQ_SE.R"
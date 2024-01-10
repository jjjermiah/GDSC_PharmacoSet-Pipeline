from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

################################################################################
## RNA-SEQ
# 
#' note: all, sanger, and broad are included here but only all will be used

rule download_RNASEQ:
    input:
        all = HTTP.remote(molecularProfiles['rnaseq']['processed']['all']['url']),
        sanger = HTTP.remote(molecularProfiles['rnaseq']['processed']['sanger']['url']),
        broad = HTTP.remote(molecularProfiles['rnaseq']['processed']['broad']['url'])
    output:
        all = "rawdata/rnaseq/rnaseq_all_20220624.zip",
        sanger = "rawdata/rnaseq/rnaseq_sanger_20210316.zip",
        broad = "rawdata/rnaseq/rnaseq_broad_20210317.zip"
    shell:
        """
        mv {input.all} {output.all} && \
        mv {input.sanger} {output.sanger} && \
        mv {input.broad} {output.broad}
        """

rule preprocess_RNASEQ:
    input:
        all = "rawdata/rnaseq/rnaseq_all_20220624.zip",
        metadata = "procdata/metadata.qs"
    output:
        preprocessed = "procdata/rnaseq/preprocessed_rnaseq.qs",
    script:
        "../scripts/process_RNASEQ.R"
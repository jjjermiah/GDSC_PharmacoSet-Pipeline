from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule download_MUTATION:
    input:
        WGS_VCF = HTTP.remote(molecularProfiles['mutation']['WGS_VCF']['url']),
        WES_VCF = HTTP.remote(molecularProfiles['mutation']['WES_VCF']['url']),
        SUMMARY = HTTP.remote(molecularProfiles['mutation']['SUMMARY']['url']),
        Genes_Metadata = HTTP.remote(molecularProfiles['mutation']['Genes_Metadata']['url']),
    output:
        WGS_VCF = "rawdata/mutation/mutations_wgs_vcf_20221123.zip",
        WES_VCF = "rawdata/mutation/mutations_wes_vcf_20221010.zip",
        SUMMARY = "rawdata/mutation/mutations_all_20230202.zip",
        Genes_Metadata = "rawdata/mutation/driver_mutations_20221208.csv",
    shell:
        """
        mv {input.WGS_VCF} {output.WGS_VCF} && \
        mv {input.WES_VCF} {output.WES_VCF} && \
        mv {input.SUMMARY} {output.SUMMARY} && \
        mv {input.Genes_Metadata} {output.Genes_Metadata}
        """

# rule preprocess_MUTATION:
#     input:
#         WGS_VCF = "rawdata/mutation/mutations_wgs_vcf_20221123.zip",
#         WES_VCF = "rawdata/mutation/mutations_wes_vcf_20221010.zip",
#         SUMMARY = "rawdata/mutation/mutations_all_20230202.zip",
#         metadata = "procdata/metadata.qs"
#     output:
#         preprocessed = "procdata/mutation/preprocessed_mutation.qs",
#     script:
#         "../scripts/preprocessMUTATION.R"

rule make_MUTATION_SE:
    input:
        preprocessed = "rawdata/mutation/mutations_all_20230202.zip",
        metadata = "procdata/metadata.qs",
        Genes_Metadata = "rawdata/mutation/driver_mutations_20221208.csv",
    output:
        mutation_SE = "procdata/mutation/mutation_SE.qs",
    script:
        "../scripts/makeMUTATION_SE.R"

rule preprocessCNV:
    input:
        WES_genes = "rawdata/cnv/WES_pureCN_CNV_genes_20221213.csv",
        WES_category = "rawdata/cnv/WES_pureCN_CNV_genes_cn_category_20221213.csv",
        WES_total_cnv = "rawdata/cnv/WES_pureCN_CNV_genes_total_copy_number_20221213.csv",
        WGS_genes = "rawdata/cnv/WGS_purple_CNV_genes_20230303.csv",
        WGS_category = "rawdata/cnv/WGS_purple_genes_cn_category_20230303.csv",
        WGS_total_cnv = "rawdata/cnv/WGS_purple_genes_total_copy_number_20230303.csv",
        metadata = "procdata/metadata.qs",
    output:
        preprocessedCNV = "procdata/cnv/preprocessedCNV.qs",
    threads:
        10
    script:
        "../scripts/preprocessCNV.R"


rule make_CNV_SE:
    input


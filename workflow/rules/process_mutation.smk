
rule download_MUTATION_processed:
    input:
        all_mutations = HTTP.remote(molecularProfiles['mutation']['SUMMARY']['url']),
        Genes_Metadata = HTTP.remote(molecularProfiles['mutation']['Genes_Metadata']['url']),
    output:
        all_mutations = "rawdata/mutation/mutations_all_20230202.zip",
        Genes_Metadata = "rawdata/mutation/driver_mutations_20221208.csv",
    shell:
        """
        mv {input.all_mutations} {output.all_mutations} && \
        mv {input.Genes_Metadata} {output.Genes_Metadata}
        """

rule preprocess_MUTATION:
    input:
        all_mutations = rules.download_MUTATION_processed.output.all_mutations, 
        Genes_Metadata = rules.download_MUTATION_processed.output.Genes_Metadata,
        metadata = rules.preprocess_METADATA.output.metadata,
    output:
        preprocessed = "procdata/mutation/preprocessed_mutation.qs",
    log:
        "logs/mutation/preprocess_MUTATION.log",
    script:
        "../scripts/mutation/preprocess_MUTATION.R"

rule make_MUTATION_SE:
    input:
        preprocessed = rules.preprocess_MUTATION.output.preprocessed,
    output:
        mutation_SE = "procdata/mutation/mutation_SE.qs",
    script:
        "../scripts/mutation/make_MUTATION_SE.R"

# Until figure out how to process VCF files for mutations, use the processed data
# rule download_MUTATION_VCF:
#     input:
#         WGS_VCF = HTTP.remote(molecularProfiles['mutation']['WGS_VCF']['url']),
#         WES_VCF = HTTP.remote(molecularProfiles['mutation']['WES_VCF']['url']),
#     output:
#         WGS_VCF = "rawdata/mutation/mutations_wgs_vcf_20221123.zip",
#         WES_VCF = "rawdata/mutation/mutations_wes_vcf_20221010.zip",
#     shell:
#         """
#         mv {input.WGS_VCF} {output.WGS_VCF} && \
#         mv {input.WES_VCF} {output.WES_VCF}
#         """


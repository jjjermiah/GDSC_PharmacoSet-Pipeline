cnv_conda_env = "../envs/cnv.yaml"
cnv = molecularProfiles['cnv']


rule download_CNV_WESData:
    input:
        WES = HTTP.remote(cnv['WES_CNV']['url'])
    output:
        WES_genes = "rawdata/cnv/WES_pureCN_CNV_genes_20221213.csv",
        WES_category = "rawdata/cnv/WES_pureCN_CNV_genes_cn_category_20221213.csv",
        WES_total_cnv = "rawdata/cnv/WES_pureCN_CNV_genes_total_copy_number_20221213.csv",
    log:
        "logs/cnv/downloadCNV_WESData.log"
    shell:
        """
        unzip -d $(dirname {output.WES_genes}) {input.WES} && \
        rm {input.WES} > {log} 2>&1
        """

rule preprocess_CNV:
    input:
        WES_genes = "rawdata/cnv/WES_pureCN_CNV_genes_20221213.csv",
        WES_category = "rawdata/cnv/WES_pureCN_CNV_genes_cn_category_20221213.csv",
        WES_total_cnv = "rawdata/cnv/WES_pureCN_CNV_genes_total_copy_number_20221213.csv",
        metadata = rules.preprocess_METADATA.output.metadata
    output:
        preprocessedCNV = "procdata/cnv/preprocessedCNV.qs",
    log:
        "logs/cnv/preprocessCNV.log",
    conda:
        cnv_conda_env,
    threads:
        6
    script:
        "../scripts/cnv/preprocess_CNV.R"

rule make_CNV_SE:
    input:
        preprocessedCNV = rules.preprocess_CNV.output.preprocessedCNV,
    output:
        CNV_se = "results/data/cnv/CNV_SE.qs",
    log:
        "logs/cnv/make_CNV_SE.log",
    conda:
        cnv_conda_env,
    script:
        "../scripts/cnv/make_CNV_SE.R"


# rule downloadCNV_WGSData:
#     input:
#         WGS = HTTP.remote(molecularProfiles['cnv']['WGS_CNV']['url']),
#     output:
#         WGS_genes = "rawdata/cnv/WGS_purple_CNV_genes_20230303.csv",
#         WGS_category = "rawdata/cnv/WGS_purple_genes_cn_category_20230303.csv",
#         WGS_total_cnv = "rawdata/cnv/WGS_purple_genes_total_copy_number_20230303.csv",
#     shell:
#         "unzip -d $(dirname {output.WGS_genes}) {input.WGS}; rm {input.WGS}"
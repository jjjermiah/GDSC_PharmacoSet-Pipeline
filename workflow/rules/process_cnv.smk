
# rule downloadCNV_WGSData:
#     input:
#         WGS = HTTP.remote(molecularProfiles['cnv']['WGS_CNV']['url']),
#     output:
#         WGS_genes = "rawdata/cnv/WGS_purple_CNV_genes_20230303.csv",
#         WGS_category = "rawdata/cnv/WGS_purple_genes_cn_category_20230303.csv",
#         WGS_total_cnv = "rawdata/cnv/WGS_purple_genes_total_copy_number_20230303.csv",
#     shell:
#         "unzip -d $(dirname {output.WGS_genes}) {input.WGS}; rm {input.WGS}"

rule downloadCNV_WESData:
    input:
        WES = HTTP.remote(molecularProfiles['cnv']['WES_CNV']['url'])
    output:
        WES_genes = "rawdata/cnv/WES_pureCN_CNV_genes_20221213.csv",
        WES_category = "rawdata/cnv/WES_pureCN_CNV_genes_cn_category_20221213.csv",
        WES_total_cnv = "rawdata/cnv/WES_pureCN_CNV_genes_total_copy_number_20221213.csv",
    shell:
        "unzip -d $(dirname {output.WES_genes}) {input.WES}; rm {input.WES}"

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
    threads:
        6
    script:
        "../scripts/cnv/preprocess_CNV.R"

rule make_CNV_SE:
    input:
        preprocessedCNV = rules.preprocess_CNV.output.preprocessedCNV,
    output:
        CNV_se = "procdata/cnv/CNV_SE.qs",
    log:
        "logs/cnv/make_CNV_SE.log",
    script:
        "../scripts/cnv/make_CNV_SE.R"
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import json
import pandas as pd
HTTP = HTTPRemoteProvider()

molecularProfiles = config['molecularProfiles']



# rule downloadMUTATION_WGSData:
#     input:
#         WGS = HTTP.remote(molecularProfiles['mutation']['WGS_VCF']['url']),
#     output:
#         WGS = directory("rawdata/mutation/WGS_VCF"),
#     shell:
#         "unzip -d {output.WGS} {input.WGS}" 
# TODO:: work with the zip file in R scripts instead of unzipping it all the time

# rule downloadMUTATION_WESData:
#     input:
#         WES = HTTP.remote(molecularProfiles['mutation']['WES_VCF']['url']),
#     output:
#         WES = directory("rawdata/mutation/WES_VCF"),
#     shell:
#         "unzip -d {output.WES} {input.WES}"

# rule downloadMUTATION_summary:
#     input:
#         summary = HTTP.remote(molecularProfiles['mutation']['SUMMARY']),
#     output:
#         summary = "rawdata/mutation/mutations_all_20230202.zip",
#     shell:
#         "unzip -d {output.summary} {input.summary}"

# rule downloadCNV_WGSData:
#     input:
#         WGS = HTTP.remote(molecularProfiles['cnv']['WGS_CNV']['url']),
#     output:
#         WGS_genes = "rawdata/cnv/WGS_purple_CNV_genes_20230303.csv",
#         WGS_category = "rawdata/cnv/WGS_purple_genes_cn_category_20230303.csv",
#         WGS_total_cnv = "rawdata/cnv/WGS_purple_genes_total_copy_number_20230303.csv",
#     shell:
#         "unzip -d $(dirname {output.WGS_genes}) {input.WGS}; rm {input.WGS}"

# rule downloadCNV_WESData:
#     input:
#         WES = HTTP.remote(molecularProfiles['cnv']['WES_CNV'])
#     output:
#         WES_genes = "rawdata/cnv/WES_pureCN_CNV_genes_20221213.csv",
#         WES_category = "rawdata/cnv/WES_pureCN_CNV_genes_cn_category_20221213.csv",
#         WES_total_cnv = "rawdata/cnv/WES_pureCN_CNV_genes_total_copy_number_20221213.csv",
#     shell:
#         "unzip -d $(dirname {output.WES_genes}) {input.WES}; rm {input.WES}"



# rule downloadMethylationData:
#     input:
#         processed = HTTP.remote(molecularProfiles['methylation']['processed']),
#         signalIntensity = HTTP.remote(molecularProfiles['methylation']['signalIntensity']),
#     output:
#         processed = "rawdata/methylation/GSE68379_Matrix.processed.txt",
#         signalIntensity = "rawdata/methylation/GSE68379_Matrix.Signal.Intensities.txt",
#     shell:
#         "gunzip -c {input.processed} > {output.processed} && gunzip -c {input.signalIntensity} > {output.signalIntensity}"

# rule downloadMETHYLATION:
#     input:
#         meth_cell_data = HTTP.remote(molecularProfiles['methylation']['meth_cell_data'])
#     output:
#         meth_cell_data = "rawdata/methylation/METH_CELL_DATA.txt",
#     shell:
#         "unzip -p {input.meth_cell_data} F2_METH_CELL_DATA.txt > {output.meth_cell_data} && rm {input.meth_cell_data}" # unzip, extract the file and remove the zip file



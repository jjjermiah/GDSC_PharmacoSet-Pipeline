from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import json
import pandas as pd
HTTP = HTTPRemoteProvider()






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



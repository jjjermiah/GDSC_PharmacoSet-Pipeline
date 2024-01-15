
methylation = molecularProfiles['methylation']

rule download_METHYLATION:
    input: 
        meth_cell_data = HTTP.remote(methylation['meth_cell_data']['url']),
        sample_mapping_metadata = HTTP.remote(methylation['sample_mapping_metadata']['url']),
    output:
        meth_cell_data = "rawdata/methylation/METH_CELL_DATA.txt.zip",
        sample_mapping_metadata = "rawdata/methylation/methSampleId_2_cosmicIds.xlsx",
    shell:
        """
        mv {input.meth_cell_data} {output.meth_cell_data} && \
        mv {input.sample_mapping_metadata} {output.sample_mapping_metadata}
        """

rule preprocess_METHYLATION:
    input:
        meth_cell_data = rules.download_METHYLATION.output.meth_cell_data,
        sample_mapping_metadata = rules.download_METHYLATION.output.sample_mapping_metadata,
        metadata = rules.preprocess_METADATA.output.metadata
    output:
        preprocessed = "procdata/methylation/preprocessed_methylation.qs"
    log:
        "logs/methylation/preprocess_METHYLATION.log"
    script:
        "../scripts/methylation/preprocess_METHYLATION.R"

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


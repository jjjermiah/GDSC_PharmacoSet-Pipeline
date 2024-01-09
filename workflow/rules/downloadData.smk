from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import json
import pandas as pd
HTTP = HTTPRemoteProvider()

molecularProfiles = config['molecularProfiles']

################################################################################
## MICROARRAY 
# todo: this is a pretty exhaustive method to download the data. 
# think of a better solution to download the data. 
rule download_MicroArrayMetadata:
    input:
        srdf = HTTP.remote(molecularProfiles['microarray']['metadata_srdf']['url']),
        fileList = HTTP.remote(molecularProfiles['microarray']['metadata_json']['url']),
        # tsv = HTTP.remote(molecularProfiles['microarray']['metadata_tsv']['url'])
    output:
        srdf = "metadata/microarray/E-MTAB-3610.sdrf.txt",
        filelist = "metadata/microarray/E-MTAB-3610_filelist.json",
        # tsv = "metadata/microarray/E-MTAB-3610.tsv",
    shell:
        """
        wget -O {output.srdf} {input.srdf} && \
        wget -O {output.filelist} {input.fileList} 
        """

# This checkpoint rule is used to parse the metadata file and save it as a json file
# The json file is used to download the expression data in downloadExpressionData
checkpoint load_MicroArrayMetadata:
    input:
        filelist = "metadata/microarray/E-MTAB-3610_filelist.json",
    output:
        microarrayFiles = "metadata/microarray/E-MTAB-3610_expressionFiles.json",
    run:
        import json
        with open(input.filelist) as f:
            filelist = json.load(f)
        
        microarrayFiles = dict()

        for i, file in enumerate(filelist):
            samplename = file['attributes'][0]['value']
            microarrayFiles[samplename] = dict()

            microarrayFiles[samplename]['filename'] = file['path']
            microarrayFiles[samplename]['size'] = file['size']
            microarrayFiles[samplename]['description'] = file['attributes'][1]['value']
       
        # save to json file
        with open(output.microarrayFiles, 'w') as f:
            json.dump(microarrayFiles, f, indent=4)

# input function to get the files from the FTP server and save them locally
def getMicroArrayFiles(wildcards):
    basepath = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/610/E-MTAB-3610/Files/"
    
    with checkpoints.load_MicroArrayMetadata.get().output[0].open() as f:
        arrayFiles = json.load(f)
    
    # get the first 10 keys of the dictionary
    files = dict(list(arrayFiles.items()))

    # append filename to basepath for each sample and return 
    ftpFilePaths = [basepath + arrayFiles[sample]['filename'] for sample in files]

    # return HTTP.remote(ftpFilePaths)
    return ["rawdata/microarray/" + arrayFiles[sample]['filename'] for sample in files]

rule preprocess_MicroArray:
    input: 
        CELfiles = getMicroArrayFiles,
        CEL_metadata = "metadata/microarray/E-MTAB-3610.sdrf.txt",
        CEL_FileList = "metadata/microarray/E-MTAB-3610_expressionFiles.json",
    output:
        preprocessedMicroArray = "procdata/expression/preprocessedMicroArray.qs",
    threads:
        10
    script:
        "../scripts/preprocess_MICROARRAY.R"

rule download_MicroArrayFILE: 
    output:
        "rawdata/microarray/{sample}.cel"
    log:
        "logs/microarray/Download_{sample}.log"
    retries: 5
    shell:
        """
        ftpFilePath="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/610/E-MTAB-3610/Files"
        wget -O {output} $ftpFilePath/{wildcards.sample}.cel > {log} 2>&1
        """

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

rule process_RNASEQ:
    input:
        all = "rawdata/rnaseq/rnaseq_all_20220624.zip",
        sampleMetadata = "metadata/sampleAnnotation.xlsx",
        cmp_Metadata = "metadata/cellModelPassportsAnnotation.csv",
        cmp_GeneAnnotation = "metadata/cellmodelpassportsGeneAnnotation.csv",
    output:
        processed = "procdata/expression/rnaseq_processed.qs",
    script:
        "../scripts/process_RNASEQ.R"

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
#         summary = "rawdata/mutation/mutations_all_20230202.csv",
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



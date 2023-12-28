
rule downloadAllData:
    input:
        "rawdata/mutation/WGS_VCF",
        "rawdata/mutation/WES_VCF",
        "rawdata/mutation/mutations_all_20230202.csv",
        "rawdata/cnv/WGS_purple_CNV_genes_20230303.csv",
        "rawdata/cnv/WGS_purple_genes_cn_category_20230303.csv",
        "rawdata/cnv/WGS_purple_genes_total_copy_number_20230303.csv",
        "rawdata/cnv/WES_pureCN_CNV_genes_20221213.csv",
        "rawdata/cnv/WES_pureCN_CNV_genes_cn_category_20221213.csv",
        "rawdata/cnv/WES_pureCN_CNV_genes_total_copy_number_20221213.csv",
        "rawdata/methylation/GSE68379_Matrix.processed.txt",
        "rawdata/methylation/GSE68379_Matrix.Signal.Intensities.txt",
        "rawdata/methylation/METH_CELL_DATA.txt",
        "metadata/expression/E-MTAB-3610.sdrf.txt",
        "metadata/expression/E-MTAB-3610_filelist.json",



rule downloadMUTATION_WGSData:
    input:
        WGS = HTTP.remote(config['rawdata']['mutation']['WGS_VCF']),
    output:
        WGS = directory("rawdata/mutation/WGS_VCF"),
    shell:
        "unzip -d {output.WGS} {input.WGS}" # TODO:: work with the zip file in R scripts instead of unzipping it all the time

rule downloadMUTATION_WESData:
    input:
        WES = HTTP.remote(config['rawdata']['mutation']['WES_VCF']),
    output:
        WES = directory("rawdata/mutation/WES_VCF"),
    shell:
        "unzip -d {output.WES} {input.WES}"

rule downloadMUTATION_summary:
    input:
        summary = HTTP.remote(config['rawdata']['mutation']['SUMMARY']),
    output:
        summary = "rawdata/mutation/mutations_all_20230202.csv",
    shell:
        "unzip -d {output.summary} {input.summary}"

rule downloadCNV_WGSData:
    input:
        WGS = HTTP.remote(config['rawdata']['cnv']['WGS_CNV']),
    output:
        WGS_genes = "rawdata/cnv/WGS_purple_CNV_genes_20230303.csv",
        WGS_category = "rawdata/cnv/WGS_purple_genes_cn_category_20230303.csv",
        WGS_total_cnv = "rawdata/cnv/WGS_purple_genes_total_copy_number_20230303.csv",
    shell:
        "unzip -d $(dirname {output.WGS_genes}) {input.WGS}; rm {input.WGS}"

rule downloadCNV_WESData:
    input:
        WES = HTTP.remote(config['rawdata']['cnv']['WES_CNV'])
    output:
        WES_genes = "rawdata/cnv/WES_pureCN_CNV_genes_20221213.csv",
        WES_category = "rawdata/cnv/WES_pureCN_CNV_genes_cn_category_20221213.csv",
        WES_total_cnv = "rawdata/cnv/WES_pureCN_CNV_genes_total_copy_number_20221213.csv",
    shell:
        "unzip -d $(dirname {output.WES_genes}) {input.WES}; rm {input.WES}"



rule downloadMethylationData:
    input:
        processed = HTTP.remote(config['rawdata']['methylation']['processed']),
        signalIntensity = HTTP.remote(config['rawdata']['methylation']['signalIntensity']),
    output:
        processed = "rawdata/methylation/GSE68379_Matrix.processed.txt",
        signalIntensity = "rawdata/methylation/GSE68379_Matrix.Signal.Intensities.txt",
    shell:
        "gunzip -c {input.processed} > {output.processed} && gunzip -c {input.signalIntensity} > {output.signalIntensity}"

rule downloadMETHYLATION:
    input:
        meth_cell_data = HTTP.remote(config['rawdata']['methylation']['meth_cell_data'])
    output:
        meth_cell_data = "rawdata/methylation/METH_CELL_DATA.txt",
    shell:
        "unzip -p {input.meth_cell_data} F2_METH_CELL_DATA.txt > {output.meth_cell_data} && rm {input.meth_cell_data}" # unzip, extract the file and remove the zip file


rule loadExpressionMetadata:
    input:
        filelist = "metadata/expression/E-MTAB-3610_filelist.json",
    output:
        expressionFiles = "metadata/expression/E-MTAB-3610_expressionFiles.json",
    run:
        import json
        with open(input.filelist) as f:
            filelist = json.load(f)
        
        expressionFiles = dict()

        for i, file in enumerate(filelist):
            samplename = file['attributes'][0]['value']
            expressionFiles[samplename] = dict()

            expressionFiles[samplename]['filename'] = file['path']
            expressionFiles[samplename]['size'] = file['size']
            expressionFiles[samplename]['description'] = file['attributes'][1]['value']

            if i == 0:
                print(expressionFiles)
        
        # save to json file
        with open(output.expressionFiles, 'w') as f:
            json.dump(expressionFiles, f, indent=4)



# input function to get the files from the FTP server and save them locally
def getExpressionFiles(wildcards):
    basepath = "https://ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/610/E-MTAB-3610/Files/"
    
    with open("metadata/expression/E-MTAB-3610_expressionFiles.json") as f:
        expressionFiles = json.load(f)
    
    # get the first 10 keys of the dictionary
    files = dict(list(expressionFiles.items())[:10])

    # append filename to basepath for each sample and return 
    ftpFilePaths = [basepath + expressionFiles[sample]['filename'] for sample in files]

    return HTTP.remote(ftpFilePaths)

rule downloadExpressionData:
    input:
        expressionFile = getExpressionFiles,
        expressionFileList = "metadata/expression/E-MTAB-3610_expressionFiles.json",
    output:
        CELfiles = directory("rawdata/expression/"),
    shell:
        # wget into the directory CELfiles
        "wget -P {output.CELfiles} {input.expressionFile}"


rule downloadExpressionMetadata:
    input:
        srdf = HTTP.remote("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/610/E-MTAB-3610/Files/E-MTAB-3610.sdrf.txt"),
        fileList = HTTP.remote("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/610/E-MTAB-3610/Files/raw-data_filelist.json"),
    output:
        SRDF = "metadata/expression/E-MTAB-3610.sdrf.txt",
        filelist = "metadata/expression/E-MTAB-3610_filelist.json",
    shell:
        "wget -O {output.SRDF} {input.srdf} && wget -O {output.filelist} {input.fileList}"
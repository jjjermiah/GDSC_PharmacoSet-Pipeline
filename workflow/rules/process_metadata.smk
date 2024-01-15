from pathlib import Path

metadata = config['metadata']
metadata_conda_env = "../envs/metadata.yaml"

rule downloadSampleMetadata:
    input:
        sampleMetadata = HTTP.remote(metadata['sampleAnnotation'])
    output:
        sampleMetadata = "metadata/sampleAnnotation.xlsx",
    shell:
        "mv {input.sampleMetadata} {output.sampleMetadata}"

rule downloadTreatmentMetadata:
    input:
        treatmentMetadata = HTTP.remote(metadata['treatmentAnnotation'])
    output:
        treatmentMetadata = "metadata/treatmentAnnotation.csv",
    shell:
        "mv {input.treatmentMetadata} {output.treatmentMetadata}"

rule downloadCellModelPassportsMetadata:
    input:
        cellModelPassportsMetadata = HTTP.remote(metadata['cellmodelpassportsAnnotation']),
        cellmodelpassportsGeneAnnotation = HTTP.remote(metadata['cellmodelpassportsGeneAnnotation'])
    output:
        cellModelPassportsMetadata = "metadata/cellModelPassportsAnnotation.csv",
        cellmodelpassportsGeneAnnotation = "metadata/cellmodelpassportsGeneAnnotation.csv",
    shell:
        """
        mv {input.cellModelPassportsMetadata} {output.cellModelPassportsMetadata} && \
        mv {input.cellmodelpassportsGeneAnnotation} {output.cellmodelpassportsGeneAnnotation}
        """

module gencode:
    snakefile:
        # "https://github.com/jjjermiah/PharmacoSet-Pipelines/blob/main/workflow/rules/downloadGencode.smk"
        github("jjjermiah/PharmacoSet-Pipelines", path = "workflow/rules/downloadGencode.smk", branch = "main")
    
use rule * from gencode as gencode_*

rule preprocess_METADATA:
    input:
        sampleMetadata = "metadata/sampleAnnotation.xlsx",
        treatmentMetadata = "metadata/treatmentAnnotation.csv",
        cellModelPassportsMetadata = "metadata/cellModelPassportsAnnotation.csv",
        cellmodelpassportsGeneAnnotation = "metadata/cellmodelpassportsGeneAnnotation.csv",
        gencode_annotation_file = gencode.gencodeAnnotation(
            dirPath = "metadata/references",
            ref_build = config["referenceGenome"]["Gencode"]["build"],
            gencode_ver = config["referenceGenome"]["Gencode"]["version"]
        ) 
    output:
        metadata = "procdata/metadata/metadata.qs"
    log:
        "logs/metadata/preprocess_METADATA.log"
    conda:
        metadata_conda_env
    script:
        "../scripts/metadata/preprocessMetadata.R"




rule downloadSampleMetadata:
    input:
        sampleMetadata = HTTP.remote(config['metadata']['sampleAnnotation'])
    output:
        sampleMetadata = "metadata/sampleAnnotation.xlsx",
    shell:
        "mv {input.sampleMetadata} {output.sampleMetadata}"

rule downloadTreatmentMetadata:
    input:
        treatmentMetadata = HTTP.remote(config['metadata']['treatmentAnnotation'])
    output:
        treatmentMetadata = "metadata/treatmentAnnotation.csv",
    shell:
        "mv {input.treatmentMetadata} {output.treatmentMetadata}"

rule downloadCellModelPassportsMetadata:
    input:
        cellModelPassportsMetadata = HTTP.remote(config['metadata']['cellmodelpassportsAnnootation']),
        cellmodelpassportsGeneAnnotation = HTTP.remote(config['metadata']['cellmodelpassportsGeneAnnotation'])
    output:
        cellModelPassportsMetadata = "metadata/cellModelPassportsAnnotation.csv",
        cellmodelpassportsGeneAnnotation = "metadata/cellmodelpassportsGeneAnnotation.csv",
    shell:
        "mv {input.cellModelPassportsMetadata} {output.cellModelPassportsMetadata} && mv {input.cellmodelpassportsGeneAnnotation} {output.cellmodelpassportsGeneAnnotation}"

rule preprocessMetadata:
    input:
        sampleMetadata = "metadata/sampleAnnotation.xlsx",
        treatmentMetadata = "metadata/treatmentAnnotation.csv",
        cellModelPassportsMetadata = "metadata/cellModelPassportsAnnotation.csv",
        cellmodelpassportsGeneAnnotation = "metadata/cellmodelpassportsGeneAnnotation.csv",
    output:
        metadata = "procdata/metadata.qs"
    script:
        "../scripts/preprocessMetadata.R"

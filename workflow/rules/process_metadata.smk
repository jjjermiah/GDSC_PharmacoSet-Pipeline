
metadata = config['metadata']

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

rule preprocess_METADATA:
    input:
        sampleMetadata = "metadata/sampleAnnotation.xlsx",
        treatmentMetadata = "metadata/treatmentAnnotation.csv",
        cellModelPassportsMetadata = "metadata/cellModelPassportsAnnotation.csv",
        cellmodelpassportsGeneAnnotation = "metadata/cellmodelpassportsGeneAnnotation.csv",
    output:
        metadata = "procdata/metadata/metadata.qs"
    script:
        "../scripts/preprocessMetadata.R"



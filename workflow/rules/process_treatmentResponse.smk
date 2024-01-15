
treatmentResponse = config["treatmentResponse"]

dataset = "GDSC2"
release = 8.5

def getTreatmentResponseURLs(wildcards):
    rawdata = HTTP.remote(treatmentResponse[float(release)][dataset]["rawdata"]["url"])
    processed = HTTP.remote(treatmentResponse[float(release)][dataset]["processed"]["url"])
    return {"rawdata": rawdata, "processed": processed}

# rule download_treatmentResponse:
#     input:
#         # rawdata = HTTP.remote(treatmentResponse[{release}][{dataset}]["rawdata"]["url"]),
#         # processed = HTTP.remote(treatmentResponse[{release}][{dataset}]["processed"]["url"])
#         unpack(getTreatmentResponseURLs),
#     output:
#         rawdata = "rawdata/treatmentResponse/release_{release}/{dataset}_public_raw_data_{date}.csv.zip",
#         processed = "rawdata/treatmentResponse/release_{release}/{dataset}_fitted_dose_response_{date}.xlsx"
#     # log:
#         # "logs/treatmentResponse/{release}/{dataset_}_{date}/download_treatmentResponse.log"
#     shell:
#         """
#         mv {input.rawdata} {output.rawdata} && \
#         mv {input.processed} {output.processed}
#         """


# def treatmentResponseFilePaths(wildcards):
#     rawdata = "rawdata/treatmentResponse/release_{release}/{dataset}_public_raw_data_{date}.csv.zip"
#     processed = "rawdata/treatmentResponse/release_{release}/{dataset}_fitted_dose_response_{date}.xlsx"
#     return {"rawdata": rawdata, "processed": processed}

rule build_treatmentResponseExperiment:
    input:
        unpack(getTreatmentResponseURLs),
        metadata = rules.preprocess_METADATA.output.metadata
    output:
        tre = "results/data/treatmentResponse/treatmentResponseExperiment.qs"
    log:
        "logs/treatmentResponse/build_treatmentResponseExperiment.log"
    script:
        "../scripts/treatmentResponse/build_treatmentResponseExperiment.R"

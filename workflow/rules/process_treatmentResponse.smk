
treatmentResponse = config["treatmentResponse"]
treatmentResponse_conda_env = "../envs/treatmentResponse.yaml"



# This function is a Input Function (snakemake feature) that uses the wildcards
# to get the URLs of treatmentResponse data based on the config file 
# (config["treatmentResponse"])
# The input for download_treatmentResponse will then have named inputs of :
# `rawdata and `processed
def getTreatmentResponseURLs(wildcards):
    release = float(wildcards.release)
    dataset = wildcards.dataset
    rawdata = treatmentResponse[release][dataset]["rawdata"]["url"]
    processed = treatmentResponse[release][dataset]["processed"]["url"]
    return {"rawdata": HTTP.remote(rawdata), "processed": HTTP.remote(processed)}

rule download_treatmentResponse:
    input:
        unpack(getTreatmentResponseURLs),
    output:
        rawdata = "rawdata/treatmentResponse/release_{release}/{dataset}_public_raw_data.csv",
        processed = "rawdata/treatmentResponse/release_{release}/{dataset}_fitted_dose_response.xlsx"
    log:
        "logs/treatmentResponse/{dataset}_{release}/download_treatmentResponse.log"
    shell:
        """
        mv {input.rawdata} {output.rawdata} && \
        mv {input.processed} {output.processed} > {log} 2>&1
        """

rule preprocess_treatmentResponse:
    input:
        rawdata = "rawdata/treatmentResponse/release_{release}/{dataset}_public_raw_data.csv",
        processed = "rawdata/treatmentResponse/release_{release}/{dataset}_fitted_dose_response.xlsx",
        metadata = "procdata/metadata/metadata.qs"
    output:
        preprocessed = "procdata/treatmentResponse/{dataset}_{release}_treatmentResponse_preprocessed.qs"
    log:
        "logs/treatmentResponse/{dataset}_{release}/preprocess_treatmentResponse.log"
    conda:
        treatmentResponse_conda_env
    threads:
        4
    script:
        "../scripts/treatmentResponse/preprocess_treatmentResponse.R"

rule build_treatmentResponseExperiment:
    input:
        preprocessed = "procdata/treatmentResponse/{dataset}_{release}_treatmentResponse_preprocessed.qs"
    output:
        tre = "results/data/treatmentResponse/{dataset}_{release}_treatmentResponseExperiment.qs"
    log:
        "logs/treatmentResponse/build_treatmentResponseExperiment_{dataset}_{release}.log"
    conda:
        treatmentResponse_conda_env
    threads:
        4
    script:
        "../scripts/treatmentResponse/build_treatmentResponseExperiment.R"


rule fit_treatmentResponseExperiment:
    input:
        # preprocessed = "procdata/treatmentResponse/{dataset}_{release}_treatmentResponse_preprocessed.qs",
        tre = rules.build_treatmentResponseExperiment.output.tre
    output:
        tre = "results/data/treatmentResponse/{dataset}_{release}_treatmentResponseExperiment_fit.qs"
    conda:
        treatmentResponse_conda_env
    log:
        "logs/treatmentResponse/fit_treatmentResponseExperiment_{dataset}_{release}.log"
    threads:
        50
    resources:
        mem_mb = 60000
    retries: 2
    script:
        "../scripts/treatmentResponse/fit_treatmentResponseExperiment.R"

rule build_all_treatmentResponseExperiment:
    input:
        tre = expand(
            "results/data/treatmentResponse/{dataset}_{release}_treatmentResponseExperiment_fit.qs",
            # dataset = ["GDSC2"],
            # release = [8.5]
            dataset =  ["GDSC1", "GDSC2"],
            release = [8.2, 8.4, 8.5]
            # release = [8.5] # [8.2, 8.4, 8.5]
        )
    localrule: True
    output:
        tre = "results/data/treatmentResponse/all_treatmentResponseExperiment.qs"
#' RULE: build_PharmacoSet
#' AUTHOR: Jermiah Joseph
#' DATE: 01-15-2024
#' This script takes in the following files:
#' - INPUT$metadata
#' - INPUT$summarizedExperiments
#' and outputs the following files:
#' - OUTPUT$pset
#' 
#' Libraries Used:
#' - MultiAssayExperiment
#' - log4r
#' - BiocParallel
#' - qs
#' - SummarizedExperiment
#' - CoreGx
#' - PharmacoGx

## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}

suppressPackageStartupMessages(library(PharmacoGx))

# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting build_PharmacoSet.R\n")


# 0.2 Read in the metadata
# --------------------------
info(paste("Loading: ", INPUT$metadata, sep = "\n\t"))
metadata <- qs::qread(INPUT$metadata, nthreads = THREADS)


# 0.3 Read in the summarized experiments
# --------------------------------------
# Read the summarized experiments
info(paste("Loading: ", INPUT$summarizedExperiments, sep = "\n\t"))
se_list <- unlist(lapply(INPUT$summarizedExperiments, qs::qread, nthreads = THREADS))

# 0.4 Read in treatmentResponseExperiment
# ----------------------------------------
info(paste("Loading: ", INPUT$tre, sep = "\n\t"))
tre <- qs::qread(INPUT$tre, nthreads = THREADS)
treatment <- metadata$treatment
treatment <- treatment[(DRUG_NAME %in% rowData(tre)$DRUG_NAME) & !duplicated(treatment$DRUG_NAME), ]
treatment <- data.frame(treatment, row.names = treatment$DRUG_NAME)

# 1.0 Build MultiAssayExperiment
# ------------------------------
# Extract unique sample IDs from the summarized experiments
sampleids_all <- lapply(se_list, colnames)

info(paste(
    "Number of samples in each experiment:\n", 
    paste(capture.output(lapply(sampleids_all, length)), collapse = "\n"),
    sep = ""))

sampleid <- unique(unlist(sampleids_all))
sample <- metadata$sample
# Check that all samples are in the metadata
stopifnot(all(sampleid %in% sample$sampleid))
# Subset the metadata to only include samples in the summarized experiments
sample <- sample[sample$sampleid %in% sampleid, ]

# Remove duplicate sample IDs
sample <- sample[!duplicated(sample$sampleid), ]

# convert sample into a data frame with the rownames being the sample IDs
# and ordered by the sample IDs
sample <- data.frame(sample, row.names = sample$sampleid)
sample <- sample[order(rownames(sample)), ]

info(sprintf("Total number of samples across all experiments: %d", nrow(sample)))

# Create a data frame for the column data, including sample IDs and batch IDs
colData <- data.frame(
    sampleid = sampleid,
    batchid = rep(NA, length(sampleid)),
    row.names = sampleid
)
info(sprintf("Column data has %d rows and %d columns", nrow(colData), ncol(colData)))

# Create an ExperimentList object from the filtered summarized experiments
ExpList <- MultiAssayExperiment::ExperimentList(se_list)
info(paste("ExperimentList:\n", capture.output(show(ExpList)), sep = ""))

# Create a sample map for each experiment in the ExperimentList
sampleMapList <- lapply(se_list, function(se){
    data.frame(
        primary = colnames(se),
        colname = colnames(se),
        stringsAsFactors = FALSE
    )
})
names(sampleMapList) <- names(ExpList)
info(paste("Sample map list:\n", capture.output(str(sampleMapList)), sep = ""))

# Convert the sample map list to a single sample map
mae_sampleMap <- MultiAssayExperiment::listToMap(sampleMapList)
info(paste("Sample map:\n", capture.output(str(mae_sampleMap)), sep = ""))

# Create a MultiAssayExperiment object with the ExperimentList, column data, and sample map
mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = ExpList,
    colData = colData,
    sampleMap = MultiAssayExperiment::listToMap(sampleMapList)
)
info(paste("MultiAssayExperiment:\n", capture.output(show(mae)), sep = ""))

pset <- PharmacoGx::PharmacoSet2(
    name = "GDSC",
    treatment = treatment,
    sample = sample,
    molecularProfiles = mae,
    treatmentResponse = tre,
    perturbation = list(),
    curation = list(sample = data.frame(), treatment = data.frame(), tissue = data.frame()),
    datasetType = "sensitivity"
)


# OUTPUT
# ------


qs::qsave(pset, file = OUTPUT$pset, nthreads = THREADS)



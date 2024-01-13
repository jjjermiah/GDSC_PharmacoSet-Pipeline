## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    # Extract input, output, wildcards, threads, and log file information from the snakemake object
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(CoreGx))
# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting build_PharmacoSet.R\n")

# Read the summarized experiments
info(paste("Loading: ", INPUT$summarizedExperiments, sep = "\n\t"))
se_list <- unlist(lapply(INPUT$summarizedExperiments, qs::qread, nthreads = THREADS))

# Extract unique sample IDs from the summarized experiments
sampleids_all <- lapply(se_list, colnames)


info(paste(
    "Number of samples in each experiment:\n", 
    paste(capture.output(lapply(sampleids_all, length)), collapse = "\n"),
    sep = ""))

sampleid <- unique(unlist(sampleids_all))
sample <- data.frame(sampleid, row.names = sampleid)
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

# PharmacoSet2(
#        name = "emptySet",
#        treatment = data.frame(),
#        sample = data.frame(),
#        molecularProfiles = MultiAssayExperiment(),
#        treatmentResponse = TreatmentResponseExperiment(),
#        perturbation = list(),
#        curation = list(sample = data.frame(), treatment = data.frame(), tissue = data.frame()),
#        datasetType = "sensitivity"
#      )
# save.image()

# make fake tre
filePath <- system.file('extdata', 'merckLongTable.csv', package='CoreGx',
  mustWork=TRUE)
merckDT <- data.table::fread(filePath, na.strings=c('NULL'))

# Our guesses of how we may identify rows, columns and assays
groups <- list(
  justDrugs=c('drug1id', 'drug2id'),
  drugsAndDoses=c('drug1id', 'drug2id', 'drug1dose', 'drug2dose'),
  justCells=c('cellid'),
  cellsAndBatches=c('cellid', 'batchid'),
  assays1=c('drug1id', 'drug2id', 'cellid'),
  assays2=c('drug1id', 'drug2id', 'drug1dose', 'drug2dose', 'cellid', 'batchid')
)

# Decide if we want to subset out mapped columns after each group
subsets <- c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)

# First we put our data in the `TRE`
TREdataMapper <- CoreGx::TREDataMapper(rawdata=merckDT)

# Then we can test our hypotheses, subset=FALSE means we don't remove mapped
#   columns after each group is mapped
guess <- CoreGx::guessMapping(TREdataMapper, groups=groups, subset=subsets)

CoreGx::rowDataMap(TREdataMapper) <- guess$drugsAndDose
CoreGx::colDataMap(TREdataMapper) <- guess$justCells

CoreGx::assayMap(TREdataMapper) <- list(
  sensitivity=list(
    guess$assays2[[1]],
    guess$assays2[[2]][seq_len(4)]
  ),
  profiles=list(
    guess$assays2[[1]],
    guess$assays2[[2]][c(5, 6)]
  )
)

tre <- CoreGx::metaConstruct(TREdataMapper)


pset <- PharmacoGx::PharmacoSet2(
    name = "GDSC",
    treatment = data.frame(),
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



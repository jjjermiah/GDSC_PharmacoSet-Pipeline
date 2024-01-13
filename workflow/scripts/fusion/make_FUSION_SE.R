## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}
library(data.table)
suppressPackageStartupMessages(library(SummarizedExperiment))
# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting make_FUSION_SE.R\n")

# 0.2 Read in the input files
# ---------------------------
input <- qs::qread(INPUT$preprocessed_fusions, nthreads = THREADS)

assays <- input$assays
rowData <- input$rowData
metadata <- input$metadata


info(sprintf(
    "\nTotal number of samples: %d\nTotal number of genes: %d", 
    ncol(assays), nrow(assays)))
info(paste(
    "Matrices:", 
    capture.output(str(assays)),
    collapse = "\n"))

info("Creating SummarizedExperiment object")
se <- list(
    "fusion" = SummarizedExperiment::SummarizedExperiment(
        assays = list(fusions = assays),
        colData = data.frame(
            sampleid = colnames(assays),
            batchid = rep(NA, ncol(assays))
        ),
        rowData = rowData,
        metadata = metadata))

info(paste("SummarizedExperiment:", capture.output(se),collapse = "\n"))
info(sprintf("Saving SummarizedExperiment object to %s", OUTPUT$fusion_se))
dir.create(dirname(OUTPUT$fusion_se), recursive = TRUE, showWarnings = FALSE)
qs::qsave(se, OUTPUT$fusion_se, nthreads = THREADS)

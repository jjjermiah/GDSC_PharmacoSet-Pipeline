#' RULE: make_MUTATION_SE
#' AUTHOR: Jermiah Joseph
#' DATE: 01-15-2024
#' This script takes in the following files:
#' - INPUT$preprocessed
#' and outputs the following files:
#' - OUTPUT$mutation_se
#' 
#' Libraries Used:
#' - data.table
#' - SummarizedExperiment
#' - log4r
#' - BiocParallel
#' - qs
#' 

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
info("Starting make_MUTATION_SE.R\n")

# 0.2 Read in the input files
# ---------------------------
info(paste("Loading", INPUT$preprocessed))
input <- qs::qread(INPUT$preprocessed, nthreads = THREADS)

assays <- input$assays
rowRanges <- input$GRanges
metadata <- input$metadata

assays <- lapply(assays, function(matrix){
  # only include rownames of matrix that are in rowRanges$symbol
  matrix[rownames(matrix) %in% rowRanges$symbol,]
})

info(paste(
    "Matrices:", 
    capture.output(str(assays)),
    collapse = "\n"))

sampleids <- unique(unlist(lapply(assays, colnames)))
geneids <- unique(unlist(lapply(assays, rownames)))

info(sprintf(
    "\nTotal number of samples: %d\nTotal number of genes: %d", 
    length(sampleids), length(geneids)))

rowRanges <- unique(rowRanges[rowRanges$symbol %in% geneids])


colData <- data.frame(
      sampleid = sampleids,
      batchid = rep(NA, length(sampleids))
)

rse <- list(
  "mutation" = SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowRanges = rowRanges,
    colData = colData,
    metadata = metadata))

info(paste("SummarizedExperiment:", capture.output(rse),collapse = "\n"))

dir.create(dirname(OUTPUT$mutation_se), recursive = TRUE, showWarnings = FALSE)
qs::qsave(rse, file = OUTPUT$mutation_se)

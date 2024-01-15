#' RULE: make_CNV_SE
#' AUTHOR: Jermiah Joseph
#' DATE: 01-15-2024
#' This script takes in the following files:
#' - INPUT$preprocessedCNV
#' and outputs the following files:
#' - OUTPUT$CNV_se
#' 
#' Libraries Used:
#' - data.table
#' - SummarizedExperiment
#' - log4r
#' - BiocParallel
#' - qs
#' - GenomicRanges

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
suppressPackageStartupMessages(library(GenomicRanges))
# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting make_CNV_SE.R\n")

# 0.2 Read in the input files
# ---------------------------
info(paste("Loading", INPUT$preprocessedCNV, " "))
preproc <- qs::qread(INPUT$preprocessedCNV, nthreads = THREADS)

assays_ <- preproc$assays
rowRanges_ <- preproc$GRanges
metadata_ <- preproc$metadata

info(paste(
    "Matrices:", 
    capture.output(str(assays_)),
    collapse = "\n"))


sampleids <- unique(unlist(lapply(assays_, colnames)))
geneids <- unique(unlist(lapply(assays_, rownames)))

info(sprintf(
    "\nTotal number of samples: %d\nTotal number of genes: %d", 
    length(sampleids), length(geneids)))

colData <- data.frame(
      sampleid = sampleids,
      batchid = rep(NA, length(sampleids))
)

rowRanges <- rowRanges_[rowRanges_$symbol %in% geneids,]

# rowRanges <- rowRanges[!duplicated(rowRanges$symbol),]

rse <- SummarizedExperiment::SummarizedExperiment(
    assays = assays_,
    rowRanges = rowRanges,
    colData = colData,
    metadata = metadata_
)
rse <- list(
    "cnv" = rse
)
info(paste("SummarizedExperiment:", capture.output(rse),collapse = "\n"))

dir.create(dirname(OUTPUT$CNV_se), recursive = TRUE, showWarnings = FALSE)
qs::qsave(rse, file = OUTPUT$CNV_se)

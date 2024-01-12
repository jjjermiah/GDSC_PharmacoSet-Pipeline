
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    # save.image()
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
info(paste("Loading", INPUT$metadata, " "))
metadata <- qs::qread(INPUT$metadata, nthreads = THREADS)
sample <- metadata$sample
geneAnnot <- metadata$GRanges

info(paste("Loading", INPUT$preprocessedCNV, " "))
matrices <- qs::qread(INPUT$preprocessedCNV, nthreads = THREADS)

# info(paste(
#     "Matrices:", 
#     capture.output(str(matrices)),
#     collapse = "\n"))


sampleids <- unique(unlist(lapply(matrices, colnames)))
geneids <- unique(unlist(lapply(matrices, rownames)))

info(sprintf(
    "\nTotal number of samples: %d\nTotal number of genes: %d", 
    length(sampleids), length(geneids)))

colData <- data.frame(
      sampleid = sampleids,
      batchid = rep(NA, length(sampleids))
)

rowRanges <- geneAnnot[geneAnnot$symbol %in% geneids,]

rse <- SummarizedExperiment::SummarizedExperiment(
    assays = matrices,
    rowRanges = rowRanges,
    colData = colData,
    metadata = list(
        data_source = snakemake@config$molecularProfiles$cnv
    )
)

info(paste(
    "SummarizedExperiment:", 
    capture.output(rse),
    collapse = "\n"
))

dir.create(dirname(OUTPUT$CNV_se), recursive = TRUE, showWarnings = FALSE)
qs::qsave(rse, file = OUTPUT$CNV_se)

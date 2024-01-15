## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}

library(data.table, quietly = TRUE)
suppressPackageStartupMessages(library(GenomicRanges))


# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting make_RNASEQ_SE.R\n")

# 0.2 read in data
# ----------------
info(paste0("Reading in: ", INPUT$preprocessed))
preproc <- qs::qread(INPUT$preprocessed, nthreads = THREADS)    

# 0.3 read in rnaseq data
# -----------------------
info(
  paste0("Parsing rnaseq data from: ", paste(names(preproc), collapse = ", ")))

assays_ <- preproc$assays
rowRanges_ <- preproc$GRanges
metadata_ <- preproc$metadata

# 1. Create SummarizedExperiment objects for each rnaseq assay
# ------------------------------------------------------------
assayNames_ <- names(assays_)

# TODO:: CHANGE colnames to sampleid instead of cell model passport ID
rse_list <- lapply(
  assayNames_, 
  function(assayName_){
    
    assay <- assays_[[assayName_]]
    # remove duplicated rownames 
    assay <- assay[!duplicated(rownames(assay)),]

    sampleids <- colnames(assay)
    geneids <- rownames(assay)
    info(paste("Creating colData for", assayName_))
    colData <- data.frame(
      sampleid = sampleids,
      batchid = rep(NA, length(sampleids))
    )

    info(paste("Creating rowRanges for", assayName_))
    rowRanges <- rowRanges_[rowRanges_$symbol %in% geneids,][!duplicated(rowRanges_$symbol),]

    assays <- list(assay)
    names(assays) <- paste0("rnaseq.", assayName_)

    info(paste0("Creating SummarizedExperiment for ", assayName_))
    rse <- SummarizedExperiment::SummarizedExperiment(
      assays = assays,
      rowRanges = rowRanges,
      colData = colData,
      metadata = metadata_[[assayName_]]
    )
    rse_show <- capture.output(rse)
    info(paste0(" SummarizedExperiment for: ", assayName_))
    info(paste0("\n\t", rse_show))
    return(rse)
  }
)
names(rse_list) <- paste0("rnaseq.", assayNames_)

# 2. Save output files
# --------------------
info(paste0("Saving output files to ", OUTPUT$rnaseq_se))
dir.create(dirname(OUTPUT$rnaseq_se), recursive=T, showWarnings = FALSE)
qs::qsave(rse_list, file = OUTPUT$rnaseq_se, nthreads = THREADS)

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
library(GenomicRanges, quietly = TRUE)
# set up logging
logger <- log4r::logger(
    appenders = list(
        log4r::file_appender(LOGFILE, append = TRUE)
    )
)

# 0.1 read in data
# ----------------
log4r::info(logger, paste0("Reading in: ", INPUT$metadata))
metadata <- qs::qread(INPUT$metadata)

log4r::info(logger, paste0("Reading in: ", INPUT$preprocessed))
preproc <- qs::qread(INPUT$preprocessed, nthreads = THREADS)    

# 0.2 read in rnaseq data
# -----------------------
log4r::info(logger, 
  paste0("Parsing rnaseq data from: ", paste(names(preproc), collapse = ", ")))

assays_ <- preproc$assays
rawData <- preproc$rawData
rowRanges_ <- preproc$GRanges

# 1. Create SummarizedExperiment objects for each rnaseq assay
# ------------------------------------------------------------
assayNames_ <- names(assays_)

# TODO:: CHANGE colnames to sampleid instead of cell model passport ID
# TODO:: CHANGE rownames to ensembl gene id instead of gene symbol
rse_list <- BiocParallel::bplapply(
  assayNames_, 
  function(assayName_){
    log4r::info(logger, paste0("Creating SummarizedExperiment for ", assayName_))
    assay <- assays_[[assayName_]]
    sampleids <- unique(colnames(assay))

    colData <- data.frame(
      sampleid = sampleids,
      batchid = rep(NA, length(sampleids))
    )

    rowRanges <- rowRanges_[rowRanges_$CMP_gene_id %in% rownames(assay),]
    assays <- list(assay)
    names(assays) <- paste0("rnaseq.", assayName_)

    metadata <- list(
      data_source = snakemake@config$molecularProfiles$rnaseq$processed,
      filename = unique(preproc$rawData[[assayName_]][["file"]])
    )

    rse <- SummarizedExperiment::SummarizedExperiment(
      assays = assays,
      rowRanges = rowRanges,
      colData = colData,
      metadata = metadata
    )
    rse_show <- capture.output(rse)
    log4r::info(logger, paste0(" SummarizedExperiment for: ", names(assays)))
    log4r::info(logger, paste0("\n\t", rse_show))
  },
  BPPARAM = BiocParallel::MulticoreParam(workers = THREADS)
)
names(rse_list) <- paste0("rnaseq.", assayNames_)

# 2. Save output files
# --------------------
log4r::info(logger, paste0("Saving output files to ", OUTPUT$rnaseq_se))
dir.create(dirname(OUTPUT$rnaseq_se), recursive=T, showWarnings = FALSE)
qs::qsave(rse_list, file = OUTPUT$rnaseq_se, nthreads = THREADS)

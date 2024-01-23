# AUTHOR: Jermiah Joseph
# CREATED: 01-08-2024
# This script takes in the following files:
#  - INPUT$metadata
#  - INPUT$WES_zipped
# and outputs the following files:
#  - OUTPUT$preprocessedCNV
# 
# PACKAGE DEPENDENCIES:
#  - data.table
#  - GenomicRanges
#  - log4r
#  - BiocParallel
#  - qs
#
# NOTES:
# 1. Categorisation of Total Copy Number values
# source(https://depmap.sanger.ac.uk/documentation/datasets/copy-number/)
    # The total copy number values have been categorised (CNA Call) 
    # using the following calculation:

    # Val = round( 2 * 2^log2(C/Ploidy) )

    # if Val == 0: Category = 'Deletion'
    # if Val == 1: Category = 'Loss'
    # if Val == 2: Category = 'Neutral'
    # if Val == 3: Category = 'Gain'
    # if Val >= 4: Category = 'Amplification'

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
    appenders = list(
      log4r::file_appender(LOGFILE, append = TRUE),
      log4r::console_appender()
    )
)

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting preprocessCNV.R\n")
    
# 0.2 Read in the input files
# ---------------------------
metadata <- qs::qread(INPUT$metadata, nthreads = THREADS)
sampleMetadata <- metadata$sample
geneAnnot <- metadata$GRanges


WES_CNV_dir <- paste0(dirname(INPUT$WES_zipped), "/WES_CNV")
dir.create(WES_CNV_dir, recursive = TRUE, showWarnings = FALSE)
info(paste0("Unzipping ", INPUT$WES_zipped, " into ", WES_CNV_dir))
unzip(INPUT$WES_zipped, exdir = WES_CNV_dir)

# list.files(WES_CNV_dir)


inputFilesNames <- file.path(WES_CNV_dir, list.files(WES_CNV_dir))

# 0.3 Read in the raw CNV data
# -----------------------
# Read in large file "WES_pureCN_CNV_genes_20221213.csv" :
# script was written to handle both WES and WGS, only to find out that WGS 
# does not contain data for any samples in GDSC sample metadata
# rawdata/cnv/WES_pureCN_CNV_genes_20221213.csv

gene_files <- inputFilesNames[grepl("genes", inputFilesNames)]

# Get the largest file in the list of gene_files
gene_files <- gene_files[which.max(file.size(gene_files))]

genes_dt <- lapply(gene_files, function(file){
    info(paste("Loading", file," "))
    df <- data.table::fread(
        file, header = TRUE, showProgress = FALSE,
        sep = ",", stringsAsFactors = FALSE, nThread = THREADS)
    
    info(sprintf("Loaded %s with %d rows and %d columns", file, nrow(df), ncol(df)))

    info(sprintf("Subsetting %s to only samples in GDSC sample metadata", file))
    df <- df[model_id %in% sampleMetadata[, model_id]]

    info(sprintf("Subsetting %s to only genes in GDSC gene annotation", file))
    # drop the symbol column 
    df <- merge(
        df[gene_id %in% geneAnnot$CMP_gene_id, !c("symbol"), with = FALSE], 
        data.table::as.data.table(geneAnnot)[,.(symbol, CMP_gene_id)], 
        by.x = "gene_id", by.y = "CMP_gene_id", all.x = TRUE)

    info(sprintf("data.table now has %d rows and %d columns", nrow(df), ncol(df)))
    info(sprintf("Total number of model_ids: %d", uniqueN(df[, model_id])))
    info(sprintf("Total number of gene_ids: %d", uniqueN(df[, gene_id])))

    return(df)
})
names(genes_dt) <- basename(gene_files)


# 1. Build data structures for each datatype:
# -------------------------------------------
wes_gene_dt <- genes_dt[[1]]

setnames(wes_gene_dt, "model_name", "sampleid")

# ------------------------------
info(paste("\n",capture.output(wes_gene_dt[, .N, by = source])))
#    source        N
# 1: Sanger 16884273
# 2:  Broad  1843646

cols <- c(
    "sampleid", "symbol",
    "total_copy_number", "cn_category", 
    "seg_mean", "gene_mean", "num_snps", "gatk_mean_log2_copy_ratio", "source")
dt <- wes_gene_dt[, ..cols, with = FALSE]

# subset assay_dt to only include rows with source == "Broad" 
source_ <- "Sanger"
assay_dt <- dt[source == source_, !c("source"), with = FALSE]

# for each column that isnt sampleid or symbol, create a matrix with genes rows
# and samples as columns and set the rownames to the gene_id
# the `dcast` function to reshape the data from 
# long to wide format, with genes as rows and samples as columns.
# The resulting data frame `assay_dt.t` contains the `col` values 
# (i.e total_copy_number, cn_category, etc) for each gene-sample combination.
# If the `col` values are of type character, the function `first` 
# is used to aggregate the values, otherwise the mean is calculated.
assayNames <- cols[!cols %in% c("sampleid", "symbol", "source")]
matrices <- BiocParallel::bplapply(
    assayNames,
    function(col){
        
        info(paste("Casting ", col))
        assay_dt.t <- dcast(
            assay_dt[, c("sampleid", "symbol", col), with = FALSE],
            symbol ~ sampleid,
            value.var = col,
            fun.aggregate = if(is.character(assay_dt[[col]])) dplyr::first else mean
        )

        info(paste("Converting ", col, " to matrix"))
        mtx <- as.matrix(
            assay_dt.t[, !c("symbol"), with = FALSE],
            rownames = assay_dt.t[["symbol"]])
        
        info(sprintf(
            "Matrix %s has %d rows and %d columns", col, nrow(mtx), ncol(mtx)))
        return(mtx)
    },
    BPPARAM = BiocParallel::MulticoreParam(
        workers = THREADS, timeout=3600)
)
names(matrices) <- paste0(assayNames, ".sanger")

# Each matrix should have the same number of rows and columns
# check that the number of rows and columns are the same for each matrix
info("Checking that the number of rows and columns are the same for each matrix")
if(all.equal(lengths(lapply(matrices, nrow)),lengths(lapply(matrices, ncol)))){
    info("All matrices have the same number of rows and columns")
} else {
    log4r::error(logger, "Not all matrices have the same number of rows and columns")
    stop("preprocessCNV.R: Not all matrices have the same number of rows and columns")
}

metadata <- list(
    data_source = snakemake@config$molecularProfiles$cnv,
    filename = basename(gene_files))

# 3. Save Output
# --------------
info("Saving Output Files")
outputFiles <- list(
    "assays" = matrices,
    "GRanges" = geneAnnot,
    "metadata" = metadata)

# make output directory if it doesnt exist
dir.create(dirname(OUTPUT$preprocessedCNV), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputFiles, file = OUTPUT$preprocessedCNV, nthreads = THREADS)

#' @title Preprocess FUSION data
#' RULE: preprocess_FUSION
#' AUTHOR: Jermiah Joseph
#' DATE: 01-15-2024
#' This script takes in the following files:
#' - INPUT$metadata
#' - INPUT$gene_fusions
#' and outputs the following files:
#' - OUTPUT$preprocessed_fusions
#' 
#' Libraries Used:
#' - data.table
#' - GenomicRanges
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
suppressPackageStartupMessages(library(GenomicRanges))
# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting make_CNV_SE.R\n")


# 0.1 read metadata
# -----------------
info(paste("Loading", INPUT$metadata, " "))
metadata <- qs::qread(INPUT$metadata, nthreads = THREADS)

sample <- metadata$sample
info(paste("Sample metadata:", capture.output(str(sample)), collapse = "\n"))

geneAnnot <- metadata$GRanges
info(paste("Gene annotation:", capture.output(str(geneAnnot)), collapse = "\n"))

# 0.2 read gene fusions
# ---------------------
zipDir <- dirname(INPUT$gene_fusions)
unzipDir <- file.path(dirname(zipDir), "gene_fusions")
unzip(INPUT$gene_fusions, exdir = unzipDir)

info(paste("Loading", list.files(unzipDir), " "))
dt <- data.table::fread(file.path(unzipDir, list.files(unzipDir)))

info(paste("Gene fusions:", capture.output(str(dt)), collapse = "\n"))
cols <- c(
    "model_id", "model_name", "tissue",
    "chr_3prime", "chr_5prime", 
    "gene_id_3prime", "gene_id_5prime", 
    "gene_symbol_3prime", "gene_symbol_5prime")

# 1.0 subset gene fusions to only genes in GDSC sample metadata
# ----------------------------------------------------------------
fusion_dt <- dt[
    dt$model_id %in% sample[, model_id], ..cols]
info(paste("Before Subsetting gene fusions:", capture.output(str(fusion_dt)), collapse = "\n"))
fusion_dt <- fusion_dt[gene_id_3prime %in% geneAnnot$CMP_gene_id]
fusion_dt <- fusion_dt[gene_id_5prime %in% geneAnnot$CMP_gene_id]

info(paste("Subsetting gene fusions:", capture.output(str(fusion_dt)), collapse = "\n"))

fusion_dt$status <- "fusion"

dt_t <- data.table::dcast(
    fusion_dt[, .(gene_symbol_5prime, model_name, gene_symbol_3prime, status)], 
    gene_symbol_3prime + gene_symbol_5prime ~ model_name, 
    value.var = "status", 
    fun.aggregate = dplyr::first,
    fill = "wt",
    sep = "_"
)

# make a matrix from dt_t with a combination of 
# gene_symbol_3prime and gene_symbol_5prime as rownames and
# model_id as colnames
mtx <- as.matrix(
    dt_t[, -c("gene_symbol_3prime", "gene_symbol_5prime")],
    rownames = paste0(dt_t$gene_symbol_3prime, "_", dt_t$gene_symbol_5prime))
info(paste("Matrix:", capture.output(str(mtx)), collapse = "\n"))

metadata <- list(
    data_source = snakemake@config$molecularProfiles$fusion,
    filename = INPUT$gene_fusions
)


dt_t[, rownames := paste0(gene_symbol_3prime, "_", gene_symbol_5prime)]
rowData <- unique(merge(
    dt_t[, c("gene_symbol_3prime", "gene_symbol_5prime", "rownames")],
    fusion_dt[, 
        c("chr_3prime", "chr_5prime", 
        "gene_id_3prime", "gene_id_5prime", 
        "gene_symbol_3prime", "gene_symbol_5prime")],
    by = c("gene_symbol_3prime", "gene_symbol_5prime")
))

rowData <- data.frame(
    rowData,
    stringsAsFactors = FALSE,
    row.names = rowData$rownames)

outputFiles <- list(
    "assays" = mtx,
    "rowData" = rowData,
    "metadata" = metadata)


dir.create(dirname(OUTPUT$preprocessed_fusions), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputFiles, file = OUTPUT$preprocessed_fusions, nthreads = THREADS)

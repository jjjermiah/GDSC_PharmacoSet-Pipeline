#' RULE: preprocess_MUTATION
#' AUTHOR: Jermiah Joseph
#' DATE: 01-15-2024
#' This script takes in the following files:
#'  - INPUT$metadata
#' - INPUT$all_mutations
#' - INPUT$Genes_Metadata
#' and outputs the following files:
#' - OUTPUT$preprocessed
#' 
#' Libraries Used:
#' - data.table
#' - GenomicRanges
#' - log4r
#' - BiocParallel
#' - qs

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

# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting preprocess_MUTATION.R\n")

# 0.2 read in metadata
# -------------------- 
info(paste("Loading", INPUT$metadata))
metadata <- qs::qread(INPUT$metadata)
sample <- metadata$sample
geneAnnot <- metadata$GRanges

# 0.3 read in mutation data
# -------------------------
dir <- paste0(dirname(INPUT$all_mutations), "/all")
allDir <- paste0(dirname(INPUT$all_mutations), "/all")
dir.create(allDir, recursive = TRUE, showWarnings = FALSE)
unzip(INPUT$all_mutations, exdir = allDir)

# 0.4 read mutation gene metadata
# -------------------------------
info(paste("Loading", file.path(allDir, list.files(allDir))))
mut_data <- data.table::fread(
    file.path(allDir, list.files(allDir)), 
    header = TRUE, sep = ",", showProgress = F)

# 0.5 read mutation gene metadata
# -------------------------------
info(paste("Loading", INPUT$Genes_Metadata))
mutGenesAnnot <- data.table::fread(INPUT$Genes_Metadata, header = TRUE, sep = ",")

# 1.0 Subset mutation data to only include samples from GDSC metadata
# -------------------------------------------------------------------
info("Subsetting mutation data to only include samples from GDSC metadata")
data.table::setkey(mut_data, model_id)
mut_dt <- unique(merge(
    sample[, .(sampleid, model_id)], 
    mut_data, 
    by = "model_id"))

# 2.0 Create GRanges object from mutGenesAnnot
# --------------------------------------------
# GRanges_dt <- data.table::as.data.table(geneAnnot)
# GRanges_dt <- GRanges_dt[, 
#     .(CMP_gene_id, ensembl_gene_id, entrez_id, strand, seqnames, 
#     hgnc_id, refseq_id, uniprot_id, gene_type, source)]

# mutGenesAnnot <- merge(
#     mutGenesAnnot[, !c("strand")], 
#     GRanges_dt, 
#     by.x = "gene_id", by.y = "CMP_gene_id")

# mutGenesAnnot <- mutGenesAnnot[gene_id %in% mut_dt$gene_id]

# gr <- GenomicRanges::makeGRangesFromDataFrame(
#     df = mutGenesAnnot, keep.extra.columns=TRUE, na.rm=TRUE,
#     start.field = "chr_start", end.field = "chr_end", seqnames.field = "seqnames",)

# 3.0 Create Assays object from mut_dt
# -------------------------------------------------
# isolate only columns to use for assay
assay_cols <- c("protein_mutation", "rna_mutation",
    "cdna_mutation", "vaf", "effect")
cols <- c("sampleid", "gene_symbol", assay_cols)

# use only sanger data
assay_dt <- mut_dt[source == "Sanger", ..cols]

matrices <- BiocParallel::bplapply(assay_cols, function(x){
    info(sprintf("Subsetting df for %s", x))
    cols_ <- c("sampleid", "gene_symbol", x)
    assay_mtx <- assay_dt[, ..cols_]

    # replace "-" in the x col with "wt"
    assay_mtx[, (x) := lapply(.SD, function(x) gsub("-", "wt", x)), .SDcols = x]

    info(sprintf("casting df for %s", x))
    # dcast assay_mtx to wide format
    assay_mtx <- data.table::dcast(
        unique(assay_mtx), gene_symbol ~ sampleid, 
        value.var = x, fun.aggregate = dplyr::first, fill = "wt")

    info(sprintf("Converting %s into matrix", x))
    # create a matrix with gene_symbol as rownames and sampleid as colnames
    mtx <- as.matrix(assay_mtx[, -c("gene_symbol"), with = FALSE], 
        rownames = assay_mtx[["gene_symbol"]])
    
    info(sprintf(
            "Matrix %s has %d rows and %d columns", x, nrow(mtx), ncol(mtx)))
    mtx
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS))
names(matrices) <- assay_cols

# 4.0 Setup metadata for SummarizedExperiment object
# --------------------------------------------------
metadata <- list(
    data_source = snakemake@config$molecularProfiles$mutation$SUMMARY,
    filename.data = INPUT$all_mutations,
    filename.gene_metadata = INPUT$Genes_Metadata)

# 5.0 Save Output
# ---------------
info("Saving Output Files")
outputFiles <- list(
    "assays" = matrices,
    "GRanges" = geneAnnot,
    "metadata" = metadata)

dir.create(dirname(OUTPUT$preprocessed), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputFiles, file = OUTPUT$preprocessed, nthreads = THREADS)

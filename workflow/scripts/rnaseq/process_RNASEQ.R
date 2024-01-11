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

# set up logging
logger <- log4r::logger(
    appenders = list(
        log4r::file_appender(LOGFILE, append = TRUE)
    )
)

# 0.1 read in metadata
# --------------------
metadata <- qs::qread(INPUT$metadata)
sample <- metadata$sample

genes <- metadata$geneAnnot

# 0.2 read in rnaseq data
# -----------------------
allDir <- paste0(dirname(INPUT$all), "/all")
dir.create(allDir, recursive = TRUE, showWarnings = FALSE)
unzip(INPUT$all, exdir = allDir)

# list.files(allDir)
# [1] "rnaseq_all_data_20220624.csv"   "rnaseq_fpkm_20220624.csv"      
# [3] "rnaseq_read_count_20220624.csv" "rnaseq_tpm_20220624.csv" 

# file names without date or extension
parsedNames <- lapply(
    list.files(allDir), 
    function(x) gsub("^rnaseq_|_\\d{8}\\.csv", "", x)
)

# for each csv in allDir, read in the csv into a data.table and add a source col
rnaseq_data <- lapply(list.files(allDir), function(file){
    # message("Reading in ", file)
    log4r::info(logger, paste0("Reading in ", file))
    data <- data.table::fread(file.path(allDir, file), header = TRUE, sep = ",")
    data$file <- file
    data
})
names(rnaseq_data) <- parsedNames


# . Get Gene annotation from Gencode
# ----------------------------------
log4r::info(logger, "Getting gene annotation from Gencode")
path <- "/home/bioinf/bhklab/jermiah/psets/PharmacoSet-Pipelines/GDSC/metadata/human/GRCh38_v44/annotation.gtf"
dsGencode <- rtracklayer::import(path)

dsGencode_dt <- data.table::as.data.table(dsGencode)
dsGencode_dt$gene_id <- gsub("\\.\\d+$", "", dsGencode_dt$gene_id)

genes_dt <- merge(
    genes,
    dsGencode_dt[type == "gene"],
    by.x = "ensembl_gene_id",
    by.y = "gene_id",
    all.x=TRUE,      #
    sort=FALSE
)
genes_dt$CMP_gene_id <- genes_dt$gene_id
setnames(genes_dt, old = "ensembl_gene_id", new = "gene_id")

gene_rRanges <- GenomicRanges::makeGRangesFromDataFrame(
    genes_dt, keep.extra.columns=TRUE, na.rm = T)

# 1.0 Subset rnaseq data to only include samples from GDSC metadata
# -----------------------------------------------------------------
log4r::info(
    logger,"Subsetting rnaseq data to only include samples from GDSC metadata")

dataset_types <- c("fpkm", "read_count", "tpm")
assays <- BiocParallel::bplapply(
    dataset_types, 
    function(x){ 
        datatype_dt <- 
            rnaseq_data[[x]][
                -(1:4),
                colnames(rnaseq_data[[x]]) %in% c("model_id", sample$model_id), 
                with = FALSE]
        setnames(datatype_dt,old = "model_id", new = "gene_id")
        datatype_dt <- datatype_dt[gene_id %in% gene_rRanges$CMP_gene_id,]

        # convert to matrix and set rownames to gene_id and remove gene_id column
        as.matrix(
            datatype_dt[, !c("gene_id"), with = FALSE],
            rownames = datatype_dt[["gene_id"]]
        )
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = 3))
names(assays) <- dataset_types

log4r::info(logger, "Saving output files")
outputFiles <- list(
    "assays" = assays,
    "rawData" = rnaseq_data,
    "GRanges" = gene_rRanges
)

dir.create(dirname(OUTPUT$preprocessed), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputFiles, file = OUTPUT$preprocessed, nthreads = THREADS)

# OLD ::: KEEPING FOR NOW


# # 2.0 Cast rnaseq data into matrices
# # ----------------------------------

# # The following function takes a data.table and a datatype and returns a matrix
# # with genes as rows and samples as columns
# # in the matrix, each cell is the mean of the datatype for the gene and sample
# dcast_datatypes <- function(dt, datatype){
#     log4r::info(logger, paste0("Casting ", datatype, " into matrix"))
#     columns_to_subset <- c("sampleid", "ensembl_gene_id", datatype)
#     dt <- dt[, ..columns_to_subset]
#     matrix <- data.table::dcast(
#         data = dt, 
#         formula = ensembl_gene_id ~ sampleid, 
#         value.var = datatype,
#         fun.aggregate = mean)
#     mtx <- as.matrix(
#         matrix,
#         rownames = matrix[["ensembl_gene_id"]]
#     )
#     # remove ensembl_gene_id column and any unnamed rows
#     mtx <- mtx[,-1]
#     mtx <- mtx[!grepl("^\\s*$", rownames(mtx)),]

#     # Assuming 'mtx' is your matrix
#     mtx_numeric <- apply(mtx, 2, as.numeric)

#     # Convert the matrix back to a data.frame with row names
#     mtx_numeric_df <- as.data.frame(mtx_numeric)
#     rownames(mtx_numeric_df) <- rownames(mtx)
#     as.matrix(mtx_numeric_df)
# }

# datatypes <- c("fpkm", "read_count", "tpm")
# # run dcast_datatypes in parallel for each datatype
# log4r::info(logger, "Casting rnaseq data into matrices")
# matrices <- BiocParallel::bplapply(
#     X = datatypes, 
#     FUN = function(x) {
#         dcast_datatypes(rnaseq, x)},
#     BPPARAM = BiocParallel::MulticoreParam(workers = 3)
# )
# names(matrices) <- paste0("rnaseq.", datatypes)

# tmp <- lapply(names(matrices), function(x) {
#     dims <- dim(matrices[[x]])
#     log4r::info(
#         logger, 
#         paste0("Matrix ", x, " has ", dims[1], " rows and ", dims[2], " columns")
#     )
# })


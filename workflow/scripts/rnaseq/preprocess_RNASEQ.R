# RULE: preprocess_RNASEQ
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
info("Starting preprocess_RNASEQ.R\n")

# 0.2 read in metadata
# --------------------
metadata <- qs::qread(INPUT$metadata)
sample <- metadata$sample
geneAnnot <- metadata$GRanges

# 0.2 read in rnaseq data
# -----------------------
allDir <- paste0(dirname(INPUT$all), "/all")
dir.create(allDir, recursive = TRUE, showWarnings = FALSE)
info(paste0("Unzipping ", INPUT$all, " into ", allDir))
unzip(INPUT$all, exdir = allDir)

# list.files(allDir)
# [1] "rnaseq_all_data_20220624.csv"   "rnaseq_fpkm_20220624.csv"      
# [3] "rnaseq_read_count_20220624.csv" "rnaseq_tpm_20220624.csv" 

# file names without date or extension
parsedNames <- lapply(list.files(allDir), 
    function(x) gsub("^rnaseq_|_\\d{8}\\.csv", "", x)
)
# for each csv in allDir, read in the csv into a data.table and add source col
rnaseq_data <- lapply(file.path(allDir, list.files(allDir)), function(file){
    info(paste0("Reading in ", file))
    data <- data.table::fread(file, header = TRUE, sep = ",", showProgress = F)
    data$file <- file
    data
})
names(rnaseq_data) <- parsedNames

# 1.0 Subset rnaseq data to only include samples from GDSC metadata
# -----------------------------------------------------------------

info("Subsetting rnaseq data to only include samples from GDSC metadata")

# NOTE: we are not going to use the "all" data due to the resources required
dataset_types <- c("fpkm", "read_count", "tpm")
assays <- BiocParallel::bplapply(
    dataset_types, 
    function(x){ 
        info(sprintf("Subsetting %s", x))

        datatype_dt <- 
            rnaseq_data[[x]][
                -(1:4),
                colnames(rnaseq_data[[x]]) %in% c("model_id", sample$model_id), 
                with = FALSE]

        sample_subset <- 
            sample[model_id %in% colnames(datatype_dt), .(sampleid, model_id)]
        
        data.table::setnames(
            x = datatype_dt,
            old = c("model_id", sample_subset$model_id), 
            new = c("gene_id", sample_subset$sampleid))

        # set column order alphabetically
        ordered_cols <- c("gene_id", sort(colnames(datatype_dt)[-1]))
        datatype_dt <- datatype_dt[, ..ordered_cols]

        info(sprintf(
            "Subsetting %s for only genes in GDSC gene annotation", x))
        datatype_dt <- merge(
            datatype_dt[gene_id %in% geneAnnot$CMP_gene_id,],
            data.table::as.data.table(geneAnnot)[,.(symbol, CMP_gene_id)],
            by.x = "gene_id", by.y = "CMP_gene_id", all.x = TRUE)

        info(sprintf("Converting %s to matrix", x))
        mtx <- as.matrix(
            datatype_dt[, !c("gene_id", "symbol"), with = FALSE],
            rownames = datatype_dt[["symbol"]]
        )
        info(sprintf(
            "Matrix %s has %d rows and %d columns", x, nrow(mtx), ncol(mtx)))
        return(mtx)
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS))
names(assays) <- dataset_types

metadata <- lapply(dataset_types, function(x) {
    list(
        data_source = snakemake@config$molecularProfiles$rnaseq$processed,
        filename = unique(rnaseq_data[[x]][["file"]]))})
names(metadata) <- dataset_types

# 3. Save Output
# --------------
info("Saving Output Files")
outputFiles <- list(
    "assays" = assays,
    "GRanges" = geneAnnot,
    "metadata" = metadata)

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



## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image()
}


library(data.table)
suppressPackageStartupMessages(library(GenomicRanges))


# 0.1 read metadata
# -----------------
metadata <- qs::qread(INPUT$metadata, nthreads = THREADS)
sample <- metadata$sample
geneAnnot <- metadata$GRanges
# 0.2 read gene fusions
# ---------------------
zipDir <- dirname(INPUT$gene_fusions)
unzipDir <- file.path(dirname(zipDir), "gene_fusions")
unzip(INPUT$gene_fusions, exdir = unzipDir)
dt <- data.table::fread(file.path(unzipDir, list.files(unzipDir)))

cols <- c(
    "model_id", "model_name", "tissue",
    "chr_3prime", "chr_5prime", 
    "gene_id_3prime", "gene_id_5prime", 
    "gene_symbol_3prime", "gene_symbol_5prime")


fusion_dt <- dt[
    dt$model_id %in% sample[, model_id], ..cols]

fusion_dt <- fusion_dt[gene_id_3prime %in% geneAnnot$CMP_gene_id]
fusion_dt <- fusion_dt[gene_id_5prime %in% geneAnnot$CMP_gene_id]

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

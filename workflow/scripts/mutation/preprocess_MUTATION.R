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
suppressPackageStartupMessages(library(SummarizedExperiment))
# 0.1 read in metadata
# -------------------- 
metadata <- qs::qread(INPUT$metadata)
sample <- metadata$sample
geneAnnot <- metadata$GRanges

# 0.2 read in mutation data
# -------------------------
dir <- paste0(dirname(INPUT$all_mutations), "/all")
allDir <- paste0(dirname(INPUT$all_mutations), "/all")
dir.create(allDir, recursive = TRUE, showWarnings = FALSE)
unzip(INPUT$all_mutations, exdir = allDir)


# 0.3 read mutation gene metadata
# -------------------------------
mut_data <- data.table::fread(
    file.path(allDir, list.files(allDir)), 
    header = TRUE, sep = ",", showProgress = F)

# 0.4 read mutation gene metadata
# -------------------------------
mutGenesAnnot <- data.table::fread(INPUT$Genes_Metadata, header = TRUE, sep = ",")

# 1.0 Subset mutation data to only include samples from GDSC metadata
# -------------------------------------------------------------------
data.table::setkey(mut_dt, model_id)
mut_dt <- unique(merge(
    sample[, .(sampleid, model_id)], 
    mut_data, 
    by = "model_id"))

# 2.0 Create GRanges object from mutGenesAnnot
# --------------------------------------------
GRanges_dt <- data.table::as.data.table(geneAnnot)
GRanges_dt <- GRanges_dt[, 
    .(CMP_gene_id, ensembl_gene_id, entrez_id, strand, seqnames, 
    hgnc_id, refseq_id, uniprot_id, gene_type, source)]

mutGenesAnnot <- merge(
    mutGenesAnnot[, !c("strand")], 
    GRanges_dt, 
    by.x = "gene_id", by.y = "CMP_gene_id")

mutGenesAnnot <- mutGenesAnnot[gene_id %in% mut_dt$gene_id]

gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = mutGenesAnnot, keep.extra.columns=TRUE, na.rm=TRUE,
    start.field = "chr_start", end.field = "chr_end", seqnames.field = "seqnames",)

# 3.0 Create SummarizedExperiment object from mut_dt
# -------------------------------------------------
# isolate only columns to use for assay
assay_cols <- c("protein_mutation", "rna_mutation",
    "cdna_mutation", "vaf", "effect")
cols <- c("sampleid", "gene_symbol", assay_cols)

# use only sanger data
assay_dt <- mut_dt[source == "Sanger", ..cols]

matrices <- BiocParallel::bplapply(assay_cols, function(x){
    cols_ <- c("sampleid", "gene_symbol", x)
    assay_mtx <- assay_dt[, ..cols_]

    # replace "-" in the x col with "wt"
    assay_mtx[, (x) := lapply(.SD, function(x) gsub("-", "wt", x)), .SDcols = x]

    # dcast assay_mtx to wide format
    assay_mtx <- data.table::dcast(
        unique(assay_mtx), gene_symbol ~ sampleid, 
        value.var = x, fun.aggregate = dplyr::first, fill = "wt")

    # create a matrix with gene_symbol as rownames and sampleid as colnames
    mtx <- as.matrix(assay_mtx[, -c("gene_symbol"), with = FALSE], 
        rownames = assay_mtx[["gene_symbol"]])
    
    mtx
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS))



# # . Get Gene annotation from Gencode
# path <- "/home/bioinf/bhklab/jermiah/psets/PharmacoSet-Pipelines/GDSC/metadata/human/GRCh38_v44/annotation.gtf"
# dsGencode <- rtracklayer::import(path)

# # remove gene_id version from dsGencode 
# dsGencode$gene_id <- gsub("\\.\\d+$", "", dsGencode$gene_id)

# geneAnnot <- geneAnnot[ensembl_gene_id %in% dsGencode$gene_id]
# # genes <- merge(geneAnnot, dsGencode, by.x = "ensembl_gene_id", by.y = "gene_id")

# mut <- merge(mut_dt, geneAnnot, by.x = "gene_id", by.y = "gene_id")



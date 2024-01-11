
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image()
}

library(data.table)

# 0.1 read in metadata
# -------------------- 
metadata <- qs::qread(INPUT$metadata)
sample <- metadata$sample
geneAnnot <- metadata$geneAnnot

# 0.2 read in mutation data
# -------------------------
dir <- paste0(dirname(INPUT$all_mutations), "/all")
allDir <- paste0(dirname(INPUT$all_mutations), "/all")
dir.create(allDir, recursive = TRUE, showWarnings = FALSE)
unzip(INPUT$all_mutations, exdir = allDir)
mut_dt <- data.table::fread(file.path(allDir, list.files(allDir)), header = TRUE, sep = ",")
mut_dt <- mut_dt[model_id %in% sample$model_id]
mut_dt <- merge(sample[, .(sampleid, model_id)], mut_dt, by.x = "model_id", by.y = "model_id")

# 0.3 read mutation gene metadata
# -------------------------------
mutGenesAnnot <- data.table::fread(INPUT$Genes_Metadata, header = TRUE, sep = ",")
mutGenesAnnot <- merge(mutGenesAnnot, geneAnnot, by.x = "gene_id", by.y = "gene_id")
mutGenesAnnot <- mutGenesAnnot[gene_id %in% mut_dt$gene_id]

mutGenes_metadata <- mutGenesAnnot[, -c("symbol", "chr_start", "chr_end", "strand")]
gr <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(mutGenesAnnot$symbol),
  ranges = IRanges::IRanges(mutGenesAnnot$chr_start, mutGenesAnnot$chr_end),
  strand = S4Vectors::Rle(mutGenesAnnot$strand),
  metadata = mutGenes_metadata
)

# TODO::Need to create matrices of mutations to be the assay for the SummarizedExperiment




# # . Get Gene annotation from Gencode
# path <- "/home/bioinf/bhklab/jermiah/psets/PharmacoSet-Pipelines/GDSC/metadata/human/GRCh38_v44/annotation.gtf"
# dsGencode <- rtracklayer::import(path)

# # remove gene_id version from dsGencode 
# dsGencode$gene_id <- gsub("\\.\\d+$", "", dsGencode$gene_id)

# geneAnnot <- geneAnnot[ensembl_gene_id %in% dsGencode$gene_id]
# # genes <- merge(geneAnnot, dsGencode, by.x = "ensembl_gene_id", by.y = "gene_id")

# mut <- merge(mut_dt, geneAnnot, by.x = "gene_id", by.y = "gene_id")


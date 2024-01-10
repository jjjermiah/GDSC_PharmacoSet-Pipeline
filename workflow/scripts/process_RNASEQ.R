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
metadata <- qs::qread(INPUT$metadata)
sample <- metadata$sample

genes <- metadata$geneAnnot

################################################################################
# rnaseq data
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
    message("Reading in ", file)
    data <- data.table::fread(file.path(allDir, file), header = TRUE, sep = ",")
    data$file <- file
    data
})
names(rnaseq_data) <- parsedNames

################################################################################
## Subset data
rnaseq <- rnaseq_data$all_data[model_id %in% sample$model_id]
rnaseq <- merge(sample[, .(sampleid, model_id)], rnaseq, by.x = "model_id", by.y = "model_id")
rnaseq <- merge(rnaseq, genes, by.x = "gene_id", by.y = "gene_id")
rnaseq <- rnaseq[order(sampleid)]

# . Get Gene annotation from Gencode
path <- "/home/bioinf/bhklab/jermiah/psets/PharmacoSet-Pipelines/GDSC/metadata/human/GRCh38_v44/annotation.gtf"
dsGencode <- rtracklayer::import(path)

# remove gene_id version from dsGencode 
dsGencode$gene_id <- gsub("\\.\\d+$", "", dsGencode$gene_id)

# subset rnaseq data to only include genes in dsGencode
# lose ~ 300 genes
rnaseq <- rnaseq[ensembl_gene_id %in% dsGencode$gene_id]

outputFiles <- list(
    "rnaseq" = rnaseq,
    "GRanges" = dsGencode
)

dir.create(dirname(OUTPUT$preprocessed), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputFiles, file = OUTPUT$preprocessed)


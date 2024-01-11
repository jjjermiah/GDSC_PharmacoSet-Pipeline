
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image()
}

# 0.1 read metadata
# -----------------
metadata <- qs::qread(INPUT$metadata, nthreads = THREADS)
sample <- metadata$sample
geneAnnot <- metadata$geneAnnot
# 0.2 read gene fusions
# ---------------------
zipDir <- dirname(INPUT$gene_fusions)
unzipDir <- file.path(dirname(zipDir), "gene_fusions")
unzip(INPUT$gene_fusions, exdir = unzipDir)
dt <- data.table::fread(file.path(unzipDir, list.files(unzipDir)))

fusion_dt <- dt[model_id %in% sample[, model_id]]

# names(fusion_dt)
#  [1] "chr_3prime"         "chr_5prime"         "dataset"           
#  [4] "dataset_id"         "gene_id_3prime"     "gene_id_5prime"    
#  [7] "gene_symbol_3prime" "gene_symbol_5prime" "in_cosmic_fusions" 
# [10] "in_frame"           "in_patient"         "model_id"          
# [13] "model_name"         "tissue"   

# I have no idea what to do with this data...


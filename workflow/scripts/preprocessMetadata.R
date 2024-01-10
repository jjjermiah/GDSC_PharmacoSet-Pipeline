
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image()
}

sample <- readxl::read_excel(INPUT$sampleMetadata, sheet = 1, col_names = TRUE, na = "NA")
sample <- data.table::as.data.table(sample)

treatment <- data.table::fread(INPUT$treatmentMetadata, header = TRUE)

cellmodels <- data.table::fread(INPUT$cellModelPassportsMetadata, header = TRUE)
geneAnnot <- data.table::fread(INPUT$cellmodelpassportsGeneAnnotation, header = TRUE)

gdsc_cmp <- cellmodels[model_name %in% sample$`Sample Name`]
cols <- c(
    'model_id', 'sample_id', 'model_name', 'tissue', 'cancer_type', 
    'cancer_type_ncit_id', 'sample_site', 'COSMIC_ID', 'BROAD_ID', 'CCLE_ID',
    'gender', 'ethnicity', 'age_at_sampling', 'clinical_staging'
)

gdsc_cmp <- gdsc_cmp[, ..cols]

sampleMetadata <- merge(sample[,.(`Sample Name`)], gdsc_cmp, by.x = "Sample Name", by.y = "model_name", all.x = TRUE)
# rename Sample Name to sample_id
data.table::setnames(sampleMetadata, "Sample Name", "sampleid")

outputFiles <- list(
    "sample" = sampleMetadata,
    "treatment" = treatment,
    "geneAnnot" = geneAnnot
)

dir.create(dirname(OUTPUT$metadata), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputFiles, file = OUTPUT$metadata)

## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image()
}
exit


# metadata
sample <- readxl::read_excel(INPUT$sampleMetadata, sheet = 1, col_names = TRUE, na = "NA")
sample <- data.table::as.data.table(sample)

cellmodels <- data.table::fread(INPUT$cmp_Metadata, header = TRUE)


################################################################################
# rnaseq data
allDir <- paste0(dirname(INPUT$all), "/all")
dir.create(allDir, recursive = TRUE, showWarnings = FALSE)
unzip(INPUT$all, exdir = allDir)
list.files(allDir)

# for each csv in allDir, read in the csv into a data.table and add a source col
all_data_CSVs <- lapply(list.files(allDir), function(file){
    message("Reading in ", file)
    data <- data.table::fread(file.path(allDir, file), header = TRUE, sep = ",")
    data$file <- file
    data
})
names(all_data_CSVs) <- list.files(allDir)

################################################################################
## Subset data


gdsc_cmp <- cellmodels[model_name %in% sample$`Sample Name`]
                           

cols <- c(
    'model_id', 'sample_id', 'model_name', 'tissue', 'cancer_type', 
    'cancer_type_ncit_id', 'sample_site', 'COSMIC_ID', 'BROAD_ID', 'CCLE_ID',
    'gender', 'ethnicity', 'age_at_sampling', 'clinical_staging'
)

gdsc_cmp <- gdsc_cmp[, ..cols]

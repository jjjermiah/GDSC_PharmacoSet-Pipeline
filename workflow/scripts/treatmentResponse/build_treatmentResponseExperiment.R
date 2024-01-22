
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}
pak::pkg_install("bhklab/CoreGx")

suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(PharmacoGx))
suppressPackageStartupMessages(library(GenomicRanges))

# devtools::install_github("bhklab/CoreGx", quiet = TRUE)
# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(
      log4r::file_appender(LOGFILE, append = TRUE),
      log4r::console_appender()
    )
)

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting build_treatmentResponseExperiment.R\n")

# 0.2 Load Data
# -------------

input <- qs::qread(INPUT$preprocessed, nthreads = THREADS)

rawdata <- input$rawdata
fitted_data <- input$fitted_data
normData <- input$normData
metadata <- input$metadata

info(paste0(c("rawdata:\n", capture.output(str(rawdata))), collapse = "\n"))
info(paste0(c("fitted_data:\n", capture.output(str(fitted_data))), collapse = "\n"))
info(paste0(c("normData:\n", capture.output(str(normData))), collapse = "\n"))
info(paste0(c("metadata:\n", capture.output(str(metadata))), collapse = "\n"))

# 1.0 Subsetting Data
# -------------------
info(paste0("Subsetting normData to not include missing values"))
proc_data <- unique(normData[
  !is.na(normalized_intensity) & 
    !is.na(DRUG_NAME) & 
    !is.na(CONC) & 
    !is.na(CELL_LINE_NAME), 
  .(CELL_LINE_NAME, DRUG_NAME, CONC, Viability = normalized_intensity)])


subsetted_samples <- unique(proc_data$CELL_LINE_NAME)

info(paste0("Subsetting proc_data to only have n = ", length(subsetted_samples), " samples"))
subset_procdata <- proc_data[CELL_LINE_NAME %in% subsetted_samples,]
subset_procdata <- subset_procdata[order(CELL_LINE_NAME, DRUG_NAME, CONC)]


subset_fitted_data <- unique(
  fitted_data[
    CELL_LINE_NAME %in% subsetted_samples &
      DRUG_NAME %in% subset_procdata$DRUG_NAME, 
    .(CELL_LINE_NAME, DRUG_NAME, LN_IC50, AUC, RMSE, Z_SCORE)])


info(paste0("Loading TREDataMapper"))
TREDataMapper <- CoreGx::TREDataMapper(rawdata=subset_procdata)
{
  CoreGx::rowDataMap(TREDataMapper) <- list(
    id_columns = (c("DRUG_NAME", "CONC")),
    mapped_columns = c())

  CoreGx::colDataMap(TREDataMapper) <- list(
    id_columns = c("CELL_LINE_NAME"),
    mapped_columns = c())

  CoreGx::assayMap(TREDataMapper) <- list(
    raw = list(
      c("DRUG_NAME", "CONC", "CELL_LINE_NAME"),
      c("Viability")))

  info(paste0(capture.output(TREDataMapper), collapse = "\n"))
  # devtools::load_all("/home/bioinf/bhklab/jermiah/Bioconductor/CoreGx/")
  
  gdsc_tre <- CoreGx::metaConstruct(TREDataMapper)
  tre <- gdsc_tre
}


published_profiles <- subset_fitted_data[
  CELL_LINE_NAME %in% tre$raw$CELL_LINE_NAME & 
    DRUG_NAME %in% tre$raw$DRUG_NAME,]

info(paste0("Adding published profiles to tre"))
assay(tre, "profiles_published") <- published_profiles
info(paste0(capture.output(tre), collapse = "\n"))

info(paste0("Saving tre to: ", OUTPUT$tre))
qs::qsave(tre, OUTPUT$tre, nthreads = THREADS)



# tre_list <- lapply(unique(proc_data$CELL_LINE_NAME)[1:2], function(x){
#   subset_rawDT <- proc_data[CELL_LINE_NAME == x]
#   # # (rawControl <- rawdata[TAG.1 == control.column, .(INTENSITY)])
#   TREdataMapper <- CoreGx::TREDataMapper(rawdata=subset_rawDT)

#   CoreGx::rowDataMap(TREdataMapper) <- list(
#     id_columns = (c("DRUG_NAME", "CONC")),
#     mapped_columns = c()
#   )

#   CoreGx::colDataMap(TREdataMapper) <- list(
#     id_columns = c("CELL_LINE_NAME"),
#     mapped_columns = c()
#   )

#   CoreGx::assayMap(TREdataMapper) <- list(
#     raw = list(
#       c("DRUG_NAME", "CONC"),
#       c("Viability")
#     )
#   )

#   gdsc_tre <- CoreGx::metaConstruct(TREdataMapper)
#   gdsc_tre
# })
# tre_list

# rawdata_2 <- removeMissingDrugs(rawdata_1)

# normalized_dt <- normalizeData(rawdata_2, trim=T, neg_control="NC-1", pos_control="B")
# setConcsForNlme(normalized_dt, group_conc_ranges = F)


# rawdata[, TAG.1 := sapply(strsplit(TAG, split = "-"), function(x) {
#   return(x[[1]])
# })]

# ## Taking mean if more than 5 measurements, otherwise median
# control.column <- c("NC-0", "NC-1")

# # Calculate viability based on intensity values
# rawdata[, Viability := {
#   control <- .SD[TAG == ..control.column, ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
#   background <- .SD[TAG.1 == "B", ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
#   Viability <- (INTENSITY - background) / (control - background) * 100
#   list(Viability = Viability)
# }, .(MASTER_CELL_ID, BARCODE)]

# 



### INVESTIGATE RAWDATA ###
# [1] "RESEARCH_PROJECT" "BARCODE"          "SCAN_ID"          "DATE_CREATED"    
#  [5] "SCAN_DATE"        "CELL_ID"          "MASTER_CELL_ID"   "COSMIC_ID"       
#  [9] "CELL_LINE_NAME"   "SANGER_MODEL_ID"  "SEEDING_DENSITY"  "DRUGSET_ID"      
# [13] "ASSAY"            "DURATION"         "POSITION"         "TAG"             
# [17] "DRUG_ID"          "CONC"             "INTENSITY"        "TAG.1"           

# # get count of rows for each cell line sorted
# rawdata[, .N, by = CELL_LINE_NAME][order(-N)]

# # get count of rows for each drug sorted
# rawdata[, .N, by = DRUG_ID][order(-N)]

# # get count of each "TAG" sorted
# rawdata[, .N, by = TAG][order(-N)]

# #  sort this: rawdata[, .N, by = TAG][order(-N)]$TAG
# rawdata[, .N, by = CONC][order(-CONC)][1:10]

# Number unique "SEEDING_DENSITY" and "BARCODE" for each "DRUGSET_ID"




# tre_list <- lapply(unique(rawDT$CELL_LINE_NAME)[1:2], function(x){
#   subset_rawDT <- rawDT[CELL_LINE_NAME == x]
#   # # (rawControl <- rawdata[TAG.1 == control.column, .(INTENSITY)])
#   TREdataMapper <- CoreGx::TREDataMapper(rawdata=subset_rawDT)

#   CoreGx::rowDataMap(TREdataMapper) <- list(
#     id_columns = (c("DRUG_NAME", "BARCODE", "CONC", "SEEDING_DENSITY", "ASSAY", "DURATION")),
#     mapped_columns = c()
#   )

#   CoreGx::colDataMap(TREdataMapper) <- list(
#     id_columns = c("CELL_LINE_NAME"),
#     mapped_columns = c()
#   )

#   CoreGx::assayMap(TREdataMapper) <- list(
#     raw = list(
#       c("DRUG_NAME", "BARCODE", "CONC", "SEEDING_DENSITY", "ASSAY", "DURATION", "CELL_LINE_NAME"),
#       c("Viability")
#     )
#   )

#   gdsc_tre <- CoreGx::metaConstruct(TREdataMapper)
#   gdsc_tre
# })
# tre_list

# lapply(tre_list, function(x){
#   rowData_ <- CoreGx::rowData(x)
#   colData_ <- CoreGx::colData(x)
#   x@.intern
# })

# # guess <- CoreGx::guessMapping(TREdataMapper, groups, subsets)
# # guess




# tre_list <- lapply(unique(rawDT$CELL_LINE_NAME)[1:2], function(x){
#   subset_rawDT <- rawDT[CELL_LINE_NAME == x]
#   # # (rawControl <- rawdata[TAG.1 == control.column, .(INTENSITY)])
#   TREdataMapper <- CoreGx::TREDataMapper(rawdata=subset_rawDT)

#   CoreGx::rowDataMap(TREdataMapper) <- list(
#     id_columns = (c("DRUG_NAME", "BARCODE", "CONC", "SEEDING_DENSITY", "ASSAY", "DURATION")),
#     mapped_columns = c()
#   )

#   CoreGx::colDataMap(TREdataMapper) <- list(
#     id_columns = c("CELL_LINE_NAME"),
#     mapped_columns = c()
#   )

#   CoreGx::assayMap(TREdataMapper) <- list(
#     raw = list(
#       c("DRUG_NAME", "BARCODE", "CONC", "SEEDING_DENSITY", "ASSAY", "DURATION", "CELL_LINE_NAME"),
#       c("Viability")
#     )
#   )

#   gdsc_tre <- CoreGx::metaConstruct(TREdataMapper)
#   gdsc_tre
# })
# tre_list

# lapply(tre_list, function(x){
#   rowData_ <- CoreGx::rowData(x)
#   colData_ <- CoreGx::colData(x)
#   x@.intern
# })

# # guess <- CoreGx::guessMapping(TREdataMapper, groups, subsets)
# # guess



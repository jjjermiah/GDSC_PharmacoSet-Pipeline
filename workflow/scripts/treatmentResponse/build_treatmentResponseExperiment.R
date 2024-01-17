
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    # save.image()
}

suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(CoreGx))
suppressPackageStartupMessages(library(GenomicRanges))


# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting build_treatmentResponseExperiment.R\n")

# 0.2 read in treatmentResponse data
# ----------------------------------
info(paste0("Reading in: ", INPUT$rawdata))
rawdata <- data.table::fread(INPUT$rawdata)
info(paste0(capture.output(str(rawdata)), collapse = "\n"))
info(paste0("Number of unique cell lines: ", uniqueN(rawdata$CELL_LINE_NAME)))
info(paste0("Number of unique drugs: ", uniqueN(rawdata$DRUG_ID), "\n\n"))

info(paste0("Reading in: ", INPUT$processed))
data <- readxl::read_xlsx(INPUT$processed)
data <- data.table::as.data.table(data)
info(paste0(capture.output(str(data)), collapse = "\n"))
info(paste0("Number of unique cell lines: ", uniqueN(data$CELL_LINE_NAME)))
info(paste0("Number of unique drugs: ", uniqueN(data$DRUG_ID), "\n\n"))

# 0.3 read in metadata 
# --------------------
info(paste0("Reading in: ", INPUT$metadata))
metadata <- qs::qread(INPUT$metadata)
treatment <- metadata$treatment

# 1.0 Subsetting Data
# -------------------

rawdata[, TAG.1 := sapply(strsplit(TAG, split = "-"), function(x) {
  return(x[[1]])
})]


## Taking mean if more than 5 measurements, otherwise median
control.column <- c("NC-0", "NC-1")
# Calculate viability based on intensity values
# Parameters:
#   - rawdata: The input data table
#   - control.column: The column name for the control condition
# Returns:
#   - A modified data table with a new column "Viability" containing the calculated viability values

# Calculate viability based on intensity values
rawdata[, Viability := {
  control <- .SD[TAG == ..control.column, ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
  background <- .SD[TAG.1 == "B", ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
  Viability <- (INTENSITY - background) / (control - background) * 100
  list(Viability = Viability)
}, .(MASTER_CELL_ID, BARCODE)]

rawdata <- merge(rawdata, treatment[,.(DRUG_ID, DRUG_NAME)], by = "DRUG_ID", all.x = T)

# FOR NOW, REMOVE ALL ROWS WHERE DRUG_NAME HAS A PUNCTUATION OR SPACE
rawdata <- rawdata[!grepl("[[:punct:]]", DRUG_NAME) & !grepl(" ", DRUG_NAME)]
# TODO:: WHEN MAPPING TREATMENTS, CLEAN NAMES!!


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

###############################################
# 2.0 Constructing the Experiment
###############################################
rawDT <- rawdata[!is.na(DRUG_NAME), .(
  CELL_LINE_NAME, DRUG_NAME, 
  BARCODE, SEEDING_DENSITY, ASSAY, DURATION, CONC, Viability)]

# create a new column using SEEDING_DENSITY, BARCODE, DURATION, ASSAY
rawDT[, EXPERIMENT := paste("seed", SEEDING_DENSITY, "barcode", BARCODE, "dur", DURATION, "assay", ASSAY, sep = "_")]
# drop the constituent columns
rawDT[, c("SEEDING_DENSITY", "BARCODE", "DURATION", "ASSAY") := NULL]

runtre <- function(){
  subset_rawDT <- 
    rawDT[CELL_LINE_NAME %in% unique(rawDT$CELL_LINE_NAME)[1:250],]
  # # (rawControl <- rawdata[TAG.1 == control.column, .(INTENSITY)])
  TREdataMapper <- CoreGx::TREDataMapper(rawdata=subset_rawDT)

  CoreGx::rowDataMap(TREdataMapper) <- list(
    # id_columns = c("DRUG_NAME",  "CONC", "SEEDING_DENSITY", "BARCODE", "ASSAY", "DURATION"),
    id_columns = c("DRUG_NAME", "EXPERIMENT", "CONC"),
    mapped_columns = c()
  )

  CoreGx::colDataMap(TREdataMapper) <- list(
    id_columns = c("CELL_LINE_NAME"),
    mapped_columns = c()
  )

  CoreGx::assayMap(TREdataMapper) <- list(
    raw = list(
      # c("DRUG_NAME",  "CELL_LINE_NAME", "CONC", "SEEDING_DENSITY", "BARCODE", "ASSAY", "DURATION"),
      c("DRUG_NAME", "EXPERIMENT", "CELL_LINE_NAME", "CONC"),
      c("Viability")
    )
  )

  gdsc_tre <- CoreGx::metaConstruct(TREdataMapper)
  gdsc_tre
}
devtools::load_all("/home/bioinf/bhklab/jermiah/Bioconductor/CoreGx/")
tre <- runtre()

qs::qsave(tre, OUTPUT$tre, nthreads = THREADS)


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


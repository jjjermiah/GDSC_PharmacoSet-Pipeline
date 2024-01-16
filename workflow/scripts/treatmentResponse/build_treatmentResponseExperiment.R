
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
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

rawdata

## Taking mean if more than 5 measurements, otherwise median
control.column <- "NC-1"
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

# NEED TO DO TRE DATA MAPPING! 


# groups <- list(
#   justDrugs=c('drug1id', 'drug2id'),
#   drugsAndDoses=c('drug1id', 'drug2id', 'drug1dose', 'drug2dose'),
#   justCells=c('cellid'),
#   cellsAndBatches=c('cellid', 'batchid'),
#   assays1=c('drug1id', 'drug2id', 'cellid'),
#   assays2=c('drug1id', 'drug2id', 'drug1dose', 'drug2dose', 'cellid', 'batchid')
# )


colnames(rawdata)
# [1] "RESEARCH_PROJECT" "BARCODE"          "SCAN_ID"          "DATE_CREATED"    
#  [5] "SCAN_DATE"        "CELL_ID"          "MASTER_CELL_ID"   "COSMIC_ID"       
#  [9] "CELL_LINE_NAME"   "SANGER_MODEL_ID"  "SEEDING_DENSITY"  "DRUGSET_ID"      
# [13] "ASSAY"            "DURATION"         "POSITION"         "TAG"             
# [17] "DRUG_ID"          "CONC"             "INTENSITY"        "TAG.1"           

# get count of rows for each cell line sorted
rawdata[, .N, by = CELL_LINE_NAME][order(-N)]

# get count of rows for each drug sorted
rawdata[, .N, by = DRUG_ID][order(-N)]

# get count of each "TAG" sorted
rawdata[, .N, by = TAG][order(-N)]

rawDT <- rawdata[, .(
  CELL_LINE_NAME, SANGER_MODEL_ID, 
  DRUG_ID, SEEDING_DENSITY, ASSAY, DURATION,
  TAG.1, CONC, INTENSITY)]
rawDT
unique(rawDT)
groups <- list(
  justDrugs = c("DRUG_ID"),
)


# (rawControl <- rawdata[TAG.1 == control.column, .(INTENSITY)])
